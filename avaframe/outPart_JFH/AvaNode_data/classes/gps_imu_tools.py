# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 14:25:16 2021

@author: neuhauser
"""

import numpy as np
import pandas as pd


def transform_gps_to_local(path):
    
    from geopy import distance
    
    # For Emlid Data
    if 'Emlid' in path:
        df = pd.read_csv(path, header=9, sep='  ')
        df = df.rename({'GPST': 'latitude[deg]'}, axis=1)
        df = df.rename({'Unnamed: 2': 'longitude[deg]'}, axis=1)
        df = df.rename({'Unnamed: 3': 'height[m]'}, axis=1)
        dt = 0.2
    #For Node Data
    else:
        df = pd.read_csv(path)
        dt = 0.1
        if (df['longitude[deg]'] > 90).any():
            df['longitude[deg]'] = df['longitude[deg]'] / 10000000
        if (df['latitude[deg]'] > 180).any():
            df['latitude[deg]'] = df['latitude[deg]'] / 10000000
            print("Converted Lat / Long ... ")
        df['height[m]'] = df['height[mm]'] / 1000
    print(df.head()) #7 columns, including the Date. 

    # Calculating the distances between the GPS coordinates
    x = [0]
    y = [0]
    z = [0]
    
    # lat = North/south // long = East/West
    for i in range(1, df.shape[0]):
        delta_x_0 = (df["latitude[deg]"][0], df["longitude[deg]"][0])
        delta_x_1 = (df["latitude[deg]"][i], df["longitude[deg]"][0])
        dx = distance.distance(delta_x_0, delta_x_1).m
        if df["latitude[deg]"][i] > df["latitude[deg]"][0]:
            x.append(dx)
        else:
            x.append(-dx)

        delta_y_0 = (df["latitude[deg]"][0], df["longitude[deg]"][0])
        delta_y_1 = (df["latitude[deg]"][0], df["longitude[deg]"][i])
        dy = distance.distance(delta_y_0, delta_y_1).m
        if df["longitude[deg]"][i] > df["longitude[deg]"][0]:
            y.append(dy)
        else:
            y.append(-dy)
        
        df['height[m]'] = df['height[mm]'] / 1000
        delta_z = df['height[m]'][0] - df['height[m]'][i]
        z.append(delta_z)
        
    gps_vx = np.zeros(len(x)-1)
    gps_vy = np.zeros(len(x)-1)
    gps_vz = np.zeros(len(x)-1)

    for i in range(len(x)-1):
        gps_vx[i] = (x[i+1] - x[i]) / dt
        gps_vy[i] = (y[i+1] - y[i]) / dt
        gps_vz[i] = (z[i+1] - z[i]) / dt
        
    gps_v_tot = np.sqrt(gps_vx ** 2 + gps_vy ** 2 + gps_vz ** 2)
    time = np.ones(len(gps_vx))
    for idx in range(len(time)):
        time[idx] = dt * idx
        
    time = np.ones(len(x))
    for idx in range(len(time)):
        time[idx] = dt * idx
    return time, x, y, z, gps_v_tot

def get_coord_transform(source_epsg, target_epsg):
    
    from osgeo import osr
    
    '''
    Creates an OGR-framework coordinate transformation for use in projecting
    coordinates to a new coordinate reference system (CRS). Used as, e.g.:
        transform = get_coord_transform(source_epsg, target_epsg)
        transform.TransformPoint(x, y)
    Arguments:
        source_epsg     The EPSG code for the source CRS
        target_epsg     The EPSG code for the target CRS
    '''
    # Develop a coordinate transformation, if desired
    transform = None
    source_ref = osr.SpatialReference()
    target_ref = osr.SpatialReference()
    source_ref.ImportFromEPSG(source_epsg)
    target_ref.ImportFromEPSG(target_epsg)
    return osr.CoordinateTransformation(source_ref, target_ref) 

def transform_imu_to_global(path, lat, long, altitude, date):
    
    from magnetic_field_calculator import MagneticFieldCalculator #Python API for British Geological Survey magnetic field calculator
    
    '''
    This functions is used to determine the magnetic declination and rotate the 
    given path with this angle.
    Input: input_path_csv_with_x_and_y, lat, long, altitude[km], date(yyyy-mm-dd)
    output: csv with new Fields, North and East are rotated and in the right projection.
    '''

    # Get the declination of the magnetic north pole to the geographic noth pole
    calculator = MagneticFieldCalculator()

    calculator = MagneticFieldCalculator(
        model='igrf',
        revision='current',
        custom_url='http://geomag.bgs.ac.uk/'
    )

    result = calculator.calculate(
        latitude=lat,
        longitude=long,
        altitude=altitude,
        date=date
    )

    declination = result['field-value']['declination']['value']
    print("Deklination = {} °".format(declination))

    theta = np.radians(declination)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))

    transform = get_coord_transform(4326, 31287) #transform in cartesian coordinates to add the path
    north, east, z = transform.TransformPoint(lat, long)

    col = ['Time', 'x', 'y', 'z']
    inp = pd.read_csv(path, delimiter = ",", header=None, skiprows=1, names=col)

    inp['x_rot'] = inp['x']
    inp['y_rot'] = inp['y']

    for i in range(len(inp['x'])):
        vec = np.array((inp['x'].iloc[i], inp['y'].iloc[i]))
        x_rot, y_rot = R @ vec
        inp['x_rot'].iloc[i] = x_rot
        inp['y_rot'].iloc[i] = y_rot

    inp['North'], inp['East'] = north + inp['x_rot'], east + inp['y_rot']

    inp.to_csv(path[:-4] + '_wgs84_magdeclination.csv')
    
def gps_to_mercator(gps_class):
    
    north = np.zeros_like(gps_class.data['latitude[deg]'])
    east = np.zeros_like(north)
    z = np.zeros_like(north)
    
    transform = get_coord_transform(4326, 31254) #transform in cartesian coordinates to add the path
    for i in range(len(gps_class.data['latitude[deg]'])):
        north[i], east[i], z[i] = transform.TransformPoint(gps_class.data['latitude[deg]'].values[i], gps_class.data['longitude[deg]'].values[i], gps_class.data['height[m]'].values[i])
    return north, east, z
