# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 17:00:48 2021

@author: neuhauser
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as dates

from classes.gps_imu_tools import transform_gps_to_local, gps_to_mercator


class GPSData(object):
    
    def __init__(self):

        self.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data'
        self.file = '220222_C10_avalanche_GPS.txt'
        self.work_dir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data'
        self.data = 0
        self.data_pos = 0
        self.is_synced = False
        self.sync = 0
        self.sync_counter = []
        self.gravi = 9.80665
        self.temp_bool = False
        self.vel_bool = False
        self.calc_pos_vel_bool = False
        self.pDop_bool = False


    def read_data(self):
        
        self.work_dir, self.file = os.path.split(self.path)
        
        print("Reading GPS File...")
        #Check how long header is
        skipheader = -1
        header = True
        with open(self.path) as f:
            while header:
                if f.readline()[0] == '#':
                    skipheader += 1
                else:
                    header = False       
        
        data = pd.read_csv(self.path, index_col=None, low_memory=False, skiprows=skipheader)
        print("Done!")
        print("Searching for Sync Keys...")
        # Search for sync keys:
        print("Sync Pulse...")
        if data.columns.str.contains('rising edge sync signal').any():
            sync_key = data.columns.str.contains('sync')  # Search header for sync
        elif data.columns.str.contains('gpsFlag').any():
            sync_key = data.columns.str.contains('gpsFlag')  
        elif data.columns.str.contains('sync signal').any():
            sync_key = data.columns.str.contains('sync')  # Search header for sync
        sync = data.columns[np.where(sync_key)[0]]
        sync = sync[0]
        self.sync = data[sync]
        print("OK")
        
        # Search for iTow
        print("iTow...")
        if data.columns.str.contains('ms').any():
            itow_key = data.columns.str.contains('ms')
        itow = data.columns[np.where(itow_key)[0]]
        itow = itow[0]
        self.itow = itow
        print("OK")
        
        #Search for UTC key:
        print("UTC...")
        if data.columns.str.contains('utc').any():
            utc_key = data.columns.str.contains('utc')
        if data.columns.str.contains('UTC ').any():
            utc_key = data.columns.str.contains('UTC')
        utc = data.columns[np.where(utc_key)[0]]
        utc = utc[0]
        self.utc = utc
        print("OK")
        
        # Start when we have fix Type 3 == 3D Position
        # Search for fix Type
        print("Fix Type...")
        if data.columns.str.contains('Type').any():
            type_key = data.columns.str.contains('Type')
            typ = data.columns[np.where(type_key)]
            typ = typ[0]
        print("OK")
        
        #Searching for pDOP
        print("pDOP...")
        if data.columns.str.contains('pDop').any():
            self.pDop_bool = True
        print("OK")
        
        #Searching for temp data
        print("Temperature...")
        if data.columns.str.contains('ambient').any():
            self.temp_bool = True
        print("OK")
        
        print("Velocity...")
        if data.columns.str.contains('velE').any():
            #v_tot = np.sqrt((self.data['velN[mm/s]']/1000)**2 + (self.data['velE[mm/s]']/1000)**2 + (self.data['velD[mm/s]']/1000)**2)
            #self.data_pos['v_total'] = v_tot
            print("GNSS File contains velocities.")
            self.vel_bool = True
        else:
            print("No velocities in this file!")
        
        #Add microeconds to date
        a = data[self.itow] % 1000
        a = a.astype(str)

        newdate = data[self.utc] + ":" + a
        data["date[ms]"] = newdate
        # Convert date column to date object
        #data[self.utc] = pd.to_datetime(data[utc], format='%Y:%m:%d:%H:%M:%S')
        data['timestamp'] = pd.to_datetime(data["date[ms]"], format='%Y:%m:%d:%H:%M:%S:%f')

        
        # Convert GPS Data
        if (data['longitude[deg]'] > 90).any():
            data['longitude[deg]'] = data['longitude[deg]'] / 10000000
        if (data['latitude[deg]'] > 180).any():
            data['latitude[deg]'] = data['latitude[deg]'] / 10000000
            print("Converted Lat / Long ... ")
        data['height[m]'] = data['height[mm]'] / 1000
               
        
        # Save sorted dataframe to preprocessed file:
        #pos.to_csv(self.path[:-4] + "_preprocessed.csv", index=False)
        try:
            pos = data.loc[data[typ] == 3]
            self.first_position = pos.iloc[0, :]
            self.data_pos = pos
        except:
            print("GNSS File contains no valuable data!")
        self.data = data
        if self.vel_bool:
            v_tot = np.sqrt((self.data['velN[mm/s]']/1000)**2 + (self.data['velE[mm/s]']/1000)**2 + (self.data['velD[mm/s]']/1000)**2)
            self.data['v_total'] = v_tot
            self.data_pos['v_total'] = v_tot
        
        if self.pDop_bool:
            self.data['pDop[]'] *= 0.01
            self.data_pos['pDop[]'] *= 0.01
        print("GPS finished!")
        
    def sync_gps(self):
        
        startvalue = self.sync[0]
        gps_peaks = 0
        gps_sync_counter = []

        counter_even = 0
        counter_odd = 1

        for value in self.sync:
            if value == startvalue:
                if value == 0:
                    gps_sync_counter.append(counter_even)
                if value == 1:
                    gps_sync_counter.append(counter_odd)
            if value != startvalue:
                gps_peaks += 1
                if value == 0:
                    counter_even += 2
                    gps_sync_counter.append(counter_even)
                if value == 1:
                    if gps_peaks > 1:
                        counter_odd += 2
                    gps_sync_counter.append(counter_odd)
                startvalue = value
        self.sync_counter = gps_sync_counter
        self.is_synced = True
        print("GPS has {} peaks!".format(gps_peaks))
        
    def calc_pos_vel(self):
        
        from geopy import distance
        
        # Calculating the distances between the GPS coordinates
        x = [0]
        y = [0]
        z = [0]
        
        dt = 0.1 # 10 Hz
        
        # lat = North/south // long = East/West
        for i in range(1, self.data_pos.shape[0]):
            delta_x_0 = (self.data_pos["latitude[deg]"].iloc[0], self.data_pos["longitude[deg]"].iloc[0])
            delta_x_1 = (self.data_pos["latitude[deg]"].iloc[i], self.data_pos["longitude[deg]"].iloc[0])
            dx = distance.distance(delta_x_0, delta_x_1).m
            if self.data_pos["latitude[deg]"].iloc[i] > self.data_pos["latitude[deg]"].iloc[0]:
                x.append(dx)
            else:
                x.append(-dx)

            delta_y_0 = (self.data_pos["latitude[deg]"].iloc[0], self.data_pos["longitude[deg]"].iloc[0])
            delta_y_1 = (self.data_pos["latitude[deg]"].iloc[0], self.data_pos["longitude[deg]"].iloc[i])
            dy = distance.distance(delta_y_0, delta_y_1).m
            if self.data_pos["longitude[deg]"].iloc[i] > self.data_pos["longitude[deg]"].iloc[0]:
                y.append(dy)
            else:
                y.append(-dy)
            
            self.data_pos['height[m]'] = self.data_pos['height[mm]'] / 1000
            delta_z = self.data_pos['height[m]'].iloc[0] - self.data_pos['height[m]'].iloc[i]
            z.append(delta_z)
            
        gps_vx = np.zeros(len(x))
        gps_vy = np.zeros(len(x))
        gps_vz = np.zeros(len(x))

        for i in range(len(x)-1):
            gps_vx[i+1] = (x[i+1] - x[i]) / dt
            gps_vy[i+1] = (y[i+1] - y[i]) / dt
            gps_vz[i+1] = (z[i+1] - z[i]) / dt
            
        gps_v_tot = np.sqrt(gps_vx ** 2 + gps_vy ** 2 + gps_vz ** 2)
        self.data_pos['v_total_pos'] = gps_v_tot
        self.data_pos['x_dist'] = x
        self.data_pos['y_dist'] = y
        self.data_pos['z_dist'] = z
        self.data_pos['tot_dist'] = np.linalg.norm((x,y,z), axis=0)
                
        
    def plot_vel(self, saveplot=False, plotindex=False, plot_pos_vel=False, starttime=0, endtime=0):
        
        if not self.calc_pos_vel_bool:
            self.calc_pos_vel()
            self.calc_pos_vel_bool = True
        
        date = dates.date2num(self.data_pos['timestamp'])
        # Define the date format
        date_form = dates.DateFormatter("%H:%M:%S")
        
        if starttime == 0:
            starttime = self.data_pos['timestamp'].iloc[0]
        if endtime == 0:
            endtime = self.data_pos['timestamp'].iloc[-1]
        
        linewidth = 0.5
        fig, ax = plt.subplots()
        ax.xaxis.set_major_formatter(date_form)
        fig.autofmt_xdate(rotation=45)
        
        #v_tot = np.sqrt((self.data['velN[mm/s]']/1000)**2 + (self.data['velE[mm/s]']/1000)**2 + (self.data['velD[mm/s]']/1000)**2)
        if self.vel_bool:
            #v_tot = np.linalg.norm([self.data_pos["velD[mm/s]"].values, self.data_pos["velE[mm/s]"].values, self.data_pos["velN[mm/s]"].values], axis=0)/1000
            #self.data_pos['v_total'] = v_tot

            if plotindex:
                ax.plot(self.data_pos['v_total'], color='b', linewidth=linewidth, label='Doppler Effect')
                ax.set_xlabel('Index')
            else:
                ax.plot(date, self.data_pos['v_total'], color='b', linewidth=linewidth, label='Doppler Effect')
                ax.set_xlabel('Time [UTC]')
        
        if plot_pos_vel and plotindex:
            ax.plot(self.data_pos['v_total_pos'], color='r', linewidth=linewidth, label='From Pos. Data')

        if plot_pos_vel:
            ax.plot(date, self.data_pos['v_total_pos'], color='r', linewidth=linewidth, label='From Pos. Data')
        
        if self.vel_bool and plot_pos_vel:
             ax.legend()
        ax.set_xlim([starttime, endtime])
        ax.set_ylim([-1, 30])
        ax.set_ylabel('Velocity [m/s]')
        ax.set_title(self.file[:-4] + " GNSS Velocities")
        ax.grid()
        
        if saveplot:
            fig.savefig(self.work_dir + '/' + self.file[:-4] + '_GNSS_vel.png', dpi=300)
            
    def export_in_radar_cs(self):
        
        from osgeo import gdal
        from geopy import distance
        
        # Calc v_tot
        if self.vel_bool:
            v_tot = np.linalg.norm([self.data_pos["velD[mm/s]"].values, self.data_pos["velE[mm/s]"].values, self.data_pos["velN[mm/s]"].values], axis=0)/1000
            self.data_pos['v_total'] = v_tot
        # Calc Distance to radar
        radar_pos = (47.306584, 11.380406) #lat, long, z
        radar_z = 1900.881
        # Export for Radar Comparison
        filepath = r"C:\Users\neuhauser\OneDrive - Bundesforschungszentrum fuer Wald\20070-AvaRange-data\GIS\Nordkette\Nordkette_wgs84.tif"
        # Open the file:
        nordkette_raster = gdal.Open(filepath)

        # Check type of the variable 'raster'
        type(nordkette_raster)
        # Projection
        nordkette_raster.GetProjection()
        band = nordkette_raster.GetRasterBand(1)

        cols = nordkette_raster.RasterXSize
        rows = nordkette_raster.RasterYSize

        transform = nordkette_raster.GetGeoTransform()

        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = -transform[5]

        data = band.ReadAsArray(0, 0, cols, rows)
        point_list = []
        for i in range(len(self.data["latitude[deg]"])):
            point_list.append([self.data["longitude[deg]"].iloc[i], self.data["latitude[deg]"].iloc[i]])

        z = []
        for point in point_list:
            col = int((point[0] - xOrigin) / pixelWidth)
            row = int((yOrigin - point[1] ) / pixelHeight)
            z.append(data[row][col])
            #print(row, col, data[row][col])
        self.data["z_raster"] = z
        
        dist = []
        for i in range(len(self.data["latitude[deg]"])):
            pos = (self.data["latitude[deg]"].iloc[i], self.data["longitude[deg]"].iloc[i])
            horizontal_distance = distance.distance(pos, radar_pos).m
            #print(horizontal_distance)
            vertical_distance = self.data["z_raster"].iloc[i] - radar_z
            dist.append(np.linalg.norm([horizontal_distance, vertical_distance]))
        #print(dist)
        # Export to csv
        self.export_df = pd.DataFrame(list(zip(self.data_pos["timestamp"].to_list(), dist, self.data_pos['v_total'].to_list())), columns=["Time", "Distance_to_radar[m]", "velocity[m/s]"])
        self.export_df.to_csv("220222_Cxx_radarCS.csv", index=None) #ToDo: Auto for date and export!!!
        
        
    def check_velocity(self):
        """Check if the integration of the velocity corresponds to the length from the position data"""
        
        time, x, y, z, gps_v_tot = transform_gps_to_local(self.path)
        #use x,y,z to get path length...
        s = 0
        for i in range(len(x)-1):
            ds = np.linalg.norm([x[i+1] - x[i], y[i+1] - y[i], z[i+1] - z[i]])
            s += ds
            #ToDo: just use values where velocity is higher then ...
            
        length_vel = np.trapz(self.data_pos['v_total'], time)
        print("Length from velocity Integration = {} meters".format(length_vel))
        print("Length from Position data = {} meters".format(s))
        
        
        
    def plot_temp(self, real_temp=1000, saveplot=False):
        
        #self.get_object_temp()
        linewidth = 0.5
        fig, ax = plt.subplots()
        
        if self.temp_bool:
            
            date = dates.date2num(self.data_pos["timestamp"])
            # Define the date format
            date_form = dates.DateFormatter("%H:%M:%S")
            ax.xaxis.set_major_formatter(date_form)
            ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=(0, 30)))
            ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=(15, 45)))
            
            #matplotlib.pyplot.plot_date(dates, y_values)
            
            ax.plot(date, self.data_pos['ambientTemp[C]'], color='r', label='Ambient Temp.', linewidth=linewidth)
            ax.plot(date, self.data_pos['objectTemp[C]'], color='b', label='Object Temp.', linewidth=linewidth)
            #ax.plot(date, self.data_pos['corr.objectTemp[C]'], color='magenta', label='Object Temp.; E=0.9', linewidth=linewidth)
            if real_temp < 1000:
                ax.axhline(real_temp, 0, max(date), color='g', linewidth=linewidth+0.5, label='Measured temp.')
                ax.set_ylim([real_temp -2, max(self.data_pos['ambientTemp[C]'])+2])
            ax.set_xlabel('Time [UTC +1]')
            ax.set_ylabel('Temperature [°C]')
            ax.set_title('Temperature development')
            ax.grid()
            ax.legend()
            if saveplot:
                fig.savefig(self.work_dir + '/' + self.file[:-4] + '_temperatures.png', dpi=300)
            
        else:
            print("No temperature data found!")
            
    def calc_E(self, real_temp, saveplot=False):
        E = ((self.data['objectTemp[C]'] + 273.15)**4 - (self.data['ambientTemp[C]'] + 273.15)**4) / ((real_temp + 273.15)**4 - (self.data['ambientTemp[C]'] + 273.15)**4)
        print("Mean Value of E = ",  np.mean(E))
        
        # Aproximate line through points
        predict_ambient = np.poly1d(np.polyfit(self.data['ambientTemp[C]'], E, 1))
        predict_object = np.poly1d(np.polyfit(self.data['objectTemp[C]'], E, 1))
        
        #print("Polynoms: ", predict)
        e_amb = predict_ambient(self.data['ambientTemp[C]'])
        e_obj = predict_object(self.data['objectTemp[C]'])
        
        real_t_amb = (((self.data['objectTemp[C]'] + 273.15)**4 - (self.data['ambientTemp[C]'] + 273.15)**4) / e_amb + (self.data['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        real_t_obj = (((self.data['objectTemp[C]'] + 273.15)**4 - (self.data['ambientTemp[C]'] + 273.15)**4) / e_obj + (self.data['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        
        c = (self.data['objectTemp[C]'] - real_temp) / self.data['ambientTemp[C]'] # Error / Ambient = c
        predict_c = np.poly1d(np.polyfit(self.data['ambientTemp[C]'], c, 2))
        print("Mean value c = ", np.mean(c))
        
        #Plot correltaiton Ambient Temp - E
        fig, ax = plt.subplots()
        ax.scatter(self.data['ambientTemp[C]'], E, marker='x', s=0.5, color='b', label='Emessivity')
        ax.plot(self.data['ambientTemp[C]'], e_amb, color='r')
        ax.set_title('Correlation Abmient Temp. - Emissivity')
        ax.set_xlabel('Ambient Temperature [°C]')
        ax.set_ylabel('Emissivity')
        if saveplot:
            fig.savefig(self.work_dir + '/' + self.file[:-4] + '_ambient_emissivity.png', dpi=300)
            
        #Plot correltaiton Object Temp - E
        fig, ax = plt.subplots()
        ax.scatter(self.data['objectTemp[C]'], E, marker='x', s=0.5, color='b')
        ax.plot(self.data['objectTemp[C]'], e_obj, color='r')
        ax.set_title('Correlation object Temp. - Emissivity')
        ax.set_xlabel('Object Temperature [°C]')
        ax.set_ylabel('Emissivity')
        if saveplot:
            fig.savefig(self.work_dir + '/' + self.file[:-4] + '_object_emissivity.png', dpi=300)
            
        #Plot correltaiton Error - Ambient
        fig, ax = plt.subplots()
        ax.scatter(self.data['ambientTemp[C]'], (self.data['objectTemp[C]'] - real_temp), marker='o', s=0.5, color='g', label='c') #ToDO: make plot error over ambient
        ax.set_title('Correlation object temp. error - ambient temp.')
        ax.set_xlabel('Ambient Temperature [°C]')
        ax.set_ylabel('Temp. Error [°C]')
        if saveplot:
            fig.savefig(self.work_dir + '/' + self.file[:-4] + '_object_emissivity.png', dpi=300)
        
        # Plot processed Temperature, mean value and measured value
        fig, ax = plt.subplots()
        
        date = dates.date2num(self.data[self.utc])
        # Define the date format
        date_form = dates.DateFormatter("%H:%M:%S")
        ax.xaxis.set_major_formatter(date_form)
        ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=(0, 30)))
        ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=(15, 45)))
        
        ax.plot(date, real_t_amb, label='IR Object Temp with ambient corr.', color='r')
        ax.plot(date, real_t_obj, label='IR Object Temp with object corr.', color='b')
        ax.axhline(np.mean(real_t_amb),0,100000, color='r', label='mean value amb corr')
        ax.axhline(np.mean(real_t_obj),0,100000, color='g', label='mean value obj corr')
        ax.axhline(np.mean(real_temp),0,100000, color='black', label='mean temp. measured', linestyle='--')
        ax.set_xlabel=('Index')
        ax.set_ylabel=('Temperature [°C]')
        ax.legend()
        if saveplot:
            fig.savefig(self.work_dir + '/' + self.file[:-4] + '_IR_temperatures.png', dpi=300)

        
        
        return predict_ambient, predict_object, predict_c
    
    def calc_temp(self):
        
        emissivity_snow = 1
        #e_amb = predict_amb(self.data_pos['ambientTemp[C]'])
        #e_obj = predict_obj(self.data_pos['objectTemp[C]'])
        #e_amb = predict(self.data_pos['objectTemp[C]'] + self.data_pos['ambientTemp[C]']) / (self.data_pos['objectTemp[C]'] - self.data_pos['ambientTemp[C]'])
        #real_t_amb = (((self.data_pos['objectTemp[C]'] + 273.15)**4 - (self.data_pos['ambientTemp[C]'] + 273.15)**4) / e_amb + (self.data_pos['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        #real_t_obj = (((self.data_pos['objectTemp[C]'] + 273.15)**4 - (self.data_pos['ambientTemp[C]'] + 273.15)**4) / e_obj + (self.data_pos['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        temp_corr = (((self.data_pos['objectTemp[C]'] + 273.15)**4 - (self.data_pos['ambientTemp[C]'] + 273.15)**4) / emissivity_snow + (self.data_pos['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        
        
        #t_c = self.data_pos['objectTemp[C]'] - 0.5573355725160639 * self.data_pos['ambientTemp[C]']
        #emm_folie = (temp_corr**4 - self.data_pos['ambientTemp[C]']**4)
        
        #c_predicted = predict_c(self.data_pos['ambientTemp[C]'])
        #t_c_predict = self.data_pos['objectTemp[C]'] - c_predicted * self.data_pos['ambientTemp[C]']
        #with emissivity = 0.9
        #t_c_em = temp_corr - 0.5573355725160639 * self.data_pos['ambientTemp[C]']
        #t_c_predict_em = real_t - c_predicted * self.data_pos['ambientTemp[C]']
        
# =============================================================================
#         # Measured data from 23.01.2022
#         self.data_pos_copy = self.data_pos
#         self.data_pos_copy['measuredTemp'] = np.nan
#         idx_start = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:27:00')[0][0]
#         idx_end = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:30:00')[0][0]
#         self.data_pos_copy['measuredTemp'].iloc[idx_start:idx_end] = -5.4
#         idx_start = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:30:00')[0][0]
#         idx_end = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:36:00')[0][0]
#         self.data_pos_copy['measuredTemp'].iloc[idx_start:idx_end] = -6.4
#         idx_start = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:36:00')[0][0]
#         idx_end = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:39:00')[0][0]
#         self.data_pos_copy['measuredTemp'].iloc[idx_start:idx_end] = -5.6
#         idx_start = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:39:00')[0][0]
#         idx_end = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:48:00')[0][0]
#         self.data_pos_copy['measuredTemp'].iloc[idx_start:idx_end] = -5.1
#         idx_start = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 07:48:00')[0][0]
#         idx_end = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 08:04:00')[0][0]
#         self.data_pos_copy['measuredTemp'].iloc[idx_start:idx_end] = -5.0
#         idx_start = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 08:04:00')[0][0]
#         idx_end = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 08:10:00')[0][0]
#         self.data_pos_copy['measuredTemp'].iloc[idx_start:idx_end] = -4.3
#         idx_start = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 08:10:00')[0][0]
#         idx_end = np.where(self.data_pos_copy.iloc[:, 2] == '2022-01-23 08:24:00')[0][0]
#         self.data_pos_copy['measuredTemp'].iloc[idx_start:idx_end] = 0
# =============================================================================
        
        
        #self.data_pos_copy.between_time('08:00:00', '08:05:00')
        
# =============================================================================
#         fig, ax = plt.subplots()
#         
#         date = dates.date2num(self.data_pos[self.utc])
#         # Define the date format
#         date_form = dates.DateFormatter("%H:%M")
#         #date_form = dates.DateFormatter("%m-%d")
#         ax.xaxis.set_major_formatter(date_form)
#         ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=(0, 10, 20, 30, 40, 50)))
#         ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=(5, 15, 25, 35, 45, 55)))
#         
#         ax.plot(date, real_t_amb, label='emessivity amb correct.')
#         #ax.plot(date, real_t_obj, label='Temperatures from IR - Object')
#         #ax.plot(date, t_c, label='c mean')
#         ax.plot(date, t_c_predict, label='c_predicted')
#         #ax.plot(date, t_c_em, label='c mean; E=0.9')
#         ax.plot(date, t_c_predict_em, label='c_predicted; E=0.9')
#         ax.plot(date, self.data_pos_copy['measuredTemp'], label='Measured Temp.', color='black')
#         #ax.axhline(np.mean(real_t_amb),min(date),max(date), color='r', label='mean value')
#         #ax.axhline(np.mean(real_temp),0,100000, color='black', label='mean IR measured', linestyle='--')
#         ax.set_xlabel=('Time [UTC+1]')
#         ax.set_ylabel=('Temperature [°C]') # ToDo: Add measured temperatures
#         
#         ax.legend()
#         ax.grid(b=True, which='major', linestyle='-')
#         ax.grid(b=True, which='minor', linestyle='--')
#         plt.show()
# =============================================================================
        
    def get_snow_temp(self, e_snow=1, e_folie=0.45):
        self.data_pos['corr.objectTemp[C]'] = (((self.data_pos['objectTemp[C]'] + 273.15)**4 - (self.data_pos['ambientTemp[C]'] + 273.15)**4) / e_snow + (self.data_pos['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        self.data['corr.objectTemp[C]'] = (((self.data['objectTemp[C]'] + 273.15)**4 - (self.data['ambientTemp[C]'] + 273.15)**4) / e_snow + (self.data['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15

        self.data['corr.snowTemp[C]'] = (((self.data['corr.objectTemp[C]'] + 273.15)**4 - (self.data['ambientTemp[C]'] + 273.15)**4) / e_folie + (self.data['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        self.data_pos['corr.snowTemp[C]'] = (((self.data_pos['corr.objectTemp[C]'] + 273.15)**4 - (self.data_pos['ambientTemp[C]'] + 273.15)**4) / e_folie + (self.data_pos['ambientTemp[C]'] + 273.15)**4)**(1/4) - 273.15
        
        self.data_pos['corr.snowTemp2[C]'] = ((self.data_pos['objectTemp[C]']+ 273.15) - ((self.data_pos['ambientTemp[C]']+ 273.15)* e_folie)) / ((1 - e_folie) * e_snow) - 273.15

        #self.data['corr.objectTemp[C]'] = (((self.data['objectTemp[C]'] + 273.15)**4 - (1-e_snow) * (self.data['ambientTemp[C]'] + 273.15)**4) / e_snow )**(1/4) - 273.15
        #self.data['corr.objectTemp[C]'] = self.data['objectTemp[C]'] - 0.5573355725160639 * self.data['ambientTemp[C]']


    # def save_to_vtk(self, Node=str):
        
    #     from pyevtk.hl import pointsToVTK
    #     from osgeo import gdal
        
    #     # Get Altitude from DHM because GPS is really wrong!
    #     filepath = r"C:\Users\neuhauser\OneDrive - Bundesforschungszentrum fuer Wald\20070-AvaRange-data\GIS\Nordkette\Nordkette_wgs84.tif"
    #     # Open the file:
    #     nordkette_raster = gdal.Open(filepath)

    #     # Check type of the variable 'raster'
    #     type(nordkette_raster)
    #     # Projection
    #     nordkette_raster.GetProjection()
    #     band = nordkette_raster.GetRasterBand(1)

    #     cols = nordkette_raster.RasterXSize
    #     rows = nordkette_raster.RasterYSize

    #     transform = nordkette_raster.GetGeoTransform()

    #     xOrigin = transform[0]
    #     yOrigin = transform[3]
    #     pixelWidth = transform[1]
    #     pixelHeight = -transform[5]

    #     data = band.ReadAsArray(0, 0, cols, rows)
    #     point_list = []
    #     for i in range(len(self.data["latitude[deg]"])):
    #         point_list.append([self.data["longitude[deg]"].iloc[i], self.data["latitude[deg]"].iloc[i]])

    #     z = []
    #     for point in point_list:
    #         col = int((point[0] - xOrigin) / pixelWidth)
    #         row = int((yOrigin - point[1] ) / pixelHeight)
    #         z.append(data[row][col])
    #         #print(row, col, data[row][col])
    #     self.data["z_raster"] = z
        
    #     north, east, height = gps_to_mercator(self)
    #     for i in range(len(north)):
    #         x = np.array(east[i])
    #         y = np.array(north[i])
    #         z = np.array(self.data["z_raster"].iloc[i], dtype=np.float64)
    #         vel = self.data_pos["v_total"].iloc[i]
    #         vel = np.array(vel, ndmin=1)
    #         pointsToVTK(r"C:\Users\neuhauser\OneDrive - Bundesforschungszentrum fuer Wald\20070-AvaRange-data\ParaView\220222_Nordkette_GNSS\C{}\GNSS_points_{}".format(Node, i), x, y, z, data = {"u" : vel})




