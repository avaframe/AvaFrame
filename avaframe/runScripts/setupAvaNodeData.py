# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:28:56 2023

@author: neuhauser

This script is used to implement AvaNode datasets into the AvaFrame Plots.
For this the datasets are ...

"""
# Python imports 
import numpy as np
import os 
# Local imports 
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
from classes.AvaNode_GNSS_Class import AvaNode_GNSS

#%% Class measurement

class measurement:
    
    def __init__(self, name):
        self.name = name
        self.gps = AvaNode_GNSS()
        self.start = 0 # Starttime in secons
        self.end = 0 # Endtime in seconds
        
        
    def define_paths(self, gps_path):
        self.gps.path = gps_path

      
    def define_times(self, start, end):
        self.start = start
        self.end = end
        
    def read_data(self):
        self.gps.read_data()
        self.gps.calc_pos_vel()
        
    def process_data(self):
        
        if self.name[:2] == '21':
            self.gps_v = self.gps.data_pos['v_total_pos']
        else:
            self.gps_v = self.gps.data_pos['v_total']

        
        #Smooth GNSS Velocity
        self.gps_v_smooth = np.convolve(self.gps_v, self.gps_kernel, mode='same')
        self.gps_v_smooth_acc = np.convolve(self.gps_v, self.gps_acc_kernel, mode='same')
        
        # Get GNSS Time
        if self.name[:2] == '21' or self.name == '220123':
            self.gps_time = self.gps.data_pos.iloc[:,0]/1000 - self.gps.data_pos.iloc[0,0]/1000
        else:
            self.gps_time = self.gps.data_pos.iloc[:,1]/1000 - self.gps.data_pos.iloc[0,1]/1000
        
        self.gps_time_null = (self.gps_time - self.start)
        
        # Normalize Imu Time
        self.imu_time_normalized = (self.imu.time - self.start)/ (self.end - self.start)
                
        # Calculate GNSS Accelerations
        self.gnss_acc = np.diff(self.gps_v_smooth_acc) / np.diff(self.gps_time)
        


#%% Set Data

measurment_name_list = ['210315', '210316-C01', '210316-C03', '220123', '220203', 
                        '220222-C07', '220222-C09', '220222-C10', '230203-C06', '230204-C07', '230315']

gps_path_list = [r'..\data\2021-03-15_Nordkette_avalanche\C01\GPS\ava210315_C01.txt',
                 r'..\data\2021-03-16_Nordkette_avalanche\C01\GPS\ava210316_C01_GPS.txt',
                 r'..\data\2021-03-16_Nordkette_avalanche\C03\GPS\ava210316_C03_GPS.txt',
                 r'..\data\2022-01-23_Nordkette_avalanche\C01\GPS\ava220123_C01.txt',
                 r'..\data\2022-02-03_Nordkette_avalanche\C10\GPS\ava220203_seilbahn_C10_gnss.txt',
                 r'..\data\2022-02-22_Nordkette_avalanche\C07\GPS\220222_C07_avalanche_GPS.txt',
                 r'..\data\2022-02-22_Nordkette_avalanche\C09\GPS\220222_C09_avalanche_GPS.txt',
                 r'..\data\2022-02-22_Nordkette_avalanche\C10\GPS\220222_C10_avalanche_GPS.txt',
                 r'..\data\2023-02-03_Nordkette_avalanche\C06\GPS\ava230203_C06_seilbahn_gnss.txt',
                 r'..\data\2023-02-04_Nordkette_avalanche\C07\GPS\ava230204_C07_gnss.txt',
                 r'..\data\2023-03-15_Nordkette_avalanche\C10\GPS\ava230315_C10_seilbahn_gnss.txt']

start_time_list = [37.5, 39.16, 71.05, 49-0.25, 137.55, 46-.585, 28+0.1152, 61-0.461, 20+1.32, 144-2.07, 148+1.72]
end_time_list = [80+2.55, 82-1.071, 112-5.11, 98-3.94, 189-4.373, 85-4.95, 70-4.82, 102-4.56, 60+6.77, 190+0.35, 188+6.88]

diff_times = np.array(start_time_list) - np.array(end_time_list)

# Colors
color_list = ['brown', 'gray', 'pink', '#0400fa', '#fa9f18', '#18fa2e', '#fa0005', '#811df3', '#00f1ee',
              '#ff04ff', '#ffff00']

#%%Read Data
# Init Experiments in list with name
exp = []
for meas in measurment_name_list:
    exp.append(measurement(meas))

# Define path and Start/Endtime
for (meas, gps_path) in zip(exp, gps_path_list):
    meas.define_paths(gps_path)
# Read Data
for meas in exp:
    meas.read_data()
#%% Define Setting

for (meas, start, end) in zip(exp, start_time_list, end_time_list):
    meas.define_times(start, end)
    
#%% Process Data

for meas in exp:
    meas.process_data()

#%% Make data dict

# Preparing the output Dictionary
AvaNodeDict = {}



#%% function to prepare a dictionary with all the nodes information 
def produceAvaNodesDict(avalancheDir):
    if avalancheDir == 'data/avaSeilbahn': 
        dictNodes = {} 
        #Calculating the Nodes velocity
        NumberNodes = [7, 9, 10]  # Numbers of the nodes you want to plot
        velNodes = findNodeVelocity(NumberNodes, avalancheDir)
        dictNodes['temporalDomain'] = velNodes
        #Checking that the extraction of experimental data has been succesful 
        if len(velNodes) != len(NumberNodes):
            print("Not all experimental data have been found!")
        
        #Changing nodes data coordinates
        e07, n07, z07, e09, n09, z09, e10, n10, z10 = changeNodeCoordinate()
        dictNodes['Coordinates'] = {}
        dictNodes['Coordinates']['e07'] = e07 
        dictNodes['Coordinates']['n07'] = n07 
        dictNodes['Coordinates']['z07'] = z07 
        dictNodes['Coordinates']['e09'] = e09 
        dictNodes['Coordinates']['n09'] = n09 
        dictNodes['Coordinates']['z09'] = z09
        dictNodes['Coordinates']['e10'] = e10 
        dictNodes['Coordinates']['n10'] = n10 
        dictNodes['Coordinates']['z10'] = z10
        
        # Calculating the velocity envelope along the thalweg for the AvaNodes 
        dictNodThalweg = findNodeThalweg(NumberNodes,avalancheDir)
        dictNodes['thalwegDomain'] = dictNodThalweg
        # Calculating the travel length in XYZ for the AvaNodes 
        dictNodes = travelNodesXYZ(dictNodes)
        
        return dictNodes
    else: 
        print('No dictionary made for the AvaNodes, as the topography is not adapted')
        
#%% function to extract the Nodes velocity 
def findNodeVelocity(NumberNodes,avalancheDir):
    """ Reads the node text file and return velocity over time data for the Nodes 

    Parameters
    ----------
    NumberNodes: list
        number of the nodes (e.x [7, 9, 10])
    avalancheDir : str or pathlib object
        path to avalanche directory

    Returns
    -------
    DictNod : list
        dictionary that contains information on the velocity over time data for each node in the inputs
        
    """
    
    # Preparing the output Dictionary
    DictNodVel = {}
    
    for j in range(0,len(NumberNodes)): 

        # Reading the experiment files for the AvaNode 
        if NumberNodes[j] // 10 == 1 :
            node = "C"+str(NumberNodes[j])     
        else:
            node = "C0"+str(NumberNodes[j])   
        file = avalancheDir + "/AvaNode_data/220222_"+node+"_avalanche_GPS.txt"
        isExist = os.path.exists(file)
        if isExist== False:
            print("No experiment data found: the AvaNode data should be in "+avalancheDir+"/AvaNode_data/220222_"+node+"_avalanche_GPS.txt")
        else:        
            Time_experiment = [x.split(',')[0] for x in open(file).readlines()]
            velN = [x.split(',')[7] for x in open(file).readlines()]
            velE = [x.split(',')[8]for x in open(file).readlines()]
            velD = [x.split(',')[9] for x in open(file).readlines()]
            
            # Calculating the velocity norm
            vel = [None]*(len(velN)-1)
            time = [None]*(len(velN)-1)
            for i in range(1, len(velN)):
                vel[i-1] = np.sqrt(int(velN[i])**2 + int(velE[i])**2 + int(velD[i])**2) * 10**-3
                time[i-1] = int(Time_experiment[i])*10**-6
            
            # Deleting the null velocity before the blasting  
            threshold = 0.5           # Warning: arbitrary 
            index = vel.index(next((i for i in vel if int(i)>=threshold),None))
            vel_experimental = vel[index:]
            time_experimental = [i-time[index] for i in time[index:]] 

            DictNodVel[node] = {"Velocity":vel,"VelocityAfterBlasting":vel_experimental, "Time":time ,"TimeAfterBlasting":time_experimental, "index":index}

    return DictNodVel
    