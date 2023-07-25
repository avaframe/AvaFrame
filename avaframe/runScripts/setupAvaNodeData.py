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
from classes.AvaNode_GNSS import AvaNode_GNSS

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
    