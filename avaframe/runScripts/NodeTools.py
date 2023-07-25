# -*- coding: utf-8 -*-
"""
Created on Oct Mon 12 2022

@author: dicko

Plot velocity envelopes for different friction model parameters
Warning : can only be used to compare simulations run with the same friction 
model  
If the friction model used for the simulations is Coulomb, compares plots for different mu values.
If the friction model is Voellmy, compares plots for different mu ans xsi values.
If the friction model is SamosAT, compares plots for different mu and tau0. 

"""
# Python imports 
import numpy as np
import os 


# Local imports 
from classes.AvaNode_GNSS_Class import AvaNode_GNSS 

import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools


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
        e10, n10, z10 = changeNodeCoordinate()
        dictNodes['Coordinates'] = {}

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
    

#%% function to generate the range time data for the nodes 
def findNodeRangeTime(NumberNodes,avalancheDir):
    """ Reads the node csv file and return range time data for the Nodes 

    Parameters
    ----------
    NumberNodes: list
        number of the nodes (e.x [7, 9, 10])
    avalancheDir : str or pathlib object
        path to avalanche directory

    Returns
    -------
    DictNodRangeTime : list
        dictionary that contains information on the range time data for each node in the inputs
        
    """
    
    # Preparing the output Dictionary
    DictNodRangeTime = {}
    
    for j in range(0,len(NumberNodes)): 

        # Reading the experiment files for the AvaNode 
        if NumberNodes[j] // 10 == 1 :
            node = "C"+str(NumberNodes[j])     
        else:
            node = "C0"+str(NumberNodes[j])   
        file = avalancheDir + "/AvaNode_data/220222_"+node+"_radarCS.csv"
        isExist = os.path.exists(file)
        if isExist== False:
            print("No experiment data found: the AvaNode data should be in "+avalancheDir+"/AvaNode_data/220222_"+node+"_radarCS.csv")
        else:        
            Time_experiment = [x.split(',')[0] for x in open(file).readlines()]
            Distance_to_radar = [x.split(',')[1] for x in open(file).readlines()]
            velocity = [x.split(',')[2] for x in open(file).readlines()]

        time = [i*10**-1 for i in range(0, len(Distance_to_radar[1:]))]
        distance_to_radar = [float(i) for i in Distance_to_radar[1:]]

        # Selecting the correct time  start for the experimental data
        threshold = 0.5             # Warning: arbitrary  
        indexes = velocity.index(next((i for i in velocity[1:] if float(i)>threshold),None))
        time_treated = [i-time[indexes] for i in time[indexes:]]  
        distance_to_radar_treated = distance_to_radar[indexes:]    
        DictNodRangeTime[node] = {"RangeTime":distance_to_radar_treated, "Time":time_treated}
        
    return DictNodRangeTime
        


#%% function to generate the thalweg distance of the nodes 
def findNodeThalweg(NumberNodes,avalancheDir):
    """ Reads the node csv file and return thalweg altitude time data for the Nodes 

    Parameters
    ----------
    NumberNodes: list
        number of the nodes (e.x [7, 9, 10])
    avalancheDir : str or pathlib object
        path to avalanche directory
    
    Returns
    -------
    DictNodThalweg : list
        dictionary that contains information on the thalweg altitude time data for each node in the inputs
        
    """
    
    # Preparing the output Dictionary
    DictNodThalweg = {}
    
    X10, Y10, Z10 = changeNodeCoordinate()
    

    sC10 = np.array(X10)
    
    #cfgSetup
    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    anaMod = cfgSetup['anaMod']
    # define reference simulation
    inputsDF, resTypeList = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfg) 
    
    # pathDict
    refSimRowHash, refSimName, inputsDF, colorParameter, valRef = aimecTools.fetchReferenceSimNo(avalancheDir, inputsDF, anaMod,cfg)
    pathDict = {'refSimRowHash': refSimRowHash, 'refSimName': refSimName, 'compType': ['singleModule', anaMod],
            'colorParameter': colorParameter, 'resTypeList': resTypeList, 'valRef': valRef}
    pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)
    pathDict = aimecTools.checkAIMECinputs(cfgSetup, inputsDF, pathDict)
    
    # refHeader
    refSimRowHash = pathDict['refSimRowHash']
    refResultSource = inputsDF.loc[refSimRowHash, cfgSetup['runoutResType']]
    refRaster = IOf.readRaster(refResultSource)
    refHeader = refRaster['header']
    
    # dem
    demSource = pathDict['demSource']
    dem = IOf.readRaster(demSource)
    
    rasterTransfo = aT.makeDomainTransfo(pathDict, dem, refHeader['cellsize'], cfgSetup)
        

        
    for i in range(0,np.shape(X10)[0]):
        distance = np.sqrt((X10[i]-rasterTransfo['gridx'])**2 + (Y10[i]-rasterTransfo['gridy'])**2)
        (sIndex, lIndex) = np.unravel_index(np.argmin(distance, axis=None), distance.shape)
        sC10[i] = rasterTransfo['s'][sIndex]
        
    # sC10[0] = 0
    # for i in range(1,np.shape(X10)[0]):
    #     D = np.sqrt((X10[i] - X10[i-1])**2 + (Y10[i] - Y10[i-1])**2)
    #     sC10[i] = sC10[i-1] + D 
        
    DictNodVel = findNodeVelocity(NumberNodes,avalancheDir)
            

    DictNodThalweg['sC10AfterBlasting'] = sC10[DictNodVel['C10']['index']:]
    
    DictNodThalweg['sC10'] = sC10
    
    return DictNodThalweg



#%% function to change the coordinate system of the nodes and put them in the same coordinate system than the simulations

def changeNodeCoordinate():
    # Load information on AvAnodes data and changing coordinate system
    gps_c10 = AvaNode_GNSS()
    #gps_c09 = GPS_Class.GPSData() 
    #gps_c07 = GPS_Class.GPSData() 
    gps_c10.path = r"C:\git_rep\AvaFrame\avaframe\data\avaSeilbahnrinne\Inputs\NODES\2022-02-22_Nordkette_avalanche\C10\GPS/220222_C10_avalanche_GPS.txt"
    #gps_c09.path = "/home/dick/Documents/AvaFrame/avaframe/data/avaSeilbahn/AvaNode_data/220222_C09_avalanche_GPS.txt"
    #gps_c07.path = "/home/dick/Documents/AvaFrame/avaframe/data/avaSeilbahn/AvaNode_data/220222_C07_avalanche_GPS.txt"
    gps_c10.read_data() 
    #gps_c09.read_data() 
    #gps_c07.read_data() 
    n10,e10,z10 = gps_c10.gnss_to_mercator()
    #n09,e09,z09 = git.gps_to_mercator(gps_c09)
    #n07,e07,z07 = git.gps_to_mercator(gps_c07)
    return e10, n10, z10 

    
#%% function to calculate the average velocity of different nodes 

def averageNodesVelocity(NumberNodes,avalancheDir):
    """ calculate the average velocity of different nodes 

    Parameters
    ----------
    NumberNodes: list
        number of the nodes (e.x [7, 9, 10])
    avalancheDir : str or pathlib object
        path to avalanche directory

    Returns
    -------
    meanNodeVelocity : float 
        average value of the velocity of NumberNodes 
        
    """
    

    DictNodVel = findNodeVelocity(NumberNodes,avalancheDir)

    min_node_length = min(len(DictNodVel['C07']['VelocityAfterBlasting'][:]), 
                      len(DictNodVel['C09']['VelocityAfterBlasting'][:]),
                      len(DictNodVel['C09']['VelocityAfterBlasting'][:]))
    meanNodeVelocity = [None]*min_node_length
    for i in range (0,min_node_length):
        meanNodeVelocity[i] = (DictNodVel['C07']['VelocityAfterBlasting'][i]+
                                  DictNodVel['C09']['VelocityAfterBlasting'][i]+
                                  DictNodVel['C10']['VelocityAfterBlasting'][i])/3
        
    return meanNodeVelocity


#%% function to calculate the travel length in the xyz space of the AvaNodes

def travelNodesXYZ(dictNodes):
    

    sC10xyz =  dictNodes['Coordinates']['e10'].copy()
    sC10xyz[0] = 0 

  
    for i in range(0,len(dictNodes['Coordinates']['e10'])-1):
           sC10xyz[i+1] = sC10xyz[i] + np.sqrt( (dictNodes['Coordinates']['e10'][i+1] - dictNodes['Coordinates']['e10'][i])**2 + (dictNodes['Coordinates']['n10'][i+1] - dictNodes['Coordinates']['n10'][i])**2 +(dictNodes['Coordinates']['z10'][i+1] - dictNodes['Coordinates']['z10'][i])**2 )     
        

    sC10xyz[-1] =  sC10xyz[len(dictNodes['Coordinates']['e10'])-1]
    

    dictNodes['thalwegDomain']['sC10xyz'] = sC10xyz

    # After blasting

    dictNodes['thalwegDomain']['sC10xyzAfterBlasting'] = sC10xyz[dictNodes['temporalDomain']['C10']['index']:]
    
    return dictNodes

if __name__ == "__main__":
    # Example usage of this script
   produceAvaNodesDict('data/avaSeilbahn')


