# -*- coding: utf-8 -*-
"""
Created on Mon Wed 12 2022

@author: dicko

Tools to extract information on the avalanche simulations run in the Output files  

"""
# Python imports 
import numpy as np
import statistics as stat
import pathlib
from matplotlib import pyplot as plt

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
import avaframe.in1Data.getInput as gI
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
from avaframe.Tools import PlotTools


#%% Function to extract the different avalanche simulations using postProcess 

# def FindAvaSimu(): 
#     # Load avalanche directory from general configuration file
#     cfgMain = cfgUtils.getGeneralConfig()
#     avalancheDir = cfgMain['MAIN']['avalancheDir']
    
#     # Extracting the different avalanche simulations of the output file 
#     F = cfgpostProcess(avalancheDir)
    
#     return F 


#%% function to generate the tracked particles dictionary 

def trackedParticlesDictionary(cfg, particlesList):
    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']
    
    # fetch dem data
    dem = gI.readDEM(avalancheDir)
    dem['originalHeader'] = {} 
    dem['originalHeader']['xllcenter'] = dem['header']['xllcenter']
    dem['originalHeader']['yllcenter'] = dem['header']['yllcenter']
    
    cfgTrackPart = cfg['TRACKPARTICLES']
    
    particlesList, trackedPartProp, track = com1DFA.trackParticles(cfgTrackPart, dem, particlesList)

    return trackedPartProp


#%% function to calculate the velocity magnitude 

def velocityMagnitude(avaDict):
    # Preparing  magnitude list  
    number_ava = len(avaDict)
    Velocity = [None]*number_ava
    
    for i in range(0, number_ava):
        number_time_steps = len(avaDict[i])
        # Reading and calculating the velocity magnitude for avalanche i
        df = [None]*number_time_steps
        for k in range(0, number_time_steps):
            df[k] = np.sqrt(np.array(avaDict[i][k]['ux']**2 + avaDict[i]
                            [k]['uy']**2 + avaDict[i][k]['uz']**2))
    
        # Preparing velocity magnitude for avalanche i
        Velocity[i] = [None]*number_time_steps
        for j in range(0, number_time_steps):
            Velocity[i][j] = df[j] 
            
    return Velocity


#%% function to generate the velocity envelope from simulated flows

def velocityEnvelope(avaDict):
    # calculating velocity magnitude for every simulations 
    Velocity = velocityMagnitude(avaDict)
    # Preparing time, maximum values, minimum values, mean, median
    number_ava = len(avaDict)
    Time = [None]*number_ava
    Max = [None]*number_ava
    Min = [None]*number_ava
    Mean = [None]*number_ava
    Median = [None]*number_ava
    SxyzMax = [None]*number_ava
    SxyzMin = [None]*number_ava
    
    for i in range(0, number_ava):
        number_time_steps = len(avaDict[i])

        # Preparing time, maximum values, minimum values, mean, median for avalanche i
        Time[i] = [None]*number_time_steps
        Max[i] = [None]*number_time_steps
        Min[i] = [None]*number_time_steps
        Mean[i] = [None]*number_time_steps
        Median[i] = [None]*number_time_steps
        SxyzMax[i] = [None]*number_time_steps
        SxyzMin[i] = [None]*number_time_steps
    
        for j in range(0, number_time_steps):
            Time[i][j] = avaDict[i][j]['t']
            Max[i][j] = max(Velocity[i][j])
            Min[i][j] = min(Velocity[i][j])
            Mean[i][j] = stat.mean(Velocity[i][j])
            Median[i][j] = stat.median(Velocity[i][j])
            SxyzMax[i][j] = max(avaDict[i][j]['sCor'])
            SxyzMin[i][j] = min(avaDict[i][j]['sCor'])
        
    dictVelEnvelope = {} 
    dictVelEnvelope['Velocity'] = Velocity
    dictVelEnvelope['Mean'] = Mean
    dictVelEnvelope['Median'] = Median
    dictVelEnvelope['Min'] = Min
    dictVelEnvelope['Max'] = Max
    dictVelEnvelope['Time'] = Time
    dictVelEnvelope['SxyzMax'] = SxyzMax
    dictVelEnvelope['SxyzMin'] = SxyzMin
            
    return dictVelEnvelope


#%% function to generate the velocity envelope along the thalweg from simulated flows

def velocityEnvelopeThalweg(avaDict):
    # calculating velocity magnitude for every simulations 
    simVelocity = velocityMagnitude(avaDict)
    # Preparing time, maximum values, minimum values, mean, median
    number_ava = len(avaDict)
    time = [None]*number_ava
    maxZ = [None]*number_ava
    minZ = [None]*number_ava
    meanZ = [None]*number_ava
    medianZ = [None]*number_ava
    sProj = [None]*number_ava
    s = [None]*number_ava
    z = [None]*number_ava
    sBetaPoint = [None]*number_ava
    sSortedDeduplicated = [None]*number_ava
    maxVelocity = [None]*number_ava
    minVelocity = [None]*number_ava
    medianVelocity = [None]*number_ava
    meanVelocity = [None]*number_ava
    listZsorted = [None]*number_ava
    maxZ = [None]*number_ava
    minZ = [None]*number_ava
    Velocity = [None]*number_ava
    maxSxyz = [None]*number_ava
    minSxyz = [None]*number_ava
    
    for i in range(0, number_ava):
    
        # Preparing time, maximum values, minimum values, mean, median for avalanche i
        number_time_steps = len(avaDict[i])
        time[i] = [None]*number_time_steps
        maxZ[i] = [None]*number_time_steps
        minZ[i] = [None]*number_time_steps
        meanZ[i] = [None]*number_time_steps
        medianZ[i] = [None]*number_time_steps
        sProj[i] = [None]*number_time_steps
        s[i] = [None]*number_time_steps
        z[i] = [None]*number_time_steps
        sBetaPoint[i] = [None]*number_time_steps
    
        for j in range(0, number_time_steps):
            time[i][j] = avaDict[i][j]['t']
            maxZ[i][j] = max(avaDict[i][j]['z'])
            minZ[i][j] = min(avaDict[i][j]['z'])
            meanZ[i][j] = stat.mean(avaDict[i][j]['z'])
            medianZ[i][j] = stat.median(avaDict[i][j]['z'])
            sProj[i][j] = avaDict[i][j]['sAimec']
            s[i][j] = avaDict[i][j]['s']
            z[i][j] = avaDict[i][j]['z']
            sBetaPoint[i][j] = avaDict[i][j]['sBetaPoint']
        # Changing the format of the data 
        Z = [None]*len(z[i])
        SProj = [None]*len(z[i])
        SXYZ = [None]*len(z[i])
        Vel = [None]*len(z[i])
        for k in range(0, len(z[i])):
            Z[k] = list(z[i][k])
            SProj[k] = list(sProj[i][k])
            SXYZ[k] = list(s[i][k])
            Vel[k] = list(simVelocity[i][k]) 
        listZ = sum(Z, [])
        listSProj = sum(SProj, [])
        listS = sum(SXYZ, [])
        listVel = sum(Vel, [])
        
        # Reading data for each s value 
        sSortedDeduplicated[i] = sorted(set(listSProj))
        maxVelocity[i] = [None]*len(sSortedDeduplicated[i])
        minVelocity[i] = [None]*len(sSortedDeduplicated[i])
        medianVelocity[i] = [None]*len(sSortedDeduplicated[i])
        meanVelocity[i] = [None]*len(sSortedDeduplicated[i])
        listZsorted[i] = [None]*len(sSortedDeduplicated[i])
        maxZ[i] = [None]*len(sSortedDeduplicated[i])
        minZ[i] = [None]*len(sSortedDeduplicated[i])
        Velocity[i] = [None]*len(sSortedDeduplicated[i])
        maxSxyz[i] = [None]*len(sSortedDeduplicated[i])
        minSxyz[i] = [None]*len(sSortedDeduplicated[i])
        # Loop on the s values 
        for k in range(0, len(sSortedDeduplicated[i])):
            listIndex = np.where(listSProj == sSortedDeduplicated[i][k])
            V = []
            L = []
            S = [] 
            for j in listIndex[0]:
                V.append(listVel[j])
                L.append(listZ[j])
                S.append(listS[j])
            listZsorted[i][k] = L
            maxVelocity[i][k] = max(V)
            minVelocity[i][k] = min(V)
            meanVelocity[i][k] = np.mean(V)
            medianVelocity[i][k] = np.median(V)
            maxZ[i][k] = max(L)
            minZ[i][k] = min(L)
            Velocity[i][k] = V    
            maxSxyz[i][k] = max(S)
            minSxyz[i][k] = min(S)
        
    dictVelAltThalweg = {}
    dictVelAltThalweg['sXYPart'] = sProj
    dictVelAltThalweg['maxSxyz'] = maxSxyz
    dictVelAltThalweg['minSxyz'] = minSxyz
    dictVelAltThalweg['sXYThalweg'] = sSortedDeduplicated
    dictVelAltThalweg['sBetaPoint'] = sBetaPoint
    #dictVelAltThalweg['Velocity'] = Velocity
    dictVelAltThalweg['meanVelocity'] = meanVelocity
    dictVelAltThalweg['medianVelocity'] = medianVelocity
    dictVelAltThalweg['minVelocity'] = minVelocity
    dictVelAltThalweg['maxVelocity'] = maxVelocity
    dictVelAltThalweg['listZsorted'] = listZsorted
    dictVelAltThalweg['minZ'] = minZ
    dictVelAltThalweg['maxZ'] = maxZ
    
    return dictVelAltThalweg
    

#%% function to adat the dictionary from the tracked particles and put it in the form of the avalanche particles dictionary

def dictChangeTrackedPart(trackedPartProp):
    trackedPartPropAdapted = {} 
    
    number_ava = len(trackedPartProp)
    
    for i in range(number_ava):
        number_time_steps = len(trackedPartProp[i]['time'])
        trackedPartPropAdapted[i] = {}
        for j in range(number_time_steps):
            trackedPartPropAdapted[i][j] = {}
            trackedPartPropAdapted[i][j]['h'] = trackedPartProp[i]['h'][j,:]
            trackedPartPropAdapted[i][j]['m'] = trackedPartProp[i]['m'][j,:]
            trackedPartPropAdapted[i][j]['t'] = trackedPartProp[i]['time'][j]
            trackedPartPropAdapted[i][j]['ux'] = trackedPartProp[i]['ux'][j,:]
            trackedPartPropAdapted[i][j]['uy'] = trackedPartProp[i]['uy'][j,:]
            trackedPartPropAdapted[i][j]['uz'] = trackedPartProp[i]['uz'][j,:]
            trackedPartPropAdapted[i][j]['x'] = trackedPartProp[i]['x'][j,:]
            trackedPartPropAdapted[i][j]['y'] = trackedPartProp[i]['y'][j,:]
            trackedPartPropAdapted[i][j]['z'] = trackedPartProp[i]['z'][j,:]
    
    return trackedPartPropAdapted


#%% function to calculate sAimec and lAimec for the tracked particles  

def sAimeclAimecCalculation(trackedPartPropAdapted,avalancheDir):
    
    rasterTransfo, resAnalysisDF, plotDict, newRasters = aimecPostProcess(avalancheDir)
    number_ava = len(trackedPartPropAdapted)
    # fetch dem data
    dem = gI.readDEM(avalancheDir)
    xllcenter = dem['header']['xllcenter']
    yllcenter = dem['header']['yllcenter']
    
    for i in range(number_ava):
        number_time_steps = len(trackedPartPropAdapted[i])
        for j in range(number_time_steps):
            lList = []  # mettre en array 
            sList  = [] 
            for x, y in zip(trackedPartPropAdapted[i][j]['x'], trackedPartPropAdapted[i][j]['y']):      
                # calculating the distance between the particle position and the grid points 
                distance = np.sqrt((x+xllcenter-rasterTransfo['gridx'])**2 + (y+yllcenter-rasterTransfo['gridy'])**2)  
                (sIndex, lIndex) = np.unravel_index(np.argmin(distance, axis=None), distance.shape)
                lList.append(rasterTransfo['l'][lIndex]) 
                sList.append(rasterTransfo['s'][sIndex]) # call the already existing function        
            trackedPartPropAdapted[i][j]['lAimec'] = lList
            trackedPartPropAdapted[i][j]['sAimec'] = sList
                
            #distance = np.sqrt((trackedPartPropAdapted[i][j]['x']-rasterTransfo['gridx'])**2 + (trackedPartPropAdapted[i][j]['y']-rasterTransfo['gridy'])**2)
            #(sIndex, lIndex) = np.unravel_index(np.argmin(distance, axis=None), distance.shape)
            #trackedPartPropAdapted[i][j]['sAimec'] = rasterTransfo['s'][sIndex]
            trackedPartPropAdapted[i][j]['sBetaPoint'] = rasterTransfo['s'][rasterTransfo['indStartOfRunout']]
    return(trackedPartPropAdapted)


#%% function to calculate the travel length of particles 

def sCalculation(trackedPartPropAdapted):

    number_ava = len(trackedPartPropAdapted)
        
    for i in range(number_ava):
        number_time_steps = len(trackedPartPropAdapted[i])
        sList  = [0]*len(trackedPartPropAdapted[i][0]['x'])
        trackedPartPropAdapted[i][0]['s'] = sList
        for j in range(number_time_steps-1):
            distance = np.sqrt((trackedPartPropAdapted[i][j+1]['x']-trackedPartPropAdapted[i][j]['x'])**2+(trackedPartPropAdapted[i][j+1]['y']-trackedPartPropAdapted[i][j]['y'])**2)      
            sList = sList + distance 
            trackedPartPropAdapted[i][j+1]['s'] = sList
    return(trackedPartPropAdapted)


#%% function to do the post process of Aimec results 

def aimecPostProcess(avalancheDir):
    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    anaMod = cfgSetup['anaMod']

    # Setup input from computational module
    inputsDF, resTypeList = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfg)
    # define reference simulation
    refSimRowHash, refSimName, inputsDF, colorParameter, valRef = aimecTools.fetchReferenceSimNo(avalancheDir, inputsDF, anaMod,
                                                                                         cfg)
    pathDict = {'refSimRowHash': refSimRowHash, 'refSimName': refSimName, 'compType': ['singleModule', anaMod],
                'colorParameter': colorParameter, 'resTypeList': resTypeList, 'valRef': valRef}
    pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)
    pathDict = aimecTools.checkAIMECinputs(cfgSetup, inputsDF, pathDict)

    # Run AIMEC postprocessing
    rasterTransfo, resAnalysisDF, plotDict, newRasters = ana3AIMEC.mainAIMEC(pathDict, inputsDF, cfg)
    plt.close()
    
    return rasterTransfo, resAnalysisDF, plotDict, newRasters

 
#%% function to gather the information for the energy line test plot

from avaframe.ana1Tests import energyLineTest 
# import analysis modules
import avaframe.ana5Utils.DFAPathGeneration as DFAPath

def energyLinePostProcessing():
    
    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']   
    energyLineTestCfg = cfgUtils.getModuleConfig(energyLineTest) 
    com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
                                          onlyDefault=energyLineTestCfg['com1DFA_override']['defaultConfig'])
    
    # Load configuration info of all com1DFA simulations
    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)
    simName = simDF.index[0] 
    
    # dem
    dem = gI.readDEM(avalancheDir)
    
    pathFromPart = energyLineTestCfg['energyLineTest'].getboolean('pathFromPart')
    avaProfileMass, particlesIni = DFAPath.generateAveragePath(avalancheDir, pathFromPart, simName, dem,
                                                               addVelocityInfo=True)
    
    # extend path profile and find intersection between the alpha line and the profile
    mu = com1DFACfg['GENERAL'].getfloat('mu')
    csz = dem['header']['cellsize']
    slopeExt, sIntersection, zIntersection, coefExt = energyLineTest.getAlphaProfileIntersection(energyLineTestCfg, avaProfileMass,
                                                                                  mu, csz)
    # compute errors on runout and veloctity altitude
    g = com1DFACfg['GENERAL'].getfloat('gravAcc')
    alphaRad = np.arctan(mu)
    alphaDeg = np.rad2deg(alphaRad)
    # compute simulation run out angle
    runOutAngleRad, runOutAngleDeg = energyLineTest.getRunOutAngle(avaProfileMass)
    zEne, u2Path, sGeomL, zGeomL, resultEnergyTest = energyLineTest.getEnergyInfo(avaProfileMass, g, mu, sIntersection, zIntersection,
                                                                   runOutAngleDeg, alphaDeg)
    
    energyLineDict = {} 
    energyLineDict['avaProfileMass'] = avaProfileMass
    energyLineDict['particlesIni'] = particlesIni
    energyLineDict['zEne'] = zEne
    energyLineDict['runOutAngleDeg'] = runOutAngleDeg
    energyLineDict['u2Path'] = u2Path
    energyLineDict['sGeomL'] = sGeomL
    energyLineDict['zGeomL'] = zGeomL
    energyLineDict['alphaDeg'] = alphaDeg
    energyLineDict['energyLineTestCfg'] = energyLineTestCfg
    energyLineDict['mu'] = mu
    
    slopeExt, sIntersection, zIntersection, coefExt = energyLineTest.getAlphaProfileIntersection(energyLineTestCfg, avaProfileMass,
                                                                                  mu, csz)
    
    energyLineDict['coefExt'] = coefExt
    energyLineDict['slopeExt'] = slopeExt
    energyLineDict['sIntersection'] = sIntersection
    energyLineDict['zIntersection'] = zIntersection
    
    return energyLineDict


# %% Calculating the max and min pfv pft and ppr envelopes from the raster files 

def rasterVelField(simNum,avalancheDir,pfvMinThreshold=1,pftMinThreshold=0.1,pprMinThreshold=10):

    # Aimec post processing
    rasterTransfo, resAnalysisDF, plotDict, newRasters = aimecPostProcess(avalancheDir)
    
    # getting the max peak flow velocity 
    maxPeakFlowVelocity = resAnalysisDF['pfvCrossMax']
    #minPeakFlowVelocity = resAnalysisDF['pfvCrossMin']
    meanPeakFlowVelocity = resAnalysisDF['pfvCrossMean']
    maxPeakFlowThickness = resAnalysisDF['pftCrossMax']
    #minPeakFlowThickness = resAnalysisDF['pftCrossMin']
    meanPeakFlowThickness = resAnalysisDF['pftCrossMean']
    maxPeakFlowPressure = resAnalysisDF['pprCrossMax']
    #minPeakFlowPressure = resAnalysisDF['pprCrossMin']
    meanPeakFlowPressure = resAnalysisDF['pprCrossMean']
     
    #avaDir = pathlib.Path(avalancheDir)
    #plotDict = PlotTools.PeakFields(avaDir, peakFilesDF, simNum, demData='')
    
    # Calculating the masked arraw and extracting the max, min and mean of altitude as well as the min of peak flow thickness
    # peak flow velocity and peak flow pressure 
    demMasked = np.ma.masked_where(newRasters['newRasterPFV'] == 0.0, newRasters['newRasterDEM'])
    pfvMasked = np.ma.masked_less(newRasters['newRasterPFV'],pfvMinThreshold) 
    pftMasked = np.ma.masked_less(newRasters['newRasterPFT'],pftMinThreshold) 
    pprMasked = np.ma.masked_less(newRasters['newRasterPPR'],pprMinThreshold) 
    AltitudeCrossMax = np.nanmax(demMasked, 1)
    AltitudeCrossMin = np.nanmin(demMasked, 1)
    AltitudeCrossMean = np.nanmean(demMasked, 1)
    minPeakFlowVelocity = np.nanmin(pfvMasked,1)
    minPeakFlowThickness = np.nanmin(pftMasked,1)
    minPeakFlowPressure = np.nanmin(pprMasked,1)
    
    # preparing the output dictionary 
    dictRaster = {}
    dictRaster['maxPeakFlowVelocity'] = maxPeakFlowVelocity
    dictRaster['minPeakFlowVelocity'] = minPeakFlowVelocity
    dictRaster['meanPeakFlowVelocity'] = meanPeakFlowVelocity
    dictRaster['maxPeakFlowThickness'] = maxPeakFlowThickness
    dictRaster['minPeakFlowThickness'] = minPeakFlowThickness
    dictRaster['meanPeakFlowThickness'] = meanPeakFlowThickness
    dictRaster['maxPeakFlowPressure'] = maxPeakFlowPressure
    dictRaster['minPeakFlowPressure'] = minPeakFlowPressure
    dictRaster['meanPeakFlowPressure'] = meanPeakFlowPressure
    dictRaster['AltitudeCrossMax'] = AltitudeCrossMax
    dictRaster['AltitudeCrossMin'] = AltitudeCrossMin
    dictRaster['AltitudeCrossMean'] = AltitudeCrossMean

    return dictRaster







