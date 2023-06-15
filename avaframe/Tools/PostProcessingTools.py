# -*- coding: utf-8 -*-
"""
Created on Mon Wed 12 2022

@author: dicko

Tools to extract information on the avalanche simulations run in the Output files

modified by AvaFrame

"""
# Python imports
import numpy as np
import statistics as stat
import pathlib

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
import avaframe.in1Data.getInput as gI
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
from avaframe.Tools import PlotTools


def reshapeParticlesDicts(particlesList, propertyList):
    """ reshape particlesList from one dict per time step with all particle properties for each particle,
        to one dict with an array of property values for all time steps for each particle
        shape: nx x ny; nx time steps, ny number of particles

        Parameters
        -----------
        particlesList: list
            list of particle dicts, one dict per time step
        propertyList: list
            list of property names that shall be reshaped and saved to particlesTimeArrays

        Returns
        --------
        particlesTimeArrays: dict
            dict with time series of properties of particles
            key: property, item: timeSteps x particlesID array of property values
    """

    # create particlesTimeArrays
    particlesTimeArrays = {}

    for props in propertyList:
        if isinstance(particlesList[0][props], int) or isinstance(particlesList[0][props], float):
            particlesTimeArrays[props] = np.zeros(len(particlesList))
            particlesTimeArrays[props][:] = np.asarray([p[props] for p in particlesList])
        else:
            particlesTimeArrays[props] = np.zeros((len(particlesList), particlesList[0]['nPart']))
            for idx in particlesList[0]['ID']:
                particlesTimeArrays[props][:,idx] = np.asarray([p[props][idx] for p in particlesList])

    return particlesTimeArrays


def velocityEnvelope(particlesTimeArrays):
    """ compute the velocity envelope from particle values

        Parameters
        -----------
        particlesTimeArrays: dict
            dict with time series for particle properties

        Returns
        --------
        dictVelEnvelope: dict
            max, mean, min values of velocity over all particles for each time step
            velocity, acceleration, min and max of trajectoryLengthXYZ of particles for each time step
    """

    # Preparing time, maximum values, minimum values, mean, median
    Max = np.zeros(len(particlesTimeArrays['t']))
    Min = np.zeros(len(particlesTimeArrays['t']))
    Mean = np.zeros(len(particlesTimeArrays['t']))
    Median = np.zeros(len(particlesTimeArrays['t']))
    SxyzMax = np.zeros(len(particlesTimeArrays['t']))
    SxyzMin = np.zeros(len(particlesTimeArrays['t']))

    # loop over time steps and compute max, min, ...
    for j in range(len(particlesTimeArrays['t'])):
        Max[j] = np.amax(particlesTimeArrays['velocityMag'][j,:])
        Min[j] = np.amin(particlesTimeArrays['velocityMag'][j,:])
        Mean[j] = np.mean(particlesTimeArrays['velocityMag'][j,:])
        Median[j] = np.median(particlesTimeArrays['velocityMag'][j,:])
        SxyzMax[j] = np.nanmax(particlesTimeArrays['trajectoryLengthXYZ'][j,:])
        SxyzMin[j] = np.nanmin(particlesTimeArrays['trajectoryLengthXYZ'][j,:])

    dictVelEnvelope = {}
    dictVelEnvelope['Velocity'] = particlesTimeArrays['velocityMag']
    dictVelEnvelope['Acc'] = particlesTimeArrays['uAcc']
    dictVelEnvelope['Mean'] = Mean
    dictVelEnvelope['Median'] = Median
    dictVelEnvelope['Min'] = Min
    dictVelEnvelope['Max'] = Max
    dictVelEnvelope['Time'] = particlesTimeArrays['t']
    dictVelEnvelope['SxyzMax'] = SxyzMax
    dictVelEnvelope['SxyzMin'] = SxyzMin

    return dictVelEnvelope


#%% function to generate the velocity envelope along the thalweg from simulated flows
def velocityEnvelopeThalweg(particlesTimeArrays):
    """ function to generate the velocity envelope along the thalweg from simulated flows

        Parameters
        -----------
        particlesTimeArrays: dict
            dict with time series of particle properties

        Returns
        --------
        dictVelAltThalweg: dict
            dict with min, max, mean, ... values of particles elevation, velocity and acceleration
            for each thalweg coordinate

    """

    time = particlesTimeArrays['t']
    sBetaPoint = particlesTimeArrays['sBetaPoint']

    elevation = particlesTimeArrays['z'].flatten()
    sXY = particlesTimeArrays['sAimec'].flatten()

    velocityMag = particlesTimeArrays['velocityMag'].flatten()
    uAcc = particlesTimeArrays['uAcc'].flatten()
    trajectoryLengthXYZ = particlesTimeArrays['trajectoryLengthXYZ'].flatten()

    # sort particle coordinates in thalweg system and only keep unique values
    sXYSorted = np.unique(sXY)

    # get max, min, mean, median values of elevation, velocity mag and trajectorylengthXYZ along thalweg
    maxVelocity = np.zeros(len(sXYSorted))
    minVelocity = np.zeros(len(sXYSorted))
    medianVelocity = np.zeros(len(sXYSorted))
    meanVelocity = np.zeros(len(sXYSorted))
    maxAcc = np.zeros(len(sXYSorted))
    minAcc = np.zeros(len(sXYSorted))
    medianAcc = np.zeros(len(sXYSorted))
    meanAcc = np.zeros(len(sXYSorted))
    maxZ = np.zeros(len(sXYSorted))
    minZ = np.zeros(len(sXYSorted))
    maxSxyz = np.zeros(len(sXYSorted))
    minSxyz = np.zeros(len(sXYSorted))

    # create array of properties along S_xy (thalweg/aimec) coordinate
    for indS, sXYCoor in enumerate(sXYSorted):
        sXYInd = np.where(sXY == sXYCoor)
        V = velocityMag[sXYInd[0]]
        A = uAcc[sXYInd[0]]
        S = trajectoryLengthXYZ[sXYInd[0]]
        maxVelocity[indS] = np.nanmax(velocityMag[sXYInd[0]])
        minVelocity[indS] = np.nanmin(velocityMag[sXYInd[0]])
        meanVelocity[indS] = np.nanmean(velocityMag[sXYInd[0]])
        medianVelocity[indS] = np.nanmedian(velocityMag[sXYInd[0]])
        maxZ[indS] = np.nanmax(elevation[sXYInd[0]])
        minZ[indS] = np.nanmin(elevation[sXYInd[0]])
        maxSxyz[indS] = np.nanmax(trajectoryLengthXYZ[sXYInd[0]])
        minSxyz[indS] = np.nanmin(trajectoryLengthXYZ[sXYInd[0]])
        maxAcc[indS] = np.nanmax(uAcc[sXYInd[0]])
        minAcc[indS] = np.nanmin(uAcc[sXYInd[0]])
        meanAcc[indS] = np.nanmean(uAcc[sXYInd[0]])
        medianAcc[indS] = np.nanmedian(uAcc[sXYInd[0]])

    dictVelAltThalweg = {}
    dictVelAltThalweg['sXYPart'] = particlesTimeArrays['sAimec']
    dictVelAltThalweg['maxSxyz'] = maxSxyz
    dictVelAltThalweg['minSxyz'] = minSxyz
    dictVelAltThalweg['sXYThalweg'] = sXYSorted
    dictVelAltThalweg['sBetaPoint'] = sBetaPoint
    dictVelAltThalweg['meanVelocity'] = meanVelocity
    dictVelAltThalweg['medianVelocity'] = medianVelocity
    dictVelAltThalweg['minVelocity'] = minVelocity
    dictVelAltThalweg['maxVelocity'] = maxVelocity
    dictVelAltThalweg['minZ'] = minZ
    dictVelAltThalweg['maxZ'] = maxZ
    dictVelAltThalweg['meanAcc'] = meanAcc
    dictVelAltThalweg['medianAcc'] = medianAcc
    dictVelAltThalweg['minAcc'] = minAcc
    dictVelAltThalweg['maxAcc'] = maxAcc

    return dictVelAltThalweg


#%% function to do the post process of Aimec results
def aimecPostProcess(avalancheDir, cfg):

    # TODO: this could be reduced by calling individual functions

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
    pathDict = aimecTools.checkAIMECinputs(cfgSetup, pathDict)

    # Run AIMEC postprocessing
    rasterTransfo, resAnalysisDF, plotDict, newRasters = ana3AIMEC.mainAIMEC(pathDict, inputsDF, cfg)

    return rasterTransfo, resAnalysisDF, plotDict, newRasters


# %% Calculating the max and min pfv pft and ppr envelopes from the raster files
def rasterVelField(simNum,avalancheDir,cfgAimec,pfvMinThreshold=1,pftMinThreshold=0.1,pprMinThreshold=10):

    # Aimec post processing
    rasterTransfo, resAnalysisDF, plotDict, newRasters = aimecPostProcess(avalancheDir, cfgAimec)

    # getting the max peak flow velocity
    maxPeakFlowVelocity = resAnalysisDF['pfvCrossMax']
    meanPeakFlowVelocity = resAnalysisDF['pfvCrossMean']
    maxPeakFlowThickness = resAnalysisDF['pftCrossMax']
    meanPeakFlowThickness = resAnalysisDF['pftCrossMean']
    maxPeakFlowPressure = resAnalysisDF['pprCrossMax']
    meanPeakFlowPressure = resAnalysisDF['pprCrossMean']

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
