# -*- coding: utf-8 -*-
"""
Created on Oct Mon 12 2022

@author: Oscar Dick

Plot velocity envelopes and several comparison tools for different simulations.
One specific simulation can also be investigated.
The present script has been used to produce the plots in my master thesis as well.

Last change: 10/05/2023

modified by AvaFrame

"""
# Python imports
import pathlib
import logging
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.com1DFA import particleTools
from avaframe.ana3AIMEC import ana3AIMEC
from avaframe.Tools import PostProcessingTools
from avaframe.Tools import PlotFunctions
from avaframe.com1DFA import com1DFA
import avaframe.in1Data.getInput as gI
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU


###############################################################################################################################################################################################################################################################################################################################################################################
##### README ##################################################################################################################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################################################################################################################
# The present script calls each of the plot functions made so far.

# The plotting functions are in PlotFunctions. Useful tools for the plotting are in PlotTools. Useful tools to post-process the
# data and prepare them for the plots are in PostProcessingTools. NodeTools contain useful tools to process the AvaNodes data

# For all: the three most interesting plots are produced by:
#   * Peak flow quantities and thalweg-altitude diagram: PlotFunctions.plotPeakVelVelThalwegEnvelope
#   * Thalweg-Altitude energy line, peak flow quantities and thalweg time diagram: PlotFunctions.plotPeakQuantThalTimeEnergyLine
#   * Comparative plots (vel,acc,s_{xyz}..): PlotFunctions.plotPeakQuantTrackedPartVel

#%%############################################################################################################################################################################################################################################################################################################################################################################
##### SECTION TO PREPARE THE DATA FOR THE PLOTS ###############################################################################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################################################################################################################
# WARNING: THE BELOW SECTON IS SPECIFICALLY DEDICATED TO EXTRACT POST-PROCESSING INFORMATION. IT IS LIKELY TO BE IMPROVED


# TO run this script:
# run com1DFA with particles

# create local logger
log = logging.getLogger(__name__)

# -----------Required settings-----------------
# log file name; leave empty to use default runLog.log
logName = 'runFrictionTest'
resTypePlots = ['ppr', 'pfv', 'pft']
pfvMinThreshold = 1
pftMinThreshold = 0.1
pprMinThreshold = 10
# -----------------------

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Load configuration info of all com1DFA simulations and extract the different avalanche simulations of the output file
log.info('Extracting the different avalanche simulations from the ini file')
SimDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# Some general and useful parameters
number_ava = len(SimDF.index)
avaDir = pathlib.Path(avalancheDir)
modName = 'com1DFA'

# get the configuration of aimec
cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC)
# fetch anaMod from aimec settings
# TODO: shall we reset to modName or only use anaMod as modName - however run script targeted towards com1DFA
anaMod = cfgAimec['AIMECSETUP']['anaMod']

# set input directory to load particles dicts from sims
inputDir = pathlib.Path(avalancheDir,'Outputs', modName, 'particles')

# loop over all sims
for i, simIndex in enumerate(SimDF.index):

    # first reset all dicts, lists,.. that are created for each sim in SimDF
    dictVelEnvelope = ''
    dictVelAltThalweg = ''
    dictVelAltThalwegPart = ''
    trackedPartProp = ''
    trackedPartPropAdapted = ''
    dictRaster = ''
    particlesTimeArrays = ''
    rasterTransfo = ''
    particlesList = ''
    demSim = ''

    # fetch name of simulation
    simName = SimDF['simName'].loc[simIndex]

    # fetch particle dicts from sim
    particlesList, _ = particleTools.readPartFromPickle(inputDir,
        simName=SimDF['simName'].loc[simIndex], flagAvaDir=False, comModule=modName)
    #
    # add aimec (thalweg) s, l coordinates to particle dicts and save to pickle
    particlesList, rasterTransfo, demSim = ana3AIMEC.addSLToParticles(avalancheDir,
        cfgAimec, SimDF['DEM'].loc[simIndex], particlesList, saveToPickle=True)

    # recreate the tracked particles dictionary using the particles
    # and config info of each SimDF.index[i]
    log.info('Extracting information on the tracked particles for simulation %s' % simIndex)
    # create path of sim config file and create config object for that sim
    simConfigFile = avaDir / 'Outputs' / modName / 'configurationFiles' / (SimDF['simName'].loc[simIndex] + '.ini')
    cfgParticles = cfgUtils.readCfgFile(avaDir, module='', fileName=simConfigFile)
    cfgParticles['TRACKPARTICLES']['particleProperties'] = 'ID|x|y|z|ux|uy|uz|m|h|uAcc|velocityMag|trajectoryLengthXY|trajectoryLengthXYZ'
    # reinitialize tracked particles using sim config
    _, trackedPartProp, _ = com1DFA.trackParticles(cfgParticles['TRACKPARTICLES'], demSim, particlesList)

    # reshape particle dicts from one dict per time step to time series for each particle property
    particlesTimeArrays = PostProcessingTools.reshapeParticlesDicts(particlesList,
        ['ux', 'uy', 'uz', 'uAcc', 'velocityMag', 'trajectoryLengthXY', 'trajectoryLengthXYZ',
        'x', 'y', 'z', 'm', 't', 'sAimec', 'sBetaPoint', 'beta'])

    # Calculating the velocity envelope of the simulations (time, max, min, mean, median)
    log.info('Calculating velocity time envelope...')
    dictVelEnvelope = PostProcessingTools.velocityEnvelope(particlesTimeArrays)

    # Calculating the velocity envelope in the thalweg coordinate system (velocity-altitude-thalweg envelope)
    log.info('Calculating velocity altitude thalweg envelopes...')
    dictVelAltThalweg = PostProcessingTools.velocityEnvelopeThalweg(particlesTimeArrays)

    # Calculating the velocity-altitude-thalweg information for the tracked particles
    log.info('Calculating velocity altitude thalweg envelopes for the tracked numerical particles...')
    # calculating and adding for each particle sAimec and lAimec (s and l along the thalweg in the projected xy plan)
    trackedPartPropAdapted = ana3AIMEC.aimecTransform(rasterTransfo, trackedPartProp, rasterTransfo['demHeader'], timeSeries=True)
    # Calculating the velocity envelope in the thalweg coordinate system (velocity-altitude-thalweg envelope) for the tracked particles
    dictVelAltThalwegPart = PostProcessingTools.velocityEnvelopeThalweg(trackedPartPropAdapted)

    # The thresholds pfvMinThreshold pftMinThreshold pprMinThreshold are used to define the masked array
    cfgAimec['FILTER'] = {'simName': SimDF['simName'].loc[simIndex]}
    log.info('Filter for: simName %s' % SimDF['simName'].loc[simIndex])
    dictRaster = PostProcessingTools.rasterVelField(i, avalancheDir, cfgAimec,
        pfvMinThreshold=pfvMinThreshold, pftMinThreshold=pftMinThreshold,
        pprMinThreshold=pprMinThreshold)

    # %% Plotting peak flow quantities and the velocity thalweg diagram
    PlotFunctions.plotPeakVelVelThalwegEnvelope(avalancheDir, simIndex, SimDF, rasterTransfo,
        dictVelAltThalweg, resTypePlots, anaMod, demData=demSim)

    PlotFunctions.plotPeakQuantThalTimeEnergyLine(avalancheDir, simIndex, SimDF,
        rasterTransfo, dictRaster, modName, demSim)

    PlotFunctions.plotPeakQuantTrackedPartVel(avalancheDir, simName, dictVelAltThalweg,
        dictVelAltThalwegPart, trackedPartProp, dictVelEnvelope, demSim, modName, rasterTransfo)
