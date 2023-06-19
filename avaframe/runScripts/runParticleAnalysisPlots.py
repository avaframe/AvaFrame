# -*- coding: utf-8 -*-
"""

Plot velocity envelopes and several comparison tools for different simulations.
One specific simulation can also be investigated.
The present script has been used to produce the plots in my master thesis as well.


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
from avaframe.out3Plot import outParticlesAnalysis as oPartAna
from avaframe.com1DFA import com1DFA
import avaframe.in1Data.getInput as gI
from avaframe.in3Utils import logUtils
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.ana5Utils.distanceTimeAnalysis as dtAna
from avaframe.in3Utils import cfgHandling


##### README #######

# The present script calls each of the plot functions made so far.

# The plotting functions are in PlotFunctions. Useful tools for the plotting are in PlotTools. Useful tools to post-process the
# data and prepare them for the plots are in PostProcessingTools. NodeTools contain useful tools to process the AvaNodes data

# For all: the three most interesting plots are produced by:
#   * Peak flow quantities and thalweg-altitude diagram: PlotFunctions.plotPeakVelVelThalwegEnvelope
#   * Thalweg-Altitude energy line, peak flow quantities and thalweg time diagram: PlotFunctions.plotPeakQuantThalTimeEnergyLine
#   * Comparative plots (vel,acc,s_{xyz}..): PlotFunctions.plotPeakQuantTrackedPartVel

##### SECTION TO PREPARE THE DATA FOR THE PLOTS ######
# WARNING: THE BELOW SECTON IS SPECIFICALLY DEDICATED TO EXTRACT POST-PROCESSING INFORMATION. IT IS LIKELY TO BE IMPROVED

# -----------Required settings-----------------
# log file name; leave empty to use default runLog.log
logName = 'runFrictionTest'
modName = 'com1DFA'
# -----------------------

# load configuration of particle analysis using com1DFA as computational module
cfgPartAna = cfgUtils.getModuleConfig(oPartAna)
resTypePlots = fU.splitIniValueToArraySteps(cfgPartAna['GENERAL']['resTypePlots'])
pfvMinThreshold = cfgPartAna['GENERAL'].getfloat('pfvMinThreshold')
pftMinThreshold = cfgPartAna['GENERAL'].getfloat('pftMinThreshold')
pprMinThreshold = cfgPartAna['GENERAL'].getfloat('pprMinThreshold')

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']
avaDir = pathlib.Path(avalancheDir)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# check if running com1DFA sim or using already exisiting outputs
if cfgPartAna['GENERAL'].getboolean('runCom1DFA'):
    log.info('Perform com1DFA runs using the override section in the outParticlesAnalysis ini file')
    # get the configuration of com1DFA using overrides
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
        onlyDefault=cfgPartAna['com1DFA_override'].getboolean('defaultConfig'))
    cfgCom1DFA, cfgPartAna = cfgHandling.applyCfgOverride(cfgCom1DFA, cfgPartAna, com1DFA, addModValues=False)

    # call com1DFA and perform simulations
    _, _, _, SimDF = com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgCom1DFA)
else:
    # Load configuration info of all com1DFA simulations and extract the different avalanche simulations of the output file
    log.info('Fetching the different avalanche simulations from %s' % str(avalancheDir))
    SimDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# get the configuration of aimec using overrides
cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC, fileOverride='', modInfo=False, toPrint=False,
    onlyDefault=cfgPartAna['ana3AIMEC_override'].getboolean('defaultConfig'))
cfgAimec, cfgPartAna = cfgHandling.applyCfgOverride(cfgAimec, cfgPartAna, ana3AIMEC, addModValues=False)
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
    resAnalysisDF = ''
    newRasters = ''

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
    particlesTimeArrays = particleTools.reshapeParticlesDicts(particlesList,
        ['ux', 'uy', 'uz', 'uAcc', 'velocityMag', 'trajectoryLengthXY', 'trajectoryLengthXYZ',
        'x', 'y', 'z', 'm', 't', 'sAimec', 'sBetaPoint', 'beta'])

    # Calculating the velocity envelope of the simulations (time, max, min, mean, median)
    log.info('Calculating velocity time envelope...')
    dictVelEnvelope = oPartAna.velocityEnvelope(particlesTimeArrays)

    # Calculating the velocity envelope in the thalweg coordinate system (velocity-altitude-thalweg envelope)
    log.info('Calculating velocity altitude thalweg envelopes...')
    dictVelAltThalweg = oPartAna.velocityEnvelopeThalweg(particlesTimeArrays)

    # Calculating the velocity-altitude-thalweg information for the tracked particles
    log.info('Calculating velocity altitude thalweg envelopes for the tracked numerical particles...')
    # calculating and adding for each particle sAimec and lAimec (s and l along the thalweg in the projected xy plan)
    trackedPartPropAdapted = ana3AIMEC.aimecTransform(rasterTransfo, trackedPartProp, rasterTransfo['demHeader'], timeSeries=True)
    # Calculating the velocity envelope in the thalweg coordinate system (velocity-altitude-thalweg envelope) for the tracked particles
    dictVelAltThalwegPart = oPartAna.velocityEnvelopeThalweg(trackedPartPropAdapted)

    # The thresholds pfvMinThreshold pftMinThreshold pprMinThreshold are used to define the masked array
    cfgAimec['FILTER'] = {'simName': SimDF['simName'].loc[simIndex]}
    log.info('Filter for: simName %s' % SimDF['simName'].loc[simIndex])

    # call full aimec analysis to get resAnalysisDF and newRasters of result fields
    _, resAnalysisDF, _, newRasters, _ = ana3AIMEC.fullAimecAnalysis(avalancheDir, cfgAimec)
    # create dicRaster
    dictRaster = oPartAna.rasterVelField(resAnalysisDF, newRasters,
        pfvMinThreshold=pfvMinThreshold, pftMinThreshold=pftMinThreshold,
        pprMinThreshold=pprMinThreshold)

    # create mtiInfo dicts for tt-diagram
    cfgRangeTime = cfgUtils.getModuleConfig(dtAna, fileOverride='', modInfo=False, toPrint=False,
            onlyDefault=cfgPartAna['distanceTimeAnalysis_override'].getboolean('defaultConfig'))
    cfgRangeTime, cfgPartAna = cfgHandling.applyCfgOverride(cfgRangeTime, cfgPartAna, dtAna, addModValues=False)
    cfgRangeTime['GENERAL']['avalancheDir'] = avalancheDir
    cfgRangeTime, mtiInfo = dtAna.createThalwegTimeInfoFromSimResults(avalancheDir,
        cfgRangeTime, 'com1DFA', simIndex, SimDF, demSim)


    # PLOTTING
    # %% Plotting peak flow quantities and the velocity thalweg diagram
    oPartAna.plotPeakVelVelThalwegEnvelope(avalancheDir, simIndex, SimDF, rasterTransfo,
        dictVelAltThalweg, resTypePlots, anaMod, demSim)

    oPartAna.plotPeakQuantThalTimeEnergyLine(avalancheDir, simIndex, SimDF,
        rasterTransfo, dictRaster, modName, demSim, mtiInfo, cfgRangeTime)

    oPartAna.plotPeakQuantTrackedPartVel(avalancheDir, simName, dictVelAltThalweg,
        dictVelAltThalwegPart, trackedPartProp, dictVelEnvelope, demSim, modName, rasterTransfo)
