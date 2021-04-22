"""
    Run script for running python DFA kernel and similarity solution test
"""

import os
import time
import numpy as np

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA
from avaframe.com1DFAPy import runCom1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana1Tests.simiSol as simiSol


# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# Load general configuration
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaSimilaritySol'
# set module name, reqiured as long we are in dev phase
# - because need to create e.g. Output folder for com1DFAPy to distinguish from
# current com1DFA
modName = 'com1DFAPy'

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for similarity solution test
simiSolCfg = os.path.join(avalancheDir, 'Inputs', 'simiSol_com1DFACfg.ini')
cfg = cfgUtils.getModuleConfig(com1DFA, simiSolCfg)

# for timing the sims
startTime = time.time()

# create output directory for test result plots
outDirTest = os.path.join(avalancheDir, 'Outputs', 'ana1Tests')
fU.makeADir(outDirTest)

# Define release thickness distribution
demFile, relFiles, entFiles, resFile, flagEntRes = gI.getInputData(avalancheDir, cfg['FLAGS'])
relDict = simiSol.getReleaseThickness(avalancheDir, cfg, demFile)
relTh = relDict['relTh']

# call com1DFAPy to perform simulation - provide configuration file and release thickness function
Particles, Fields, Tsave, dem, plotDict, reportDictList = runCom1DFA.runCom1DFAPy(avaDir=avalancheDir, cfgFile=simiSolCfg, relThField=relTh)

# compute similartiy solution
log.info('Computing similarity solution')
solSimi = simiSol.runSimilarity()

# +++++++++POSTPROCESS++++++++++++++++++++++++
# -------------------------------
if cfgMain['FLAGS'].getboolean('showPlot'):
    simiSol.plotContoursSimiSol(Particles, Fields, solSimi, relDict, cfg, outDirTest)


# TODO here is still user interaction
# option for user interaction
if cfg['SIMISOL'].getboolean('flagInteraction'):
    value = input("give time step to plot (float in s):\n")
else:
    value = cfg['SIMISOL'].getfloat('dtSol')

try:
    value = float(value)
except ValueError:
    value = 'n'
while isinstance(value, float):

    # determine index for time step
    ind_t = np.searchsorted(Tsave, value)
    ind_time = np.searchsorted(solSimi['Time'], value)

    # get similartiy solution h, u at reuired time step
    simiDict = simiSol.getSimiSolParameters(solSimi, relDict, ind_time, cfg)

    # get particle parameters
    comSol = simiSol.prepareParticlesFieldscom1DFAPy(Fields, Particles, ind_t, relDict, simiDict, 'xaxis')
    comSol['outDirTest'] = outDirTest
    comSol['showPlot'] = cfgMain['FLAGS'].getboolean('showPlot')
    comSol['Tsave'] = Tsave[ind_t]

    # make plot
    simiSol.plotProfilesSimiSol(ind_time, relDict, comSol, simiDict, solSimi, 'xaxis')

    # get particle parameters
    comSol = simiSol.prepareParticlesFieldscom1DFAPy(Fields, Particles, ind_t, relDict, simiDict, 'yaxis')
    comSol['outDirTest'] = outDirTest
    comSol['showPlot'] = cfgMain['FLAGS'].getboolean('showPlot')
    comSol['Tsave'] = Tsave[ind_t]

    # make plot
    simiSol.plotProfilesSimiSol(ind_time, relDict, comSol, simiDict, solSimi, 'yaxis')

    # # option for user interaction
    if cfg['SIMISOL'].getboolean('flagInteraction'):
        value = input("give time step to plot (float in s):\n")
        try:
            value = float(value)
        except ValueError:
            value = 'n'
    else:
        value = 'n'
