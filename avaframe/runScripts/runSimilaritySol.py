"""
    Run script for running python DFA kernel and similarity solution test
    This script calculates the similarity solution for a gliding avalanche on
    a inclined plane according to similarity solution from :
    Hutter, K., Siegel, M., Savage, S.B. et al.
    Two-dimensional spreading of a granular avalanche down an inclined plane
    Part I. theory. Acta Mechanica 100, 37–68 (1993).
    https://doi.org/10.1007/BF01176861
    and compares it to the DFA kernel com1DFA
"""

import os
import time
import numpy as np

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in1Data import getInput as gI
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana1Tests.simiSol as simiSol


# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# Load general configuration
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaSimilaritySol'

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

# call com1DFA to perform simulation - provide configuration file and release thickness function
particlesList, fieldsList, Tsave, dem, plotDict, reportDictList = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=simiSolCfg,
relThField=relTh)
relDict['dem'] = dem

# compute similartiy solution
log.info('Computing similarity solution')
solSimi = simiSol.mainSimilaritySol()

hErrorL2Array, vErrorL2Array = simiSol.analyzeResults(particlesList, fieldsList, solSimi, relDict, cfg, outDirTest)

# +++++++++POSTPROCESS++++++++++++++++++++++++
# -------------------------------
if cfgMain['FLAGS'].getboolean('showPlot'):
    simiSol.plotContoursSimiSol(particlesList, fieldsList, solSimi, relDict, cfg, outDirTest)


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
    ind_t = min(np.searchsorted(Tsave, value), len(Tsave)-1)
    ind_time = np.searchsorted(solSimi['Time'], Tsave[ind_t])

    # get similartiy solution h, u at reuired time step
    simiDict = simiSol.getSimiSolParameters(solSimi, relDict, ind_time, cfg)

    # get particle parameters
    comSol = simiSol.prepareParticlesFieldscom1DFA(fieldsList, particlesList, ind_t, relDict, simiDict, 'xaxis')
    comSol['outDirTest'] = outDirTest
    comSol['showPlot'] = cfgMain['FLAGS'].getboolean('showPlot')
    comSol['Tsave'] = Tsave[ind_t]

    # make plot
    simiSol.plotProfilesSimiSol(ind_time, relDict, comSol, simiDict, solSimi, 'xaxis')

    # get particle parameters
    comSol = simiSol.prepareParticlesFieldscom1DFA(fieldsList, particlesList, ind_t, relDict, simiDict, 'yaxis')
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
