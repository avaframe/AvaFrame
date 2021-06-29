"""
    Run script for performing an avalanche simulation with parameter variation and performing a probability analysis with the simulation results
"""

# Load modules
import os
import time
import glob
import pathlib

# Local imports
from avaframe.runCom1DFA import runCom1DFA
from avaframe.out3Plot import statsPlots as sP
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.ana4Stats import probAna
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.com1DFA.com1DFA as com1DFA


# log file name; leave empty to use default runLog.log
logName = 'runCom1DFAandProbAna'
modName = 'com1DFA'

# Load general configuration filee
cfgMain = cfgUtils.getGeneralConfig()

# Define avalanche directories for tests
avalancheDirectories = ['data/avaHockeyChannel']

# run simulations sequentially
for avaDir in avalancheDirectories:

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avaDir)

    # Load input parameters from configuration file
    # write config to log file
    avaName = os.path.basename(avaDir)
    probSimCfg = os.path.join('..', 'benchmarks', '%sStatsTest' % avaName, '%sProbAna_com1DFACfg.ini' % avaName)
    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # Run Standalone DFA
    particlesList, fieldsList, Tsave, dem, plotDict, reportDictList = runCom1DFA(avaDir=avaDir, cfgFile=probSimCfg, relThField='', variationDict='')

    # Load input parameters from configuration file
    cfgProb = cfgUtils.getModuleConfig(probAna)

    # provide optional filter criteria for simulations
    parametersDict = fU.getFilterDict(cfgProb, 'FILTER')

    # perform probability analysis
    probAna.probAnalysis(avaDir, cfgProb, com1DFA, parametersDict=parametersDict)

    # make a plot of the map
    inputDir = pathlib.Path(avaDir, 'Outputs', 'ana4Stats')
    sP.plotProbMap(avaDir, inputDir, cfgProb)
