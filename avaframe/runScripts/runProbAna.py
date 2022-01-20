"""
    Run script for performing an avalanche simulation with parameter variation and performing a probability analysis with the simulation results
"""

# Load modules
import time
import glob
import pathlib

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.out3Plot import statsPlots as sP
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.ana4Stats import probAna
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.fileHandlerUtils as fU


# log file name; leave empty to use default runLog.log
logName = 'runCom1DFAandProbAna'

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
    avaDir = pathlib.Path(avaDir)
    avaName = avaDir.name
    probSimCfg = pathlib.Path('..', 'benchmarks', '%sStatsTest' % avaName, '%sProbAna_com1DFACfg.ini' % avaName)
    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # Run Standalone DFA
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avaDir, cfgMain, cfgFile=probSimCfg)

    # Load input parameters from configuration file
    cfgProb = cfgUtils.getModuleConfig(probAna)

    # provide optional filter criteria for simulations
    parametersDict = fU.getFilterDict(cfgProb, 'FILTER')

    # perform probability analysis
    probAna.probAnalysis(avaDir, cfgProb, com1DFA, parametersDict=parametersDict)

    # make a plot of the map
    inputDir = pathlib.Path(avaDir, 'Outputs', 'ana4Stats')
    sP.plotProbMap(avaDir, inputDir, cfgProb)
