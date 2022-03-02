"""
    Run script for running the dam break problem on a sloping bed and compare to simulation results
"""

# Load modules
import numpy as np
import pathlib
import pandas as pd

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.ana1Tests import damBreak
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.out3Plot.outAna1Plots as outAna1Plots



# log file name; leave empty to use default runLog.log
logName = 'runAnalyzeDamBreakProblem'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaDamBreak'
cfgMain['MAIN']['avalancheDir'] = avalancheDir

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# create output directory for test result plots
outDirTest = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests', 'damBreak')
fU.makeADir(outDirTest)

# Load configuration
damBreakCfg = pathlib.Path(avalancheDir, 'Inputs', 'damBreak_com1DFACfg.ini')
cfg = cfgUtils.getModuleConfig(com1DFA, damBreakCfg)
cfgGen = cfg['GENERAL']

# Load flow depth from analytical solution
solDam = damBreak.damBreakSol(avalancheDir, cfgMain, cfg, outDirTest)

# create dataFrame of results
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

pathToResults = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests', 'damBreak', 'results10.p')
if pathToResults.is_file():
    simDF = pd.read_pickle(pathToResults)
# simDF = damBreak.postProcessSimiSol(avalancheDir, cfgMain, cfg['DAMBREAK'], simDF, solDam, outDirTest)

simDF = simDF[simDF['iniStep']==True]

# make convergence plot
fig1, ax1, ax2, slopeU, slopeH = outAna1Plots.plotErrorConvergence(simDF, outDirTest, cfg['DAMBREAK'], 'nPart', ['hErrorL2', 'vhErrorL2'],
                          'aPPK', 'nPPK0', logScale=True, fit=True)
