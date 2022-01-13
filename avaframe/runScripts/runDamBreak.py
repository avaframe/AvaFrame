"""
    Run script for running the dam break problem on a sloping bed and compare to simulation results
"""

# Load modules
import os
import numpy as np
import pathlib

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.ana1Tests import damBreak
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runDamBreakProblem'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avaDir = 'data/avaDamBreak'
cfgMain['MAIN']['avalancheDir'] = avaDir

# Clean input directory(ies) of old work and output files
initProj.cleanModuleFiles(avaDir, com1DFA)

# Start logging
log = logUtils.initiateLogger(avaDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avaDir)


# Load configuration
damBreakCfg = os.path.join(avaDir, 'Inputs', 'damBreak_com1DFACfg.ini')
cfg = cfgUtils.getModuleConfig(com1DFA, damBreakCfg)
cfgGen = cfg['GENERAL']

# Load flow depth from analytical solution
hL, hR, uR, phi, xR = damBreak.damBreakSol(avaDir, cfgMain, cfg)
xR = xR * np.cos(phi)  # projected on the horizontal plane
dtAnalysis = cfg['DAMBREAK'].getfloat('dtStep')

# call com1DFAPy to perform simulation - provide configuration file and release thickness function
dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avaDir, cfgMain, cfgFile=damBreakCfg)

# create dataFrame of results
inputDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'peakFiles', 'timeSteps')
dataComSolDF = fU.makeSimDF(inputDir, avaDir=avaDir)

# make comparison plots
damBreak.plotComparison(dataComSolDF, hL, xR, hR, uR, dtAnalysis, cfgMain)
