"""
    Run script for running the dam break problem on a sloping bed and compare to simulation results
"""

# Load modules
import os
import numpy as np

# Local imports
import avaframe.com1DFAPy.com1DFA as com1DFAPy
from avaframe.com1DFAPy import runCom1DFA
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
initProj.cleanSingleAvaDir(avaDir, keep=logName)

# Start logging
log = logUtils.initiateLogger(avaDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avaDir)


# Load configuration
damBreakCfg = os.path.join(avaDir, 'Inputs', 'damBreak_com1DFACfg.ini')
cfg = cfgUtils.getModuleConfig(com1DFAPy, damBreakCfg)
cfgGen = cfg['GENERAL']

# Load flow depth from analytical solution
hL, hR, uR, phi, xR = damBreak.damBreakSol(avaDir, cfgMain, cfg)
xR = xR * np.cos(phi)  # projected on the horizontal plane
dtAnalysis = cfg['DAMBREAK'].getfloat('dtStep')

# call com1DFAPy to perform simulation - provide configuration file and release thickness function
Particles, Fields, Tsave = runCom1DFA.runCom1DFAPy(avaDir=avaDir, cfgFile=damBreakCfg, flagAnalysis=True)

# create simDict of results
inputDir = 'data/avaDamBreak/Outputs/com1DFAPy/peakFiles/'
dataComSol = fU.makeSimDict(inputDir, avaDir=avaDir)

# make comparison plots
damBreak.plotComparison(dataComSol, hL, xR, hR, uR, dtAnalysis, cfgMain)
