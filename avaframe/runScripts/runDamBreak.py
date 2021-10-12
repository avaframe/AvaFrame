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
initProj.cleanModuleFiles(avaDir, com1DFA, 'com1DFA')

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

# call com1DFAPy to perform simulation - provide configuration file and release thickness function
particlesList, fieldsList, Tsave, dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avaDir, cfgMain, cfgFile=damBreakCfg)

# create dataFrame of results
inputDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'peakFiles', 'timeSteps')
dataComSolDF = fU.makeSimDF(inputDir, avaDir=avaDir)


fieldsList, fieldHeader = com1DFA.readFields(avaDir, ['FD', 'FV'], simName='', flagAvaDir=True, comModule='com1DFA')

# user interaction?
if cfg['DAMBREAK'].getboolean('flagInteraction'):
    dtAnalysis = input("give time step to plot (float in s):\n")
else:
    dtAnalysis = cfg['DAMBREAK'].getfloat('dtStep')

try:
    dtAnalysis = float(dtAnalysis)
except ValueError:
    dtAnalysis = 'n'
while isinstance(dtAnalysis, float):

    # determine index for time step
    ind_t = min(np.searchsorted(Tsave, dtAnalysis), min(len(Tsave)-1, len(fieldsList)-1))
    dtAnalysis = Tsave[ind_t]
    # ind_time = np.absolute(solSimi['Time'] - Tsave[ind_t]).argmin()
    # make comparison plots
    damBreak.plotComparison(fieldsList[0], fieldsList[ind_t], fieldHeader, hL, xR, hR, uR, dtAnalysis, cfgMain)
    # # option for user interaction
    if cfg['DAMBREAK'].getboolean('flagInteraction'):
        dtAnalysis = input("give time step to plot (float in s):\n")
        try:
            dtAnalysis = float(dtAnalysis)
        except ValueError:
            dtAnalysis = 'n'
    else:
        dtAnalysis = 'n'
