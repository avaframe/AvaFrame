"""
    Run script for plotting two matrix datasets
    This file is part of Avaframe.
"""

# Load modules
import os

# Local imports
from avaframe.out3Plot import outQuickPlot
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'runQuickPlot'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)


# Load input parameters from configuration file
cfg = cfgUtils.getModuleConfig(outQuickPlot)


# REQUIRED+++++++++++++++++++
# Which parameter to filter data, e.g. varPar = 'simType', values = ['null'] or varPar = 'Mu', values = ['0.055', '0.155']
 # values need to be given as list, also if only one value
outputVariable = cfg['PLOT']['outputVariable'].split('_')
values = cfg['PLOT']['values'].split('_')
parameter = cfg['PLOT']['parameter']
#++++++++++++++++++++++++++++

# Plot data comparison for all output variables defined in suffix
for val in values:
    for var in outputVariable:
        outQuickPlot.quickPlot(avalancheDir, var, val, parameter, cfgMain, cfg)
