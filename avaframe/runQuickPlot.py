"""
    Run script for plotting two matrix datasets
    This file is part of Avaframe.
"""

# Load modules
import os

# Local imports
from avaframe.out3SimpPlot import outQuickPlot
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

# REQUIRED+++++++++++++++++++
com1DFAOutput = 'FullOutput_mu_'
outputVariable = ['pfd', 'ppr']
simName = 'entres'
#++++++++++++++++++++++++++++

# Plot data comparison for all output variables defined in suffix
for var in outputVariable:
    # Run Standalone DFA
    outQuickPlot.quickPlot(avalancheDir, var, com1DFAOutput, simName)
