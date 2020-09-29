"""
    Run script for plotting two matrix datasets
    This file is part of Avaframe.
"""

# Load modules
import os

# Local imports
from avaframe.out3Plot import outQuickPlotSimple
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'runQuickPlotSimple'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Set directory where input files is located
inputDir = os.path.join(avalancheDir, 'Work', 'simplePlot')
outQuickPlotSimple.quickPlotSimple(avalancheDir, inputDir)
