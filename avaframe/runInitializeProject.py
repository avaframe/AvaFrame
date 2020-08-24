"""
    Run script for Initialization
    This file is part of Avaframe.
"""

# Load modules
import os
import shutil

# Local imports
from avaframe.in3Utils import initializeProject
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'initializeProject'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger('.', logName)
log.info('MAIN SCRIPT')
log.info('Initializing Project: %s', avalancheDir)

# Initialize project
initializeProject.initializeFolderStruct(avalancheDir)

logOrigin = os.path.join('.', logName + '.log')

logDest = os.path.join(avalancheDir, logName + '.log')

shutil.move(logOrigin, logDest)
