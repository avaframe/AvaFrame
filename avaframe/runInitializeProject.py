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


#----------Rewquired Settings---------------

#TODO move this to a main cfg file!

# Avalanche directory; see doc.avaframe.org for setup
# TODO: full path needed?
avalancheDir = 'data/TestAva'

# log file name; leave empty to use default runLog.log
logName = 'initializeProject'
#------------------------------------------


# Start logging
log = logUtils.initiateLogger('./', logName)
log.info('MAIN SCRIPT')
log.info('Initializing Project: %s', avalancheDir)

# Initialize project
initializeProject.initializeFolderStruct(avalancheDir)
logOrigin = os.path.join('./', logName + '.log')
logDest = os.path.join(avalancheDir, logName)
shutil.move(logOrigin, logDest)
