"""
    Run script for getting rename dataframe 
"""

# Load modules
import logging
import pathlib

# Local imports
#  from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runRenaming'

#  Load general configuration file
#  cfgMain = cfgUtils.getGeneralConfig()

avaDir = 'data/avaAlr'

log = logUtils.initiateLogger(avaDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avaDir)

avaDir = pathlib.Path(avaDir)
csvString = 'mu,tau0'

renameDF = cfgHandling.addInfoToSimName(avaDir,csvString)
print(renameDF)
