"""
    Run script for getting rename dataframe, no actuall renameing is done!
"""

# Load modules
import logging
import pathlib

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import logUtils

# -------------Required settings ------
# which parameters to add to simulation name
csvString = 'mu,tau0'
# -------------Required settings ------

# log file name; leave empty to use default runLog.log
logName = 'runRenaming'

#  Load general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avaDir = cfgMain['MAIN']['avalancheDir']


log = logUtils.initiateLogger(avaDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avaDir)

avaDir = pathlib.Path(avaDir)

renameDF = cfgHandling.addInfoToSimName(avaDir, csvString)

print(renameDF)
