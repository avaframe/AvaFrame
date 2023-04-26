"""
    Run script for creating a release area info file
"""

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in1Data.getInput as gI


# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runReleaseInfo'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration file for com1DFA module
cfg = cfgUtils.getModuleConfig(com1DFA)

# fetch input data and create release info dataframe and csv file
relDFDict = gI.createReleaseStats(avalancheDir, cfg)
