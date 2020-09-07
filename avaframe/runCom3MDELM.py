"""Run script for module com3MDELM
"""

# Local imports
from avaframe.com3MDELM import com3MDELM
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runCom3MDELM'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(com3MDELM)

cfgPath = cfg['PATH']
cfgPath['saveOutPath'] = avalancheDir
cfgSetup = cfg['MDELMSETUP']
cfgFlags = cfg['FLAGS']

# Dump config to log file
logUtils.writeCfg2Log(cfg,'com3MDELM')

# Calculate MDELM
com3MDELM.com3MDELMMain(cfgPath, cfgSetup)
