"""Run script for module com4FlowPy
"""

# Local imports
from avaframe.com4FlowPy import com4FlowPy
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runcom4FlowPy'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
# write config to log file
cfg = cfgUtils.getModuleConfig(com4FlowPy)

cfgSetup = cfg['SETUP']
cfgFlags = cfg['FLAGS']

# Extract input file locations
cfgPath = com4FlowPy.readFlowPyinputs(avalancheDir)

log.info("Running com4FlowPyMain model on DEM \n \t %s \n \t with profile \n \t %s ",
         cfgPath['demSource'], cfgPath['releasePath'])

com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)
