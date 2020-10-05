"""
    Run script for running Standalone DFA
    This file is part of Avaframe.
"""

# Load modules
import os

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import time

# log file name; leave empty to use default runLog.log
logName = 'runCom1DFA'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load input parameters from configuration file
cfg = cfgUtils.getModuleConfig(com1DFA)

# write config to log file
logUtils.writeCfg2Log(cfg, 'com1DFA')

startTime = time.time()
# Run Standalone DFA
com1DFA.runSamos(cfg, avalancheDir)

endTime = time.time()

log.info(('Took %s seconds to calculate.' % (endTime - startTime)))
