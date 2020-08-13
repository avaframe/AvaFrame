"""Run script for module template"""

# Local imports
from avaframe.tmp1Ex import tmp1Ex
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# -----------Required settings-----------------

# Avalanche directory; see doc.avaframe.org for setup
avalancheDir = './'

# log file name; leave empty to use default runLog.log
logName = 'runTmp1Ex'

# ---------------------------------------------

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(tmp1Ex)

# Dump config to log file
logUtils.writeCfg2Log(cfg,'tmp1Ex')

# Different ways to call functions
tmp1Ex.tmp1ExMain(cfg)

