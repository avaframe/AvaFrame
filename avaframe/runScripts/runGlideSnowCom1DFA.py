"""
Run the glide snow tool of com1DFA
"""

import pathlib

# Local imports
# import config and init tools
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initializeProject as initProj
# import computation modules
from avaframe.com1DFA import com1DFA
from avaframe.com1DFA import glideSnow

# +++++++++REQUIRED+++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runGlideSnowCom1DFA'
# ++++++++++++++++++++++++++++++

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Clean input directory of old work and output files from module
initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=False)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# load glide snow tool config
glideSnowCfg = cfgUtils.getModuleConfig(glideSnow)

# perform com1DFA simulation with glide snow settings
glideSnow.runGlideSnowTool(avalancheDir, cfgMain, glideSnowCfg)
