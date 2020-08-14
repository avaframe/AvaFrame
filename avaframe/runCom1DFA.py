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


#----------Rewquired Settings---------------

#TODO move this to a main cfg file!

# Avalanche directory; see doc.avaframe.org for setup
# TODO: full path needed?
avalancheDir = 'data/avaSlide'

# log file name; leave empty to use default runLog.log
logName = 'runCom1DFA'
#------------------------------------------


# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Get path of module
modPath = os.path.dirname(com1DFA.__file__)

# Load input parameters from configuration file
cfg = cfgUtils.getModuleConfig(com1DFA)


# Run Standalone DFA
com1DFA.runSamos(cfg, avalancheDir, modPath)
