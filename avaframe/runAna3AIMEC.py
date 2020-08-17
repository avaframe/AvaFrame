#! /usr/bin/python
"""
    Run AIMEC

    This file is part of Avaframe.

"""

# Load modules
import logging
import time

# Local imports
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# -----------Required settings-----------------

# TODO: move this to a main cfg file!

# Avalanche directory; see doc.avaframe.org for setup
# TODO: full path needed?
avalancheDir = 'data/avaSlide'


# log file name; leave empty to use default runLog.log
logName = 'runAna3AIMEC'

# ---------------------------------------------

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(ana3AIMEC)

cfgSetup = cfg['AIMECSETUP']
cfgDFA = cfg['DFA']

# Dump config to log file
logUtils.writeCfg2Log(cfg, 'ana3AIMEC')

# Setup input from com1DFA
dfa2Aimec.makeAimecDirs(avalancheDir)
dfa2Aimec.getDFAData(avalancheDir, cfgDFA)
dfa2Aimec.writeAimecPathsFile(cfgSetup, avalancheDir)

saveOutPath = avalancheDir + '/Outputs'
# Extract input file locations
cfgPath = ana3AIMEC.readAIMECinputs(avalancheDir)

# write config to log file
logUtils.writeCfg2Log(cfg, 'ana3AIMEC')

startTime = time.time()

log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
         cfgPath['demSource'], cfgPath['profileLayer'])

ana3AIMEC.mainAIMEC(cfgPath, cfg)

endTime = time.time()

log.info(('Took %s seconds to calculate.' % (endTime - startTime)))
