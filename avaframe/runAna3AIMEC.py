#! /usr/bin/python
"""
    Run ana3AIMEC

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
# log file name; leave empty to use default runLog.log
logName = 'runAna3AIMEC'

# ---------------------------------------------
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
cfg = cfgUtils.getModuleConfig(ana3AIMEC)

cfgSetup = cfg['AIMECSETUP']

# Setup input from com1DFA
dfa2Aimec.mainDfa2Aimec(avalancheDir, cfgSetup)

saveOutPath = avalancheDir + '/Outputs'
# Extract input file locations
cfgPath = ana3AIMEC.readAIMECinputs(avalancheDir, dirName='com1DFA')

startTime = time.time()

log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
         cfgPath['demSource'], cfgPath['profileLayer'])
# Run AIMEC postprocessing
_, _, _ = ana3AIMEC.mainAIMEC(cfgPath, cfg)

endTime = time.time()

log.info(('Took %s seconds to calculate.' % (endTime - startTime)))
