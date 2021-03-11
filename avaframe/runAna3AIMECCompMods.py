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
from avaframe.in3Utils import initializeProject as iP
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# -----------Required settings-----------------
# log file name; leave empty to use default runLog.log
logName = 'runAna3AIMECCompMods'

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

iP.cleanModuleFiles(avalancheDir, ana3AIMEC)

# write configuration to file
cfgUtils.writeCfgFile(avalancheDir, ana3AIMEC, cfg)

cfgSetup = cfg['AIMECSETUP']

# Setup input from com1DFA
pathDict = dfa2Aimec.dfaComp2Aimec(avalancheDir, cfgSetup)

for simN in pathDict:

    pDict = pathDict[simN]
    
    # TODO rename these
    cfgPath = {}
    cfgPath['pressurefileList'] = pathDict['ppr']
    cfgPath['depthfileList'] = pathDict['pfd']
    cfgPath['massfileList'] = pathDict['mb']
    cfgPath['speedfileList'] = pathDict['pfv']
    cfgPath['numSim'] = len(cfgPath['pressurefileList'])

    # TODO define referenceFile
    cfgPath['referenceFile'] = 0

    # Extract input file locations
    cfgPath = ana3AIMEC.readAIMECinputs(avalancheDir, cfgPath, dirName=simN)

    startTime = time.time()

    log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
             cfgPath['demSource'], cfgPath['profileLayer'])
    # Run AIMEC postprocessing
    ana3AIMEC.mainAIMEC(cfgPath, cfg)

    endTime = time.time()

    log.info(('Took %s seconds to calculate.' % (endTime - startTime)))
