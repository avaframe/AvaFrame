"""
    Run ana3AIMEC
"""

# Load modules
import time

# Local imports
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
from avaframe.in3Utils import initializeProject as iP
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def runAna3AIMEC(avalancheDir=''):
    """ run script for AIMEC analysis
    reads the avalancheDir from the configuration file or the one given in input
    proceeds to AIMEC analysis and produces plots and reports
    """
    # -----------Required settings-----------------
    # log file name; leave empty to use default runLog.log
    logName = 'runAna3AIMEC'

    # ---------------------------------------------
    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
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
    anaMod = cfgSetup['anaMod']

    # Setup input from com1DFA
    pathDict = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod)

# Setup input from com1DFA
pathDict = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfgSetup)

# TODO: define referenceFile
pathDict['numSim'] = len(pathDict['ppr'])

# define reference simulation
pathDict = aimecTools.fetchReferenceSimNo(pathDict, cfgSetup)

    log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
             pathDict['demSource'], pathDict['profileLayer'])
    # Run AIMEC postprocessing
    ana3AIMEC.mainAIMEC(pathDict, cfg)

    endTime = time.time()

log.info("Running ana3AIMEC model on test case DEM: \n %s \n with profile: \n %s ",
         pathDict['demSource'], pathDict['profileLayer'])
# Run AIMEC postprocessing
ana3AIMEC.mainAIMEC(pathDict, cfg)


if __name__ == '__main__':
    runAna3AIMEC()
