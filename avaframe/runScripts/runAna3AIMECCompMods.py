"""
    Run ana3AIMEC

    This file is part of Avaframe.

"""

# Load modules
import time

# Local imports
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
from avaframe.in3Utils import initializeProject as iP
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def runAna3AIMECCompMods(avalancheDir=''):
    """ run script for AIMEC module comparison analysis
    reads the avalancheDir from the configuration file or the one given in input
    proceeds to AIMEC analysis and produces plots and reports
    """
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

    # compList = [['relAlr', 'null'], ['relKot', 'ent']]
    compList = [['relAlr', 'ent']]

    for comp in compList:

        # get config
        rel = comp[0]
        simType = comp[1]

        # Setup input from com1DFA
        pathDict = dfa2Aimec.dfaComp2Aimec(avalancheDir, cfg, rel, simType)

        # TODO: define referenceFile
        pathDict['numSim'] = len(pathDict['ppr'])
        pathDict['referenceFile'] = 0

        # Extract input file locations
        pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=rel+'_'+simType)

        startTime = time.time()

        log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
                 pathDict['demSource'], pathDict['profileLayer'])
        # Run AIMEC postprocessing
        ana3AIMEC.mainAIMEC(pathDict, cfg)

        endTime = time.time()

        log.info(('Took %s seconds to calculate.' % (endTime - startTime)))


if __name__ == '__main__':
    runAna3AIMECCompMods()
