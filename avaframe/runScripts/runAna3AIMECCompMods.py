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

    simTypeList = ['null']

    for simType in simTypeList:

        # Setup input from com1DFA
        inputsDF, pathDict = dfa2Aimec.dfaBench2Aimec(avalancheDir, cfg, simType, simType)

        # TODO: define referenceFile
        pathDict['refSimulation'] = inputsDF.index[0]

        # Extract input file locations
        pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName='aimec_'+simType)

        startTime = time.time()

        log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
                 pathDict['demSource'], pathDict['profileLayer'])
        # Run AIMEC postprocessing
        rasterTransfo, resAnalysisDF = ana3AIMEC.mainAIMEC(pathDict, inputsDF, cfg)

        endTime = time.time()

        log.info(('Took %s seconds to calculate.' % (endTime - startTime)))


if __name__ == '__main__':
    runAna3AIMECCompMods()
