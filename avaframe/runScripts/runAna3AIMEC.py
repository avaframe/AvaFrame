"""
    Run ana3AIMEC
"""
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
    # Load avalanche directory from general configuration file or if available from inputs
    if avalancheDir != '':
        cfgMain = cfgUtils.getGeneralConfig()
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

    # Setup input from computational module
    inputsDF, resTypeList = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfg)
    # define reference simulation
    refSimHash, refSimName, inputsDF, colorParameter = aimecTools.fetchReferenceSimNo(avalancheDir, inputsDF, anaMod,
                                                                                      cfgSetup)
    pathDict = {'refSimHash': refSimHash, 'refSimName': refSimName, 'compType': ['singleModule', anaMod],
                'colorParameter': colorParameter}
    pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)
    pathDict['resTypeList'] = resTypeList
    pathDict = aimecTools.checkAIMECinputs(cfgSetup, inputsDF, pathDict)
    log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
             pathDict['demSource'], pathDict['profileLayer'])
    # Run AIMEC postprocessing
    rasterTransfo, resAnalysisDF, plotDict = ana3AIMEC.mainAIMEC(pathDict, inputsDF, cfg)
    return pathDict, rasterTransfo, resAnalysisDF


if __name__ == '__main__':
    runAna3AIMEC()
