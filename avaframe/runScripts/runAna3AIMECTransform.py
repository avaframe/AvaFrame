"""
    Run ana3AIMEC tranform to add the projected s and l in the aimec coordinate system
"""
# Load modules
import logging
import copy

# Local imports
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
from avaframe.in3Utils import initializeProject as iP
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.com1DFA import particleTools as paT
from avaframe.com1DFA.com1DFA import savePartToPickle
import avaframe.in2Trans.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)

# -----------Required settings-----------------
# log file name; leave empty to use default runLog.log
logName = 'runAna3AIMECTransform'

# ---------------------------------------------
import pathlib
import avaframe.in3Utils.fileHandlerUtils as fU


def runAna3AIMECTransform(avalancheDir, particlesProperties, cfg):
    """ run script for AIMEC analysis
    proceeds to AIMEC analysis and add projected r and l to the particles dictionary
    Also adds the angle beta and the s of the beta point to the particle file
    """

    cfgSetup = cfg['AIMECSETUP']
    anaMod = cfgSetup['anaMod']

    # Setup input from computational module
    inputsDF, resTypeList = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfg)
    # define reference simulation
    refSimRowHash, refSimName, inputsDF, colorParameter, valRef = aimecTools.fetchReferenceSimNo(avalancheDir, inputsDF, anaMod,
                                                                                         cfg)
    pathDict = {'refSimRowHash': refSimRowHash, 'refSimName': refSimName, 'compType': ['singleModule', anaMod],
                'colorParameter': colorParameter, 'resTypeList': resTypeList, 'valRef': valRef}
    pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)
    pathDict = aimecTools.checkAIMECinputs(cfgSetup, pathDict)
    log.info("Running ana3AIMEC model on test case DEM \n %s \n with profile \n %s ",
             pathDict['demSource'], pathDict['profileLayer'])
    # Run AIMEC postprocessing
    rasterTransfo, resAnalysisDF, plotDict, newRasters = ana3AIMEC.mainAIMEC(pathDict, inputsDF, cfg)

    #DEM
    demSource = pathDict['demSource']
    dem = IOf.readRaster(demSource)

    # fetch all flow parameter result fields
    modName = 'com1DFA'
    configDir = pathlib.Path(avalancheDir, 'Outputs', modName, 'configurationFiles')
    if (configDir.is_dir() is False) or ( len(list(configDir.glob('*.ini'))) == 0):
        fU.fileNotFoundMessage(('No configuration files found in %s - consider first running avalanche simulations' %
            configDir))

    # load info on particles
    inputDir = pathlib.Path(avalancheDir,'Outputs',modName, 'particles')
    particleList, timeStepInfo = paT.readPartFromPickle(inputDir, simName='', flagAvaDir=False, comModule='com1DFA')

    # Create new particles list by copying the already existing one
    particlesListDict  = copy.copy(particleList)

    # Prepare output direction
    outDir = pathlib.Path(avalancheDir, 'Outputs', 'Aimec')
    outDirPicData = outDir / 'particles'
    fU.makeADir(outDirPicData)

    log.info('Perform particles projection')
    for particle in particlesListDict:
        particle = ana3AIMEC.aimecTransform(rasterTransfo, particle, dem['header'])
        particle['sBetaPoint'] = rasterTransfo['s'][rasterTransfo['indStartOfRunout']]
        particle['beta'] = rasterTransfo['startOfRunoutAreaAngle']
        simName = particle['simName']
        # save pickle file
        log.info('Saving pickle file')
        savePartToPickle(particle, outDirPicData, simName)

    # save csv files
    log.info('Saving csv files')
    paT.savePartToCsv(particlesProperties[0], particlesListDict, outDir)

    return particlesListDict


if __name__ == '__main__':
    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Load configuration info of all com1DFA simulations
    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)
    particlesProperties = simDF['particleProperties']

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
    runAna3AIMECTransform(avalancheDir, particlesProperties, cfg)
