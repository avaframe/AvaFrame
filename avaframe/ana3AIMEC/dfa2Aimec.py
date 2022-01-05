"""
    Helper function to export required data from com1DFA to be used by Aimec.
"""

# Load modules
import logging
import pathlib
import numpy as np

# local modules
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.out3Plot.plotUtils as pU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def extractCom1DFAMBInfo(avaDir, inputsDF, simNameInput=''):
    """ Extract the mass balance info from the log file """

    if simNameInput != '':
        simNames = [simNameInput]
        nameDir = simNameInput
    else:
        # Get info from ExpLog
        nameDir = 'com1DFAOrig'
        logLoc = pathlib.Path(avaDir, 'Outputs', 'com1DFAOrig')
        logName = logLoc / 'ExpLog.txt'
        logDictExp = fU.readLogFile(logName)
        names = logDictExp['fullName']
        simNames = sorted(set(names), key=lambda s: (s.split("_")[0], s.split("_")[1], s.split("_")[3]))

    # Read mass data from log and save to file for each simulation run
    for simName in simNames:
        log.debug('This is the simulation name %s for mod com1DFAOrig ' % (simName))

        # Read log file and extract mass and time info
        locFiles = pathlib.Path(avaDir, 'Outputs', 'com1DFAOrig')
        fileName = locFiles / ('start%s.log' % (simName))
        fU.checkIfFileExists(fileName, fileType='mass log info')
        logDict = fU.extractLogInfo(fileName)
        indRun = logDict['indRun']

        # Write mass balance info files
        for k in range(len(indRun)-1):
            saveName = locFiles / ('mass_%s.txt' % (simName))
            with open(saveName, 'w') as MBFile:
                MBFile.write('time, current, entrained\n')
                for m in range(indRun[k], indRun[k] + indRun[k+1] - indRun[k]-1):
                    MBFile.write('%.02f,    %.06f,    %.06f\n' %
                                 (logDict['time'][m], logDict['mass'][m], logDict['entrMass'][m]))
            inputsDF.loc[simName, 'massBal'] = saveName
            log.debug('Mass file saved to %s ' % (saveName))
            log.info('Added to pathDict[massBal] %s ' % (saveName))

    return inputsDF


def getMBInfo(avaDir, inputsDF, comMod, simName=''):
    """ Get mass balance info """

    # Get info from ExpLog
    if simName != '':
        mbFile = pathlib.Path(avaDir, 'Outputs', comMod, 'mass_%s.txt' % simName)
        fU.checkIfFileExists(mbFile, fileType='mass')
        inputsDF.loc[simName, 'massBal'] = mbFile
        log.info('Added to pathDict[massBal] %s' % (mbFile))

    else:
        dir = pathlib.Path(avaDir, 'Outputs', comMod)
        mbFiles = list(dir.glob('mass*.txt'))
        if len(mbFiles) == 0:
            message = 'No mass log file found in directory %s' % (str(dir))
            log.error(message)
            raise FileNotFoundError(message)
        mbNames = sorted(set(mbFiles), key=lambda s: (str(s).split("_")[1], str(s).split("_")[2], str(s).split("_")[4]))

        for mFile in mbNames:
            name = mFile.stem
            nameParts = name.split('_')
            simName = ('_'.join(nameParts[1:]))
            inputsDF.loc[simName, 'massBal'] = mFile
            log.debug('Added to pathDict[massBal] %s' % (mFile))
    return inputsDF


def getRefMB(testName, pathDict, simName):
    """ Get mass balance info """

    # Get info from ExpLog
    mbFile = pathlib.Path('..', 'benchmarks', testName, 'mass_%s.txt' % simName)
    pathDict['massBal'].append(mbFile)
    log.info('Added to pathDict[massBal] %s' % (mbFile))

    return pathDict


def dfaComp2Aimec(avaDir, cfg, rel, simType):
    """ Create a pathDict where the paths to all the files required by aimec are saved for two modules -
        in order to compare always two simulations at a time

        for now matching simulations are identified via releaseScenario and simType

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        cfgSetup: configParser object
            configuration for aimec
        rel: str
            releaseScenario
        simType: str
            simulation type (null, ent, entres, ..)

        Returns
        --------
        pathDict: dict
            dictionary with paths to result and optionally mass files for matching simulations from
            two modules - matching in terms of releaseScenario and simType
    """

    cfgSetup = cfg['AIMECSETUP']
    comModules = cfg['AIMECSETUP']['comModules'].split('|')
    # get directories where simulation results can be found for both modules
    inputDirRef, inputDirComp, pathDict = getCompDirs(avaDir, cfgSetup)

    # Load all infos on reference simulations
    refData = fU.makeSimDF2(avaDir, None, inputDir=inputDirRef)
    refData = refData[(refData['releaseArea']==rel) & (refData['simType']==simType)]
    # Load all infos on comparison module simulations
    compData = fU.makeSimDF2(avaDir, None, inputDir=inputDirComp)
    compData = compData[(compData['releaseArea']==rel) & (compData['simType']==simType)]

    # merge dataFrames and locate the rows with common releaseArea ans simType
    dfMerged = refData[['simName', 'releaseArea','simType']].merge(compData[['simName', 'releaseArea','simType']], on=['releaseArea','simType'],
                   how='inner', indicator=True, suffixes=('Ref', 'Comp'))
    dfMerged = dfMerged[dfMerged['_merge'] == 'both']
    if dfMerged.shape[0] > 0:
        refSimName = dfMerged.loc[0, 'simNameRef']
        compSimName = dfMerged.loc[0, 'simNameComp']
        log.info('Reference simulation: %s and to comparison simulation: %s ' % (refSimName, compSimName))
    else:
        message = 'No matching simulations found for reference and comparison simulation \
                   for releaseScenario: %s and simType: %s' % (rel, simType)
        log.error(message)
        raise FileNotFoundError(message)

    # build input dataFrame
    inputsDF = refData[refData['simName'] == refSimName]
    inputsDF = inputsDF.append(compData[compData['simName'] == compSimName])

    # if desired set path to mass log files
    comModules = pathDict['compType']
    sims = {comModules[1]: refSimName, comModules[2]: compSimName}
    if cfg['FLAGS'].getboolean('flagMass'):
        for comMod, sim in sims.items():
            log.info('mass file for comMod: %s and sim: %s' % (comMod, sim))
            if comMod == 'benchmarkReference':
                inputsDF = getRefMB(cfg['AIMECSETUP']['testName'], inputsDF, sim)
            elif comMod == 'com1DFAOrig':
                inputsDF = extractCom1DFAMBInfo(avaDir, inputsDF, simNameInput=sim)
            else:
                inputsDF = getMBInfo(avaDir, inputsDF, comMod, simName=sim)
    return inputsDF, pathDict


def dfaBench2Aimec(avaDir, cfg, simNameRef, simNameComp):
    """ Exports the required data from com1DFA to be used by Aimec

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        cfg: confiParser object
            configuration settings for aimec
        simNameRef: str
            name of reference simulation results
        simNameComp: str
            name of comparison simulation results

        Returns
        --------
        pathDict: dict
            dictionary with paths to simulation results

    """

    cfgSetup = cfg['AIMECSETUP']
    comModules = cfg['AIMECSETUP']['comModules'].split('|')
    # get directories where simulation results can be found for both modules
    inputDirRef, inputDirComp, pathDict = getCompDirs(avaDir, cfgSetup)

    # Load all infos on reference simulations
    refData = fU.makeSimDF2(avaDir, None, inputDir=inputDirRef, simName=simNameRef)

    # Load all infos on comparison module simulations
    compData = fU.makeSimDF2(avaDir, None, inputDir=inputDirComp, simName=simNameComp)

    # build input dataFrame
    inputsDF = refData.append(compData)

    # if desired set path to mass log files
    comModules = pathDict['compType']
    sims = {comModules[1]: simNameRef, comModules[2]: simNameComp}
    if cfg['FLAGS'].getboolean('flagMass'):
        for comMod, sim in sims.items():
            log.info('mass file for comMod: %s and sim: %s' % (comMod, sim))
            if comMod == 'benchmarkReference':
                inputsDF = getRefMB(cfg['AIMECSETUP']['testName'], inputsDF, sim)
            elif comMod == 'com1DFAOrig':
                inputsDF = extractCom1DFAMBInfo(avaDir, inputsDF, simNameInput=sim)
            else:
                inputsDF = getMBInfo(avaDir, inputsDF, comMod, simName=sim)
    return inputsDF, pathDict


def getCompDirs(avaDir, cfgSetup):
    """ Determine dictionaries where simulation results can be found

        Parameters
        ----------
        avaDir: str
            path to avalanche directory
        cfgSetup: configParser object
            configuration settings for aimec

        Returns
        --------
        inputDirRef: str
            path to reference simulation results
        inputDirComp: str
            path to comparison simulation results
        pathDict: dict
            dictionary with paths for aimec
        refModule: str
            name of reference module
    """

    # look for matching simulations
    comModules = cfgSetup['comModules'].split('|')
    refModule = comModules[0]
    compModule = comModules[1]
    log.info('Reference data is from module: %s' % refModule)
    log.info('Comparison data is from module: %s' % compModule)

    # Lead all infos on refernce simulations
    if refModule == 'benchmarkReference':
        testName = cfgSetup['testName']
        inputDirRef = pathlib.Path('..', 'benchmarks', testName)
    else:
        inputDirRef = pathlib.Path(avaDir, 'Outputs', refModule, 'peakFiles')

    inputDirComp = pathlib.Path(avaDir, 'Outputs', compModule, 'peakFiles')
    pathDict = {}
    pathDict['compType'] = ['comModules', refModule, compModule]
    # info about colourmap
    pathDict['contCmap'] = True

    return inputDirRef, inputDirComp, pathDict


def mainDfa2Aimec(avaDir, comModule, cfg):
    """ Fetch available raster results path from com1DFA to be used by Aimec

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        comModule: str
            computational module name that has been used to produce simulation results
        cfg: configParser object
            configuration for aimec

        Returns
        --------
        inputsDF: dataFrame
            with a line for each simulation available and the corresponding simulation results: ppr, pfd, pfv
    """

    inputsDF = fU.makeSimDF2(avaDir, comModule)

    if cfg['FLAGS'].getboolean('flagMass'):
        # Extract mb info
        if comModule == 'com1DFAOrig':
            # path dictionary for Aimec
            inputsDF = extractCom1DFAMBInfo(avaDir, inputsDF)
        else:
            inputsDF = getMBInfo(avaDir, inputsDF, comModule)
    return inputsDF
