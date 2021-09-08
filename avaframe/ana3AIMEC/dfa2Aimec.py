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


def extractCom1DFAMBInfo(avaDir, pathDict, simNameInput=''):
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
    countFile = 0
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
            if simNameInput != '':
                pathDict['massBal'].append(saveName)
            else:
                pathDict['massBal'].append(saveName)
            log.debug('Mass file saved to %s ' % (saveName))
            log.info('Added to pathDict[massBal] %s ' % (saveName))
            countFile = countFile + 1

    return pathDict


def getMBInfo(avaDir, pathDict, comMod, simName=''):
    """ Get mass balance info """

    # Get info from ExpLog
    if simName != '':
        mbFile = pathlib.Path(avaDir, 'Outputs', comMod, 'mass_%s.txt' % simName)
        fU.checkIfFileExists(mbFile, fileType='mass')
        pathDict['massBal'].append(mbFile)
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
            pathDict['massBal'].append(mFile)
            log.debug('Added to pathDict[massBal] %s' % (mFile))

    return pathDict


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

    # initialise empty pathDict for all required files
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'massBal': []}
    # get directories where simulation results can be found for both modules
    inputDirRef, inputDirComp, pathDict, refModule = getCompDirs(avaDir, cfgSetup, pathDict)

    # Load all infos on reference simulations
    refData = fU.makeSimDF(inputDirRef)

    # Load all infos on comparison module simulations
    compData = fU.makeSimDF(inputDirComp)

    # identify names of simulations that match criteria for both cases
    simSearch = True
    for countRef, simNameShort in enumerate(refData['simName']):
        for countComp, simNameComp in enumerate(compData['simName']):
            if simSearch:
                if (refData['releaseArea'][countRef] == compData['releaseArea'][countComp] == rel and
                    refData['simType'][countRef] == compData['simType'][countComp] == simType):
                    refSimName = refData['simName'][countRef]
                    compSimName = compData['simName'][countComp]
                    log.info('Reference simulation: %s and to comparison simulation: %s ' % (refSimName, compSimName))
                    simSearch = False
    if simSearch == True:
        message = 'No matching simulations found for reference and comparison simulation \
                   for releaseScenario: %s and simType: %s' % (rel, simType)
        log.error(message)
        raise FileNotFoundError(message)

    # fill pathDict
    pathDict = getPathsFromSimName(pathDict, avaDir, cfg, inputDirRef, refSimName, inputDirComp, compSimName)

    return pathDict


def getPathsFromSimName(pathDict, avaDir, cfg, inputDirRef, simNameRef, inputDirComp, simNameComp):
    """ Set paths of reference and comparison files

        Parameters
        -----------
        pathDict: dict
            dictionary where paths should be filled
        avaDir: str
            path to avalanche directory
        cfg: configparser object
            configuration for aimec
        inputDirRef: str
            path to reference simulation results
        simNameRef: str
            name of reference simulation results
        inputDirComp: str
            path to comparison simulation results
        simNameComp: str
            name of comparison simulation results

        Returns
        -------
        pathDict: dictionary
            dictionary with paths to simulation results

        """

    comModules = cfg['AIMECSETUP']['comModules'].split('|')

    # set path to peak field results
    suffix = ['pfd', 'ppr', 'pfv']
    for suf in suffix:
        refFile = inputDirRef / (simNameRef + '_' + suf + '.asc')
        fU.checkIfFileExists(refFile, fileType='reference simulation')
        pathDict[suf].append(refFile)
        log.info('Added to pathDict[%s] %s ' % (suf, refFile))
        compFile = inputDirComp / (simNameComp + '_' + suf + '.asc')
        fU.checkIfFileExists(compFile, fileType='comparison simulation')
        pathDict[suf].append(compFile)
        log.info('Added to pathDict[%s] %s ' % (suf, compFile))

    # if desired set path to mass log files
    sims = {comModules[0]: simNameRef, comModules[1]: simNameComp}
    if cfg['FLAGS'].getboolean('flagMass'):
        for comMod, sim in sims.items():
            log.info('mass file for comMod: %s and sim: %s' % (comMod, sim))
            if comMod == 'benchmarkReference':
                pathDict = getRefMB(cfg['AIMECSETUP']['testName'], pathDict, sim)
            elif comMod == 'com1DFAOrig':
                pathDict = extractCom1DFAMBInfo(avaDir, pathDict, simNameInput=sim)
            else:
                pathDict = getMBInfo(avaDir, pathDict, comMod, simName=sim)

    return pathDict


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

    # initialise empty pathDict for all required files
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'massBal': []}
    # get directories where simulation results can be found for both modules
    inputDirRef, inputDirComp, pathDict, refModule = getCompDirs(avaDir, cfgSetup, pathDict)

    # fill pathDict
    pathDict = getPathsFromSimName(pathDict, avaDir, cfg, inputDirRef, simNameRef, inputDirComp, simNameComp)

    return pathDict


def getCompDirs(avaDir, cfgSetup, pathDict):
    """ Determine dictionaries where simulation results can be found

        Parameters
        ----------
        avaDir: str
            path to avalanche directory
        cfgSetup: configParser object
            configuration settings for aimec
        pathDict: dict
            dictionary for paths for aimec

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

    pathDict['compType'] = ['comModules', refModule, compModule]
    pathDict['referenceFile'] = 0
    # info about colourmap
    pathDict['contCmap'] = True

    return inputDirRef, inputDirComp, pathDict, refModule


def mainDfa2Aimec(avaDir, comModule, cfg):
    """ Exports the required data from com1DFA to be used by Aimec

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
        pathDict: dict
            dictionary with path to simulation results for result types - keys: ppr, pfd, pfv
            and if configuration for ordering is provided, key: colorParameter with values for color coding results
    """

    # path dictionary for Aimec
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'massBal': [], 'colorParameter': []}

    # Setup input from com1DFA and save file paths to dictionary
    suffix = ['pfd', 'ppr', 'pfv']
    cfgSetup = cfg['AIMECSETUP']
    pathDict = fU.getDFADataPaths(avaDir, pathDict, cfgSetup, suffix, comModule)

    if cfg['FLAGS'].getboolean('flagMass'):
        # Extract mb info
        if comModule == 'com1DFAOrig':
            pathDict = extractCom1DFAMBInfo(avaDir, pathDict)
        else:
            pathDict = getMBInfo(avaDir, pathDict, comModule)

    pathDict['compType'] = ['singleModule', comModule]

    # info about colormap
    pathDict['contCmap'] = pU.contCmap

    return pathDict


def indiDfa2Aimec(avaDir, suffix, cfg, inputDir=''):
    """ Exports the required data from com1DFA to be used by Aimec for only one results parameter
        and with an option to specify the input directory

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        cfg: configParser object
            configuration for aimec
        inputDir: str
            optional - path to an input directory where simulation results can be found

        Returns
        --------
        pathDict: dict
            dictionary with path to simulation results for result types - keys: suffix
            and if configuration for ordering is provided, key: colorParameter with values for color coding results

     """

    # path dictionary for Aimec
    pathDict = {suffix: [], 'colorParameter': []}

    # Setup input from com1DFA and save file paths to dictionary
    comModule = cfg['anaMod']
    pathDict = fU.getDFADataPaths(avaDir, pathDict, cfg, suffix, comModule, inputDir=inputDir)

    pathDict['compType'] = ['singleModule', comModule]

    # info about colormap
    pathDict['contCmap'] = pU.contCmap

    return pathDict
