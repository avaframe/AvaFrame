'''
    Functions to initialize project, i.e create folder structure

    This file is part of Avaframe.
'''

import os
import logging
import shutil

log = logging.getLogger(__name__)


def _checkForFolderAndDelete(baseDir, folderName):
    '''Helper function ONLY USE WITH CHECKED INPUT'''
    fPath = os.path.join(baseDir, folderName)
    try:
        shutil.rmtree(fPath)
    except FileNotFoundError:
        log.debug("No %s folder found.", folderName)


def _checkAvaDirVariable(avaDir):
    '''Check for empty or nonString avaDir variable'''
    # check for empty avaDir name, abort if empty
    if not avaDir:
        log.warning("AvaDir variable is undefined, returning. ")
        return 'AvaDir is empty'

    # make sure avaDir is a string
    isString = isinstance(avaDir, str)
    if not isString:
        log.warning("AvaDir is not a string, returning")
        return 'AvaDir is NOT a string'

    return 'SUCCESS'


def cleanModuleFiles(avaDir, module, alternativeName=''):
    '''Cleans all generated files from the provided module in the outputs
    folder and work folder

    Parameters
    ----------
    avaDir : path/string
        Avalanche directory path
    module : module object
        The module do delete.
        e.g. from avaframe.com2AB import com2AB
        leads to cleanModuleFiles(com2AB)
        whereas
        from avaframe.com2AB import com2AB as c2
        leads to cleanModuleFiles(c2)
    '''

    # check for empty or non string variable
    result = _checkAvaDirVariable(avaDir)
    if not 'SUCCESS' in result:
        return result

    # get filename of module
    # TODO: remove this option when com1DFAPy replaces com1DFA
    if alternativeName:
        modName = alternativeName
    else:
        modName = os.path.basename(module.__file__)
        modName = os.path.splitext(modName)[0]

    # output dir
    outDir =  os.path.join(avaDir, 'Outputs')
    workDir =  os.path.join(avaDir, 'Work')

    log.info("Cleaning module %s in folder: %s ", modName, outDir)
    _checkForFolderAndDelete(outDir, modName)
    _checkForFolderAndDelete(workDir, modName)

    return 'SUCCESS'


def cleanSingleAvaDir(avaDir, keep=None, deleteOutput=True):
    '''
    Clean a single avalanche directory of the work and output directories
    Expects a avalanche directory name as string
    and optional:
    a log name to keep (and not delete)
    Boolean to be able to avoid deletion of Outputs (true by default)
    '''

    # check for empty or non string variable
    result = _checkAvaDirVariable(avaDir)
    if not 'SUCCESS' in result:
        return result

    # Info to user
    log.debug("Cleaning folder: %s ", avaDir)

    # Try to remove OUTPUTS folder, only pass FileNotFoundError, i.e. folder
    # does not exist
    if deleteOutput:
        _checkForFolderAndDelete(avaDir, 'Outputs')

    # Try to remove Work folder, only pass FileNotFoundError, i.e. folder
    # does not exist
    _checkForFolderAndDelete(avaDir, 'Work')

    # check for *.log files, go to remove only if exists
    allFiles = os.listdir(avaDir)
    logFiles = [fi for fi in allFiles if fi.endswith(".log")]

    # keep the log file specified in keep
    if keep:
        logFiles = [item for item in logFiles if keep not in item]

    for fi in logFiles:
        filePath = os.path.join(avaDir, fi)
        os.remove(filePath)

    log.debug("Directory is squeaky clean")
    return 'Cleaned directory'


def initializeFolderStruct(pathAvaName, removeExisting=False):
    ''' Initialize the standard folder structure. If removeExisting is true,
    deletes any existing folders! BEWARE!

    Parameters
    ----------
    pathAvaName: str
        string with base avalanche path

    removeExisting: boolean, optional, default: False
        remove existing directory, use this to completely clean out the
        directory
    '''

    avaName = os.path.basename(os.path.normpath(pathAvaName))

    log.info("Running initialization sequence for avalanche %s", avaName)

    if os.path.exists(pathAvaName):
        if removeExisting:
            log.warning('Will delete %s !', pathAvaName)
            # remove and remake
            shutil.rmtree(pathAvaName)
            os.makedirs(pathAvaName)
        else:
            log.warning("Folder %s already exists. Continuing", avaName)

    else:
        try:
            os.makedirs(pathAvaName)
        except IOError as e:
            print('I/O error({0}): {1}'.format(e.errno, e.strerror))
            return

    createFolderStruct(pathAvaName)

    log.info("Done initializing avalanche %s", avaName)


def createFolderStruct(pathAvaName):
    '''creates the folder structure with avalanche base path'''

    Inputs = checkMakeDir(pathAvaName, 'Inputs')

    inputsSubDirs = ['RES', 'REL', 'SECREL', 'ENT',
                     'POINTS', 'LINES']

    for cuDir in inputsSubDirs:
        checkMakeDir(Inputs, cuDir)

    checkMakeDir(pathAvaName, 'Outputs')

    checkMakeDir(pathAvaName, 'Work')


def checkMakeDir(base, dirName):
    '''combines base and dir, checks whether it exists, if not creates it
    '''

    pathName = os.path.join(base, dirName)
    if not os.path.exists(pathName):
        os.makedirs(pathName)
        log.debug('Creating folder: %s: ', pathName)
    else:
        log.info('Folder exists: %s, Skipping', pathName)
    return pathName
