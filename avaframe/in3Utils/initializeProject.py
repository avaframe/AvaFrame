'''
    Functions to initialize project, i.e create folder structure

    This file is part of Avaframe.
'''

import os
import logging
import shutil
import pathlib

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
        log.warning("avaDir variable is undefined, returning. ")
        return 'avaDir is empty'

    # make sure avaDir is a string
    isString = isinstance(avaDir, str)
    isPurePath = isinstance(avaDir, pathlib.PurePath)
    if not ((isString) or (isPurePath)):
        log.warning("avaDir is not a string or a PurePath, returning")
        return 'avaDir is NOT a string or PurePath'

    return 'SUCCESS'


def cleanModuleFiles(avaDir, module, alternativeName='', deleteOutput=True):
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
        Boolean to be able to avoid deletion of Outputs (true by default)
    '''

    # check for empty or non string variable
    result = _checkAvaDirVariable(avaDir)
    if 'SUCCESS' not in result:
        return result

    # get filename of module
    # TODO: remove this option when com1DFAPy replaces com1DFA
    if alternativeName:
        modName = alternativeName
    else:
        modName = os.path.basename(module.__file__)
        modName = os.path.splitext(modName)[0]

    # output dir
    outDir = os.path.join(avaDir, 'Outputs')
    workDir = os.path.join(avaDir, 'Work')

    if deleteOutput:
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
    if 'SUCCESS' not in result:
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


def cleanDEMremeshedDir(avaDir):
    """ clean DEMremeshed folder in avaDir Inputs

        Parameters
        ------------
        avaDir: str or pathlib patch
            path to avalanche directory
    """

    avaDirInputs = pathlib.Path(avaDir, 'Inputs')
    avaDirInputsString = str(avaDirInputs)

    # check for empty or non string variable
    result = _checkAvaDirVariable(avaDirInputsString)
    if 'SUCCESS' not in result:
        return result

    # clean directory
    folderName = 'DEMremeshed'
    _checkForFolderAndDelete(avaDirInputsString, folderName)
    log.info('Cleaned %s/%s directory' % (avaDirInputsString,folderName))

    return 'SUCCESS'


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
    
    checkMakeDir(pathAvaName, 'Avaflow_Input')


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
