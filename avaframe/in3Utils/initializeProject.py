'''Functions to initialize project, i.e create folder structurec'''

import os
import logging

log = logging.getLogger(__name__)


def initializeFolderStruct(pathAvaName):
    ''' Initialize the standard folder structure

    str: pathAvaName : string with base avalanche path

    '''

    avaName = os.path.basename(os.path.normpath(pathAvaName))

    log.info("Running initialization sequence for avalanche %s", avaName)
    if os.path.exists(pathAvaName):
        log.warning("Folder %s already exists. Continuing", avaName)
        createFolderStruct(pathAvaName)

    else:
        try:
            os.makedirs(pathAvaName)
        except IOError as e:
            print('I/O error({0}): {1}'.format(e.errno, e.strerror))
            return

        log.info("Running initialization sequence for avalanche %s", avaName)

        createFolderStruct(pathAvaName)

    log.info("Done initializing avalanche %s", avaName)


def checkMakeDir(base, dirName):
    '''combines base and dir, checks whether it exists, if not creates it
    '''

    pathName = os.path.join(base, dirName)
    if not os.path.exists(pathName):
        os.makedirs(pathName)
    return pathName


def createFolderStruct(pathAvaName):
    '''creates the folder structure with avalanche base path'''

    Inputs = checkMakeDir(pathAvaName, 'Inputs')

    inputsSubDirs = ['RES', 'REL', 'ENT',
                     'POINTS', 'LINES']

    for cuDir in inputsSubDirs:
        checkMakeDir(Inputs, cuDir)

    checkMakeDir(pathAvaName, 'Outputs')
