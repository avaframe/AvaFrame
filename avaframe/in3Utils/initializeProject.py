import sys
import os
import logging
log = logging.getLogger(__name__)


def initializeFolderStruct(pathAvaName):
    if os.path.exists(pathAvaName):
        log.info("Running initialization sequence for avalanche %s", str(pathAvaName))
        log.warning("Folder %s already exists. Checking folder structure and adding missing files and folders.", str(
            pathAvaName).split('/')[-1])
        createFolderStruct(pathAvaName)
    else:
        try:
            os.makedirs(pathAvaName)
        except IOError as e:
            print('I/O error({0}): {1}'.format(e.errno, e.strerror))
            return
        log.info("Running initialization sequence for avalanche %s", str(pathAvaName))
        createFolderStruct(pathAvaName)
        log.info("Done initializing avalanche %s", str(pathAvaName).split('/')[-1])


def createFolderStruct(pathAvaName):

    Inputs = os.path.join(pathAvaName, 'Inputs/')
    if not os.path.exists(Inputs):
        os.makedirs(Inputs)
    Res = os.path.join(Inputs, 'RES/')
    if not os.path.exists(Res):
        os.makedirs(Res)
    Rel = os.path.join(Inputs, 'REL/')
    if not os.path.exists(Rel):
        os.makedirs(Rel)
    Ent = os.path.join(Inputs, 'ENT/')
    if not os.path.exists(Ent):
        os.makedirs(Ent)
    Points = os.path.join(Inputs, 'POINTS/')
    if not os.path.exists(Points):
        os.makedirs(Points)
    Lines = os.path.join(Inputs, 'LINES/')
    if not os.path.exists(Lines):
        os.makedirs(Lines)

    Outputs = os.path.join(pathAvaName, 'Outputs/')
    if not os.path.exists(Outputs):
        os.makedirs(Outputs)
