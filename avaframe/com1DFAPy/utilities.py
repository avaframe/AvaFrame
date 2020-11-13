"""
    utility functions for com1DFA

    This file is part of Avaframe.
"""

# Load modules
import os
import glob
import logging
import numpy as np

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getInputData(avaDir, cfg):
    """ Fetch input datasets required for simulation

    Input:
            Directory

    Output:
            release area file
            entrainment area file
            resistance area file
    """

    # Set directories for inputs, outputs and current work
    inputDir = os.path.join(avaDir, 'Inputs')
    outputDir = os.path.join(avaDir, 'Outputs', 'com1DFAPy')
    fU.makeADir(outputDir)
    workDir = os.path.join(avaDir, 'Work', 'com1DFAPy')
    fU.makeADir(workDir)

    # Set flag if there is an entrainment or resistance area
    flagEntRes = False

    # Initialise release areas, default is to look for shapefiles
    relFiles = glob.glob(inputDir+os.sep + 'devREL'+os.sep + '*.shp')
    log.info('Release area files are: %s' % relFiles)

    # Initialise resistance areas
    if cfg['FLAGS'].getboolean('flagRes'):
        resFiles = glob.glob(inputDir+os.sep + 'RES' + os.sep+'*.shp')
        if len(resFiles) < 1:
            log.warning('No resistance file')
            resFiles.append('')  # Kept this for future enhancements
        else:
            try:
                message = 'There shouldn\'t be more than one resistance .shp file in ' + inputDir + '/RES/'
                assert len(resFiles) < 2, message
            except AssertionError:
                raise
            flagEntRes = True
    else:
        resFiles = []
        resFiles.append('')

    # Initialise entrainment areas
    if cfg['FLAGS'].getboolean('flagEnt'):
        entFiles = glob.glob(inputDir+os.sep + 'ENT' + os.sep+'*.shp')
        if len(entFiles) < 1:
            log.warning('No entrainment file')
            entFiles.append('')  # Kept this for future enhancements
        else:
            try:
                message = 'There shouldn\'t be more than one entrainment .shp file in ' + inputDir + '/ENT/'
                assert len(entFiles) < 2, message
            except AssertionError:
                raise
            flagEntRes = True
    else:
        entFiles = []
        entFiles.append('')

    # Initialise DEM
    demFile = glob.glob(inputDir+os.sep+'*.asc')
    try:
        assert len(demFile) == 1, 'There should be exactly one topography .asc file in ' + inputDir
    except AssertionError:
        raise

    # return DEM, first item of release, entrainment and resistance areas
    return demFile[0], relFiles, entFiles[0], resFiles[0], flagEntRes
