# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 11:49:44 2025

@author: Domi
"""

import argparse
import pathlib
from com6RockAvalanche.variableVoellmyShapeToRaster import generateMuXsiRasters
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import initializeProject as initProj

def runMuXsiWorkflow(configPath=''):
    """
    Run the workflow to generate \u03bc and \u03be rasters.

    Parameters
    ----------
    configPath : str
        Path to the configuration file.

    Returns
    -------
    None
    """
    logName = 'runMuXsi'

    # Load general configuration file
    cfgMain = cfgUtils.getGeneralConfig()

    # Load configuration file path from general config if not provided
    if configPath:
        cfgMain['MAIN']['configFile'] = configPath
    else:
        configPath = cfgMain['MAIN']['configFile']

    configPath = pathlib.Path(configPath)

    # Start logging
    log = logUtils.initiateLogger(configPath.parent, logName)
    log.info('MAIN SCRIPT')
    log.info('Using configuration file: %s', configPath)

    # Clean input directory(ies) of old work files
    initProj.cleanSingleAvaDir(configPath.parent, deleteOutput=False)

    # Run the raster generation process
    generateMuXsiRasters(str(configPath))

    log.info('Workflow completed successfully.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run \u03bc and \u03be raster generation workflow')
    parser.add_argument('configPath', metavar='c', type=str, nargs='?', default='',
                        help='Path to the configuration file')

    args = parser.parse_args()
    runMuXsiWorkflow(str(args.configPath))
