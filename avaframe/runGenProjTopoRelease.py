"""
    Run script for generateTopo in module in3Utils
    This file is part of Avaframe.
"""

import sys
import os
import logging
import shutil

# Local imports
from avaframe.in3Utils import initializeProject
from avaframe.in3Utils import getReleaseArea as gR
from avaframe.in3Utils import generateTopo as gT
from avaframe.in3Utils import cfgUtils, logUtils
from avaframe.out3SimpPlot import outGenerateTopo as oT

#----------Rewquired Settings---------------

#TODO move this to a main cfg file!

# Avalanche directory; see doc.avaframe.org for setup
# TODO: full path needed?
avalancheDir = 'data/avaInlcinedPlane'

# log file name; leave empty to use default runLog.log
logName = 'genProjTopoRelease'
#------------------------------------------


# Start logging
log = logUtils.initiateLogger('.', logName)
log.info('MAIN SCRIPT')
log.info('Initializing Project: %s', avalancheDir)


# Initialize project
initializeProject.initializeFolderStruct(avalancheDir)

logOrigin = os.path.join('.', logName + '.log')

logDest = os.path.join(avalancheDir, logName + '.log')

shutil.move(logOrigin, logDest)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')

# Load input parameters from configuration file
cfgT = cfgUtils.getModuleConfig(gT)
cfgR = cfgUtils.getModuleConfig(gR)

# Call main function to generate DEMs
[z, name_ext, outDir] = gT.generateTopo(cfgT, avalancheDir)

# Plot new topogrpahy
oT.plotDEM(z, name_ext, cfgT, outDir)

# Initialise DEM_type
DEM_type = cfgT['TOPO']['DEM_type']

# Initialise output dir
outDir = gT.makeOutputDir(avalancheDir)

if DEM_type == 'HX' or DEM_type == 'BL':
    log.warning('There is no release area available for this DEM type')
else:
    # Make release area
    [xv, yv, xyPoints] = gR.getReleaseArea(cfgT, cfgR, outDir)

    # Plot Release Area
    if cfgR.getboolean('GENERAL', 'showplot'):
        oT.plotReleasePoints(xv, yv, xyPoints, DEM_type)
