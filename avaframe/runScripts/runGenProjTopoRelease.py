"""
    Run script for generateTopo in module in3Utils
    This file is part of Avaframe.
"""

import os
import logging
import shutil

# Local imports
from avaframe.in3Utils import initializeProject
from avaframe.in3Utils import getReleaseArea as gR
from avaframe.in3Utils import generateTopo as gT
from avaframe.in3Utils import cfgUtils, logUtils
from avaframe.out3Plot import outTopo as oT


# log file name; leave empty to use default runLog.log
logName = 'genProjTopoRelease'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger('.', logName)
log.info('MAIN SCRIPT')
log.info('Initializing Project: %s', avalancheDir)


# Initialize project
initializeProject.initializeFolderStruct(avalancheDir)

logOrigin = log.handlers[-1].baseFilename 

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
oT.plotGeneratedDEM(z, name_ext, cfgT, outDir, cfgMain)

# Initialise demType
demType = cfgT['TOPO']['demType']

if demType == 'HX' or demType == 'BL':
    log.warning('There is no release area available for this DEM type')
else:
    # Make release area
    [xv, yv, xyPoints] = gR.getReleaseArea(cfgT, cfgR, avalancheDir)

# Plot Release Area
if cfgMain['FLAGS'].getboolean('showPlot'):
    oT.plotReleasePoints(xv, yv, xyPoints, demType)
