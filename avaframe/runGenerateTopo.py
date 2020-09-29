"""
    Run script for generateTopo in module in3Utils
    This file is part of Avaframe.
"""

import sys
import logging

# Local imports
from avaframe.in3Utils import generateTopo as gT
from avaframe.in3Utils import cfgUtils, logUtils
from avaframe.out3Plot import outGenerateTopo as oT


# log file name; leave empty to use default runLog.log
logName = 'generateTopo'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')

# Load input parameters from configuration file
cfg = cfgUtils.getModuleConfig(gT)

# Call main function to generate DEMs
[z, name_ext, outDir] = gT.generateTopo(cfg, avalancheDir)

# Plot new topogrpahy
oT.plotDEM(z, name_ext, cfg, outDir)
