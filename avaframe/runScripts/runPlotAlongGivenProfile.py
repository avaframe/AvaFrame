"""
    Run script for plotting a profile along the middle of the y-axis for peak files within an avalancheDir,
    all peakfields need to have same spatial extent 
"""
# Load modules
# importing general python modules
import pathlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.out3Plot import plotUtils as pU
import avaframe.in2Trans.rasterUtils as IOf
import avaframe.out3Plot.outProfilePlots as oT


################USER Input#############
# color indicator - options: ['files', 'names', 'resType', 'simType', 'isDefault', 'frictCalib', 'simName', 'modelType',
# 'releaseArea', 'cellSize', 'simID', 'timeStep']
cInd = "names"
resType = "pfv"
modName = 'com1DFA'
inputDir = ''
############################################################

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain["MAIN"]["avalancheDir"]
outDir = pathlib.Path(avalancheDir, 'Outputs', '../out3Plot')

# Start logging
logName = "plotProfiles_%s_%s" % (resType, cInd)
# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info("MAIN SCRIPT")
log.info("Current avalanche: %s", avalancheDir)

if inputDir == '':
    inputDir = pathlib.Path(avalancheDir, 'Outputs', modName, 'peakFiles')
elif isinstance(inputDir, pathlib.Path) is False:
    message = 'Provided path to inputDir: %s is not valid' % str(inputDir)
    log.error(message)
    raise FileNotFoundError(message)

# load dataFrame with info on peakFiles
peakFilesDF = fU.makeSimDF(inputDir)

# create plots
oT.plotAlongProfile(avalancheDir, peakFilesDF, resType='pft', cInd=cInd)
