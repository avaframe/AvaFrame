"""
    Fetch point values of simulation result fields using bilinear interpolation
"""

import numpy as np

# Local imports
from avaframe.in3Utils import cfgUtils, logUtils
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.transfromFields as tF
import avaframe.out3Plot.outTransformPlots as oT

# log file name; leave empty to use default runLog.log
logName = "FetchPointValues"

resType = 'pfv'
xyPoint = {'x': np.asarray([3710.]), 'y': np.asarray([-4000.])}

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain["MAIN"]["avalancheDir"]

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info("MAIN SCRIPT")


# load dataFrame for all configurations
simDF = cfgUtils.createConfigurationInfo(avalancheDir)

# create data frame that lists all available simulations and path to their result type result files
inputsDF, resTypeList = fU.makeSimFromResDF(avalancheDir, 'com1DFA')
# merge  parameters as columns to dataDF for matching simNames
dataDF = inputsDF.merge(simDF, left_on="simName", right_on="simName")

# derive values of resType fields at xyPint location
dataDF = tF.fetchPointValuesFromField(dataDF, xyPoint, resType, interpMethod='bilinear')


# create a plot
outName = 'pointValues_%s' % resType
oT.plotFetchedValues(avalancheDir, dataDF, ('pointValues_%s' % resType), 'relTh0', outName, scenario='simTypeActual')

