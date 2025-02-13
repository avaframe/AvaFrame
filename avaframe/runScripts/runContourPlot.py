"""
Run the contour plots of com1DFA
"""

import pathlib

# Local imports
# import config and init tools
from avaframe.in3Utils import cfgUtils
from avaframe.in2Trans import shpConversion as shpConv
from avaframe.out3Plot import outContours as oC
# import computation modules

# +++++++++REQUIRED+++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runContourPlotsfromShp'
levels = [0.5, 1., 1.5, 2., 3., 4.]
plotTitle = 'Contour line comparisons'
resType = 'pft'
unit = 'm'
# ++++++++++++++++++++++++++++++

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# read external shp file data for comparison
inDirRef = pathlib.Path(avalancheDir, 'Work', 'shpFiles')
inFileRef = list(inDirRef.glob('*.shp'))[0]
shpData = shpConv.SHP2Array(inFileRef, defname=None)
contourDictRef = oC.readShpLines(shpData, layerName='layerName')

# read simulation result field and extract contour lines)
inFileSimList = list(inDirRef.glob('*%s.asc' % resType)) + list(inDirRef.glob('*%s.tif' % resType))
inFileSim = inFileSimList[0]
contourDictSim = oC.createRasterContourDict(inFileSim, levels)

# setup pathDict where to save plots
pathDict = {'avaDir': avalancheDir, 'pathResult': pathlib.Path(avalancheDir, 'Outputs', 'contourPlots'),
    'title': plotTitle, 'parameter': resType, 'unit': unit}
oC.plotContoursFromDict(contourDictRef, contourDictSim, pathDict, levels, multiplePlots=True)
