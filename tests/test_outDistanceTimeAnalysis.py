''' Tests for module distanceTimeAnalysis  '''
import numpy as np
import numpy.ma as ma
import pandas as pd
import pathlib
import configparser
import matplotlib.pyplot as plt

# Local imports
import avaframe.ana5Utils.distanceTimeAnalysis as dtAna
import avaframe.out3Plot.outDistanceTimeAnalysis as oAna
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


def test_plotRangeTime(tmp_path):
    """ test plot range time """

    # setup required inputs
    avaDir = pathlib.Path(tmp_path)
    mti = np.repeat([np.arange(50)], 45, axis=0)
    mtiInfo = {'mti': mti, 'rangeList': np.arange(50), 'timeList': np.arange(50),
        'rangeGates': np.arange(45) -25., 'plotTitle': 'testTitle', 'referencePointName': 'betaPoint',
        'type': 'tt', 'betaPointAngle': 9.987}
    cfgRangeTime = configparser.ConfigParser()
    cfgRangeTime['GENERAL'] = {'avalancheDir': avaDir, 'rangeTimeResType': 'FT', 'simHash': 'simDI',
        'minVelTimeStep': 2.}
        # 'minVelTimeStep': 2., 'maxOrMean': 'max'}
    cfgRangeTime['PLOTS'] = {'width': 0.25, 'height': 0.25, 'lw': 0.25, 'textsize': 7}

    # call function to be tested
    oAna.plotRangeTime(mtiInfo, cfgRangeTime)


def test_radarFieldOfViewPlot(tmp_path):
    """ test plotting radar field of view and dem """

    # setup required inputs
    avaDir = pathlib.Path(tmp_path)
    radarFov = np.asarray([[0.,10.], [5., 5.]])
    radarRange = np.repeat([np.arange(50)], 45, axis=0)
    cfgRangeTime = configparser.ConfigParser()
    cfgRangeTime['GENERAL'] = {'avalancheDir': avaDir, 'aperture': 5., 'simHash': 'simDI'}
    cfgRangeTime['PLOTS'] = {'gateContours': 5}


    rangeGates = np.arange(45)
    demData = np.zeros((45, 50))
    dem = {'header': {'nrows': 45, 'ncols': 50, 'cellsize': 5., 'xllcenter': 0.0, 'yllcenter': 0.0}, 'rasterData': demData}

    oAna.radarFieldOfViewPlot(radarFov, radarRange, cfgRangeTime, rangeGates, dem)


def test_animationPlot(tmp_path):
    """ test plot range time """

    # setup required inputs
    avaDir = pathlib.Path(tmp_path)
    mti = np.ones((543, 12))
    demData = {'rasterData': np.ones((401, 1001)), 'header': {'cellsize': 5., 'xllcenter': 0.0,
        'yllcenter': 0.0, 'noDataValue': 9999, 'nrows': 401, 'ncols': 1001}}
    data = np.ones((401, 1001))
    mtiInfo = {'mti': mti, 'rangeList': np.arange(12), 'timeList': np.arange(12),
        'rangeGates': np.arange(543) -25., 'plotTitle': 'testTitle', 'referencePointName': 'betaPoint',
        'type': 'tt', 'betaPointAngle': 9.987, 'cLower': 2, 'slRaster': np.ones((543,121)),
        'rangeRaster': np.ones((543, 121)), 'sType': 'projected'}
    mtiInfo['rasterTransfo'] = {'x': np.arange(543), 'y': np.ones(543), 'cellSizeSL': 2., 'DBXl': np.asarray([0.0, 10]),
        'DBXr': np.asarray([2., 20]), 'DBYl': np.asarray([0, 9]), 'DBYr': np.asarray([0, 2.]), 'indStartOfRunout': 4, 'gridx': np.ones((543, 121)),
        'gridy': np.ones((543, 121)), 'l': np.arange(121), 's': np.arange(543)}

    cfgRangeTime = configparser.ConfigParser()
    cfgRangeTime['GENERAL'] = {'avalancheDir': avaDir, 'rangeTimeResType': 'FT', 'simHash': 'simDI',
        'minVelTimeStep': 2., 'thresholdResult': 1.e-2, 'minVelTimeStep': 2, 'maxOrMean': 'max'}
    cfgRangeTime['PLOTS'] = {'width': 0.25, 'height': 0.25, 'lw': 0.25, 'textsize': 7}
    cfgRangeTime['ANIMATE'] = {'resMin': 0.0, 'resMax': 1.9, 'xyEastNorth': True, 'panelTitles': True,
        'xMin': 0.0, 'xMax': 160, 'yMin': 1840, 'yMax': 870}
    cellSize = 5.
    resType = 'FT'
    timeStep = 1.

    # call function to be tested
    oAna.animationPlot(demData, data, cellSize, resType, cfgRangeTime, mtiInfo, timeStep)

    cfgRangeTime['ANIMATE']['panelTitles'] = 'False'
    cfgRangeTime['ANIMATE']['xyEastNorth'] = 'False'

    # call function to be tested
    oAna.animationPlot(demData, data, cellSize, resType, cfgRangeTime, mtiInfo, timeStep)
