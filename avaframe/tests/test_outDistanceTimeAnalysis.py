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
    mti = np.repeat([np.arange(50)], 45, axis=0)
    demData = {'rasterData': np.zeros((10, 15)), 'header': {'cellsize': 5., 'xllcenter': 0.0,
        'yllcenter': 0.0, 'noDataValue': 9999, 'nrows': 10, 'ncols': 15}}
    data = np.ones((10, 15))
    mtiInfo = {'mti': mti, 'rangeList': np.arange(50), 'timeList': np.arange(50),
        'rangeGates': np.arange(45) -25., 'plotTitle': 'testTitle', 'referencePointName': 'betaPoint',
        'type': 'tt', 'betaPointAngle': 9.987, 'cLower': 3, 'slRaster': np.ones((45,45)),
        'rangeRaster': np.ones((50, 45)), 'sType': 'projected'}
    mtiInfo['rasterTransfo'] = {'x': np.arange(45), 'y': np.ones(45), 'cellSizeSL': 2., 'DBXl': 0.0,
        'DBXr': 2., 'DBYl': 0, 'DBYr': 2., 'indStartOfRunout': 4, 'gridx': np.ones((10, 15)),
        'gridy': np.ones((10, 15)), 'l': np.arange(45), 's': np.arange(45)}

    cfgRangeTime = configparser.ConfigParser()
    cfgRangeTime['GENERAL'] = {'avalancheDir': avaDir, 'rangeTimeResType': 'FT', 'simHash': 'simDI',
        'minVelTimeStep': 2., 'thresholdResult': 1.e-2, 'minVelTimeStep': 2, 'maxOrMean': 'max'}
    cfgRangeTime['PLOTS'] = {'width': 0.25, 'height': 0.25, 'lw': 0.25, 'textsize': 7}
    cfgRangeTime['ANIMATE'] = {'resMin': 0.0, 'resMax': 10., 'xyEastNorth': True, 'panelTitles': True,
        'xMin': 0.0, 'xMax': 45, 'yMin': 0.0, 'yMax': 45.}
    cellSize = 5.
    resType = 'FT'
    timeStep = 1.

    # call function to be tested
    oAna.animationPlot(demData, data, cellSize, resType, cfgRangeTime, mtiInfo, timeStep)

    cfgRangeTime['ANIMATE']['panelTitles'] = 'False'
    cfgRangeTime['ANIMATE']['xyEastNorth'] = 'False'

    # call function to be tested
    oAna.animationPlot(demData, data, cellSize, resType, cfgRangeTime, mtiInfo, timeStep)
