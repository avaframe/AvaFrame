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

    # call function to be tested
    oAna.plotRangeTime(mtiInfo, cfgRangeTime['GENERAL'])


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
