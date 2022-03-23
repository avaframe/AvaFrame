''' Tests for module distanceTimeAnalysis  '''
import numpy as np
import numpy.ma as ma
import pandas as pd
import pathlib
import configparser
import matplotlib.pyplot as plt

# Local imports
import avaframe.ana5Utils.distanceTimeAnalysis as dtAna
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


def test_getRadarLocation():
    """ test get radar location """

    # setup required inputs
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'radarLocation': '10.|25.|40.|57.'}

    # call function to be tested
    radarFov = dtAna.getRadarLocation(cfg)

    assert radarFov[0][0] == 10.
    assert radarFov[1][1] == 57.
    assert np.array_equal(radarFov, np.asarray([[10., 25.],[40., 57.]]))


def test_setDEMOrigin():
    """ test setting of DEM origin back to original """

    # setup required inputs
    headerDEM = {'xllcenter': 0.0, 'yllcenter': 0.0, 'cellsize': 2., 'ncols': 8, 'nrows': 11}
    demSims = {'header': headerDEM, 'originalHeader': {'xllcenter': 1.0, 'yllcenter': 5.,
        'cellsize': 2., 'ncols': 8, 'nrows': 11},
        'rasterData': np.zeros((11, 8))}

    # call function to be tested
    demOriginal = dtAna.setDemOrigin(demSims)

    assert demOriginal['header']['xllcenter'] == 1.0
    assert demOriginal['header']['yllcenter'] == 5.0
    assert demOriginal['header']['ncols'] == 8
    assert demOriginal['header']['nrows'] == 11
    assert demOriginal['header']['cellsize'] == 2.


def test_radarMask(tmp_path):
    """ test creating radar mask array """

    # setup required input
    testAvaDir = pathlib.Path(tmp_path, 'testAva')
    headerDEM = {'xllcenter': 0.0, 'yllcenter': 0.0, 'cellsize': 1., 'ncols': 11, 'nrows': 11}
    demOriginal = {'header': headerDEM, 'rasterData': np.zeros((11, 11))}
    radarFov = [[1., 10.0], [5., 5.]]
    aperture = 40.5
    cfgRangeTime = configparser.ConfigParser()
    cfgRangeTime['GENERAL'] = {'rgWidth': 2., 'avalancheDir': testAvaDir, 'simHash': 'test123',
        'gateContours': 20}

    # call function to be tested
    radarRange, rangeGates = dtAna.radarMask(demOriginal, radarFov, aperture, cfgRangeTime)

    print('randarRange', radarRange)
    print('rangeGates', rangeGates)


    assert np.array_equal(rangeGates, np.asarray([2., 4., 6., 8., 10., 12.]))


def test_avalancheMask():
    """ test creating avalanche mask """

    # setup required input parameters
    resData = np.zeros((10, 12))
    resData[4,4:6] = 4.5
    resData[5,4:6] = 4.2
    threshold = 4.19

    # call function to be tested
    maskResType = np.ma.masked_where(resData < threshold, resData)

    assert np.all(maskResType.mask[4,4:6]) == False
    assert np.all(maskResType.mask[5,4:6]) == False
    assert np.all(maskResType.mask[0,:])


def test_minRangeSimulation():
    """ test if min range is found in simulation results """

    # setup required inputs
    flowF = np.zeros((10, 12))
    flowF[4,4:6] = 4.5
    flowF[5,4:6] = 4.2
    threshold = 4.19
    rangeMasked = np.arange(12)
    rangeMasked = np.repeat([rangeMasked], 10, axis=0)

    # call function to be tested
    losRange = dtAna.minRangeSimulation(flowF, rangeMasked, threshold)

    assert losRange == 4.
