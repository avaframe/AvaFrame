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

    assert radarFov[0] == [10., 25.]
    assert radarFov[1] == [40., 57.]
    assert radarFov == [[10., 25.],[40., 57.]]


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


def test_rotate():
    """ test rotate a vector by an angle """

    # setup required data
    locationPoints = [[0., 0.], [0., 10.]]
    theta = 45.
    deg = True

    # call function to be tested
    rotatedLine = dtAna.rotate(locationPoints, theta, deg=deg)

    sL = np.sqrt(50.)
    assert rotatedLine[0][0] == 0.
    assert -sL - 1.e-12 < rotatedLine[0][1] < -sL+1.e-12
    assert rotatedLine[1][0] == 0.
    assert sL - 1.e-12 < rotatedLine[1][1] < sL+1.e-12
    assert isinstance(rotatedLine, list)
    assert len(rotatedLine) == 2
    assert len(rotatedLine[0]) == 2


def test_cartToSpherical():
    """ test converting to spherical coordinates """

    # setup required inputs
    X = 10.
    Y = 10.
    Z = np.sqrt(10.**2 + 10.**2)

    # call function to be tested
    r, phi, theta = dtAna.cartToSpherical(X, Y, Z)

    assert r == np.sqrt(10.**2 + 10.**2 + Z**2)
    assert phi == 45.
    assert theta == 45.


def test_radarMask(tmp_path):
    """ test creating radar mask array """

    # setup required input
    testAvaDir = pathlib.Path(tmp_path, 'testAva')
    headerDEM = {'xllcenter': 0.0, 'yllcenter': 0.0, 'cellsize': 1., 'ncols': 11, 'nrows': 11}
    demOriginal = {'header': headerDEM, 'rasterData': np.zeros((11, 11))}
    radarFov = [[1., 10.0], [5., 5.]]
    aperture = 45.
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
    maskResType = dtAna.avalancheMask(resData, threshold)

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


def test_fetchFlowFields():
    """ test fetching fields in a folder """

    # setup required input
    # get input data
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannelPytest'
    flowFieldsDir = dirPath / '..'/ '..' / 'benchmarks' / avaName / 'Outputs' / 'com1DFA' / 'peakFiles'
    suffix = 'ppr'

    # call function to be tested
    flowFields = dtAna.fetchFlowFields(flowFieldsDir, suffix=suffix)
    flowFields = sorted(flowFields)
    print('flowFields', flowFields, sorted(flowFields))

    assert flowFields[0].stem == 'release1HS_ent_dfa_d10bdc1e81_ppr'
    assert flowFields[1].stem == 'release2HS_ent_dfa_e2145362b7_ppr'
    assert len(flowFields) == 2

    # call function to be tested
    flowFields = dtAna.fetchFlowFields(flowFieldsDir)
    flowFields = sorted(flowFields)
    print('flowFields', flowFields, sorted(flowFields))

    assert flowFields[0].stem == 'release1HS_ent_dfa_d10bdc1e81_pfd'
    assert len(flowFields) == 6
