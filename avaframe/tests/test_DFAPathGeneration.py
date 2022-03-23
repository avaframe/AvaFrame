"""Tests for module DFAPathGeneration"""
import numpy as np
import math
import configparser
import pytest

# Local imports
import avaframe.ana5Utils.DFAPathGeneration as DFAPathGeneration
import avaframe.in3Utils.geoTrans as gT

def test_appendAverageStd():
    values = np.array([1, 2, 3, 4, 5])
    weights = np.array([2, 1, 2, 1, 2])
    average, std = DFAPathGeneration.weightedAvgAndStd(values, weights)
    assert average == 3
    assert std == 1.5
    proList = ['x', 'y', 'u2', 'ekin']
    particles = {'x': values, 'y': values, 'u2': 2*values, 'ekin': 10*2*values}
    avaProfile = {'x': np.empty((0, 1)), 'y': np.empty((0, 1)), 'xstd': np.empty((0, 1)), 'ystd': np.empty((0, 1)),
                  'u2': np.empty((0, 1)), 'u2std': np.empty((0, 1)), 'ekin': np.empty((0, 1)),
                  'ekinstd': np.empty((0, 1)), 'totEKin': np.empty((0, 1))}
    avaProfile = DFAPathGeneration.appendAverageStd(proList, avaProfile, particles, weights)
    print(avaProfile)
    assert avaProfile['x'] == 3
    assert avaProfile['xstd'] == 1.5
    assert avaProfile['y'] == 3
    assert avaProfile['ystd'] == 1.5
    assert avaProfile['u2'] == 6
    assert avaProfile['u2std'] == 3
    assert avaProfile['ekin'] == 60
    assert avaProfile['ekinstd'] == 30


def test_getDFAPathFromPart():
    values = np.array([1, 2, 3, 4, 5])
    weights = np.array([2, 1, 2, 1, 2])
    average, std = DFAPathGeneration.weightedAvgAndStd(values, weights)
    assert average == 3
    assert std == 1.5
    particlesList = [{'nPart': 5, 'm': weights, 'x': values, 'y': values, 'z': values, 's': values, 'sCor': values,
                      'ux': values, 'uy': values, 'uz': values}]
    avaProfile = DFAPathGeneration.getDFAPathFromPart(particlesList, addVelocityInfo=False)
    print(avaProfile)
    for prop in ['x', 'y', 'z', 's', 'sCor']:
        assert avaProfile[prop] == 3
        assert avaProfile[prop + 'std'] == 1.5

    avaProfile = DFAPathGeneration.getDFAPathFromPart(particlesList, addVelocityInfo=True)
    print(avaProfile)
    for prop in ['x', 'y', 'z', 's', 'sCor']:
        assert avaProfile[prop] == 3
        assert avaProfile[prop + 'std'] == 1.5
    assert avaProfile['u2'] == 33.75
    assert avaProfile['u2std'] == pytest.approx(27.52612396, abs=1e-6)
    assert avaProfile['ekin'] == pytest.approx(30, abs=1e-6)
    assert avaProfile['ekinstd'] == pytest.approx(27.69927797, abs=1e-6)


def test_extendDFAPathKernel():
    """"""
    # setup required inputs
    cfg = configparser.ConfigParser()
    cfg['PATH'] = {'nCellsResample': '5', 'extTopOption': '1', 'nCellsMinExtend': '1',
                      'nCellsMaxExtend': '2', 'factBottomExt': 0.2}
    avaProfile = {'x': np.array([10, 20, 30]), 'y': np.array([10, 20, 30]), 'z': np.array([40, 30, 20]),
                  's': np.array([0, math.sqrt(200), math.sqrt(800)])}
    particlesIni = {'x': np.array([7, 6.9]),
                    'y': np.array([10, 20])}
    dem = {'header': {'xllcenter': 0, 'yllcenter': 0, 'cellsize': 10, 'nrows': 10, 'ncols': 11},
           'rasterData': np.array([[50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0]])}

    particlesIni, _ =  gT.projectOnRaster(dem, particlesIni, interp='bilinear')

    # using the longest runout method
    avaProfileExt = DFAPathGeneration.extendDFAPathKernel(cfg['PATH'], avaProfile, dem, particlesIni)
    print(avaProfileExt)
    atol = 1e-10
    assert avaProfileExt['x'][0] == 7.0
    assert avaProfileExt['x'][-1] == pytest.approx(34.42426406871193, abs=1e-6)
    assert avaProfileExt['y'][0] == 10
    assert avaProfileExt['y'][-1] == pytest.approx(34.42426406871193, abs=1e-6)
    assert avaProfileExt['z'][0] == pytest.approx(43, abs=1e-6)
    assert avaProfileExt['z'][-1] == pytest.approx(15.575735931288072, abs=1e-6)

    # now use the highest point method
    cfg = configparser.ConfigParser()
    cfg['PATH'] = {'nCellsResample': '5', 'extTopOption': '0', 'nCellsMinExtend': '2',
                      'nCellsMaxExtend': '30', 'factBottomExt': 0.2}
    avaProfile = {'x': np.array([10, 20, 30]), 'y': np.array([10, 20, 30]), 'z': np.array([40, 30, 20])}

    avaProfileExt = DFAPathGeneration.extendDFAPathKernel(cfg['PATH'], avaProfile, dem, particlesIni)
    print(avaProfileExt)
    atol = 1e-10
    assert np.allclose(avaProfileExt['x'][0], 6.9, atol=atol)
    assert avaProfileExt['x'][-1] == pytest.approx(36.41019804034152, abs=1e-6)
    assert np.allclose(avaProfileExt['y'][0], 20.0, atol=atol)
    assert avaProfileExt['y'][-1] == pytest.approx(34.35700456911144, abs=1e-6)
    assert avaProfileExt['z'][0] == pytest.approx(43.09999999999, abs=1e-6)
    assert avaProfileExt['z'][-1] == pytest.approx(13.589801959658478, abs=1e-6)
