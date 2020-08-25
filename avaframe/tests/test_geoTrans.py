"""Tests for module geoTrans"""
import numpy as np
import pytest

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in3Utils.ascUtils as IOf


def test_projectOnRaster(capfd):
    '''projectOnRaster'''
    dem = {}
    Points = {}
    header = IOf.cASCheader()
    header.xllcorner = -5
    header.yllcorner = 2
    header.cellsize = 5

    rasterdata = np.array(([0, 1, 2], [2, 0.5, 5]))
    Points['x'] = np.array((-5, -0.000000001, -5, 2, 0, -2.5, -2, 2))
    Points['y'] = np.array((2, 6.99999999999, 6.999999999999, 6, 10, 4.5, np.nan, 8))

    dem['header'] = header
    dem['rasterData'] = rasterdata

    Points = geoTrans.projectOnRaster(dem, Points)
    tol = 1e-8
    assert (Points['z'][0] == pytest.approx(0, abs=tol))
    assert (Points['z'][1] == pytest.approx(0.5, abs=tol))
    assert (Points['z'][2] == pytest.approx(2, abs=tol))
    assert (Points['z'][3] == pytest.approx(2.12, abs=tol))
    assert np.isnan(Points['z'][4])
    assert (Points['z'][5] == pytest.approx(0.875, abs=tol))
    assert np.isnan(Points['z'][6])
    assert np.isnan(Points['z'][7])


def test_projectOnRasterVect(capfd):
    '''projectOnRasterVect'''
    dem = {}
    Points = {}
    header = IOf.cASCheader()
    header.xllcorner = 1
    header.yllcorner = 2
    header.cellsize = 1
    header.ncols = 4
    header.nrows = 3

    rasterdata = np.array(([0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]))
    Points['x'] = np.array(([1.0000000000001, 1.00000000001],
                            [-4, 2.67], [3.9999999999, 3], [9, 2]))
    Points['y'] = np.array(([2.0000000000001, 3], [4, 2.6], [3.4, 7], [3.6, 7]))

    dem['header'] = header
    dem['rasterData'] = rasterdata

    Points, itot, ioob = geoTrans.projectOnRasterVect(dem, Points, interp='bilinear')
    print(Points)
    atol = 1e-8
    # Compare result to reference solution
    zSol = np.array([[4.99600361e-13, 4.00000000e+00],
                     [np.nan, 4.07000000e+00],
                     [8.60000000e+00,            np.nan],
                     [np.nan,            np.nan]])

    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=atol)
    assert (testRes == True)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=atol)
    assert (testRes == True)

    Points, itot, ioob = geoTrans.projectOnRasterVect(dem, Points, interp='nearest')
    atol = 1e-8
    # Compare result to reference solution
    zSol = np.array([[0.,  4.], [np.nan,  6.], [7., np.nan], [np.nan, np.nan]])
    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=atol)
    assert (testRes == True)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=atol)
    assert (testRes == True)


def test_checkProfile(capfd):
    '''checkProfile'''

    AvaProfile = {}
    AvaProfile['x'] = np.array((0, 1, 2, 3, 4))
    AvaProfile['y'] = np.array((1, 2, 3, 4, 5))
    AvaProfile['z'] = np.array((0, 1, 2, 3, 4))

    projSplitPoint, AvaProfile = geoTrans.checkProfile(AvaProfile, projSplitPoint=None)
    assert AvaProfile['z'][0] == 4
    assert AvaProfile['z'][-1] == 0
    assert AvaProfile['x'][0] == 4
    assert AvaProfile['x'][-1] == 0
    assert AvaProfile['y'][0] == 5
    assert AvaProfile['y'][-1] == 1
    assert AvaProfile['y'][1] == 4

    projSplitPoint = {}
    projSplitPoint['indSplit'] = 1
    AvaProfile = {}
    AvaProfile['x'] = np.array((0, 1, 2, 3, 4))
    AvaProfile['y'] = np.array((1, 2, 3, 4, 5))
    AvaProfile['z'] = np.array((0, 1, 2, 3, 4))
    AvaProfile['s'] = np.array((0, 2, 4, 6, 8))
    projSplitPoint, AvaProfile = geoTrans.checkProfile(AvaProfile, projSplitPoint=projSplitPoint)
    assert AvaProfile['z'][0] == 4
    assert AvaProfile['z'][-1] == 0
    assert AvaProfile['x'][0] == 4
    assert AvaProfile['x'][-1] == 0
    assert AvaProfile['y'][0] == 5
    assert AvaProfile['y'][-1] == 1
    assert AvaProfile['y'][1] == 4
    assert AvaProfile['s'][0] == 0
    assert AvaProfile['indSplit'] == 3


def test_findAngleProfile(capfd):
    '''findAngleProfile'''
    s = np.linspace(0, 400, 41)
    angle = np.linspace(40, 0, 41)
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    deltaInd = 3
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 30

    deltaInd = 1
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 30

    angle[10] = 8
    angle[11] = 8
    angle[12] = 8
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    deltaInd = 3
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 30

    angle[13] = 8
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 9


def test_path2domain(capfd):
    '''test_path2domain'''
    header = IOf.cASCheader()
    header.xllcenter = 1
    header.yllcenter = 2
    header.cellsize = 1
    xyPath = {}
    xyPath['x'] = np.array((0, 10, 20, 30, 40))
    xyPath['y'] = np.array((10, 20, 30, 40, 50))
    w = 10
    DB = geoTrans.path2domain(xyPath, w, header)
    print(DB)

    atol = 1e-8
    # Compare result to reference solution
    zSol = np.array([-4.53553391,  5.46446609, 15.46446609, 25.46446609, 35.46446609])
    testRes = np.allclose(DB['DBXl'], zSol, atol=atol)
    assert (testRes == True)
    zSol = np.array([2.53553391, 12.53553391, 22.53553391, 32.53553391, 42.53553391])
    testRes = np.allclose(DB['DBXr'], zSol, atol=atol)
    assert (testRes == True)
    zSol = np.array([11.53553391, 21.53553391, 31.53553391, 41.53553391, 51.53553391])
    testRes = np.allclose(DB['DBYl'], zSol, atol=atol)
    assert (testRes == True)
    zSol = np.array([4.46446609, 14.46446609, 24.46446609, 34.46446609, 44.46446609])
    testRes = np.allclose(DB['DBYr'], zSol, atol=atol)
    assert (testRes == True)
