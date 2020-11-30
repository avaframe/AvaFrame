"""Tests for module geoTrans"""
import numpy as np
import math
import pytest

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.ascUtils as IOf


def test_projectOnRaster(capfd):
    '''projectOnRaster'''
    dem = {}
    Points = {}
    header = IOf.cASCheader()
    header.xllcenter = 0
    header.yllcenter = 0
    header.cellsize = 1

    rasterdata = np.array(([0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]))
    Points['x'] = np.array((0.4, 0, 0.5, 1.6, 2.4, 3.4, 2.4, -1, 4, 0))
    Points['y'] = np.array((0, 1.4, 0.5, 0.6, -0.4, 1.4, 2.4, -1, 3, np.nan))

    dem['header'] = header
    dem['rasterData'] = rasterdata
    Points = geoTrans.projectOnRaster(dem, Points, interp='nearest')
    zSol = np.array([0, 1, 0, 3, 2, 4, 4, np.nan, np.nan, np.nan])
    print(Points['z'])
    tol = 1e-8
    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=tol)
    assert (testRes)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=tol)
    assert (testRes)

    Points = geoTrans.projectOnRaster(dem, Points, interp='bilinear')
    zSol = np.array([0.4, 1.4, 1, 2.2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    print(Points['z'])
    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=tol)
    assert (testRes)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=tol)
    assert (testRes)


def test_projectOnRasterVect(capfd):
    '''projectOnRasterVect'''
    dem = {}
    Points = {}
    header = IOf.cASCheader()
    header.xllcenter = 0
    header.yllcenter = 0
    header.cellsize = 1

    rasterdata = np.array(([0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]))
    Points['x'] = np.array(([0.4, 0], [0.5, 1.6], [2.4, 3.4], [2.4, -1], [4, 0]))
    Points['y'] = np.array(([0, 1.4], [0.5, 0.6], [-0.4, 1.4], [2.4, -1], [3, np.nan]))

    dem['header'] = header
    dem['rasterData'] = rasterdata

    Points, _ = geoTrans.projectOnRasterVect(dem, Points, interp='nearest')
    zSol = np.array([[0, 1], [0, 3], [2, 4], [4, np.nan], [np.nan, np.nan]])
    print(Points['z'])
    tol = 1e-8
    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=tol)
    assert (testRes)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=tol)
    assert (testRes)

    Points, _ = geoTrans.projectOnRasterVect(dem, Points, interp='bilinear')
    zSol = np.array([[0.4, 1.4], [1, 2.2], [np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]])
    print(Points['z'])
    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=tol)
    assert (testRes)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=tol)
    assert (testRes)


def test_pointsToRaster(capfd):
    '''pointsToRaster'''

    rasterdata = np.array(([0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]))
    x = np.array((0.4, 0, 0.5, 1.6, 2.4, 2.4, 1.4))
    y = np.array((0, 1.4, 0.5, 0.6, 0.4, 1.4, 1.4))
    z = np.array((1., 1., 1., 1., 1., 1., 1.))

    Z = geoTrans.pointsToRaster(x, y, z, rasterdata, csz=1, xllc=0, yllc=0, interp='bilinear')
    zSol = np.array(([0.85, 0.81, 0.6, 0.24], [0.85, 0.85, 1.2, 0.4], [0.4, 0.24, 0.4, 0.16]))
    print(Z)
    tol = 1e-8
    testRes = np.allclose(Z, zSol, atol=tol)
    assert (testRes)

    Z = geoTrans.pointsToRaster(x, y, z, rasterdata, csz=1, xllc=0, yllc=0, interp='nearest')
    zSol = np.array(([2, 0, 1, 0], [1, 1, 2, 0], [0, 0, 0, 0]))
    print(Z)
    tol = 1e-8
    testRes = np.allclose(Z, zSol, atol=tol)
    assert (testRes)


def test_findAngleProfile(capfd):
    '''findAngleProfile'''
    s = np.linspace(0, 400, 41)
    angle = np.linspace(40, 0, 41)
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    deltaInd = 3
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 31

    deltaInd = 1
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 31

    angle[10] = 8
    angle[11] = 8
    angle[12] = 8
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    deltaInd = 3
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 31

    angle[13] = 8
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    assert ids10Point == 10


def test_prepareAngleProfile(capfd):
    '''prepareAngleProfile'''
    AvaProfile = {}
    s = np.linspace(0, 400, 41)
    AvaProfile['s'] = s
    AvaProfile['z'] = s*s/400 - 2*s + 400
    theta = -np.arctan(2*s/400 - 2)*180/math.pi
    AvaProfile['indSplit'] = 12
    beta = 10
    angle, tmp, deltaInd = geoTrans.prepareAngleProfile(beta, AvaProfile)
    print(theta)
    print(angle)
    print(tmp)
    print(deltaInd)
    index = np.array([37, 38, 39, 40])
    tol = 1e-8
    testRes = np.allclose(tmp, index, atol=tol)
    assert testRes
    assert deltaInd == 3
    smaller = np.where(angle[1:39] > theta[0:38])
    print(smaller)
    assert np.size(smaller) == 0
    greater = np.where(angle[1:40] < theta[1:40])
    print(greater)
    assert np.size(greater) == 0


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


def test_path2domain(capfd):
    '''test_path2domain'''
    xllc = 1
    yllc = 2
    csz = 2
    w = 10
    xyPath = {}
    xyPath['x'] = np.array((0, 10, 20, 30, 40))
    xyPath['y'] = np.array((10, 20, 30, 40, 50))
    rasterTransfo = {}
    rasterTransfo['xllc'] = xllc
    rasterTransfo['yllc'] = yllc
    rasterTransfo['cellsize'] = csz
    rasterTransfo['domainWidth'] = w

    DB = geoTrans.path2domain(xyPath, rasterTransfo)

    atol = 1e-8
    # Compare result to reference solution
    zSol = np.array([-2.26776695,  2.73223305,  7.73223305, 12.73223305,
                    17.73223305])
    testRes = np.allclose(DB['DBXl'], zSol, atol=atol)
    assert (testRes == True)
    zSol = np.array([1.26776695,  6.26776695, 11.26776695, 16.26776695,
                    21.26776695])
    testRes = np.allclose(DB['DBXr'], zSol, atol=atol)
    assert (testRes == True)
    zSol = np.array([5.76776695, 10.76776695, 15.76776695, 20.76776695,
                    25.76776695])
    testRes = np.allclose(DB['DBYl'], zSol, atol=atol)
    assert (testRes == True)
    zSol = np.array([2.23223305,  7.23223305, 12.23223305, 17.23223305,
                    22.23223305])
    testRes = np.allclose(DB['DBYr'], zSol, atol=atol)
    assert (testRes == True)
