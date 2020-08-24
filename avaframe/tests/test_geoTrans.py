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

    rasterdata = np.array(([0,1,2],[2, 0.5, 5]))
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
    header.xllcorner = 0
    header.yllcorner = 0
    header.cellsize = 1
    header.ncols = 4
    header.nrows = 3

    rasterdata = np.array(([0, 1, 2, 3],[4, 5, 6, 7],[8, 9, 10, 11]))
    Points['x'] = np.array(([0.0000000000001,0.00000000001],[-5, 1.67],[2.9999999999, 1],[8, 1]))
    Points['y'] = np.array(([0.0000000000001,1],[2, 0.6],[1.4, 5],[1.6, 5]))

    dem['header'] = header
    dem['rasterData'] = rasterdata

    Points, itot, ioob = geoTrans.projectOnRasterVect(dem, Points, interp='bilinear')
    tol = 1e-8
    assert (Points['z'][0,0] == pytest.approx(0, abs=tol))
    assert (Points['z'][0,1] == pytest.approx(4, abs=tol))
    assert np.isnan(Points['z'][1,0])
    assert (Points['z'][1,1] == pytest.approx(4.07, abs=tol))
    assert (Points['z'][2,0] == pytest.approx(8.6, abs=tol))
    assert np.isnan(Points['z'][2,1])
    assert np.isnan(Points['z'][3,0])
    assert np.isnan(Points['z'][3,1])

    Points, itot, ioob = geoTrans.projectOnRasterVect(dem, Points, interp='nearest')
    tol = 1e-8
    assert (Points['z'][0,0] == pytest.approx(0, abs=tol))
    assert (Points['z'][0,1] == pytest.approx(4, abs=tol))
    assert np.isnan(Points['z'][1,0])
    assert (Points['z'][1,1] == pytest.approx(6, abs=tol))
    assert (Points['z'][2,0] == pytest.approx(7, abs=tol))
    assert np.isnan(Points['z'][2,1])
    assert np.isnan(Points['z'][3,0])
    assert np.isnan(Points['z'][3,1])


def test_checkProfile(capfd):
    '''checkProfile'''

    AvaProfile = {}
    AvaProfile['x'] = np.array((0, 1 , 2, 3, 4))
    AvaProfile['y'] = np.array((1, 2 , 3, 4, 5))
    AvaProfile['z'] = np.array((0, 1 , 2, 3, 4))

    projSplitPoint, AvaProfile = geoTrans.checkProfile(AvaProfile, projSplitPoint=None)
    assert AvaProfile['z'][0]==4
    assert AvaProfile['z'][-1]==0
    assert AvaProfile['x'][0]==4
    assert AvaProfile['x'][-1]==0
    assert AvaProfile['y'][0]==5
    assert AvaProfile['y'][-1]==1
    assert AvaProfile['y'][1]==4

    projSplitPoint = {}
    projSplitPoint['indSplit'] = 1
    AvaProfile = {}
    AvaProfile['x'] = np.array((0, 1 , 2, 3, 4))
    AvaProfile['y'] = np.array((1, 2 , 3, 4, 5))
    AvaProfile['z'] = np.array((0, 1 , 2, 3, 4))
    AvaProfile['s'] = np.array((0, 2 , 4, 6, 8))
    projSplitPoint, AvaProfile = geoTrans.checkProfile(AvaProfile, projSplitPoint=projSplitPoint)
    assert AvaProfile['z'][0]==4
    assert AvaProfile['z'][-1]==0
    assert AvaProfile['x'][0]==4
    assert AvaProfile['x'][-1]==0
    assert AvaProfile['y'][0]==5
    assert AvaProfile['y'][-1]==1
    assert AvaProfile['y'][1]==4
    assert AvaProfile['s'][0]==0
    assert AvaProfile['indSplit']==3

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
    # DB = path2domain(xyPath, w, header)
