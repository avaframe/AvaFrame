"""Tests for module geoTrans"""
import numpy as np
import math
import pytest
import logging
import matplotlib.path as mpltPath
import pathlib
import shutil
import os
import configparser

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU


log = logging.getLogger(__name__)


def test_projectOnRaster(capfd):
    '''projectOnRaster'''
    dem = {}
    Points = {}
    header = {}
    header['xllcenter'] = 0
    header['yllcenter'] = 0
    header['cellsize'] = 1
    # make sure it works for 2D arrays too
    rasterdata = np.array(([0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]))
    Points['x'] = np.array(([0.4, 0], [0.5, 1.6], [2.4, 3.4], [2.4, -1], [4, 0]))
    Points['y'] = np.array(([0, 1.4], [0.5, 0.6], [-0.4, 1.4], [2.4, -1], [3, np.nan]))

    dem['header'] = header
    dem['rasterData'] = rasterdata

    Points, _ = geoTrans.projectOnRaster(dem, Points, interp='nearest')
    zSol = np.array([[0, 1], [0, 3], [2, 4], [4, np.nan], [np.nan, np.nan]])
    print(Points['z'])
    tol = 1e-8
    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=tol)
    assert (testRes)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=tol)
    assert (testRes)

    Points, _ = geoTrans.projectOnRaster(dem, Points, interp='bilinear')
    zSol = np.array([[0.4, 1.4], [1, 2.2], [np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]])
    print(Points['z'])
    zSolnan = np.isnan(zSol)
    testRes = np.allclose(np.isnan(Points['z']), zSolnan, atol=tol)
    assert (testRes)
    zSol = zSol[~np.isnan(zSol)]
    testRes = np.allclose(Points['z'][~np.isnan(Points['z'])], zSol, atol=tol)
    assert (testRes)


def test_resizeData(capfd):
    '''resizeData'''
    a = 2
    b = 1
    m = 10
    n = 15
    csz1 = 5
    x = np.linspace(0, m-1, m) * csz1
    y = np.linspace(0, n-1, n) * csz1
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y
    raster1 = {}
    header1 = {}
    header1['ncols'] = m
    header1['nrows'] = n
    header1['xllcenter'] = 0
    header1['yllcenter'] = 0
    header1['cellsize'] = csz1
    raster1['header'] = header1
    raster1['rasterData'] = Z

    m = 13
    n = 18
    csz2 = 4
    x = np.linspace(0, m-1, m) * csz2 + 0.8
    y = np.linspace(0, n-1, n) * csz2 + 0.8
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y
    raster2 = {}
    header2 = {}
    header2['ncols'] = m
    header2['nrows'] = n
    header2['xllcenter'] = 0.8
    header2['yllcenter'] = 0.8
    header2['cellsize'] = csz2
    raster2['header'] = header2
    raster2['rasterData'] = Z

    data1, data2 = geoTrans.resizeData(raster1, raster2)
    zSol = np.zeros(np.shape(data2))
    print(data1)
    print(data2)
    tol = 1e-10
    testRes = np.allclose(np.nan_to_num(data1-data2), zSol, atol=tol)
    assert (testRes)


def test_findAngleProfile(capfd):
    '''findAngleProfile'''
    s = np.linspace(0, 400, 41)
    angle = np.linspace(40, 0, 41)
    tmp = np.where((angle < 10.0) & (angle >= 0.0))
    tmp = np.asarray(tmp).flatten()
    ds = np.abs(s - np.roll(s, 1))
    ds = ds[tmp]
    dsMin = 30
    ids10Point = geoTrans.findAngleProfile(tmp, ds, dsMin)
    assert ids10Point == 31

    dsMin = 10
    ids10Point = geoTrans.findAngleProfile(tmp, ds, dsMin)
    assert ids10Point == 31

    angle[10] = 8
    angle[11] = 8
    angle[12] = 8
    tmp = np.where((angle < 10.0) & (angle >= 0.0))
    tmp = np.asarray(tmp).flatten()
    ds = np.abs(s - np.roll(s, 1))
    ds = ds[tmp]
    dsMin = 30
    ids10Point = geoTrans.findAngleProfile(tmp, ds, dsMin)
    assert ids10Point == 31

    angle[13] = 8
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    tmp = np.asarray(tmp).flatten()
    ids10Point = geoTrans.findAngleProfile(tmp, ds, dsMin)
    ds = np.abs(s - np.roll(s, 1))
    ds = ds[tmp]
    assert ids10Point == 10

    tmp = np.where((angle < 0.0) & (angle >= 0.0))
    tmp = np.asarray(tmp).flatten()
    ds = np.abs(s - np.roll(s, 1))
    ds = ds[tmp]
    dsMin = 10
    # do we react properly when the input line exceeds the dem?
    with pytest.raises(Exception) as e:
        assert geoTrans.findAngleProfile(tmp, ds, dsMin)
    assert str(e.value) == "No point found. Check the angle and threshold distance."

    tmp = np.where((angle < 4.0) & (angle >= 0.0))
    tmp = np.asarray(tmp).flatten()
    ds = np.abs(s - np.roll(s, 1))
    ds = ds[tmp]
    dsMin = 40
    # do we react properly when the input line exceeds the dem?
    with pytest.raises(Exception) as e:
        assert geoTrans.findAngleProfile(tmp, ds, dsMin)
    assert str(e.value) == "No point found. Check the angle and threshold distance."


def test_prepareAngleProfile(caplog):
    '''prepareAngleProfile'''
    AvaProfile = {}
    s = np.linspace(0, 400, 41)
    AvaProfile['s'] = s
    AvaProfile['z'] = s*s/400 - 2*s + 400
    theta = -np.arctan(2*s/400 - 2)*180/math.pi
    beta = 10
    with caplog.at_level(logging.WARNING):
        angle, tmp, ds = geoTrans.prepareAngleProfile(beta, AvaProfile)
    assert 'No split Point given!' in caplog.text
    AvaProfile['indSplit'] = 12
    angle, tmp, ds = geoTrans.prepareAngleProfile(beta, AvaProfile)
    print(theta)
    print(angle)
    print(tmp)
    print(ds)
    index = np.array([37, 38, 39, 40])
    tol = 1e-8
    testRes = np.allclose(tmp, index, atol=tol)
    assert testRes
    assert np.allclose(ds, 10*np.ones(4), atol=tol)
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
    csz = 2
    w = 10
    xyPath = {}
    xyPath['x'] = np.array((0, 10, 20, 30, 40))
    xyPath['y'] = np.array((10, 20, 30, 40, 50))
    rasterTransfo = {}
    rasterTransfo['cellSizeSL'] = csz
    rasterTransfo['domainWidth'] = w

    DB = geoTrans.path2domain(xyPath, rasterTransfo)

    atol = 1e-8
    # Compare result to reference solution
    zSol = xyPath['x'] - w / np.sqrt(2) / 2
    print(zSol)
    print(DB['DBXl']*csz)
    testRes = np.allclose(DB['DBXl']*csz, zSol, atol=atol)
    assert testRes
    zSol = xyPath['x'] + w / np.sqrt(2) / 2
    testRes = np.allclose(DB['DBXr']*csz, zSol, atol=atol)
    assert testRes
    zSol = xyPath['y'] + w / np.sqrt(2) / 2
    testRes = np.allclose(DB['DBYl']*csz, zSol, atol=atol)
    assert testRes
    zSol = xyPath['y'] - w / np.sqrt(2) / 2
    testRes = np.allclose(DB['DBYr']*csz, zSol, atol=atol)
    assert testRes


def test_areaPoly(capfd):
    '''test_areaPoly'''
    a = 1
    R = a/(2*math.sin(math.pi/5))
    A = 5*a*a/4*math.sqrt(1+2/math.sqrt(5))
    theta = 2*math.pi/5*np.arange(5)
    x = R*np.cos(theta)
    y = R*np.sin(theta)

    area = geoTrans.areaPoly(x, y)
    tol = 1e-14
    assert area == pytest.approx(A, rel=tol)


def test_remeshData(tmp_path):
    """ test shape of interpolated data onto new mesh """

    headerInfo = {}
    headerInfo['cellsize'] = 5
    headerInfo['ncols'] = 4
    headerInfo['nrows'] = 5
    headerInfo['xllcenter'] = 0
    headerInfo['yllcenter'] = 0
    headerInfo['noDataValue'] = -9999
    # create an inclined plane
    z0 = 10
    data = getIPZ(z0, 15, 20, 5)
    rasterDict = {'header': headerInfo, 'rasterData': data}
    # outFile = os.path.join(tmp_path, 'test.asc')
    # IOf.writeResultToAsc(headerInfo, data, outFile, flip=False)
    atol = 1.e-10
    dataNew = geoTrans.remeshData(rasterDict, 2., larger=False)
    dataRaster = dataNew['rasterData']
    indNoData = np.where(dataRaster == -9999)
    headerNew = dataNew['header']
    xExtent = (headerNew['ncols']-1) * headerNew['cellsize']
    yExtent = (headerNew['nrows']-1) * headerNew['cellsize']

    # compute solution
    dataSol = getIPZ(z0, xExtent, yExtent, headerNew['cellsize'])

    # compare solution to result from function
    testRes = np.allclose(dataRaster, dataSol, atol=atol)

    assert dataNew['rasterData'].shape[0] == 11
    assert dataNew['rasterData'].shape[1] == 8
    assert len(indNoData[0]) == 0
    assert np.isclose(dataNew['rasterData'][0, 0], 10.)
    assert testRes

    dataNew = geoTrans.remeshData(rasterDict, 2., remeshOption='interp2d', interpMethod='cubic', larger=False)
    dataRaster = dataNew['rasterData']
    indNoData = np.where(dataRaster == -9999)
    headerNew = dataNew['header']
    xExtent = (headerNew['ncols']-1) * headerNew['cellsize']
    yExtent = (headerNew['nrows']-1) * headerNew['cellsize']

    # compute solution
    dataSol = getIPZ(z0, xExtent, yExtent, headerNew['cellsize'])

    # compare solution to result from function
    testRes = np.allclose(dataRaster, dataSol, atol=atol)

    assert dataNew['rasterData'].shape[0] == 11
    assert dataNew['rasterData'].shape[1] == 8
    assert len(indNoData[0]) == 0
    assert np.isclose(dataNew['rasterData'][0, 0], 10.)
    assert testRes

    dataNew = geoTrans.remeshData(rasterDict, 2., remeshOption='RectBivariateSpline', interpMethod='cubic', larger=False)
    dataRaster = dataNew['rasterData']
    indNoData = np.where(dataRaster == -9999)
    headerNew = dataNew['header']
    xExtent = (headerNew['ncols']-1) * headerNew['cellsize']
    yExtent = (headerNew['nrows']-1) * headerNew['cellsize']

    # compute solution
    dataSol = getIPZ(z0, xExtent, yExtent, headerNew['cellsize'])

    # compare solution to result from function
    testRes = np.allclose(dataRaster, dataSol, atol=atol)

    assert dataNew['rasterData'].shape[0] == 11
    assert dataNew['rasterData'].shape[1] == 8
    assert len(indNoData[0]) == 0
    assert np.isclose(dataNew['rasterData'][0, 0], 10.)
    assert testRes


def test_remeshDEM(tmp_path):
    """ test size of interpolated data onto new mesh """

    headerInfo = {}
    headerInfo['cellsize'] = 5
    headerInfo['ncols'] = 4
    headerInfo['nrows'] = 5
    headerInfo['xllcenter'] = 0
    headerInfo['yllcenter'] = 0
    headerInfo['noDataValue'] = -9999
    # create an inclined plane
    z0 = 10
    data = getIPZ(z0, 15, 20, 5)

    avaDir = pathlib.Path(tmp_path, 'avaTest')
    fU.makeADir(avaDir)
    fU.makeADir((avaDir / 'Inputs'))
    avaDEM = avaDir / 'Inputs' / 'avaAlr.asc'
    IOf.writeResultToAsc(headerInfo, data, avaDEM, flip=True)

    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'meshCellSizeThreshold': '0.0001', 'meshCellSize': '2.', 'avalancheDir': str(avaDir)}

    # call function
    pathDem = geoTrans.remeshDEM(avaDEM, cfg)
    fullP = avaDir / 'Inputs' / pathDem
    dataNew = IOf.readRaster(fullP)

    dataRaster = dataNew['rasterData']
    indNoData = np.where(dataRaster == -9999)
    headerNew = dataNew['header']
    xExtent = (headerNew['ncols']-1) * headerNew['cellsize']
    yExtent = (headerNew['nrows']-1) * headerNew['cellsize']

    # compute solution
    dataSol = getIPZ(z0, xExtent, yExtent, 2.)
    # dataSol = getIPZ(z0, xExtent, yExtent, headerNew['cellsize'])

    # compare solution to result from function
    testRes = np.allclose(dataRaster[:-1, :-1], dataSol[:-1, :-1], atol=1.e-6)
    print(dataRaster)
    print(dataSol)

    assert dataNew['rasterData'].shape[0] == 11
    assert dataNew['rasterData'].shape[1] == 8
    assert len(indNoData[0]) == 0
    assert testRes
    assert np.isclose(dataRaster[0, 0], 10.)

    # copy data
    avaName = 'avaParabola'
    dirPath = pathlib.Path(__file__).parents[0]
    inputDir1 = dirPath / '..' / 'data' / avaName
    inputDEM = dirPath / 'data' / 'remeshedDEM8.00.asc'
    avaDir1 = pathlib.Path(tmp_path, avaName)
    avaDEM = avaDir1 / 'Inputs' / 'DEMremeshed' / 'DEM_PF_Topo_remeshedDEM8.00.asc'
    shutil.copytree(inputDir1, avaDir1)
    inputsAVA = avaDir1 / 'Inputs' / 'DEMremeshed'
    fU.makeADir(inputsAVA)
    shutil.copy(inputDEM, avaDEM)
    inputDEM1 = inputDir1 / 'Inputs' / 'DEM_PF_Topo.asc'
    avaDEM1 = avaDir1 / 'Inputs' / 'DEM_PF_Topo.asc'
    shutil.copy(inputDEM1, avaDEM1)
    cfg['GENERAL']['avalancheDir'] = str(avaDir1)
    cfg['GENERAL']['meshCellSize'] = '8.'

    # call function
    pathDem2 = geoTrans.remeshDEM(avaDEM1, cfg)
    fullP2 = avaDir1 / 'Inputs' / pathDem2
    dataNew2 = IOf.readRaster(fullP2)
    dataRaster2 = dataNew2['rasterData']
    dataSol = IOf.readRaster(inputDEM)
    # compare solution to result from function
    testRes2 = np.allclose(dataRaster2, dataSol['rasterData'], atol=1.e-6)

    assert dataNew2['rasterData'].shape[0] == dataSol['header']['nrows']
    assert dataNew2['rasterData'].shape[1] == dataSol['header']['ncols']
    assert testRes2

    dataMod = IOf.readRaster(inputDEM)
    dataMod['header']['cellsize'] = 9.0
    IOf.writeResultToAsc(dataMod['header'], dataMod['rasterData'], avaDEM, flip=True)

    with pytest.raises(FileExistsError) as e:
        assert geoTrans.remeshDEM(avaDEM1, cfg)
    assert str(e.value) == ("Name for saving remeshedDEM already used: %s" % avaDEM.name)


def test_isCounterClockWise(capfd):
    """ test isCounterClockWise """

    xCoord = np.array([0, 1, 2, 0.5])
    yCoord = np.array([0, 1, -1, -2])

    polygon = np.stack((xCoord, yCoord), axis=-1)

    path = mpltPath.Path(polygon)
    is_ccw = geoTrans.isCounterClockWise(path)

    assert is_ccw is not True


def test_checkOverlap(capfd):
    """ test checkOverlap function with both crop set to False and True
    """
    nRows = 5
    nCols = 6
    refRaster = np.zeros((nRows, nCols))
    refRaster[:3, :3] = 1
    checkedRaster1 = np.zeros((nRows, nCols))
    checkedRaster1[2:, 2:] = 1
    checkedRaster1[2, 2] = 0

    toCheckRaster = np.zeros((nRows, nCols))
    toCheckRaster[2:, 2:] = 1
    print(refRaster)
    print(toCheckRaster)
    checkedRaster2 = geoTrans.checkOverlap(toCheckRaster, refRaster, 'to Check', 'ref', crop=True)
    atol = 1e-10
    print(checkedRaster2)
    print(checkedRaster1)
    assert np.allclose(checkedRaster2, checkedRaster1, atol=atol)

    toCheckRaster = np.zeros((nRows, nCols))
    toCheckRaster[2:, 2:] = 1
    with pytest.raises(AssertionError, match=r"to Check area features overlapping with ref area - this is not allowed"):
        checkedRaster2 = geoTrans.checkOverlap(toCheckRaster, refRaster, 'to Check', 'ref', crop=False)

    toCheckRaster = np.zeros((nRows, nCols))
    toCheckRaster[3:, 3:] = 1
    checkedRaster2 = geoTrans.checkOverlap(toCheckRaster, refRaster, 'to Check', 'ref', crop=True)
    toCheckRaster = np.zeros((nRows, nCols))
    toCheckRaster[3:, 3:] = 1
    assert np.allclose(checkedRaster2, toCheckRaster, atol=atol)
    checkedRaster2 = geoTrans.checkOverlap(toCheckRaster, refRaster, 'to Check', 'ref', crop=False)
    toCheckRaster = np.zeros((nRows, nCols))
    toCheckRaster[3:, 3:] = 1
    assert np.allclose(checkedRaster2, toCheckRaster, atol=atol)


def test_computeS(capfd):
    avaPath = {}
    avaPath['x'] = np.array([1, 4, 0, 3])
    avaPath['y'] = np.array([0, 4, 1, -3])
    avaPath = geoTrans.computeS(avaPath)
    print(avaPath)
    atol = 1e-10
    assert len(avaPath['s']) == 4
    assert np.allclose(avaPath['s'], np.array([0, 5, 10, 15]), atol=atol)


def test_makeCoordinateGrid(capfd):
    xllc = 1
    yllc = -1
    csz = 2
    ncols = 3
    nrows = 4
    x, y = geoTrans.makeCoordinateGrid(xllc, yllc, csz, ncols, nrows)
    print(x)
    print(y)

    assert x[0, 0] == 1
    assert x[0, 1] == 3
    assert x[0, 2] == 5
    assert y[0, 1] == -1
    assert y[1, 1] == 1
    assert y[3, 1] == 5

    rasterHeader = {'ncols': ncols, 'nrows': nrows, 'xllcenter': xllc, 'yllcenter': yllc, 'cellsize': csz}
    x, y, ncols, nrows = geoTrans.makeCoordGridFromHeader(rasterHeader, cellSizeNew=None, larger=True)
    print(x)
    print(y)
    print(ncols)
    print(nrows)

    assert x[0, 0] == 1
    assert x[0, 1] == 3
    assert x[0, 2] == 5
    assert y[0, 1] == -1
    assert y[1, 1] == 1
    assert y[3, 1] == 5

    rasterHeader = {'ncols': ncols, 'nrows': nrows, 'xllcenter': xllc, 'yllcenter': yllc, 'cellsize': csz}
    x, y, ncols, nrows = geoTrans.makeCoordGridFromHeader(rasterHeader, cellSizeNew=1, larger=True)
    print(x)
    print(y)
    print(ncols)
    print(nrows)

    assert x[0, 0] == 1
    assert x[0, 1] == 2
    assert x[0, 2] == 3
    assert y[0, 1] == -1
    assert y[1, 1] == 0
    assert y[3, 1] == 2


def getIPZ(z0, xEnd, yEnd, dx):

    meanAlpha = 30.
    nstepsX = int((xEnd + dx) / dx)
    nstepsY = int((yEnd + dx) / dx)
    xv = np.linspace(0, xEnd, nstepsX)
    yv = np.linspace(0, yEnd, nstepsY)
    nRows = len(yv)
    nCols = len(xv)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nRows, nCols))
    # Set surface elevation from slope and max. elevation
    zv = z0 - np.tan(np.radians(meanAlpha)) * x

    return zv


def test_cartToSpherical():
    """ test converting to spherical coordinates """

    # setup required inputs
    X = 10.
    Y = 10.
    Z = np.sqrt(10.**2 + 10.**2)

    # call function to be tested
    r, phi, theta = geoTrans.cartToSpherical(X, Y, Z)

    assert r == np.sqrt(10.**2 + 10.**2 + Z**2)
    assert phi == 45.
    assert theta == 45.

    # setup required inputs
    X = np.sqrt(10.)
    Y = 3*np.sqrt(10.)
    Z = 20

    # call function to be tested
    r, phi, theta = geoTrans.cartToSpherical(X, Y, Z)

    print('r', r, 'phi', phi, 'theta', theta)

    assert r == np.sqrt(X**2 + Y**2 + Z**2)
    assert phi == np.rad2deg(np.arctan(2)) - 45
    assert theta == 90. - np.rad2deg(np.arctan(2))


def test_rotate():
    """ test rotate a vector by an angle """

    # setup required data
    locationPoints = [[0., 0.], [0., 10.]]
    theta = 45.
    deg = True

    # call function to be tested
    rotatedLine = geoTrans.rotate(locationPoints, theta, deg=deg)

    sL = np.sqrt(50.)
    assert rotatedLine[0][0] == 0.
    assert -sL - 1.e-12 < rotatedLine[0][1] < -sL+1.e-12
    assert rotatedLine[1][0] == 0.
    assert sL - 1.e-12 < rotatedLine[1][1] < sL+1.e-12
    assert isinstance(rotatedLine, list)
    assert len(rotatedLine) == 2
    assert len(rotatedLine[0]) == 2

    locationPoints = [[0., 0.], [0., 10.]]
    theta = np.rad2deg(np.arctan(2)) - 45
    deg = True

    # call function to be tested
    rotatedLine = geoTrans.rotate(locationPoints, theta, deg=deg)
    print('rotated line', rotatedLine)

    sLX = np.sqrt(10.)
    assert rotatedLine[0][0] == 0.
    assert -sLX - 1.e-12 < rotatedLine[0][1] < -sLX+1.e-12
    assert rotatedLine[1][0] == 0.
    assert 3*sLX - 1.e-12 < rotatedLine[1][1] < 3*sLX+1.e-12
    assert isinstance(rotatedLine, list)
    assert len(rotatedLine) == 2
    assert len(rotatedLine[0]) == 2

    vector = np.diff(locationPoints)
    vectorRot = np.diff(rotatedLine)

    # check if rotation worked correctly
    lenVector = np.sqrt(vector[0]**2 + vector[1]**2)
    lenVectorRot = np.sqrt(vectorRot[0]**2 + vectorRot[1]**2)

    assert np.isclose(lenVector, lenVectorRot)


def test_rotateRaster():
    rasterData = np.zeros((5, 5))
    rasterData[2, 0:2] = 1
    rasterData[2, 3:4] = -1
    header = {'xllcenter': -2, 'yllcenter': -2, 'cellsize': 1, 'ncols': 5, 'nrows': 5}
    rasterDict = {'header': header, 'rasterData': rasterData}
    theta = 90
    print(header)
    print(rasterData)
    rotatedRaster = geoTrans.rotateRaster(rasterDict, theta, deg=True)
    print(rotatedRaster['rasterData'][0:-1, 2])
    assert np.allclose(rotatedRaster['rasterData'][0:-1, 2], np.array([1, 1, 0, -1]), atol=1e-6)
