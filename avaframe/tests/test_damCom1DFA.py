"""Tests for module com1DFAtools"""
import numpy as np
import math
import configparser
import pytest

# Local imports
import avaframe.com1DFA.damCom1DFA as damCom1DFA


def test_linesIntersect(capfd):
    '''
    '''
    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (1, 0)
    (xF1, yF1) = (0, -1)
    (xF2, yF2) = (0, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 0.5

    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (1, 0)
    (xF1, yF1) = (0.5, -1)
    (xF2, yF2) = (0.5, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 0.5

    (xOld, yOld) = (-1, 0.5)
    (xNew, yNew) = (1, 0.5)
    (xF1, yF1) = (0, -1)
    (xF2, yF2) = (0, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 0.75

    (xOld, yOld) = (-1, -1)
    (xNew, yNew) = (1, 1)
    (xF1, yF1) = (1, -1)
    (xF2, yF2) = (-1, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 0.5

    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (1, 0)
    (xF1, yF1) = (0, 0)
    (xF2, yF2) = (0, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 0

    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (1, 0)
    (xF1, yF1) = (0, -1)
    (xF2, yF2) = (0, 0)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 1

    (xOld, yOld) = (0, 0)
    (xNew, yNew) = (1, 0)
    (xF1, yF1) = (0, -1)
    (xF2, yF2) = (0, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 0.5

    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (0, 0)
    (xF1, yF1) = (0, -1)
    (xF2, yF2) = (0, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 1
    assert factor == 0.5

    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (-0.5, 0)
    (xF1, yF1) = (0, -1)
    (xF2, yF2) = (0, 1)
    intersect, factor = damCom1DFA.linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    assert intersect == 0
    assert factor == 0


def test_getIntersection(capfd):
    '''
    '''
    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (1, 0)
    (xFoot, yFoot, zFoot) = (np.array([-3., -1., 1., 3.]), np.array([3., 1., -1., -3.]), np.array([0., 0., 0., 0.]))
    (xCrown, yCrown, zCrown) = (np.array([-3., -1., 1., 3.]), np.array([3., 1., -1., -3.]), np.array([1., 2., 3., 4.]))
    (xTangent, yTangent, zTangent) = (np.array([1, 1, 1, 1])/math.sqrt(2), -np.array([1, 1, 1, 1])/math.sqrt(2), np.array([0., 0., 0., 0.]))
    nDamPoints = 4
    intersection, section, xF, yF, zF, xC, yC, zC, xT, yT, zT = damCom1DFA.getIntersection(xOld, yOld, xNew, yNew,
                                                                                  xFoot, yFoot, zFoot,
                                                                                  xCrown, yCrown, zCrown,
                                                                                  xTangent, yTangent, zTangent,
                                                                                  nDamPoints)
#    print(xT, yT, zT)
    assert intersection == 1
    assert xF == 0
    assert yF == 0
    assert zF == 0
    assert xC == 0
    assert yC == 0
    assert zC == 2.5
    atol = 1e-10
    assert xT == pytest.approx(1/math.sqrt(2), abs=atol)
    assert yT == pytest.approx(-1/math.sqrt(2), abs=atol)
    assert zT == 0

    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (1, 0)
    (xFoot, yFoot, zFoot) = (np.array([-1., 1., 3.]), np.array([1., -1., -3.]), np.array([0., 0., 0.]))
    (xCrown, yCrown, zCrown) = (np.array([-1., 1., 3.]), np.array([1., -1., -3.]), np.array([2., 3., 4.]))
    (xTangent, yTangent, zTangent) = (np.array([1, 1, 1])/math.sqrt(2), -np.array([1, 1, 1])/math.sqrt(2), np.array([0., 0., 0.]))
    nDamPoints = 3
    intersection, section, xF, yF, zF, xC, yC, zC, xT, yT, zT = damCom1DFA.getIntersection(xOld, yOld, xNew, yNew,
                                                                                  xFoot, yFoot, zFoot,
                                                                                  xCrown, yCrown, zCrown,
                                                                                  xTangent, yTangent, zTangent,
                                                                                  nDamPoints)

    assert intersection == 1
    assert xF == 0
    assert yF == 0
    assert zF == 0
    assert xC == 0
    assert yC == 0
    assert zC == 2.5
    atol = 1e-10
    assert xT == pytest.approx(1/math.sqrt(2), abs=atol)
    assert yT == pytest.approx(-1/math.sqrt(2), abs=atol)
    assert zT == 0

    (xOld, yOld) = (-1, 0)
    (xNew, yNew) = (1, 0)
    (xFoot, yFoot, zFoot) = (np.array([-3., -1., 1.]), np.array([3., 1., -1.]), np.array([0., 0., 0.]))
    (xCrown, yCrown, zCrown) = (np.array([-3., -1., 1.]), np.array([3., 1., -1.]), np.array([1., 2., 3.]))
    (xTangent, yTangent, zTangent) = (np.array([1, 1, 1])/math.sqrt(2), -np.array([1, 1, 1])/math.sqrt(2), np.array([0., 0., 0.]))
    nDamPoints = 3
    intersection, section, xF, yF, zF, xC, yC, zC, xT, yT, zT = damCom1DFA.getIntersection(xOld, yOld, xNew, yNew,
                                                                                  xFoot, yFoot, zFoot,
                                                                                  xCrown, yCrown, zCrown,
                                                                                  xTangent, yTangent, zTangent,
                                                                                  nDamPoints)

    assert intersection == 1
    assert xF == 0
    assert yF == 0
    assert zF == 0
    assert xC == 0
    assert yC == 0
    assert zC == 2.5
    atol = 1e-10
    assert xT == pytest.approx(1/math.sqrt(2), abs=atol)
    assert yT == pytest.approx(-1/math.sqrt(2), abs=atol)
    assert zT == 0


def test_initializeWallLines(capfd):
    '''
    '''
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'dam': 'True', 'damSlope': '45'}
    originalHeader = {}
    originalHeader['xllcenter'] = -5
    originalHeader['yllcenter'] = -5
    originalHeader['cellsize'] = 1.0
    originalHeader['nodata_value'] = -9999
    originalHeader['nrows'] = 20
    originalHeader['ncols'] = 20
    demHeader = {}
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0
    demHeader['cellsize'] = 1.0
    demHeader['nodata_value'] = -9999
    demHeader['nrows'] = 20
    demHeader['ncols'] = 20
    dem = {'header': demHeader, 'originalHeader': originalHeader, 'rasterData': np.zeros((20, 20)),
           'Nx': np.zeros((20, 20)), 'Ny': np.zeros((20, 20)), 'Nz': np.ones((20, 20))}
    wallLineDict = {'x': np.array([0., 0., 0., 1., 2.]), 'y': np.array([2., 1., 0., 0., 0.]),
                    'z': np.array([10., 10., 10., 10., 10.]), 'fileName': 'test.shp', 'Name': 'dam', 'slope': 45}
    wallLineDict = damCom1DFA.initializeWallLines(cfg['GENERAL'], dem, wallLineDict)
#    print(wallLineDict)
#    print(wallLineDict['cellsCrossed'])
#    print(wallLineDict['xTangent'])
#    print(wallLineDict['yTangent'])
#    print(wallLineDict['zTangent'])

    assert wallLineDict['dam'] == 1
    assert wallLineDict['nPoints'] == 5
    assert wallLineDict['height'][0] == 10
    assert np.size(wallLineDict['height']) == 5
    assert wallLineDict['slope'][0] == math.pi/4
    assert np.size(wallLineDict['slope']) == 5
    assert np.allclose(wallLineDict['x'], np.array([15, 15, 12.07106781, 6, 7]), atol=1.e-8)
    assert np.allclose(wallLineDict['y'], np.array([7, 6, 12.07106781, 15, 15]), atol=1.e-8)
    assert np.allclose(wallLineDict['z'], np.zeros((5)), atol=1.e-10)
    assert np.allclose(wallLineDict['xTangent'], np.array([0, 0, 1, 1]), atol=1.e-10)
    assert np.allclose(wallLineDict['yTangent'], np.array([-1, -1, 0, 0]), atol=1.e-10)
    assert np.allclose(wallLineDict['zTangent'], np.zeros((4)), atol=1.e-10)
    assert np.allclose(wallLineDict['xCrown'], np.array([5, 5, 5, 6, 7]), atol=1.e-8)
    assert np.allclose(wallLineDict['yCrown'], np.array([7, 6, 5, 5, 5]), atol=1.e-8)
    assert np.allclose(wallLineDict['zCrown'], 10*np.ones((5)), atol=1.e-10)

    wallLineDict = damCom1DFA.initializeWallLines(cfg['GENERAL'], dem, None)
    assert wallLineDict['dam'] == 0


def test_getWallInteraction(capfd):
    '''
    '''
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'restitutionCoefficient': '1'}
    (nrows, ncols, csz) = (31, 31, 1)
    originalHeader = {}
    originalHeader['xllcenter'] = -15
    originalHeader['yllcenter'] = -15
    originalHeader['cellsize'] = csz
    originalHeader['nodata_value'] = -9999
    originalHeader['nrows'] = nrows
    originalHeader['ncols'] = ncols
    demHeader = {}
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0
    demHeader['cellsize'] = csz
    demHeader['nodata_value'] = -9999
    demHeader['nrows'] = nrows
    demHeader['ncols'] = ncols
    dem = {'header': demHeader, 'originalHeader': originalHeader, 'rasterData': np.zeros((nrows, ncols)),
           'Nx': np.zeros((nrows, ncols)), 'Ny': np.zeros((nrows, ncols)), 'Nz': np.ones((nrows, ncols))}
    wallLineDict = {'x': np.array([0., 0.]), 'y': np.array([-10, 10.]),
                    'z': np.array([10., 10.]), 'fileName': 'test.shp', 'Name': 'dam', 'slope': 90}
    wallLineDict = damCom1DFA.initializeWallLines(cfg['GENERAL'], dem, wallLineDict)
    interpOption = 2
    restitutionCoefficient = 1
    (xOld, yOld, zOld) = (-1-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (xNew, yNew, zNew) = (2-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (uxNew, uyNew, uzNew) = (10, 0, 0)
    fI, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm = damCom1DFA.getWallInteraction(
        xOld, yOld, zOld, xNew, yNew, zNew, uxNew, uyNew, uzNew, wallLineDict['nPoints'],
        wallLineDict['x'], wallLineDict['y'], wallLineDict['z'],
        wallLineDict['xCrown'], wallLineDict['yCrown'], wallLineDict['zCrown'],
        wallLineDict['xTangent'], wallLineDict['yTangent'], wallLineDict['zTangent'],
        ncols, nrows, csz, interpOption, restitutionCoefficient, 2,
        dem['Nx'], dem['Ny'], dem['Nz'], dem['rasterData'], np.zeros((nrows, ncols)))
#    print(xNew, yNew, zNew, uxNew, uyNew, uzNew)
    atol = 1e-6
    assert fI == 1
    assert xNew+originalHeader['xllcenter'] == -2
    assert yNew+originalHeader['yllcenter'] == 0
    assert zNew == 0
    assert uxNew == pytest.approx(-10, abs=atol)
    assert uyNew == 0
    assert uzNew == 0

    restitutionCoefficient = 0.5
    (xOld, yOld, zOld) = (-1-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (xNew, yNew, zNew) = (2-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (uxNew, uyNew, uzNew) = (10, 0, 0)
    fI, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm = damCom1DFA.getWallInteraction(
        xOld, yOld, zOld, xNew, yNew, zNew, uxNew, uyNew, uzNew, wallLineDict['nPoints'],
        wallLineDict['x'], wallLineDict['y'], wallLineDict['z'],
        wallLineDict['xCrown'], wallLineDict['yCrown'], wallLineDict['zCrown'],
        wallLineDict['xTangent'], wallLineDict['yTangent'], wallLineDict['zTangent'],
        ncols, nrows, csz, interpOption, restitutionCoefficient, 2,
        dem['Nx'], dem['Ny'], dem['Nz'], dem['rasterData'], np.zeros((nrows, ncols)))
#    print(xNew, yNew, zNew, uxNew, uyNew, uzNew)
    atol = 1e-6
    assert fI == 1
    assert xNew+originalHeader['xllcenter'] == -1
    assert yNew+originalHeader['yllcenter'] == 0
    assert zNew == 0
    assert uxNew == pytest.approx(-5, abs=atol)
    assert uyNew == 0
    assert uzNew == 0

    restitutionCoefficient = 0
    (xOld, yOld, zOld) = (-1-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (xNew, yNew, zNew) = (2-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (uxNew, uyNew, uzNew) = (10, 0, 0)
    fI, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm = damCom1DFA.getWallInteraction(
        xOld, yOld, zOld, xNew, yNew, zNew,
        uxNew, uyNew, uzNew,
        wallLineDict['nPoints'],
        wallLineDict['x'], wallLineDict['y'], wallLineDict['z'],
        wallLineDict['xCrown'], wallLineDict['yCrown'], wallLineDict['zCrown'],
        wallLineDict['xTangent'], wallLineDict['yTangent'], wallLineDict['zTangent'],
        ncols, nrows, csz, interpOption, restitutionCoefficient, 2,
        dem['Nx'], dem['Ny'], dem['Nz'], dem['rasterData'], np.zeros((nrows, ncols)))
#    print(xNew, yNew, zNew, uxNew, uyNew, uzNew)
    assert fI == 1
    assert xNew+originalHeader['xllcenter'] == 0
    assert yNew+originalHeader['yllcenter'] == 0
    assert zNew == 0
    assert uxNew == 0
    assert uyNew == 0
    assert uzNew == 0

    restitutionCoefficient = 1
    (xOld, yOld, zOld) = (-1-originalHeader['xllcenter'], -5-originalHeader['yllcenter'], 0)
    (xNew, yNew, zNew) = (2-originalHeader['xllcenter'], 5-originalHeader['yllcenter'], 0)
    (uxNew, uyNew, uzNew) = (10, 3, 0)
    fI, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm = damCom1DFA.getWallInteraction(
        xOld, yOld, zOld, xNew, yNew, zNew,
        uxNew, uyNew, uzNew,
        wallLineDict['nPoints'],
        wallLineDict['x'], wallLineDict['y'], wallLineDict['z'],
        wallLineDict['xCrown'], wallLineDict['yCrown'], wallLineDict['zCrown'],
        wallLineDict['xTangent'], wallLineDict['yTangent'], wallLineDict['zTangent'],
        ncols, nrows, csz, interpOption, restitutionCoefficient, 2,
        dem['Nx'], dem['Ny'], dem['Nz'], dem['rasterData'], np.zeros((nrows, ncols)))
#    print(xNew, yNew, zNew, uxNew, uyNew, uzNew)
    atol = 1e-6
    assert fI == 1
    assert xNew+originalHeader['xllcenter'] == pytest.approx(-2, abs=atol)
    assert yNew+originalHeader['yllcenter'] == pytest.approx(5, abs=atol)
    assert zNew == 0
    assert uxNew == pytest.approx(-10, abs=atol)
    assert uyNew == pytest.approx(3, abs=atol)
    assert uzNew == 0

    restitutionCoefficient = 0
    (xOld, yOld, zOld) = (-1-originalHeader['xllcenter'], -5-originalHeader['yllcenter'], 0)
    (xNew, yNew, zNew) = (2-originalHeader['xllcenter'], 5-originalHeader['yllcenter'], 0)
    (uxNew, uyNew, uzNew) = (10, 3, 0)
    fI, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm = damCom1DFA.getWallInteraction(
        xOld, yOld, zOld, xNew, yNew, zNew,
        uxNew, uyNew, uzNew,
        wallLineDict['nPoints'],
        wallLineDict['x'], wallLineDict['y'], wallLineDict['z'],
        wallLineDict['xCrown'], wallLineDict['yCrown'], wallLineDict['zCrown'],
        wallLineDict['xTangent'], wallLineDict['yTangent'], wallLineDict['zTangent'],
        ncols, nrows, csz, interpOption, restitutionCoefficient, 2,
        dem['Nx'], dem['Ny'], dem['Nz'], dem['rasterData'], np.zeros((nrows, ncols)))
#    print(xNew, yNew, zNew, uxNew, uyNew, uzNew)
    assert fI == 1
    assert xNew+originalHeader['xllcenter'] == pytest.approx(0, abs=atol)
    assert yNew+originalHeader['yllcenter'] == pytest.approx(-5/3, abs=atol)
    assert zNew == 0
    assert uxNew == 0
    assert uyNew == 0
    assert uzNew == 0

    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'restitutionCoefficient': '1'}
    (nrows, ncols, csz) = (31, 31, 1)
    originalHeader = {}
    originalHeader['xllcenter'] = -15
    originalHeader['yllcenter'] = -15
    originalHeader['cellsize'] = csz
    originalHeader['nodata_value'] = -9999
    originalHeader['nrows'] = nrows
    originalHeader['ncols'] = ncols
    demHeader = {}
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0
    demHeader['cellsize'] = csz
    demHeader['nodata_value'] = -9999
    demHeader['nrows'] = nrows
    demHeader['ncols'] = ncols
    dem = {'header': demHeader, 'originalHeader': originalHeader, 'rasterData': np.zeros((nrows, ncols)),
           'Nx': np.zeros((nrows, ncols)), 'Ny': np.zeros((nrows, ncols)), 'Nz': np.ones((nrows, ncols))}
    wallLineDict = {'x': np.array([1., 1.]), 'y': np.array([-10, 10.]),
                    'z': np.array([1., 1.]), 'fileName': 'test.shp', 'Name': 'dam', 'slope': 45}
    wallLineDict = damCom1DFA.initializeWallLines(cfg['GENERAL'], dem, wallLineDict)
    restitutionCoefficient = 1
    (xOld, yOld, zOld) = (-1-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (xNew, yNew, zNew) = (2-originalHeader['xllcenter'], 0-originalHeader['yllcenter'], 0)
    (uxNew, uyNew, uzNew) = (10, 0, 0)
    fI, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm = damCom1DFA.getWallInteraction(
        xOld, yOld, zOld, xNew, yNew, zNew,
        uxNew, uyNew, uzNew,
        wallLineDict['nPoints'],
        wallLineDict['x'], wallLineDict['y'], wallLineDict['z'],
        wallLineDict['xCrown'], wallLineDict['yCrown'], wallLineDict['zCrown'],
        wallLineDict['xTangent'], wallLineDict['yTangent'], wallLineDict['zTangent'],
        ncols, nrows, csz, interpOption, restitutionCoefficient, 2,
        dem['Nx'], dem['Ny'], dem['Nz'], dem['rasterData'], np.zeros((nrows, ncols)))
#    print(xNew, yNew, zNew, uxNew, uyNew, uzNew)
    atol = 1e-6
    assert fI == 1
    assert xNew+originalHeader['xllcenter'] == pytest.approx(0, abs=atol)
    assert yNew+originalHeader['yllcenter'] == pytest.approx(0, abs=atol)
    assert zNew == pytest.approx(0, abs=atol)  # because of the reprojection...
    assert uxNew == pytest.approx(0, abs=atol)
    assert uyNew == pytest.approx(0, abs=atol)
    assert uzNew == pytest.approx(0, abs=atol)  # because of the reprojection...

    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'restitutionCoefficient': '1'}
    (nrows, ncols, csz) = (31, 31, 1)
    originalHeader = {}
    originalHeader['xllcenter'] = -15
    originalHeader['yllcenter'] = -15
    originalHeader['cellsize'] = csz
    originalHeader['nodata_value'] = -9999
    originalHeader['nrows'] = nrows
    originalHeader['ncols'] = ncols
    demHeader = {}
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0
    demHeader['cellsize'] = csz
    demHeader['nodata_value'] = -9999
    demHeader['nrows'] = nrows
    demHeader['ncols'] = ncols
    dem = {'header': demHeader, 'originalHeader': originalHeader, 'rasterData': np.zeros((nrows, ncols)),
           'Nx': np.zeros((nrows, ncols)), 'Ny': np.zeros((nrows, ncols)), 'Nz': np.ones((nrows, ncols))}
    wallLineDict = {'x': np.array([1., 1.]), 'y': np.array([-10, 10.]),
                    'z': np.array([1., 1.]), 'fileName': 'test.shp', 'Name': 'dam', 'slope': 45}
    wallLineDict = damCom1DFA.initializeWallLines(cfg['GENERAL'], dem, wallLineDict)
    restitutionCoefficient = 1
    (xOld, yOld, zOld) = (-1-originalHeader['xllcenter'], -5-originalHeader['yllcenter'], 0)
    (xNew, yNew, zNew) = (2-originalHeader['xllcenter'], 5-originalHeader['yllcenter'], 0)
    (uxNew, uyNew, uzNew) = (10, 3, 0)
    fI, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm = damCom1DFA.getWallInteraction(
        xOld, yOld, zOld, xNew, yNew, zNew,
        uxNew, uyNew, uzNew,
        wallLineDict['nPoints'],
        wallLineDict['x'], wallLineDict['y'], wallLineDict['z'],
        wallLineDict['xCrown'], wallLineDict['yCrown'], wallLineDict['zCrown'],
        wallLineDict['xTangent'], wallLineDict['yTangent'], wallLineDict['zTangent'],
        ncols, nrows, csz, interpOption, restitutionCoefficient, 2,
        dem['Nx'], dem['Ny'], dem['Nz'], dem['rasterData'], np.zeros((nrows, ncols)))
#    print(xNew, yNew, zNew, uxNew, uyNew, uzNew)
    atol = 1e-10
    assert fI == 1
    assert xNew+originalHeader['xllcenter'] == pytest.approx(0, abs=atol)
    assert yNew+originalHeader['yllcenter'] == pytest.approx(5, abs=atol)
    assert zNew == pytest.approx(0, abs=atol)  # because of the reprojection...
    assert uxNew == pytest.approx(0, abs=atol)
    assert uyNew == pytest.approx(3, abs=atol)
    assert uzNew == pytest.approx(0, abs=atol)
