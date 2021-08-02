"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
from avaframe.com1DFA import com1DFA
import avaframe.in2Trans.ascUtils as IOf
import pytest
import configparser



def test_setDEMOriginToZero():
    """ test if origin is set to zero """

    # setup required input
    tHeader = IOf.cASCheader()
    tHeader.xllcenter = 10.
    tHeader.yllcenter = 4.0
    dem = {'header': tHeader}

    # call function to be tested
    demTest = com1DFA.setDEMoriginToZero(dem)

    assert demTest['header'].xllcenter == 0.0
    assert demTest['header'].yllcenter == 0.0


def test_placeParticles():
    """ test placing of particles """

    # setup required inputs
    massCell = 10.
    indx = 0
    indy = 1
    csz = 5
    massPerPart = 2.
    initPartDistType = 'uniform'
    rng = np.random.default_rng(12345)

    # call funciton to be tested - uniform
    xpart, ypart, mPart, nPart =  com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest = np.asarray([-1.66666666, 0.0, 1.66666666, -1.66666666, 0., 1.66666666, -1.66666666,
                        0.0, 1.66666666])
    ypartTest = np.asarray([3.33333333, 3.33333333, 3.33333333, 5.0, 5., 5., 6.66666666, 6.66666666,
                        6.66666666])

    # call funciton to be tested - uniform
    massCell = 8.
    xpart2, ypart2, mPart2, nPart2 =  com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest2 = np.asarray([-1.25, 1.25, -1.25, 1.25])
    ypartTest2 = np.asarray([3.75, 3.75, 6.25, 6.25])

    # call funciton to be tested - random
    massCell = 11.5
    initPartDistType = 'random'
    xpart3, ypart3, mPart3, nPart3 =  com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest3 = np.asarray([-0.9162083, 1.48682729, 0.88127335, -0.54445225, -0.83593036, 0.49154377])
    ypartTest3 = np.asarray([3.43367093, 5.86378022, 7.20901433, 3.74122857, 7.24440576, 5.83618727])

    assert nPart == 9.0
    assert np.isclose(mPart, 1.111111)
    assert np.allclose(xpart, xpartTest)
    assert np.allclose(ypart, ypartTest)
    assert nPart2 == 4.0
    assert mPart2 == 2.
    assert np.allclose(xpart2, xpartTest2)
    assert np.allclose(ypart2, ypartTest2)
    assert nPart3 == 6.0
    assert np.isclose(mPart3, 1.9166666666666)
    assert np.allclose(xpart3, xpartTest3)
    assert np.allclose(ypart3, ypartTest3)


def test_initializeMesh():
    """ test mesh initialization """

    # setup required input
    demHeader = IOf.cASCheader()
    demHeader.xllcenter = 101.23
    demHeader.yllcenter = 24.54
    demHeader.cellsize = 5.0
    demHeader.noDataValue = -9999
    demHeader.nrows = 5
    demHeader.ncols = 5

    # define plane with constant slope of 45°
    demData = np.asarray([[1., 2., 3., 4., 5.],
                          [1., 2., 3., 4., 5.],
                          [1., 2., 3., 4., 5.],
                          [1., 2., 3., 4., 5.],
                          [1., 2., 3., 4., 5.]])

    demOri = {'header': demHeader, 'rasterData': demData}
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'sphKernelRadius': '2.5', 'meshCellSizeThreshold': '0.0001',
                      'meshCellSize': '5.'}
    num = 1

    # setup testResults
    demNewHeader = IOf.cASCheader()
    demNewHeader.xllcenter = 0.
    demNewHeader.yllcenter = 0.
    demNewHeader.cellsize = 5.0
    demNewHeader.noDataValue = -9999
    demNewHeader.nrows = 5
    demNewHeader.ncols = 5
    demTest = {'header': demNewHeader}
    demTest['originOri'] =  {'xllcenter': 101.23, 'yllcenter': 24.54}
    demTest['outOfDEM'] = np.asarray([[False, False, False, False, True],
                                      [False, False, False, False, True],
                                      [False, False, False, False, True],
                                      [False, False, False, False, True],
                                      [False, False, False, False, True]])
    # normal vector of plane
    demTest['Nx'] = np.asarray([[0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.]]) - 25.
    demTest['Ny'] = np.asarray([[0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.]])
    demTest['Nz'] = np.asarray([[0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0.]]) + 25.
    # setup neighbour grid
    headerNeighbourGrid = IOf.cASCheader()
    headerNeighbourGrid.cellsize = 2.5
    headerNeighbourGrid.ncols = 10
    headerNeighbourGrid.nrows = 10
    headerNeighbourGrid.xllcenter = 0
    headerNeighbourGrid.yllcenter = 0
    demTest['headerNeighbourGrid'] = headerNeighbourGrid
    areaCell = np.sqrt(-25.*-25. + 0*0 + 25*25)
    demTest['areaRaster'] = np.zeros((5, 5)) + areaCell
    demTest['rasterData'] = demData

    # call function to be tested
    demOri, dem = com1DFA.initializeMesh(cfg['GENERAL'], demOri, num)

    assert dem['header'].xllcenter == demTest['header'].xllcenter
    assert dem['header'].yllcenter == demTest['header'].yllcenter
    assert dem['header'].ncols == demTest['header'].ncols
    assert dem['header'].nrows == demTest['header'].nrows
    assert dem['header'].cellsize == demTest['header'].cellsize
    assert dem['header'].yllcenter == demTest['header'].yllcenter
    assert np.array_equal(dem['rasterData'], demTest['rasterData'])
    assert dem['Nx'][4,4] == demTest['Nx'][4,4]
    assert dem['Ny'][4,4] == demTest['Ny'][4,4]
    assert dem['Nz'][4,4] == demTest['Nz'][4,4]
    assert np.array_equal(dem['areaRaster'], demTest['areaRaster'])
    assert dem['headerNeighbourGrid'].xllcenter == demTest['headerNeighbourGrid'].xllcenter
    assert dem['headerNeighbourGrid'].yllcenter == demTest['headerNeighbourGrid'].yllcenter
    assert dem['headerNeighbourGrid'].ncols == demTest['headerNeighbourGrid'].ncols
    assert dem['headerNeighbourGrid'].nrows == demTest['headerNeighbourGrid'].nrows
    assert dem['headerNeighbourGrid'].cellsize == demTest['headerNeighbourGrid'].cellsize
    assert dem['headerNeighbourGrid'].yllcenter == demTest['headerNeighbourGrid'].yllcenter


def test_polygon2Raster():
    """ test if polygon is converted to raster properly """

    # setup required inputs
    demHeader = IOf.cASCheader()
    demHeader.cellsize = 1
    demHeader.ncols = 10
    demHeader.nrows = 10
    demHeader.xllcenter = 0
    demHeader.yllcenter = 0

    Line = {'x': np.asarray([0, 1., 0.989, 0., 0.]), 'y': np.asarray([0., 0., 0.989, 1., 0.0])}
    radius = 0.0001
    th = 1.2

    # call function to be tested
    Mask = com1DFA.polygon2Raster(demHeader, Line, radius, th=th)

    # setup test output
    maskTest = np.zeros((10, 10))
    maskTest[0, 0:2] = 1.2
    maskTest[1, 0] = 1.2

    # call function to be tested
    Mask2 = com1DFA.polygon2Raster(demHeader, Line, 0.1, th=th)
    maskTest2 = maskTest.copy()
    maskTest2[1, 1] = 1.2

    assert np.array_equal(maskTest, Mask)
    assert np.array_equal(maskTest2, Mask2)


def test_getSimTypeList():
    """ test create list of simTypes """

    # setup required input
    simTypeList = ['ent', 'res', 'null', 'available', 'entres']
    inputSimFiles = {'entResInfo': {'flagEnt': 'Yes', 'flagRes': 'Yes'}}

    # call function to be tested
    simTypeList = com1DFA.getSimTypeList(simTypeList, inputSimFiles)

    # setup test result
    simTypeListTest = ['ent', 'null', 'res', 'entres']

    assert set(simTypeListTest).issubset(simTypeList)
    assert 'available' not in simTypeList
