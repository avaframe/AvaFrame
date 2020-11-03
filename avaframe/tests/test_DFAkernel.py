"""Tests for module geoTrans"""
import numpy as np
import pytest

# Local imports
import avaframe.DFAkernel.DFAtools as DFAtools
import avaframe.in3Utils.ascUtils as IOf


def test_polygon2Raster(capfd):
    demHeader = IOf.cASCheader()
    demHeader.ncols = 6
    demHeader.nrows = 4
    demHeader.xllcenter = 0
    demHeader.yllcenter = 0
    demHeader.cellsize = 1
    x = np.array([1.4, 2.9])
    y = np.array([0.4, 2.4])
    Line = {}
    Line['x'] = x
    Line['y'] = y
    mask = DFAtools.polygon2Raster(demHeader, Line)
    print(mask)
    # assert 1 == 2


def test_normalize(capfd):
    '''normalize'''
    x = np.array([1.])
    y = np.array([1.])
    z = np.array([1.])
    xn, yn, zn = DFAtools.normalize(x, y, z)
    atol = 1e-10
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1., rel=atol)
    assert xn == 1/np.sqrt(3.)
    assert yn == 1/np.sqrt(3.)
    assert zn == 1/np.sqrt(3.)

    x = np.array([1.])
    y = np.array([2.])
    z = np.array([3.])
    xn, yn, zn = DFAtools.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1., rel=atol)

    x = np.array([1.])
    y = np.array([0.])
    z = np.array([1.])
    xn, yn, zn = DFAtools.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1, rel=atol)
    assert xn == pytest.approx(1/np.sqrt(2.), rel=atol)
    assert yn == pytest.approx(0, rel=atol)
    assert zn == pytest.approx(1/np.sqrt(2.), rel=atol)


def test_getNormalVect(capfd):
    '''projectOnRasterVect'''
    a = 1
    xllcenter = 0
    yllcenter = 0
    cellsize = 1
    ncols = 15
    nrows = 10
    x = np.linspace(0, ncols-1, ncols)
    y = np.linspace(0, nrows-1, nrows)
    X, Y = np.meshgrid(x, y)
    Z = a * X * X
    Nx, Ny, Nz = DFAtools.getNormalVect(Z, cellsize)

    atol = 1e-10
    TestNX = np.allclose(Nx[1:nrows-1, 1:ncols-1], -2*a*X[1:nrows-1, 1:ncols-1]/np.sqrt(1+4*a*a*X[1:nrows-1, 1:ncols-1]*X[1:nrows-1, 1:ncols-1]), atol=atol)
    assert TestNX
    TestNY = np.allclose(Ny[1:nrows-1, 1:ncols-1], np.zeros(np.shape(X[1:nrows-1, 1:ncols-1])), atol=atol)
    assert TestNY
    TestNZ = np.allclose(Nz[1:nrows-1, 1:ncols-1], 1/np.sqrt(1+4*a*a*X[1:nrows-1, 1:ncols-1]*X[1:nrows-1, 1:ncols-1]), atol=atol)
    assert TestNZ

    Z = a * Y
    Nx, Ny, Nz = DFAtools.getNormalVect(Z, cellsize)

    atol = 1e-10
    TestNY = np.allclose(Ny[1:nrows-1, 1:ncols-1], -a/np.sqrt(1+a*a*np.ones(np.shape(Y[1:nrows-1, 1:ncols-1]))), atol=atol)
    assert TestNY
    TestNX = np.allclose(Nx[1:nrows-1, 1:ncols-1], np.zeros(np.shape(X[1:nrows-1, 1:ncols-1])), atol=atol)
    assert TestNX
    TestNZ = np.allclose(Nz[1:nrows-1, 1:ncols-1], 1/np.sqrt(1+a*a*np.ones(np.shape(Y[1:nrows-1, 1:ncols-1]))), atol=atol)
    assert TestNZ


def test_getNeighbours(capfd):
    header = IOf.cASCheader()
    header.ncols = 6
    header.nrows = 4
    header.cellsize = 1
    dem = {}
    dem['header'] = header
    particles = {}
    particles['Npart'] = 13
    # 2 part in cell 0, 3 in 5, 2 in 8, 1 in 9, 1 in 14, 2 in 15, 1 in 18, 1 in 23
    particles['x'] = np.array([-0.4, 0.2,
                               5.4, 5.1, 4.6,
                               1.6, 2.2,
                               3,
                               2.1,
                               2.9, 3.2,
                               -0.4,
                               5.4])
    particles['y'] = np.array([-0.4, 0.2,
                               -0.4, 0.2, 0,
                               1.4, 0.9,
                               1,
                               2.1,
                               2.4, 1.8,
                               3.4,
                               3.4])
    particles = DFAtools.getNeighbours(particles, dem)
    print(particles['InCell'])
    print(particles['indPartInCell'])
    print(particles['partInCell'])
    atol = 1e-10
    indPCell = np.array([2.,  # found 2 particles
                         2., 2., 2., 2.,  # nothing happens
                         5.,  # found 3 particles
                         5., 5.,  # nothing happens
                         7.,   # found 1 particles
                         8.,  # found 1 particles
                         8., 8., 8., 8.,  # nothing happens
                         9.,  # found 1 particles
                         11.,  # found 2 particles
                         11., 11.,  # nothing happens
                         12.,  # found 1 particles
                         12., 12., 12., 12.,  # nothing happens
                         13.,  # found 1 particles
                         0.])  # always an extra zero at the end
    assert np.allclose(particles['indPartInCell'], indPCell, atol=atol)
    pInC = np.array([1., 0.,
                     4., 3., 2.,
                     6., 5.,
                     7.,
                     8.,
                     10., 9.,
                     11.,
                     12.])
    assert np.allclose(particles['partInCell'], pInC, atol=atol)
