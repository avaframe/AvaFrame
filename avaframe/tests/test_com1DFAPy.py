"""Tests for module com1DFAtools"""
import numpy as np
import pytest
import copy
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Local imports
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.com1DFAPy.SPHfunctions as SPH
import avaframe.in2Trans.ascUtils as IOf


def test_polygon2Raster(capfd):
    demHeader = IOf.cASCheader()
    demHeader.ncols = 15
    demHeader.nrows = 8
    demHeader.xllcenter = 0
    demHeader.yllcenter = 0
    demHeader.cellsize = 1
    x = np.array([0, 2, 4, 6, 8, 10, 12, 14, 2, 0])
    y = np.array([1, 3, 2, 3, 0, 3, 1, 6, 6, 1])
    Line = {}
    Line['x'] = x
    Line['y'] = y
    mask = np.zeros((demHeader.nrows, demHeader.ncols))
    mask = com1DFA.polygon2Raster(demHeader, Line, mask)
    Mask = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                     [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0],
                     [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                     [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                     [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                     [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    assert np.allclose(Mask, mask,  atol=1e-10)


def test_normalize(capfd):
    '''normalize'''
    x = np.array([1.])
    y = np.array([1.])
    z = np.array([1.])
    norme = DFAtls.norm(x, y, z)
    norme2 = DFAtls.norm2(x, y, z)
    xn, yn, zn = DFAtls.normalize(x, y, z)
    atol = 1e-10
    assert norme == np.sqrt(3.)
    assert norme2 == 3.
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1., rel=atol)
    assert xn == 1/np.sqrt(3.)
    assert yn == 1/np.sqrt(3.)
    assert zn == 1/np.sqrt(3.)

    x = np.array([1.])
    y = np.array([2.])
    z = np.array([3.])
    xn, yn, zn = DFAtls.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1., rel=atol)

    x = np.array([1.])
    y = np.array([0.])
    z = np.array([1.])
    xn, yn, zn = DFAtls.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1, rel=atol)
    assert xn == pytest.approx(1/np.sqrt(2.), rel=atol)
    assert yn == pytest.approx(0, rel=atol)
    assert zn == pytest.approx(1/np.sqrt(2.), rel=atol)


def test_getNormalMesh(capfd):
    '''projectOnRasterVect'''
    a = 2
    b = 1
    xllcenter = 0
    yllcenter = 0
    cellsize = 1
    m = 10
    n = 15
    x = np.linspace(0, m-1, m)
    y = np.linspace(0, n-1, n)
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y
    Z1 = a * X * X + b * Y * Y
    for num in [1, 4, 6, 8]:
        Nx, Ny, Nz = DFAtls.getNormalMesh(Z, cellsize, num=num)
        atol = 1e-10
        TestNX = np.allclose(Nx[1:n-1, 1:m-1], (-a*np.ones(np.shape(Y)) / np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1], atol=atol)
        assert TestNX
        TestNY = np.allclose(Ny[1:n-1, 1:m-1], (-b*np.ones(np.shape(Y)) / np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1], atol=atol)
        assert TestNY
        TestNZ = np.allclose(Nz[1:n-1, 1:m-1], (np.ones(np.shape(Y)) / np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1], atol=atol)
        assert TestNZ

        Nx, Ny, Nz = DFAtls.getNormalMesh(Z1, cellsize, num=num)

        print(Nx)
        print((-2*a*X / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1])
        print(Ny)
        print((-2*b*Y / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1])
        print(Nz)
        print((1 / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1])
        atol = 1e-10
        TestNX = np.allclose(Nx[1:n-1, 1:m-1], (-2*a*X / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1], atol=atol)
        assert TestNX
        TestNY = np.allclose(Ny[1:n-1, 1:m-1], (-2*b*Y / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1], atol=atol)
        assert TestNY
        TestNZ = np.allclose(Nz[1:n-1, 1:m-1], (1 / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1], atol=atol)
        assert TestNZ


def test_getAreaMesh(capfd):
    '''projectOnRasterVect'''
    a = 0.1
    b = 0.2
    xllcenter = 0
    yllcenter = 0
    csz = 1
    m = 15
    n = 10
    x = np.linspace(0, m-1, m)
    y = np.linspace(0, n-1, n)
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y
    Z1 = a * X * X + b * Y * Y
    Nx, Ny, Nz = DFAtls.getNormalMesh(Z, csz, num=4)
    Area = DFAtls.getAreaMesh(Nx, Ny, Nz, csz)
    print(np.sqrt((1+a*a+b*b)))
    print(Area)
    atol = 1e-10
    TestArea = np.allclose(Area[1:n-1, 1:m-1], np.sqrt((1+a*a+b*b)) * np.ones(np.shape(Y[1:n-1, 1:m-1])), atol=atol)
    assert TestArea


def test_getNeighbours(capfd):
    header = IOf.cASCheader()
    header.ncols = 4
    header.nrows = 5
    header.cellsize = 1
    dem = {}
    dem['header'] = header
    particles = {}
    particles['Npart'] = 16
    particles['x'] = np.array([1, 0, 1, 2, 1, 2, 0, 1, 0, 2, 0, 2, 1, 2, 3, 3])
    particles['y'] = np.array([2, 1, 0, 1, 3, 3, 2, 1, 0, 0, 3, 2, 2, 1, 1, 4])
    particles['z'] = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    particles['m'] = particles['z']
    particlesVect = copy.deepcopy(particles)
    atol = 1e-10
    indPCell = np.array([0.,  # always an extra zero at the begining
                         1,  # found 1 particle
                         2,  # found 1 particle
                         3,  # found 1 particle
                         3,  # nothing happens
                         4,  # found 1 particle
                         5,  # found 1 particle
                         7,   # found 2 particles
                         8,  # found 1 particle
                         9,  # found 1 particle
                         11,   # found 2 particles
                         12,  # found 1 particle
                         12,  # nothing happens
                         13,  # found 1 particle
                         14,  # found 1 particle
                         15,  # found 1 particle
                         15, 15, 15, 15,  # nothing happens
                         16])  # found 1 particle
    pInC = np.array([8,
                     2,
                     9,
                     1,
                     7,
                     3, 13,
                     14,
                     6,
                     0, 12,
                     11,
                     10,
                     4,
                     5,
                     15])

    particles = SPH.getNeighbours(particles, dem)
    print(particles['InCell'])
    print(particles['indPartInCell'])
    print(particles['partInCell'])
    assert np.allclose(particles['indPartInCell'], indPCell, atol=atol)
    assert np.allclose(particles['partInCell'], pInC, atol=atol)
    particlesVect = SPH.getNeighboursVect(particlesVect, dem)
    print(particlesVect['InCell'])
    print(particlesVect['indPartInCell'])
    print(particlesVect['partInCell'])
    assert np.allclose(particlesVect['indPartInCell'], indPCell, atol=atol)
    assert np.allclose(particlesVect['partInCell'], pInC, atol=atol)


def test_calcGradHSPH(capfd):
    header = IOf.cASCheader()
    header.ncols = 4
    ncols = header.ncols
    header.nrows = 5
    nrows = header.nrows
    header.cellsize = 1
    csz = header.cellsize
    dem = {}
    dem['header'] = header
    particles = {}
    particles['Npart'] = 16
    particles['x'] = np.array([1, 0, 1, 2, 1, 2, 0, 1, 0, 2, 0, 2, 1, 2, 3, 3])
    particles['y'] = np.array([2, 1, 0, 1, 3, 3, 2, 1, 0, 0, 3, 2, 2, 1, 1, 4])
    particles['z'] = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    particles['m'] = particles['z']
    particles = SPH.getNeighboursVect(particles, dem)
    print(particles['InCell'])
    print(particles['indPartInCell'])
    print(particles['partInCell'])
    _, _, _, index = SPH.calcGradHSPH(particles, 0, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 0, ncols, nrows, csz)
    print(index)
    atol = 1e-10
    IndexTh = np.array([1, 7, 3, 13, 6, 12, 11, 10, 4, 5])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
    _, _, _, index = SPH.calcGradHSPH(particles, 7, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 7, ncols, nrows, csz)
    print(index)
    IndexTh = np.array([8, 2, 9, 1, 3, 13, 6, 0, 12, 11])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
    _, _, _, index = SPH.calcGradHSPH(particles, 8, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 8, ncols, nrows, csz)
    print(index)
    IndexTh = np.array([2, 1, 7])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
    _, _, _, index = SPH.calcGradHSPH(particles, 9, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 9, ncols, nrows, csz)
    print(index)
    IndexTh = np.array([2, 7, 3, 13, 14])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
    _, _, _, index = SPH.calcGradHSPH(particles, 11, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 11, ncols, nrows, csz)
    print(index)
    IndexTh = np.array([7, 3, 13, 14, 0, 12, 4, 5])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
    _, _, _, index = SPH.calcGradHSPH(particles, 6, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 6, ncols, nrows, csz)
    print(index)
    IndexTh = np.array([1, 7, 0, 12, 10, 4])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
    _, _, _, index = SPH.calcGradHSPH(particles, 5, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 5, ncols, nrows, csz)
    print(index)
    IndexTh = np.array([0, 12, 11, 4, 15])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
    _, _, _, index = SPH.calcGradHSPH(particles, 15, ncols, nrows, csz)
    _, _, _, indexVect = SPH.calcGradHSPHVect(particles, 15, ncols, nrows, csz)
    print(index)
    IndexTh = np.array([5])
    assert np.allclose(index, IndexTh, atol=atol)
    assert np.allclose(indexVect, IndexTh, atol=atol)
