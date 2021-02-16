"""Tests for module com1DFAtools"""
import numpy as np
import pytest
import copy
import configparser

# Local imports
import avaframe.in2Trans.ascUtils as IOf
import avaframe.com1DFAPy.DFAfunctionsCython as DFAfunC
import avaframe.com1DFAPy.DFAtools as DFAtls


def test_normalizeC(capfd):
    '''normalize'''
    x = np.array([1.])
    y = np.array([1.])
    z = np.array([1.])
    norme = DFAfunC.normpy(x, y, z)
    norme = np.asarray(norme)
    norme2 = DFAfunC.norm2py(x, y, z)
    norme2 = np.asarray(norme2)
    xn, yn, zn = DFAfunC.normalizepy(x, y, z)
    xn = np.asarray(xn)
    yn = np.asarray(yn)
    zn = np.asarray(zn)
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
    xn, yn, zn = DFAfunC.normalizepy(x, y, z)
    xn = np.asarray(xn)
    yn = np.asarray(yn)
    zn = np.asarray(zn)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1., rel=atol)

    x = np.array([1.])
    y = np.array([0.])
    z = np.array([1.])
    xn, yn, zn = DFAfunC.normalizepy(x, y, z)
    xn = np.asarray(xn)
    yn = np.asarray(yn)
    zn = np.asarray(zn)
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
    for num in [4, 6, 8]:
        Nx, Ny, Nz = DFAtls.getNormalMesh(Z, cellsize, num=num)
        atol = 1e-10
        TestNX = np.allclose(Nx[1:n-1, 1:m-1], (-a*np.ones(np.shape(Y)) /
                                                np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1], atol=atol)
        assert TestNX
        TestNY = np.allclose(Ny[1:n-1, 1:m-1], (-b*np.ones(np.shape(Y)) /
                                                np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1], atol=atol)
        assert TestNY
        TestNZ = np.allclose(Nz[1:n-1, 1:m-1], (np.ones(np.shape(Y)) /
                                                np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1], atol=atol)
        assert TestNZ

        Nx, Ny, Nz = DFAtls.getNormalMesh(Z1, cellsize, num=num)

        print(Nx)
        print((-2*a*X / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1])
        print(Ny)
        print((-2*b*Y / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1])
        print(Nz)
        print((1 / np.sqrt(1 + 4*a*a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1])
        atol = 1e-10
        TestNX = np.allclose(Nx[1:n-1, 1:m-1], (-2*a*X / np.sqrt(1 + 4*a *
                                                                 a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1], atol=atol)
        assert TestNX
        TestNY = np.allclose(Ny[1:n-1, 1:m-1], (-2*b*Y / np.sqrt(1 + 4*a *
                                                                 a*X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1], atol=atol)
        assert TestNY
        TestNZ = np.allclose(Nz[1:n-1, 1:m-1], (1 / np.sqrt(1 + 4*a*a *
                                                            X*X + 4*b*b*Y*Y))[1:n-1, 1:m-1], atol=atol)
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
    TestArea = np.allclose(Area[1:n-1, 1:m-1], np.sqrt((1+a*a+b*b)) *
                           np.ones(np.shape(Y[1:n-1, 1:m-1])), atol=atol)
    assert TestArea


def test_getNeighboursC(capfd):
    header = IOf.cASCheader()
    header.ncols = 4
    header.nrows = 5
    header.cellsize = 1
    dem = {}
    dem['header'] = header
    particles = {}
    particles['Npart'] = 16
    particles['x'] = np.array([1., 0., 1., 2., 1., 2., 0., 1., 0., 2., 0., 2., 1., 2., 3., 3.])
    particles['y'] = np.array([2., 1., 0., 1., 3., 3., 2., 1., 0., 0., 3., 2., 2., 1., 1., 4.])
    particles['z'] = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
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

    particles = DFAfunC.getNeighboursC(particles, dem)
    print(particles['InCell'])
    print(particles['indPartInCell'])
    print(particles['partInCell'])
    assert np.allclose(particles['indPartInCell'], indPCell, atol=atol)
    assert np.allclose(particles['partInCell'], pInC, atol=atol)


def test_calcGradHSPHC(capfd):
    header = IOf.cASCheader()
    header.ncols = 4
    ncols = header.ncols
    header.nrows = 5
    nrows = header.nrows
    header.cellsize = 1
    csz = header.cellsize
    dem = {}
    dem['header'] = header
    cfg = configparser.ConfigParser()
    cfg.add_section('GENERAL')
    cfg.optionxform = str
    cfg.set('GENERAL', 'rho', '200')
    cfg.set('GENERAL', 'minRKern', '0.001')
    cfg.set('GENERAL', 'gravAcc', '9.81')
    cfg = cfg['GENERAL']
    particles = {}
    particles['Npart'] = 16
    particles['x'] = np.array([1., 0., 1., 2., 1., 2., 0., 1., 0., 2., 0., 2., 1., 2., 3., 3.])
    particles['y'] = np.array([2., 1., 0., 1., 3., 3., 2., 1., 0., 0., 3., 2., 2., 1., 1., 4.])
    particles['z'] = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    particles['ux'] = particles['z']
    particles['uy'] = particles['z']
    particles['uz'] = particles['z']
    particles['m'] = particles['z']
    particles['h'] = particles['z']
    particles = DFAfunC.getNeighboursC(particles, dem)
    print(particles['InCell'])
    print(particles['indPartInCell'])
    print(particles['partInCell'])
    indX = (particles['indX']).astype('int')
    indY = (particles['indY']).astype('int')
    Nx = np.zeros((nrows, ncols))
    Ny = np.zeros((nrows, ncols))
    Nz = np.ones((nrows, ncols))
    _, _, _, L, index = DFAfunC.computeGradC(cfg, particles, header, Nx, Ny, Nz, indX, indY, 1, 0, Pytest=1)
    L = np.asarray(L)
    index = np.asarray(index)
    print(L)
    cindex = np.cumsum(index)
    print(cindex)
    atol = 1e-10
    IndexTh = np.array([1, 7, 3, 13, 6, 12, 11, 10, 4, 5])
    print(L[cindex[0]:cindex[1]])
    assert np.allclose(L[cindex[0]:cindex[1]], IndexTh, atol=atol)
    IndexTh = np.array([8, 2, 9, 1, 3, 13, 6, 0, 12, 11])
    print(L[cindex[7]:cindex[8]])
    assert np.allclose(L[cindex[7]:cindex[8]], IndexTh, atol=atol)
    IndexTh = np.array([2, 1, 7])
    print(L[cindex[8]:cindex[9]])
    assert np.allclose(L[cindex[8]:cindex[9]], IndexTh, atol=atol)
    IndexTh = np.array([2, 7, 3, 13, 14])
    print(L[cindex[9]:cindex[10]])
    assert np.allclose(L[cindex[9]:cindex[10]], IndexTh, atol=atol)
    IndexTh = np.array([7, 3, 13, 14, 0, 12, 4, 5])
    print(L[cindex[11]:cindex[12]])
    assert np.allclose(L[cindex[11]:cindex[12]], IndexTh, atol=atol)
    IndexTh = np.array([1, 7, 0, 12, 10, 4])
    print(L[cindex[6]:cindex[7]])
    assert np.allclose(L[cindex[6]:cindex[7]], IndexTh, atol=atol)
    IndexTh = np.array([0, 12, 11, 4, 15])
    print(L[cindex[5]:cindex[6]])
    assert np.allclose(L[cindex[5]:cindex[6]], IndexTh, atol=atol)
    IndexTh = np.array([5])
    print(L[cindex[15]:cindex[16]])
    assert np.allclose(L[cindex[15]:cindex[16]], IndexTh, atol=atol)
