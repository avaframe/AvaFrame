"""Tests for module com1DFAtools"""
import numpy as np
import matplotlib.pyplot as plt
import pytest

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls


def test_getNormalMesh(capfd):
    '''projectOnRaster'''
    a = 2
    b = 1
    cellsize = 1
    m = 10
    n = 15
    x = np.linspace(0, m-1, m)
    y = np.linspace(0, n-1, n)
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y
    header = {}
    header['ncols'] = m
    header['nrows'] = n
    header['cellsize'] = cellsize
    dem = {}
    dem['header'] = header
    Z1 = a * X * X + b * Y * Y
    for num in [4, 6, 8]:
        dem['rasterData'] = Z
        Nx, Ny, Nz = DFAtls.getNormalMesh(dem, num)
        Nx, Ny, Nz = DFAtls.normalize(Nx, Ny, Nz)
        print(Nx)
        print((-a*np.ones(np.shape(Y)) / np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1])
        print(Ny)
        print((-b*np.ones(np.shape(Y)) / np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1])
        print(Nz)
        print((np.ones(np.shape(Y)) / np.sqrt(1 + a*a + b*b))[1:n-1, 1:m-1])

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

        dem['rasterData'] = Z1
        Nx, Ny, Nz = DFAtls.getNormalMesh(dem, num)
        Nx, Ny, Nz = DFAtls.normalize(Nx, Ny, Nz)

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
    '''projectOnRaster'''
    a = 0.1
    b = 0.2
    csz = 1
    m = 15
    n = 10
    x = np.linspace(0, m-1, m)
    y = np.linspace(0, n-1, n)
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y
    Z1 = a * X * X + b * Y * Y
    header = {}
    header['ncols'] = m
    header['nrows'] = n
    header['cellsize'] = csz
    dem = {}
    dem['header'] = header
    dem['rasterData'] = Z
    Nx, Ny, Nz = DFAtls.getNormalMesh(dem, 4)
    Nx, Ny, Nz = DFAtls.normalize(Nx, Ny, Nz)
    Area = DFAtls.getAreaMesh(Nx, Ny, Nz, csz, 4)
    print(np.sqrt((1+a*a+b*b)))
    print(Area)
    atol = 1e-10
    TestArea = np.allclose(Area[1:n-1, 1:m-1], np.sqrt((1+a*a+b*b)) *
                           np.ones(np.shape(Y[1:n-1, 1:m-1])), atol=atol)
    assert TestArea
