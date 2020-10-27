"""Tests for module geoTrans"""
import numpy as np
import pytest

# Local imports
import avaframe.DFAkernel.tools as tools


def test_normalize(capfd):
    '''normalize'''
    x = np.array([1])
    y = np.array([1])
    z = np.array([1])
    xn, yn, zn = tools.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == 1
    assert xn == 1/np.sqrt(3)
    assert yn == 1/np.sqrt(3)
    assert zn == 1/np.sqrt(3)

    x = np.array([1])
    y = np.array([2])
    z = np.array([3])
    xn, yn, zn = tools.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == 1


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
    Nx, Ny, Nz = tools.getNormalVect(Z, cellsize)


    tol = 1e-8
    TestNX = np.allclose(Nx[1:nrows-1, 1:ncols-1], -2*a*X[1:nrows-1, 1:ncols-1]/np.sqrt(1+4*a*a*X[1:nrows-1, 1:ncols-1]*X[1:nrows-1, 1:ncols-1]), atol=tol)
    assert (TestNX == True)
    TestNY = np.allclose(Ny[1:nrows-1, 1:ncols-1], np.zeros(np.shape(X[1:nrows-1, 1:ncols-1])), atol=tol)
    assert (TestNY == True)
    TestNZ = np.allclose(Nz[1:nrows-1, 1:ncols-1], 1/np.sqrt(1+4*a*a*X[1:nrows-1, 1:ncols-1]*X[1:nrows-1, 1:ncols-1]), atol=tol)
    assert (TestNZ == True)
