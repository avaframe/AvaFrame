"""
    Pytest for generateTopo

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.in3Utils import generateTopo as gT
import pytest
import configparser

def test_getParabolaParams():
    """ Test if A and B is computed correctly """

    # input parameters
    C = 1000.
    fLens = 2250.
    meanAlpha = 0.

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'C': C, 'fLens': fLens, 'meanAlpha': meanAlpha}

    # Get parabola Parameters
    [A, B, fLen] = gT.getParabolaParams(cfg)

    ATest1 = 0.00019753086419753085
    BTest1 = -0.8888888888888888
    fLenTest1 = 2250.

    # Test if computations work with given mean_alpha
    # input parameters
    C = 1000.
    fLens = 2250.
    meanAlpha = 34.

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'C': C, 'fLens': fLens, 'meanAlpha': meanAlpha}

    # Get parabola Parameters
    [A2, B2, fLen2] = gT.getParabolaParams(cfg)

    ATest2 = 0.0004549617392929703
    BTest2 = -1.3490170336848535
    fLenTest2 = 1482.56096851274

    assert A - ATest1 == pytest.approx(0.0, abs=1.e-12)
    assert B - BTest1 == pytest.approx(0.0, abs=1.e-12)
    assert(fLen == fLenTest1)
    assert A2 - ATest2 == pytest.approx(0.0, abs=1.e-12)
    assert B2 - BTest2 == pytest.approx(0.0, abs=1.e-12)
    assert(fLen2 == fLenTest2)


def test_getGridDefs():
    """ Test generation of grid extents """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 5., 'xEnd': 10., 'yEnd': 100.}

    dx, xEnd, yEnd = gT.getGridDefs(cfg)

    assert xEnd == pytest.approx(15.0, abs=1.e-12)
    assert yEnd == pytest.approx(105.0, abs=1.e-12)
    assert yEnd != pytest.approx(100.0, abs=1.e-12)

def test_computeCoordGrid():
    """ Test if vectors that define the domain and raster sizes are
        generated correctly
    """

    xvTest = np.array([0, 1, 2, 3 ,4])
    yvTest = np.array([0, 1, 2])
    nColsTest = 5
    nRowsTest = 3

    xv, yv, zv, x, y, nRows, nCols = gT.computeCoordGrid(1., 5., 3.)

    assert nCols == 5
    assert nRows == 3

#
def test_flatplane():
    """ Test flat plane generation """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'zElev': 0.0, 'dx': 5., 'xEnd': 10., 'yEnd': 100.}

    # Call function to be tested
    x, y, zv = gT.flatplane(cfg)
    zMax = np.amax(zv)
    zMin = np.amin(zv)
    zMean = np.mean(zv)

    assert zMax == pytest.approx(0.0, abs=1.e-12)
    assert zMin == pytest.approx(0.0, abs=1.e-12)
    assert(zMean == 0.0)


def test_inclinedplane():
    """ Test inclined plane generation """

    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 5, 'xEnd': 5000., 'yEnd': 2000., 'z0': 2200.,
               'meanAlpha': 34., 'channel': 'False', 'topoconst': 'True'}
    cfg['CHANNELS'] ={'c_ff': 250, 'c_radius': 100}
    nrows = 301
    ncols = 1001

    # Call function to be tested
    x, y, z = gT.inclinedplane(cfg)

    # Load reference Solution
    zSol = np.loadtxt(os.path.join(os.getcwd(), 'avaframe/tests/refData/myDEM_IP_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    assert (testRes == True)


def test_hockeysmooth():
    """ Test hockey smooth generation """

    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 5, 'xEnd': 5000., 'yEnd': 2000., 'rCirc': 200.,
               'z0' : 2200., 'meanAlpha': 34., 'channel': 'True',
               'narrowing' : 'True', 'topoconst': 'True'}
    cfg['CHANNELS'] ={'c_ff': 250, 'c_radius': 100, 'c_init' : 250, 'c_mustart' : 0.2, 'c_muendFP' : 0.86}
    nrows = 301
    ncols = 1001

    # Call function to be tested
    x, y, z = gT.hockeysmooth(cfg)

    # Load reference Solution
    zSol = np.loadtxt(os.path.join(os.getcwd(), 'avaframe/tests/refData/myDEM_HS2_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=5.e-6)

    assert (testRes == True)


def test_hockey():
    """ Test hockey generation """

    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'C' : 1000, 'meanAlpha' : 34, 'fLens' : 0, 'dx': 5, 'xEnd': 5000.,
                  'yEnd': 2000., 'channel': 'False', 'topoconst': 'True'}
    cfg['CHANNELS'] ={'c_ff': 250, 'c_radius': 100, 'c_mustart' : 0.2, 'c_muend' : 0.6, 'c_init' : 250}
    nrows = 301
    ncols = 1001

    # Call function to be tested
    x, y, z = gT.hockey(cfg)

    # Load reference Solution
    zSol = np.loadtxt(os.path.join(os.getcwd(), 'avaframe/tests/refData/myDEM_HS_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    assert (testRes == True)


def test_helix():
    """ Test helix generation """

    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'C' : 1000, 'meanAlpha' : 34, 'rHelix' : 1250, 'fLens' : 0, 'dx': 5, 'xEnd': 5000.,
                  'yEnd': 2000.,  'channel': 'True', 'narrowing' : 'True', 'topoconst': 'True'}
    cfg['CHANNELS'] ={'c_ff': 250, 'c_radius': 100, 'c_mustart' : 0.2, 'c_muend' : 0.6, 'c_init' : 250}
    nrows = 301
    ncols = 1001

    # Call function to be tested
    x, y, z = gT.helix(cfg)

    # Load reference Solution
    zSol = np.loadtxt(os.path.join(os.getcwd(), 'avaframe/tests/refData/myDEM_HX_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    assert (testRes == True)
