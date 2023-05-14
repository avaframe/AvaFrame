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
    """ Test if parabola parameters are computed correctly """

    # Input parameters
    C = 1000.
    fLens = 2250.
    meanAlpha = 0.

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'C': C, 'fLens': fLens, 'meanAlpha': meanAlpha}

    # Get parabola Parameters
    [A, B, fLen] = gT.getParabolaParams(cfg)

    # Expected results
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

    # Expected results
    ATest2 = 0.0004549617392929703
    BTest2 = -1.3490170336848535
    fLenTest2 = 1482.56096851274

    # Test
    assert A - ATest1 == pytest.approx(0.0, abs=1.e-12)
    assert B - BTest1 == pytest.approx(0.0, abs=1.e-12)
    assert fLen == pytest.approx(fLenTest1, abs=1.e-10)
    assert A2 - ATest2 == pytest.approx(0.0, abs=1.e-12)
    assert B2 - BTest2 == pytest.approx(0.0, abs=1.e-12)
    assert fLen2 == pytest.approx(fLenTest2, abs=1.e-10)


def test_getGridDefs():
    """ Test generation of grid extents """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 5., 'xEnd': 10., 'yEnd': 100.}

    dx, xEnd, yEnd = gT.getGridDefs(cfg)

    # second Test
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 10., 'xEnd': 10., 'yEnd': 100.}
    dx2, xEnd2, yEnd2 = gT.getGridDefs(cfg)

    # Test
    assert xEnd == pytest.approx(10.0, abs=1.e-12)
    assert yEnd == pytest.approx(100.0, abs=1.e-12)
    assert xEnd2 != pytest.approx(15.0, abs=1.e-12)
    assert yEnd2 != pytest.approx(105.0, abs=1.e-12)


def test_computeCoordGrid():
    """ Test if vectors that define the domain and raster sizes are
        generated correctly
    """

    # Load function with given inputs
    xv, yv, zv, x, y, nRows, nCols = gT.computeCoordGrid(1., 5., 3.)

    # Define correct results
    xvTest = np.array([0, 1, 2, 3, 4, 5])
    yvTest = np.array([-1.5, -0.5, 0.5, 1.5])

    # Test
    assert nCols == 6
    assert nRows == 4
    assert xv[2] == xvTest[2]
    assert yv[-1] == yvTest[-1]


def test_flatplane():
    """ Test flat plane generation """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'zElev': 50.0, 'dx': 5., 'xEnd': 10., 'yEnd': 100.}

    # Call function to be tested
    x, y, zv = gT.flatplane(cfg)
    zMax = np.amax(zv)
    zMin = np.amin(zv)
    zMean = np.mean(zv)

    # Test
    assert zMax == pytest.approx(50.0, abs=1.e-12)
    assert zMin == pytest.approx(50.0, abs=1.e-12)
    assert zMean == 50.0


def test_inclinedplane():
    """ Test inclined plane generation """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 5, 'xEnd': 5000., 'yEnd': 2000., 'z0': 2200.,
                   'meanAlpha': 34., 'channel': 'False', 'topoconst': 'True'}
    cfg['CHANNELS'] = {'cff': 250, 'cRadius': 100}

    # Call function to be tested
    x, y, z = gT.inclinedplane(cfg)

    # Load reference solution
    dirPath = os.path.dirname(__file__)
    zSol = np.loadtxt(os.path.join(dirPath, '..', 'data', 'avaInclinedPlane',
                                   'Inputs', 'DEM_IP_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    # Test
    assert (testRes is True)

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 5, 'xEnd': 5000., 'yEnd': 1500., 'z0': 2200.,
                   'meanAlpha': 34., 'channel': 'True', 'topoconst': 'True'}
    cfg['CHANNELS'] = {'cff': 250, 'cRadius': 100}

    # Call function to be tested
    x, y, z = gT.inclinedplane(cfg)

    # Load reference solution
    dirPath = os.path.dirname(__file__)
    zSol = np.loadtxt(os.path.join(dirPath, 'data',
                                   'DEM_IP_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    # Test
    assert (testRes is True)


def test_hockey():
    """ Test hockey generation """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'dx': 5, 'xEnd': 5000., 'yEnd': 2000., 'rCirc': 200.,
                   'z0': 2200., 'meanAlpha': 34., 'channel': 'True',
                   'narrowing': 'True', 'topoAdd': 'True'}
    cfg['CHANNELS'] = {'cff': 250, 'cRadius': 100,
                       'cInit': 250, 'cMustart': 0.2, 'cMuendFP': 0.86}

    # Call function to be tested
    x, y, z = gT.hockey(cfg)

    # Load reference Solution
    dirPath = os.path.dirname(__file__)
    zSol = np.loadtxt(os.path.join(dirPath, '..', 'data', 'avaHockeyChannel',
                                   'Inputs', 'DEM_HS_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=5.e-6)

    # Test
    assert (testRes is True)


def test_parabola():
    """ Test parabola generation """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'C': 2200, 'meanAlpha': 34, 'fLens': 0, 'dx': 5, 'xEnd': 5000.,
                   'yEnd': 2000., 'channel': 'False', 'topoAdd': 'False', 'dam': 'False'}
    cfg['CHANNELS'] = {'cff': 250, 'cRadius': 100,
                       'cMustart': 0.2, 'cMuend': 0.6, 'cInit': 250}

    # Call function to be tested
    x, y, z = gT.parabola(cfg)

    # Load reference Solution
    dirPath = os.path.dirname(__file__)
    zSol = np.loadtxt(os.path.join(dirPath, '..', 'data', 'avaParabola',
                                   'Inputs', 'DEM_PF_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    # Test
    assert (testRes is True)


def test_parabolaRotation():
    """ Test parabola with rotations generation
    only checks if it runs"""

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'demType': 'TPF', 'fFlat': 500, 'C': 2200, 'meanAlpha': 34, 'fLens': 0, 'dx': 5, 'xEnd': 5000.,
                   'yEnd': 2000., 'channel': 'False', 'topoAdd': 'False', 'dam': 'False'}
    cfg['CHANNELS'] = {'cff': 250, 'cRadius': 100,
                       'cMustart': 0.2, 'cMuend': 0.6, 'cInit': 250}

    # Call function to be tested
    x, y, z = gT.parabolaRotation(cfg)

    # Load reference Solution
    dirPath = os.path.dirname(__file__)
    zSol = np.loadtxt(os.path.join(dirPath, '..', 'data', 'avaParabola',
                                   'Inputs', 'DEM_PF_Topo.asc'), skiprows=6)


def test_helix():
    """ Test helix generation """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'C': 1000, 'meanAlpha': 34, 'rHelix': 1250, 'fLens': 0, 'dx': 5, 'xEnd': 5000.,
                   'yEnd': 2000.,  'channel': 'True', 'narrowing': 'True', 'topoAdd': 'True'}
    cfg['CHANNELS'] = {'cff': 250, 'cRadius': 100,
                       'cMustart': 0.2, 'cMuend': 0.6, 'cInit': 250}

    # Call function to be tested
    x, y, z = gT.helix(cfg)

    # Load reference Solution
    dirPath = os.path.dirname(__file__)
    zSol = np.loadtxt(os.path.join(dirPath, '..', 'data', 'avaHelixChannel',
                                   'Inputs', 'DEM_HX_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    # Test
    assert (testRes is True)


def test_bowl():
    """ Test helix generation """

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['TOPO'] = {'rBowl': 2500, 'fLens': 0, 'dx': 5, 'xEnd': 5000.,
                   'yEnd': 5000.}

    # Call function to be tested
    x, y, z = gT.bowl(cfg)

    # Load reference Solution
    dirPath = os.path.dirname(__file__)
    zSol = np.loadtxt(os.path.join(dirPath, '..', 'data', 'avaBowl',
                                   'Inputs', 'DEM_BL_Topo.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(z, zSol, atol=1.e-6)

    # Test
    assert (testRes is True)
