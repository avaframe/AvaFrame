"""
    Pytest for getReleaseArea

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.in3Utils import getReleaseArea as gR

import pytest
import configparser


def test_makexyPoints():
    """ Test generation of corner Points """

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR['GENERAL'] = {'lenP': 4}

    # load function
    xyPoints = gR.makexyPoints(2, 4, 5, cfgR)

    # Test
    assert xyPoints[0,0] == 6.
    assert xyPoints[0,1] == -2.5
    assert xyPoints[1,0] == 6.
    assert xyPoints[1,1] == 2.5
    assert xyPoints[2,0] == 2
    assert xyPoints[2,1] == 2.5
    assert xyPoints[3,0] == 2
    assert xyPoints[3,1] == -2.5


def test_getCornersFP():
    """ Test generate Flat Plane Release Area """

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR['GENERAL'] = {'hr': 200.0, 'vol' : 100000, 'dh' : 1.0, 'xStart' : 100.0, 'lenP' : 4}
    cfgR['FP'] = {'xExtent': 200.0}

    # Load function
    xyPoints = gR.getCornersFP(cfgR)

    # Test
    assert xyPoints[0,0] == 300.0
    assert xyPoints[0,1] == -250.
    assert xyPoints[1,0] == 300.0
    assert xyPoints[1,1] == 250.0
    assert xyPoints[2,0] == 100.0
    assert xyPoints[2,1] == 250.0
    assert xyPoints[3,0] == 100.0
    assert xyPoints[3,1] == -250.0


def test_getCornersIP():
    """ Test generate Inclined Plane Release Area """

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR['GENERAL'] = {'hr': 200.0, 'vol' : 100000, 'dh' : 1.0, 'xStart' : 100.0, 'lenP' : 4}
    cfgT = configparser.ConfigParser()
    cfgT['TOPO'] = {'meanAlpha' : 34}

    # Load function
    xyPoints = gR.getCornersIP(cfgR, cfgT)

    # Test
    assert xyPoints[0,0] == pytest.approx(396.512194, abs=1.e-6)
    assert xyPoints[0,1] == pytest.approx(-139.798226, abs=1.e-6)
    assert xyPoints[1,0] == pytest.approx(396.512194, abs=1.e-6)
    assert xyPoints[1,1] == pytest.approx(139.798226, abs=1.e-6)
    assert xyPoints[2,0] == pytest.approx(100., abs=1.e-6)
    assert xyPoints[2,1] == pytest.approx(139.798226, abs=1.e-6)
    assert xyPoints[3,0] == pytest.approx(100., abs=1.e-6)
    assert xyPoints[3,1] == pytest.approx(-139.798226, abs=1.e-6)


def test_getCornersHS():
    """ Test generate Hockey Release Area """

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR['GENERAL'] = {'hr': 200.0, 'vol' : 100000, 'dh' : 1.0, 'xStart' : 100.0, 'lenP' : 4}
    cfgR['HS'] = {'AlphaStop' : 30}
    # Initialise input in correct format
    cfgT = configparser.ConfigParser()
    cfgT['TOPO'] = {'C': 1000, 'fLens': 2200, 'meanAlpha': 34}

    # Load function
    xyPoints = gR.getCornersHS(cfgR, cfgT)

    # Test
    assert xyPoints[0,0] == pytest.approx(848.056768, abs=1.e-6)
    assert xyPoints[0,1] == pytest.approx(-125., abs=1.e-6)
    assert xyPoints[1,0] == pytest.approx(848.056768, abs=1.e-6)
    assert xyPoints[1,1] == pytest.approx(125., abs=1.e-6)
    assert xyPoints[2,0] == pytest.approx(501.646607, abs=1.e-6)
    assert xyPoints[2,1] == pytest.approx(125., abs=1.e-6)
    assert xyPoints[3,0] == pytest.approx(501.646607, abs=1.e-6)
    assert xyPoints[3,1] == pytest.approx(-125., abs=1.e-6)


def test_correctOrigin():
    """ Test if setting a new origin works fine """

    # Initialise input in correct format
    cfgT = configparser.ConfigParser()
    cfgT['TOPO'] = {'dx': 5.0, 'xEnd': 4000.0, 'yEnd': 1000.0}
    cfgT['DEMDATA'] = {'xl': 100.0, 'yl': -1500.0}

    # Define xyPoints
    xyPointsTest = np.zeros((4,2))
    xyPointsTest[0,0] = 300.0
    xyPointsTest[0,1] = -250.
    xyPointsTest[1,0] = 300.0
    xyPointsTest[1,1] = 250.0
    xyPointsTest[2,0] = 100.0
    xyPointsTest[2,1] = 250.0
    xyPointsTest[3,0] = 100.0
    xyPointsTest[3,1] = -250.

    # Load function
    xv, yv, xyPoints = gR.correctOrigin(xyPointsTest, cfgT)

    # Test
    assert xyPoints[0,1] == -1250
    assert xyPoints[2,1] == -750
    assert xyPoints[1,0] == 400.0
    assert xyPoints[0,0] == 400.0
    assert xv[-1] == 4100.
    assert yv[0] == -1500.
