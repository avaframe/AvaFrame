"""
    Pytest for getReleaseArea

    This file is part of Avaframe.

 """

#  Load modules
from avaframe.in3Utils import getReleaseArea as gR
from avaframe.in3Utils.getReleaseArea import getReleaseArea, writeReleaseArea
import avaframe.in3Utils.fileHandlerUtils as fU

import pytest
import configparser
import numpy as np
import pytest
from pathlib import Path
import pathlib
import geopandas as gpd


def test_makexyPoints():
    """Test generation of corner Points"""

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR["GENERAL"] = {"lenP": 4}

    # load function
    xyPoints = gR.makexyPoints(2, 4, 5, cfgR)

    # Test
    assert xyPoints[0, 0] == 6.0
    assert xyPoints[0, 1] == -2.5
    assert xyPoints[1, 0] == 6.0
    assert xyPoints[1, 1] == 2.5
    assert xyPoints[2, 0] == 2
    assert xyPoints[2, 1] == 2.5
    assert xyPoints[3, 0] == 2
    assert xyPoints[3, 1] == -2.5


def test_getCornersFP():
    """Test generate Flat Plane Release Area"""

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR["GENERAL"] = {"hr": 200.0, "vol": 100000, "dh": 1.0, "xStart": 100.0, "lenP": 4}
    cfgR["FP"] = {"xExtent": 200.0}

    # Load function
    xyPoints = gR.getCornersFP(cfgR)

    # Test
    assert xyPoints[0, 0] == 300.0
    assert xyPoints[0, 1] == -250.0
    assert xyPoints[1, 0] == 300.0
    assert xyPoints[1, 1] == 250.0
    assert xyPoints[2, 0] == 100.0
    assert xyPoints[2, 1] == 250.0
    assert xyPoints[3, 0] == 100.0
    assert xyPoints[3, 1] == -250.0


def test_getCornersIP():
    """Test generate Inclined Plane Release Area"""

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR["GENERAL"] = {"hr": 200.0, "vol": 100000, "dh": 1.0, "xStart": 100.0, "lenP": 4}
    cfgT = configparser.ConfigParser()
    cfgT["TOPO"] = {"meanAlpha": 34}

    # Load function
    xyPoints = gR.getCornersIP(cfgR, cfgT)

    # Test
    assert xyPoints[0, 0] == pytest.approx(396.512194, abs=1.0e-6)
    assert xyPoints[0, 1] == pytest.approx(-139.798226, abs=1.0e-6)
    assert xyPoints[1, 0] == pytest.approx(396.512194, abs=1.0e-6)
    assert xyPoints[1, 1] == pytest.approx(139.798226, abs=1.0e-6)
    assert xyPoints[2, 0] == pytest.approx(100.0, abs=1.0e-6)
    assert xyPoints[2, 1] == pytest.approx(139.798226, abs=1.0e-6)
    assert xyPoints[3, 0] == pytest.approx(100.0, abs=1.0e-6)
    assert xyPoints[3, 1] == pytest.approx(-139.798226, abs=1.0e-6)


def test_getCornersHS():
    """Test generate Hockey Release Area"""

    # Initialise input in correct format
    cfgR = configparser.ConfigParser()
    cfgR["GENERAL"] = {"hr": 200.0, "vol": 100000, "dh": 1.0, "xStart": 100.0, "lenP": 4}
    cfgR["HS"] = {"AlphaStop": 30}
    # Initialise input in correct format
    cfgT = configparser.ConfigParser()
    cfgT["TOPO"] = {"C": 1000, "fLens": 2200, "meanAlpha": 34}

    # Load function
    xyPoints = gR.getCornersHS(cfgR, cfgT)

    # Test
    assert xyPoints[0, 0] == pytest.approx(848.056768, abs=1.0e-6)
    assert xyPoints[0, 1] == pytest.approx(-125.0, abs=1.0e-6)
    assert xyPoints[1, 0] == pytest.approx(848.056768, abs=1.0e-6)
    assert xyPoints[1, 1] == pytest.approx(125.0, abs=1.0e-6)
    assert xyPoints[2, 0] == pytest.approx(501.646607, abs=1.0e-6)
    assert xyPoints[2, 1] == pytest.approx(125.0, abs=1.0e-6)
    assert xyPoints[3, 0] == pytest.approx(501.646607, abs=1.0e-6)
    assert xyPoints[3, 1] == pytest.approx(-125.0, abs=1.0e-6)


def test_correctOrigin():
    """Test if setting a new origin works fine"""

    # Initialise input in correct format
    cfgT = configparser.ConfigParser()
    cfgT["TOPO"] = {"dx": 5.0, "xEnd": 4000.0, "yEnd": 1000.0}
    cfgT["DEMDATA"] = {"xl": 100.0, "yl": -1500.0}

    # Define xyPoints
    xyPointsTest = np.zeros((4, 2))
    xyPointsTest[0, 0] = 300.0
    xyPointsTest[0, 1] = -250.0
    xyPointsTest[1, 0] = 300.0
    xyPointsTest[1, 1] = 250.0
    xyPointsTest[2, 0] = 100.0
    xyPointsTest[2, 1] = 250.0
    xyPointsTest[3, 0] = 100.0
    xyPointsTest[3, 1] = -250.0

    # Load function
    xv, yv, xyPoints = gR.correctOrigin(xyPointsTest, cfgT)

    # Test
    assert xyPoints[0, 1] == -1250
    assert xyPoints[2, 1] == -750
    assert xyPoints[1, 0] == 400.0
    assert xyPoints[0, 0] == 400.0
    assert xv[-1] == 4100.0
    assert yv[0] == -1500.0


# tests from here on are based on AI


@pytest.fixture
def cfgR():
    cfg = configparser.ConfigParser()
    cfg["FILE"] = {"relNo": "1", "relName": "TestRelease"}
    cfg["GENERAL"] = {"dh": "5.0", "outputtxt": "True"}
    return cfg


@pytest.fixture
def cfgT():
    cfg = configparser.ConfigParser()
    cfg["TOPO"] = {"demType": "FP"}
    return cfg


def test_writeReleaseArea_createsNxyzFile(tmp_path, cfgR, mocker):
    xyPoints = np.array([[0, 1], [2, 3], [4, 5], [6, 7]])
    demType = "FP"
    outDir = tmp_path
    mocker.patch("avaframe.in3Utils.getReleaseArea.log.info")

    writeReleaseArea(xyPoints, demType, cfgR, outDir)

    expectedFile = outDir / "release1FP.nxyz"
    assert expectedFile.exists()

    with open(expectedFile, "r") as f:
        lines = f.readlines()
        assert lines[0] == "name=TestRelease\n"
        assert lines[1] == "d0=5.00\n"
        assert lines[2] == "rho=None\n"
        assert lines[3] == "5\n"
        dataLines = [line.strip() for line in lines[4:]]
        assert dataLines == [
            "0.000000 1.000000 0.000000",
            "2.000000 3.000000 0.000000",
            "4.000000 5.000000 0.000000",
            "6.000000 7.000000 0.000000",
            "0.000000 1.000000 0.000000",
        ]


def test_writeReleaseArea_copiesTxtFile(tmp_path, cfgR, mocker):
    mocker.patch("avaframe.in3Utils.getReleaseArea.log.info")
    mockCopy = mocker.patch("shutil.copyfile")

    xyPoints = np.array([[0, 1], [2, 3], [4, 5], [6, 7]])
    writeReleaseArea(xyPoints, "FP", cfgR, tmp_path)

    src = tmp_path / "release1FP.nxyz"
    dst = tmp_path / "release1FP.txt"
    mockCopy.assert_called_once_with(src, dst)

def test_writeReleaseArea(tmp_path):
    """ test writing a shp file for the created polygon """

    # setup required inputs
    cfgR = configparser.ConfigParser()
    cfgR['FILE'] = {'relNo': 1, 'relName': 'releaseTest'}
    cfgR['GENERAL'] = {'dh': 1.5, 'outputtxt': True}
    outDir  = pathlib.Path(tmp_path, 'testDir')
    fU.makeADir(outDir)
    demType = 'HS'
    xyPoints = np.array([[2.0, 1.0], [5.0, 2.0], [4.0, 7.0], [1.0, 5.0]])

    gR.writeReleaseArea(xyPoints, demType, cfgR, outDir)

    # read file
    gs = gpd.read_file((outDir / 'release1HS.shp'))

    print('gs', gs)
    gsCoors = gs.get_coordinates().to_numpy()
    xyPointsWithEnd = np.array([[2.0, 1.0], [5.0, 2.0], [4.0, 7.0], [1.0, 5.0], [2.0, 1.0]])


    print('gsCoors', type(gsCoors), gsCoors.shape, 'type', xyPoints.shape, type(xyPoints))
    assert np.array_equal(xyPointsWithEnd, np.asarray(gs.get_coordinates().to_numpy()))


def test_getReleaseArea_logsMissingDirectory(tmp_path, cfgT, cfgR, mocker):
    mockLog = mocker.patch("avaframe.in3Utils.getReleaseArea.log.error")
    mocker.patch("avaframe.in3Utils.getReleaseArea.getCornersFP")
    mocker.patch("avaframe.in3Utils.getReleaseArea.correctOrigin", return_value=("xv", "yv", "xyPoints"))
    mocker.patch("avaframe.in3Utils.getReleaseArea.writeReleaseArea")

    # Do not create the required directory
    getReleaseArea(cfgT, cfgR, tmp_path)

    mockLog.assert_called_once_with(
        "Required folder structure: NameOfAvalanche/Inputs missing!                     Run runInitializeProject first!"
    )


def test_getReleaseArea_demTypeHandling(mocker, cfgR):
    cfgT = configparser.ConfigParser()
    cfgT["TOPO"] = {"demType": "IP"}
    mockGetCorners = mocker.patch("avaframe.in3Utils.getReleaseArea.getCornersIP")
    mocker.patch("avaframe.in3Utils.getReleaseArea.correctOrigin", return_value=("xv", "yv", "xyPoints"))
    mocker.patch("avaframe.in3Utils.getReleaseArea.writeReleaseArea")

    getReleaseArea(cfgT, cfgR, Path("/dummy"))
    mockGetCorners.assert_called_once_with(cfgR, cfgT)


def test_getReleaseArea_invalidDemTypeLogsWarning(mocker):
    cfgT = configparser.ConfigParser()
    cfgT["TOPO"] = {"demType": "BL"}
    cfgR = configparser.ConfigParser()
    mockLog = mocker.patch("avaframe.in3Utils.getReleaseArea.log.warning")

    with pytest.raises(UnboundLocalError):
        getReleaseArea(cfgT, cfgR, Path("/dummy"))
    mockLog.assert_called_once_with("no release area available for demType: BL")
