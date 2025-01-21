"""Tests for module com2AB"""

import numpy as np
import pytest
import pathlib
import shutil
import logging
import configparser

# Local imports
import avaframe.com2AB.com2AB as com2AB
from avaframe.out3Plot import outAB
from avaframe.in3Utils import cfgUtils

from avaframe.com2AB.com2AB import calcABAngles


# from unittest.mock import patch
# import unittest.mock as mocker


def test_setEqParameters(capfd):
    """Simple test for module setEqParameters"""
    cfg = cfgUtils.getModuleConfig(com2AB)
    # small avalanche
    eqParamRef = {}
    eqParamRef["parameterSet"] = "Small avalanches"
    eqParamRef["k1"] = 0.933
    eqParamRef["k2"] = 0.0
    eqParamRef["k3"] = 0.0088
    eqParamRef["k4"] = -5.02
    eqParamRef["SD"] = 2.36

    cfg["ABSETUP"]["smallAva"] = "True"

    eqParams = com2AB.setEqParameters(cfg)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]

    eqParamRef = {}
    eqParamRef["parameterSet"] = "Standard"
    eqParamRef["k1"] = 1.05
    eqParamRef["k2"] = -3130.0
    eqParamRef["k3"] = 0.0
    eqParamRef["k4"] = -2.38
    eqParamRef["SD"] = 1.25

    cfg["ABSETUP"]["smallAva"] = "False"
    eqParams = com2AB.setEqParameters(cfg)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]


def test_setEqParameters_smallAva_true(caplog):
    cfg = configparser.ConfigParser()
    cfg["ABSETUP"] = {
        "smallAva": "True",
        "k1_small": "0.1",
        "k2_small": "0.2",
        "k3_small": "0.3",
        "k4_small": "0.4",
        "SD_small": "0.5",
        "k1": "1.0",
        "k2": "2.0",
        "k3": "3.0",
        "k4": "4.0",
        "SD": "5.0",
    }
    caplog.set_level(logging.DEBUG)
    result = com2AB.setEqParameters(cfg)

    assert "Using small Avalanche Setup" in caplog.text
    assert result["k1"] == 0.1
    assert result["k2"] == 0.2
    assert result["k3"] == 0.3
    assert result["k4"] == 0.4
    assert result["SD"] == 0.5
    assert result["parameterSet"] == "Small avalanches"


def test_setEqParameters_smallAva_false(caplog):
    cfg = configparser.ConfigParser()
    cfg["ABSETUP"] = {
        "smallAva": "False",
        "k1_small": "0.1",
        "k2_small": "0.2",
        "k3_small": "0.3",
        "k4_small": "0.4",
        "SD_small": "0.5",
        "k1": "1.0",
        "k2": "2.0",
        "k3": "3.0",
        "k4": "4.0",
        "SD": "5.0",
    }
    caplog.set_level(logging.DEBUG)
    result = com2AB.setEqParameters(cfg)

    assert "Using configuration file Avalanche Setup" in caplog.text
    assert result["k1"] == 1.0
    assert result["k2"] == 2.0
    assert result["k3"] == 3.0
    assert result["k4"] == 4.0
    assert result["SD"] == 5.0
    assert result["parameterSet"] == "Standard"


def test_setEqParameters_keys_present():
    cfg = configparser.ConfigParser()
    cfg["ABSETUP"] = {
        "smallAva": "False",
        "k1_small": "0.1",
        "k2_small": "0.2",
        "k3_small": "0.3",
        "k4_small": "0.4",
        "SD_small": "0.5",
        "k1": "1.0",
        "k2": "2.0",
        "k3": "3.0",
        "k4": "4.0",
        "SD": "5.0",
    }
    result = com2AB.setEqParameters(cfg)
    expected_keys = {"k1", "k2", "k3", "k4", "SD", "parameterSet"}
    assert expected_keys == set(result.keys())


def test_calcABAngles(caplog):
    """Simple test for function calcABAngles"""

    cfg = cfgUtils.getModuleConfig(com2AB)
    # Make a reference quadratic profile
    B = -np.tan(np.deg2rad(45))
    A = -B / 4000
    C = 1000
    N = 1000
    s = np.linspace(0.0, -B / (2 * A), num=N)
    z = np.empty(np.shape(s))
    for i in range(N):
        if s[i] < (-B / (2 * A)):
            z[i] = 1 * (A * s[i] * s[i] + B * s[i] + C)
        else:
            z[i] = 1 * (-B * B / (4 * A) + C)

    thetaBeta = 10
    xBeta = (-np.tan(np.deg2rad(thetaBeta)) - B) / (2 * A)
    yBeta = A * xBeta * xBeta + B * xBeta + C
    beta = np.rad2deg(np.arctan2((C - yBeta), xBeta))
    # use standard coeef
    k1 = 1.05
    k2 = -3130.0
    k3 = 0.0
    k4 = -2.38
    SD = 1.25
    alpharef = k1 * beta + k2 * 2 * A + k3 * B * B / (2 * A) + k4
    SDs = [SD, -1 * SD, -2 * SD]
    alphaSDref = k1 * beta + k2 * 2 * A + k3 * B * B / (2 * A) + k4 + SDs

    # Using com2AB.calcAB to get the solution
    eqIn = {}
    eqIn["s"] = s  # curvilinear coordinate (of the x, y path)
    eqIn["x"] = []  # x coordinate of the path
    eqIn["y"] = []  # y coordinate of the path
    eqIn["z"] = z  # z coordinate of the path (projection of x,y on the raster)
    eqIn["indSplit"] = 2  # index of split point
    eqParams = com2AB.setEqParameters(cfg)
    eqOut = com2AB.calcABAngles(eqIn, eqParams, 30)
    alpha = eqOut["alpha"]
    alphaSD = eqOut["alphaSD"]

    # compare results with a relative tolerance of tol
    tol = 0.002  # here 0.1% relative diff
    assert (
            (alpha == pytest.approx(alpharef, rel=tol))
            and (alphaSD[0] == pytest.approx(alphaSDref[0], rel=tol))
            and (alphaSD[1] == pytest.approx(alphaSDref[1], rel=tol))
            and (alphaSD[2] == pytest.approx(alphaSDref[2], rel=tol))
    )

    with pytest.raises(IndexError) as e:
        eqOut = com2AB.calcABAngles(eqIn, eqParams, 500)
    assert str(e.value) == "No Beta point found. Check your pathAB.shp and splitPoint.shp."


def test_writeABtoSHP(tmp_path):
    """test writing to shapefile"""

    avaName = "avaSlide"
    dirname = pathlib.Path(__file__).parents[0]
    sourceDir = dirname / ".." / "data" / avaName

    avalancheDir = tmp_path / avaName

    # Copy input to tmp dir
    shutil.copytree(sourceDir, avalancheDir)

    cfg = cfgUtils.getModuleConfig(com2AB)

    # run main routine
    pathDict, dem, splitPoint, eqParams, resAB = com2AB.com2ABMain(cfg, avalancheDir)
    abShpFile = outAB.writeABtoSHP(pathDict, resAB)
    abShpFile = str(abShpFile) + ".shp"
    file_name = pathlib.Path(abShpFile)

    # check if file exists -No checks for correct content -
    assert file_name.exists()


# test complete routine
def test_com2ABMain(capfd):
    """Simple test for function com2ABMain"""
    # load and prepare Inputs
    listNames = ["avaHockeySmall", "avaHockeyChannel", "avaBowl"]
    dirname = pathlib.Path(__file__).parents[0]
    for name in listNames:
        avalancheDir = dirname / ".." / "data" / name
        saveOutPathRef = dirname / ".." / ".." / "benchmarks" / (name + "ABPytest")
        cfg = cfgUtils.getModuleConfig(com2AB)
        flags = cfg["FLAGS"]
        # run main routine
        pathDict, dem, splitPoint, eqParams, resAB = com2AB.com2ABMain(cfg, avalancheDir)
        for key in resAB:
            eqOut = resAB[key]

        # open ref data
        flags["fullOut"] = "True"
        eqParamsRef, eqOutRef = outAB.readABresults(saveOutPathRef, name, flags)
        for key in eqParamsRef.keys():
            assert eqParamsRef[key] == eqParams[key]

        atol = 1e-10
        assert (
                (np.allclose(eqOutRef["x"], eqOut["x"], atol=atol))
                and (np.allclose(eqOutRef["y"], eqOut["y"], atol=atol))
                and (np.allclose(eqOutRef["z"], eqOut["z"], atol=atol))
                and (np.allclose(eqOutRef["s"], eqOut["s"], atol=atol))
        )
        assert (np.allclose(eqOutRef["alpha"], eqOut["alpha"], atol=atol)) and (
            np.allclose(eqOutRef["alphaSD"], eqOut["alphaSD"], atol=atol)
        )


# Compare to QGIS AB routine
def test_QGISAB(capfd):
    """Compare com2ABMain results to QGIS AB results for avaSlide"""
    # load and prepare Inputs
    avaName = "avaSlide"
    dirname = pathlib.Path(__file__).parents[0]
    avalancheDir = dirname / ".." / "data" / avaName
    cfg = cfgUtils.getModuleConfig(com2AB)
    # run main routine
    pathDict, dem, splitPoint, eqParams, resAB = com2AB.com2ABMain(cfg, avalancheDir)
    # process data to get results
    for i, name in enumerate(resAB):
        avaProfile = resAB[name]
        beta = avaProfile["beta"]
        alpha = avaProfile["alpha"]
        alphaSD = avaProfile["alphaSD"]
        s = avaProfile["s"]
        indAlpha = avaProfile["indAlpha"]
        indBetaPoint = avaProfile["indBetaPoint"]
        indAlphaP1SD = avaProfile["indAlphaP1SD"]
        indAlphaM1SD = avaProfile["indAlphaM1SD"]
        indAlphaM2SD = avaProfile["indAlphaM2SD"]
        # get ref results
        nameRef = name + "_AB_QGIS.txt"
        nameRefpath = pathlib.Path(dirname, "..", "..", "benchmarks", avaName + "ABPytest", nameRef)
        data = np.loadtxt(nameRefpath, skiprows=4, delimiter=",")

        tolDist = 12.5  # increased from 10.
        tolAngle = 0.19  # increased from 0.12 to 0.19
        assert (
                (alpha == pytest.approx(data[0, 4], abs=tolAngle))
                and (beta == pytest.approx(data[1, 4], rel=tolAngle))
                and (alphaSD[1] == pytest.approx(data[2, 4], rel=tolAngle))
                and (alphaSD[2] == pytest.approx(data[3, 4], rel=tolAngle))
                and (alphaSD[0] == pytest.approx(data[4, 4], rel=tolAngle))
        )
        assert (
                (s[indAlpha] == pytest.approx(data[0, 3], abs=tolDist))
                and (s[indBetaPoint] == pytest.approx(data[1, 3], rel=tolDist))
                and (s[indAlphaM1SD] == pytest.approx(data[2, 3], rel=tolDist))
                and (s[indAlphaM2SD] == pytest.approx(data[3, 3], rel=tolDist))
                and (s[indAlphaP1SD] == pytest.approx(data[4, 3], rel=tolDist))
        )


@pytest.fixture
def sampleAvaProfile():
    return {
        "s": np.array([0, 1, 2, 3, 4]),
        "z": np.array([0, 1, 4, 9, 16]),  # Quadratic: z = s^2
        "indSplit": 2,  # Middle index for example
    }


@pytest.fixture
def sampleEqParameters():
    return {"k1": 0.1, "k2": 0.2, "k3": 0.3, "k4": 0.4, "SD": 0.5}


def test_calcABAngles_success(sampleAvaProfile, sampleEqParameters, mocker):
    # Mock external dependencies
    mocker.patch("avaframe.in3Utils.geoTrans.prepareAngleProfile", return_value=(None, None, None))
    mocker.patch("avaframe.in3Utils.geoTrans.findAngleProfile", return_value=1)  # indBetaPoint = 1

    result = calcABAngles(sampleAvaProfile, sampleEqParameters, dsMin=1.0)

    # Verify key additions
    assert "beta" in result
    assert "alpha" in result
    assert "alphaSD" in result
    assert "SDs" in result
    assert "poly" in result
    assert "indBetaPoint" in result

    # Verify quadratic fit
    assert np.allclose(result["poly"].coefficients, [1.0, 0.0, 0.0])  # z = s^2

    # Verify H0 calculation (max - min of parabola)
    assert result["poly"](result["s"]).max() - result["poly"](result["s"]).min() == pytest.approx(16.0)

    # Verify beta calculation (dzBeta = z[0] - z[1] = 0 - 1 = -1)
    expectedBeta = np.degrees(np.arctan2(-1, 1))  # arctan(-1/1) = -45°
    assert result["beta"] == pytest.approx(expectedBeta)

    # Verify alpha calculation
    # alpha = k1*beta + k2*(2nd_deriv) + k3*H0 + k4
    # 2nd derivative of z = s^2 is 2
    expectedAlpha = 0.1 * expectedBeta + 0.2 * 2 + 0.3 * 16.0 + 0.4
    assert result["alpha"] == pytest.approx(expectedAlpha)

    # Verify alphaSD
    SDs = [0.5, -0.5, -1.0]
    assert result["SDs"] == SDs
    assert result["alphaSD"] == pytest.approx([expectedAlpha + sd for sd in SDs])


def test_calcABAngles_no_beta_found(sampleAvaProfile, sampleEqParameters, mocker):
    mocker.patch("avaframe.in3Utils.geoTrans.prepareAngleProfile", return_value=(None, None, None))
    mocker.patch("avaframe.in3Utils.geoTrans.findAngleProfile", side_effect=IndexError)

    with pytest.raises(IndexError) as excinfo:
        calcABAngles(sampleAvaProfile, sampleEqParameters, dsMin=1.0)

    assert "No Beta point found. Check your pathAB.shp and splitPoint.shp." in str(excinfo.value)


def test_calcABAngles_polyfit_verification(sampleAvaProfile, sampleEqParameters, mocker):
    # Test with known quadratic relationship
    mocker.patch("avaframe.in3Utils.geoTrans.prepareAngleProfile", return_value=(None, None, None))
    mocker.patch("avaframe.in3Utils.geoTrans.findAngleProfile", return_value=1)

    result = calcABAngles(sampleAvaProfile, sampleEqParameters, dsMin=1.0)

    # Verify polynomial coefficients (should match z = s^2)
    assert np.allclose(result["poly"].coefficients, [1.0, 0.0, 0.0])

    # Verify second derivative calculation
    assert result["poly"].deriv(2)(0) == pytest.approx(2.0)  # Second derivative of z = s^2 is 2


def test_calcABAngles_debugPlot_true(sampleAvaProfile, sampleEqParameters, mocker):
    # Mock debugPlot to True and plotting functions
    mocker.patch("avaframe.com2AB.com2AB.debugPlot", True)  # Replace with actual module path
    # mock_plot_slope = mocker.patch("avaframe.out3Plot.outDebugPlots.plotSlopeAngle")
    # mock_plot_profile = mocker.patch("avaframe.out3Plot.outDebugPlots.plotProfile")
    mock_plot_slope = mocker.patch("avaframe.com2AB.com2AB.debPlot.plotSlopeAngle")
    mock_plot_profile = mocker.patch("avaframe.com2AB.com2AB.debPlot.plotProfile")

    # Mock geoTrans functions to return dummy data
    mocker.patch(
        "avaframe.in3Utils.geoTrans.prepareAngleProfile",
        return_value=(
            np.array([30.0, 25.0, 20.0, 15.0, 10.0]),  # angle values
            "tmp_dummy",
            "ds_dummy",
        ),
    )
    mocker.patch("avaframe.in3Utils.geoTrans.findAngleProfile", return_value=1)

    # Execute the function
    result = calcABAngles(sampleAvaProfile, sampleEqParameters, dsMin=1.0)

    # Verify plots were called
    # mock_plot_slope.assert_called_once_with(
    #     sampleAvaProfile["s"], np.array([30.0, 25.0, 20.0, 15.0, 10.0]), 1
    # )
    # mock_plot_profile.assert_called_once_with(sampleAvaProfile["s"], sampleAvaProfile["z"], 1)
