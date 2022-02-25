"""Tests for module ana1Tests"""
import numpy as np
import math
import pathlib
import shutil
import pytest

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.ana1Tests.analysisTools as anaTools
import avaframe.ana1Tests.energyLineTest as energyLineTest
import avaframe.com1DFA.com1DFA as com1DFA


# ############# Test analysis tools ##########
# ############################################

def test_ErrorNorm(capfd):
    '''test L2Norm '''
    cellSize = 2
    cosAngle = 0.5
    localError = np.array([[0, 1, 2], [3, 4, 5]])
    normL2 = anaTools.L2Norm(localError, cellSize, cosAngle)
    atol = 1e-10
    assert normL2 == pytest.approx(np.sqrt(120), abs=atol)

    refArray = np.array([[2, 0, 1], [0, 1, 0]])
    errorL2, errorL2Rel, errorMax, errorMaxRel = anaTools.computeErrorAndNorm(localError, refArray, cellSize, cosAngle)

    assert errorL2 == pytest.approx(np.sqrt(120), abs=atol)
    assert errorL2Rel == pytest.approx(np.sqrt(60)/4, abs=atol)
    assert errorMax == np.sqrt(5)
    assert errorMaxRel == pytest.approx(np.sqrt(5/2), abs=atol)

    errorL2, errorL2Rel, errorMax, errorMaxRel = anaTools.normL2Scal(refArray, localError, cellSize, cosAngle)
    assert errorL2 == pytest.approx(np.sqrt(49*8), abs=atol)
    assert errorL2Rel == pytest.approx(np.sqrt(49/6), abs=atol)
    assert errorMax == 5
    assert errorMaxRel == 5/2

    localError = {'fx': localError, 'fy': localError, 'fz': localError}
    refArray = {'fx': refArray, 'fy': refArray, 'fz': refArray}
    errorL2, errorL2Rel, errorMax, errorMaxRel = anaTools.normL2Vect(refArray, localError, cellSize, cosAngle)
    assert errorL2 == pytest.approx(np.sqrt(3*49*8), abs=atol)
    assert errorL2Rel == pytest.approx(np.sqrt(49/6), abs=atol)
    assert errorMax == pytest.approx(np.sqrt(3)*5, abs=atol)
    assert errorMaxRel == pytest.approx(5/2, abs=atol)

    localError = np.array([[0, 1, 2], [3, 4, 5]])
    refArray = np.array([[0, 0, 0], [0, 0, 0]])
    errorL2, errorL2Rel, errorMax, errorMaxRel = anaTools.computeErrorAndNorm(localError, refArray, cellSize, cosAngle)

    assert errorL2 == pytest.approx(np.sqrt(120), abs=atol)
    assert errorL2Rel == pytest.approx(np.sqrt(120), abs=atol)
    assert errorMax == np.sqrt(5)
    assert errorMaxRel == np.sqrt(5)

    localError = np.array([[0, 0, 0], [0, 0, 0]])
    normL2 = anaTools.L2Norm(localError, cellSize, cosAngle)
    atol = 1e-6
    assert normL2 == 0


# ############# Test energy line test ########
# ############################################
def test_getIntersection(capfd):
    '''test find intersection '''
    dirname = pathlib.Path(__file__).parents[0]
    energyLineTestCfgFile = dirname / '..' / 'tests' / 'data' / 'testEnergyLine' / 'energyLineTest_com1DFACfg.ini'

    energyLineTestCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=str(energyLineTestCfgFile), modInfo=True)

    mu = 1
    avaProfileMass = {'s':np.array([0, 0.1, 1.1, 2.1, 3.1, 4.1, 5.1]), 'z':np.array([5.1, 0, -0.5, -1, -1.5, -2, -2.5])}
    slopeExt, sIntersection, zIntersection, coefExt = energyLineTest.getAlphaProfileIntersection(energyLineTestCfg, avaProfileMass, mu)
    print(sIntersection)
    print(zIntersection)
    print(slopeExt)
    print(coefExt)
    atol = 1e-10
    assert sIntersection == pytest.approx(10.1, abs=atol)
    assert zIntersection == pytest.approx(-5, abs=atol)
    assert slopeExt == pytest.approx(-0.5, abs=atol)
    assert coefExt == pytest.approx(1, abs=atol)

    avaProfileMass = {'s':np.array([0, 0.1, 1.1, 2.1, 3.1, 4.1, 5.1]), 'z':np.array([2.6, 0, -0.5, -1, -1.5, -2, -2.5])}
    slopeExt, sIntersection, zIntersection, coefExt = energyLineTest.getAlphaProfileIntersection(energyLineTestCfg, avaProfileMass, mu)
    print(sIntersection)
    print(zIntersection)
    print(slopeExt)
    print(coefExt)
    atol = 1e-10
    assert sIntersection == pytest.approx(5.1, abs=atol)
    assert zIntersection == pytest.approx(-2.5, abs=atol)
    assert slopeExt == pytest.approx(-0.5, abs=atol)
    assert coefExt == pytest.approx(0.05, abs=atol)

    avaProfileMass = {'s':np.array([0, 0.1, 1.1, 2.1, 3.1, 4.1, 5.1]), 'z':np.array([2.45, 0, -0.5, -1, -1.5, -2, -2.5])}
    slopeExt, sIntersection, zIntersection, coefExt = energyLineTest.getAlphaProfileIntersection(energyLineTestCfg, avaProfileMass, mu)
    print(sIntersection)
    print(zIntersection)
    print(slopeExt)
    print(coefExt)
    atol = 1e-10
    assert sIntersection == pytest.approx(4.8, abs=atol)
    assert zIntersection == pytest.approx(-2.35, abs=atol)
    assert slopeExt == pytest.approx(-0.5, abs=atol)
    assert coefExt == pytest.approx(0.05, abs=atol)


def test_getRunOutAngle(capfd):
    '''test find intersection '''
    avaProfileMass = {'s':np.array([0, 1]), 'z':np.array([0, -1])}
    runOutAngleRad, runOutAngleDeg = energyLineTest.getRunOutAngle(avaProfileMass)

    atol = 1e-10
    assert runOutAngleRad == pytest.approx(math.pi/4, abs=atol)
    assert runOutAngleDeg == pytest.approx(45, abs=atol)


    avaProfileMass = {'s':np.array([0, 1]), 'z':np.array([0, -math.sqrt(3)])}
    runOutAngleRad, runOutAngleDeg = energyLineTest.getRunOutAngle(avaProfileMass)

    atol = 1e-10
    assert runOutAngleRad == pytest.approx(math.pi/3, abs=atol)
    assert runOutAngleDeg == pytest.approx(60, abs=atol)


def test_getEnergyInfo(capfd):
    '''test find intersection '''
    dirname = pathlib.Path(__file__).parents[0]
    energyLineTestCfgFile = dirname / '..' / 'tests' / 'data' / 'testEnergyLine' / 'energyLineTest_com1DFACfg.ini'

    energyLineTestCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=str(energyLineTestCfgFile), modInfo=True)

    g = 9.81
    alphaDeg = 30
    mu = np.tan(np.radians(alphaDeg))

    avaProfileMass = {'s': np.array([0, 2, 4, 6, 8, 10]), 'z': np.array([10, 4, 3, 2, 1, 0]),
    'v2': 2*g * np.array([0, 4, 3, 2, 1, 0])}

    runOutAngleRad, runOutAngleDeg = energyLineTest.getRunOutAngle(avaProfileMass)
    slopeExt, sIntersection, zIntersection, coefExt = energyLineTest.getAlphaProfileIntersection(energyLineTestCfg, avaProfileMass, mu)

    zEne, v2Path, sGeomL, zGeomL, errorEnergyTest = energyLineTest.getEnergyInfo(avaProfileMass, g, mu, sIntersection,
                                                                                 zIntersection, runOutAngleDeg, alphaDeg)
    print(zEne)
    print(v2Path)
    print(sGeomL)
    print(zGeomL)
    print(errorEnergyTest)
    atol = 1e-10
    assert errorEnergyTest['runOutZError'] == pytest.approx(-10, abs=atol)
    assert errorEnergyTest['runOutSError'] == pytest.approx(10, abs=atol)
    assert errorEnergyTest['runOutAngleError'] == pytest.approx(15, abs=atol)
    assert errorEnergyTest['rmseVelocityElevation'] == pytest.approx(2.559271214294478, abs=atol)


def test_mainEnegyLineTest(tmp_path):
    dirname = pathlib.Path(__file__).parents[0]
    energyLineTestCfgFile = dirname / '..' / 'tests' / 'data' / 'testEnergyLine' / 'energyLineTest_com1DFACfg.ini'
    sourceDir = dirname / '..' / 'data' / 'avaParabola' / 'Inputs'
    targetDir = tmp_path / 'avaParabola' / 'Inputs'
    avalancheDir = tmp_path / 'avaParabola'

    shutil.copytree(sourceDir, targetDir)

    # log file name; leave empty to use default runLog.log
    logName = 'runenergyLineTest'

    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    cfgMain['MAIN']['avalancheDir'] = str(avalancheDir)
    # ----------------

    errorEnergyTest = energyLineTest.mainEnergyLineTest(cfgMain, energyLineTestCfgFile)
    print(errorEnergyTest)
    assert abs(errorEnergyTest['runOutSError']) < 0.01
    assert abs(errorEnergyTest['runOutZError']) < 0.01
    assert abs(errorEnergyTest['rmseVelocityElevation']) < 0.01
    assert abs(errorEnergyTest['runOutAngleError']) < 0.0001
