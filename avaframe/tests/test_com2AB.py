"""Tests for module com2AB"""
import numpy as np
import pytest
import pathlib

# Local imports
import avaframe.com2AB.com2AB as com2AB
from avaframe.out3Plot import outAB
from avaframe.in3Utils import cfgUtils


def test_setEqParameters(capfd):
    '''Simple test for module setEqParameters'''
    cfg = cfgUtils.getModuleConfig(com2AB)
    # small avalanche
    eqParamRef = {}
    eqParamRef['ParameterSet'] = 'Small avalanches'
    eqParamRef['k1'] = 0.933
    eqParamRef['k2'] = 0.0
    eqParamRef['k3'] = 0.0088
    eqParamRef['k4'] = -5.02
    eqParamRef['SD'] = 2.36

    eqParams = com2AB.setEqParameters(cfg, smallAva=True)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]

    eqParamRef = {}
    eqParamRef['ParameterSet'] = 'Standard'
    eqParamRef['k1'] = 1.05
    eqParamRef['k2'] = -3130.0
    eqParamRef['k3'] = 0.0
    eqParamRef['k4'] = -2.38
    eqParamRef['SD'] = 1.25

    eqParams = com2AB.setEqParameters(cfg, smallAva=False)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]

    cfg['ABSETUP']['customParam'] = 'True'
    cfg['ABSETUP']['k1'] = '1'
    cfg['ABSETUP']['k2'] = '2'
    cfg['ABSETUP']['k3'] = '3'
    cfg['ABSETUP']['k4'] = '4'
    cfg['ABSETUP']['SD'] = '5'

    eqParamRef = {}
    eqParamRef['ParameterSet'] = 'Custom'
    eqParamRef['k1'] = cfg.getfloat('ABSETUP', 'k1')
    eqParamRef['k2'] = cfg.getfloat('ABSETUP', 'k2')
    eqParamRef['k3'] = cfg.getfloat('ABSETUP', 'k3')
    eqParamRef['k4'] = cfg.getfloat('ABSETUP', 'k4')
    eqParamRef['SD'] = cfg.getfloat('ABSETUP', 'SD')

    eqParams = com2AB.setEqParameters(cfg, smallAva=False)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]


def test_calcAB(capfd):
    '''Simple test for function calcAB'''

    cfg = cfgUtils.getModuleConfig(com2AB)
    # Make a reference quadratic profile
    B = -np.tan(np.deg2rad(45))
    A = -B/4000
    C = 1000
    N = 1000
    s = np.linspace(0.0, -B/(2*A), num=N)
    z = np.empty(np.shape(s))
    for i in range(N):
        if (s[i] < (-B/(2*A))):
            z[i] = 1*(A * s[i]*s[i]+B * s[i]+C)
        else:
            z[i] = 1*(-B*B / (4*A) + C)

    thetaBeta = 10
    xBeta = (-np.tan(np.deg2rad(thetaBeta)) - B)/(2*A)
    yBeta = A*xBeta*xBeta + B*xBeta + C
    beta = np.rad2deg(np.arctan2((C-yBeta), xBeta))
    # use standard coeef
    k1 = 1.05
    k2 = -3130.0
    k3 = 0.0
    k4 = -2.38
    SD = 1.25
    alpharef = k1 * beta + k2 * 2*A + k3 * B*B/(2*A) + k4
    SDs = [SD, -1*SD, -2*SD]
    alphaSDref = k1 * beta + k2 * 2*A + k3 * B*B/(2*A) + k4 + SDs

    # Using com2AB.calcAB to get the solution
    eqIn = {}
    eqIn['s'] = s  # curvilinear coordinate (of the x, y path)
    eqIn['x'] = []  # x coordinate of the path
    eqIn['y'] = []  # y coordinate of the path
    eqIn['z'] = z  # z coordinate of the path (projection of x,y on the raster)
    eqIn['indSplit'] = 2  # index of split point
    eqParams = com2AB.setEqParameters(cfg, smallAva=False)
    eqOut = com2AB.calcAB(eqIn, eqParams, 30)
    alpha = eqOut['alpha']
    alphaSD = eqOut['alphaSD']

    # compare results with a relative tolerance of tol
    tol = 0.002  # here 0.1% relative diff
    assert (alpha == pytest.approx(alpharef, rel=tol)) and (
            alphaSD[0] == pytest.approx(alphaSDref[0], rel=tol)) and (
            alphaSD[1] == pytest.approx(alphaSDref[1], rel=tol)) and (
            alphaSD[2] == pytest.approx(alphaSDref[2], rel=tol))


# test complete routine
def test_com2ABMain(capfd):
    '''Simple test for function com2ABMain'''
    # load and prepare Inputs
    listNames = ['avaHockeySmall', 'avaHockeyChannel', 'avaBowl']
    dirname = pathlib.Path(__file__).parents[0]
    for name in listNames:
        avalancheDir = dirname / '..' / 'data' / name
        saveOutPathRef = dirname / '..' / '..' / 'benchmarks' / (name + 'ABPytest')
        cfg = cfgUtils.getModuleConfig(com2AB)
        flags = cfg['FLAGS']
        # run main routine
        resAB = com2AB.com2ABMain(cfg, avalancheDir)
        eqParams = resAB['eqParams']
        eqOut = resAB[name]

        # open ref data
        flags['fullOut'] = 'True'
        eqParamsRef, eqOutRef = outAB.readABresults(saveOutPathRef, name, flags)

        for key in eqParamsRef.keys():
            assert eqParamsRef[key] == eqParams[key]

        atol = 1e-10
        assert (np.allclose(eqOutRef['x'], eqOut['x'], atol=atol)) and (
                np.allclose(eqOutRef['y'], eqOut['y'], atol=atol)) and (
                np.allclose(eqOutRef['z'], eqOut['z'], atol=atol)) and (
                np.allclose(eqOutRef['s'], eqOut['s'], atol=atol))
        assert (np.allclose(eqOutRef['alpha'], eqOut['alpha'], atol=atol)
                ) and (np.allclose(eqOutRef['alphaSD'], eqOut['alphaSD'],
                       atol=atol))


# Compare to QGIS AB routine
def test_QGISAB(capfd):
    '''Compare com2ABMain results to QGIS AB results for avaSlide'''
    # load and prepare Inputs
    avaName = 'avaSlide'
    dirname = pathlib.Path(__file__).parents[0]
    avalancheDir = dirname / '..' / 'data' / avaName
    cfg = cfgUtils.getModuleConfig(com2AB)
    # run main routine
    resAB = com2AB.com2ABMain(cfg, avalancheDir)
    # process data to get results
    AvaPath = resAB['AvaPath']
    NameAva = AvaPath['Name']
    for i in range(len(NameAva)):
        name = NameAva[i]
        resAB = outAB.processABresults(resAB, name)
        beta = resAB[name]['beta']
        alpha = resAB[name]['alpha']
        alphaSD = resAB[name]['alphaSD']
        s = resAB[name]['s']
        ids_alpha = resAB[name]['ids_alpha']
        ids10Point = resAB[name]['ids10Point']
        ids_alphaP1SD = resAB[name]['ids_alphaP1SD']
        ids_alphaM1SD = resAB[name]['ids_alphaM1SD']
        ids_alphaM2SD = resAB[name]['ids_alphaM2SD']
        # get ref results
        nameRef = name + '_AB_QGIS.txt'
        nameRefpath = pathlib.Path(dirname, '..', '..', 'benchmarks', avaName + 'ABPytest',
                                   nameRef)
        data = np.loadtxt(nameRefpath, skiprows=4, delimiter=',')
        tolDist = 10
        tolAngle = 0.12
        assert (alpha == pytest.approx(data[0, 4], abs=tolAngle)) and (
                beta == pytest.approx(data[1, 4], rel=tolAngle)) and (
                alphaSD[1] == pytest.approx(data[2, 4], rel=tolAngle)) and (
                alphaSD[2] == pytest.approx(data[3, 4], rel=tolAngle)) and (
                alphaSD[0] == pytest.approx(data[4, 4], rel=tolAngle))
        assert (s[ids_alpha] == pytest.approx(data[0, 3], abs=tolDist)) and (
                s[ids10Point] == pytest.approx(data[1, 3], rel=tolDist)) and (
                s[ids_alphaM1SD] == pytest.approx(data[2, 3], rel=tolDist)) and (
                s[ids_alphaM2SD] == pytest.approx(data[3, 3], rel=tolDist)) and (
                s[ids_alphaP1SD] == pytest.approx(data[4, 3], rel=tolDist))
