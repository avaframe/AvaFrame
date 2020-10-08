"""Tests for module com2AB"""
import numpy as np
import pytest
import os

# Local imports
import avaframe.com2AB.com2AB as com2AB
from avaframe.out3Plot import outAB
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


def test_setEqParameters(capfd):
    '''Simple test for module setEqParameters'''
    # small avalanche
    eqParamRef = {}
    eqParamRef['ParameterSet'] = 'Small avalanches'
    eqParamRef['k1'] = 0.933
    eqParamRef['k2'] = 0.0
    eqParamRef['k3'] = 0.0088
    eqParamRef['k4'] = -5.02
    eqParamRef['SD'] = 2.36

    eqParams = com2AB.setEqParameters(smallAva=True, customParam=None)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]

    eqParamRef = {}
    eqParamRef['ParameterSet'] = 'Standard'
    eqParamRef['k1'] = 1.05
    eqParamRef['k2'] = -3130.0
    eqParamRef['k3'] = 0.0
    eqParamRef['k4'] = -2.38
    eqParamRef['SD'] = 1.25

    eqParams = com2AB.setEqParameters(smallAva=False, customParam=None)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]

    customParam = {}
    customParam['k1'] = 1
    customParam['k2'] = 2
    customParam['k3'] = 3
    customParam['k4'] = 4
    customParam['SD'] = 5

    eqParamRef = {}
    eqParamRef['ParameterSet'] = 'Custom'
    eqParamRef['k1'] = customParam['k1']
    eqParamRef['k2'] = customParam['k2']
    eqParamRef['k3'] = customParam['k3']
    eqParamRef['k4'] = customParam['k4']
    eqParamRef['SD'] = customParam['SD']

    eqParams = com2AB.setEqParameters(smallAva=False, customParam=customParam)
    for key in eqParamRef.keys():
        assert eqParamRef[key] == eqParams[key]

# def test_prepareLine(capfd):
#     '''Simple test for function prepareLine'''
#     header
#     header.xllcorner = 10
#     header.yllcorner = -15
#     header.cellsize = 5
#     header.ncols = 200
#     header.nrows = 100
#     x = np.linspace(header.xllcorner)
#     rasterdata
#     avapath = np.array([[], []])
#     AvaProfile, SplitPoint, indSplit = prepareLine(
#         header, rasterdata, avapath, splitPoint, distance=10)


def test_calcAB(capfd):
    '''Simple test for function calcAB'''

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
    eqParams = com2AB.setEqParameters(smallAva=False, customParam=None)
    eqOut = com2AB.calcAB(eqIn, eqParams)
    alpha = eqOut['alpha']
    alphaSD = eqOut['alphaSD']

    # compare results with a relative tolerance of tol
    tol = 0.001  # here 0.1% relative diff
    assert (alpha == pytest.approx(alpharef, rel=tol)) and (
            alphaSD[0] == pytest.approx(alphaSDref[0], rel=tol)) and (
            alphaSD[1] == pytest.approx(alphaSDref[1], rel=tol)) and (
            alphaSD[2] == pytest.approx(alphaSDref[2], rel=tol))


# test complete routine
def test_com2ABMain(capfd):
    '''Simple test for function com2ABMain'''
    # load and prepare Inputs
    dirname = os.path.dirname(__file__)
    avalancheDir = os.path.join(dirname, '..', 'data',
                                'avaHockeySmoothChannel')
    saveOutPath = os.path.join(dirname, 'data')
    saveOutPathRef = os.path.join(dirname, '..', '..', 'benchmarks',
                                  'avaHockeySmoothChannel')
    cfg = cfgUtils.getModuleConfig(com2AB)
    cfgSetup = cfg['ABSETUP']
    flags = cfg['FLAGS']

    # Extract input file locations
    cfgPath = com2AB.readABinputs(avalancheDir)

    # Read input data for ALPHABETA
    dem = IOf.readRaster(cfgPath['demSource'])
    avaPath = shpConv.readLine(cfgPath['profileLayer'], cfgPath['defaultName'],
                               dem)
    splitPoint = shpConv.readPoints(cfgPath['splitPointSource'], dem)
    # run main routine
    com2AB.com2ABMain(dem, avaPath, splitPoint,
                      saveOutPath, cfgSetup)
    eqParams, eqOut = outAB.readABresults(saveOutPath, avaPath['Name'][0],
                                          flags)
    # open ref data
    flags['fullOut'] = 'True'
    eqParamsRef, eqOutRef = outAB.readABresults(saveOutPathRef,
                                                avaPath['Name'][0], flags)

    for key in eqParamsRef.keys():
        assert eqParamsRef[key] == eqParams[key]

    atol = 1e-10
    assert (np.allclose(eqOutRef['x'], eqOutRef['x'], atol=atol)) and (
            np.allclose(eqOutRef['y'], eqOutRef['y'], atol=atol)) and (
            np.allclose(eqOutRef['z'], eqOutRef['z'], atol=atol)) and (
            np.allclose(eqOutRef['s'], eqOutRef['s'], atol=atol)) and (
            np.allclose(eqOutRef['alpha'], eqOutRef['alpha'], atol=atol)) and (
            np.allclose(eqOutRef['alphaSD'], eqOutRef['alphaSD'], atol=atol))
