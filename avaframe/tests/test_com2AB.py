"""Tests for module com2AB"""
import numpy as np
import pytest
import avaframe.com2AB.com2AB as com2AB

def test_setEqParameters(capfd):
    '''Simple test for module setEqParameters'''
    # small avalanche
    eqParameters = {}
    eqParameters['ParameterSet'] = 'Small avalanches'
    eqParameters['k1'] = 0.933
    eqParameters['k2'] = 0.0
    eqParameters['k3'] = 0.0088
    eqParameters['k4'] = -5.02
    eqParameters['SD'] = 2.36

    eqParams = com2AB.setEqParameters(smallAva=True, customParam=None)
    assert(eqParams['k1'] == pytest.approx(eqParameters['k1'])) and (
        eqParams['ParameterSet'] == (eqParameters['ParameterSet']))

    eqParameters = {}
    eqParameters['ParameterSet'] = 'Standard'
    eqParameters['k1'] = 1.05
    eqParameters['k2'] = -3130.0
    eqParameters['k3'] = 0.0
    eqParameters['k4'] = -2.38
    eqParameters['SD'] = 1.25

    eqParams = com2AB.setEqParameters(smallAva=False, customParam=None)
    assert(eqParams['k1'] == pytest.approx(eqParameters['k1'])) and (
        eqParams['ParameterSet'] == (eqParameters['ParameterSet']))

    customParam = {}
    customParam['k1'] = 1
    customParam['k2'] = 2
    customParam['k3'] = 3
    customParam['k4'] = 4
    customParam['SD'] = 5

    eqParameters = {}
    eqParameters['ParameterSet'] = 'Custom'
    eqParameters['k1'] = customParam['k1']
    eqParameters['k2'] = customParam['k2']
    eqParameters['k3'] = customParam['k3']
    eqParameters['k4'] = customParam['k4']
    eqParameters['SD'] = customParam['SD']

    eqParams = com2AB.setEqParameters(smallAva=False, customParam=customParam)
    assert(eqParams['k1'] == pytest.approx(eqParameters['k1'])) and (
        eqParams['ParameterSet'] == (eqParameters['ParameterSet']))

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

    theta_beta = 10
    x_beta = (-np.tan(np.deg2rad(theta_beta)) - B)/(2*A)
    y_beta = A*x_beta*x_beta + B*x_beta + C
    beta = np.rad2deg(np.arctan2((C-y_beta), x_beta))
    # use standard coeef
    k1 = 1.05
    k2 = -3130.0
    k3 = 0.0
    k4 = -2.38
    SD = 1.25
    alpha_ref = k1 * beta + k2 * 2*A + k3 * B*B/(2*A) + k4
    SDs = [SD, -1*SD, -2*SD]
    alphaSD_ref = k1 * beta + k2 * 2*A + k3 * B*B/(2*A) + k4 + SDs

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
    assert (alpha == pytest.approx(alpha_ref, rel=tol)) and (alphaSD[0] == pytest.approx(alphaSD_ref[0], rel=tol)) and (
        alphaSD[1] == pytest.approx(alphaSD_ref[1], rel=tol)) and (alphaSD[2] == pytest.approx(alphaSD_ref[2], rel=tol))


def test_find_10Point(capfd):
    '''find_10Point'''
    s = np.linspace(0,400,41)
    angle = np.linspace(40,0,41)
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    delta_ind = 3
    ids_10Point = com2AB.find_10Point(tmp,delta_ind)
    assert ids_10Point==30

    delta_ind = 0
    ids_10Point = com2AB.find_10Point(tmp,delta_ind)
    assert ids_10Point==30

    angle[10] = 8
    angle[11] = 8
    angle[12] = 8
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    delta_ind = 3
    ids_10Point = com2AB.find_10Point(tmp,delta_ind)
    assert ids_10Point==30

    angle[13] = 8
    tmp = np.where((angle < 10.0) & (angle > 0.0))
    ids_10Point = com2AB.find_10Point(tmp,delta_ind)
    assert ids_10Point==9
