"""Tests for module com2AB"""
import numpy as np
import math
import pytest
import avaframe.com2AB.com2AB as com2AB


def test_com2AB(capfd):
    '''Simple test for module com2AB'''

    # Make a reference quadratic profile
    B = -np.tan(np.deg2rad(45))
    A = -B/4000
    C = 1000
    N = 1000
    s =  np.linspace(0.0, -B/(2*A), num=N)
    z = np.empty(np.shape(s))
    for i in range(N):
        if (s[i] < (-B/(2*A))):
            z[i] = 1*(A * s[i]*s[i]+B * s[i]+C)
        else:
            z[i] = 1*(-B*B / (4*A) + C)

    theta_beta = 10
    x_beta = (-np.tan(np.deg2rad(theta_beta)) - B)/(2*A)
    y_beta = A*x_beta*x_beta + B*x_beta + C
    beta = np.rad2deg(np.arctan2((C-y_beta),x_beta))
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
    eqParams = com2AB.setEqParameters(smallAva=False)
    eqOut = com2AB.calcAB(eqIn, eqParams)
    alpha = eqOut['alpha']
    alphaSD = eqOut['alphaSD']

    # compare results with a relative tolerance of tol
    tol = 0.001 # here 0.1% relative diff
    assert (alpha == pytest.approx(alpha_ref,rel=tol)) and (alphaSD[0] == pytest.approx(alphaSD_ref[0],rel=tol)) and (alphaSD[1] == pytest.approx(alphaSD_ref[1],rel=tol)) and (alphaSD[2] == pytest.approx(alphaSD_ref[2],rel=tol))
