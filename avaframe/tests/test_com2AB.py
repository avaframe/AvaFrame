"""Tests for module com2AB"""
import numpy as np

import avaframe.com2AB.com2AB as com2AB


def test_com2AB(capfd):
    '''Simple test for module com2AB'''

    eqIn = {}
    s = np.linespace(0.0, 1000.0, 101)
    eqIn['s'] = s  # curvilinear coordinate (of the x, y path)
    eqIn['x'] = []  # x coordinate of the path
    eqIn['y'] = []  # y coordinate of the path
    A = 0.000158
    B = -0.821131
    C = 1061.832865
    z = np.empty(np.shape(s))
    for i in range(101):
        if (s[i] < -B/2/A):
            z[i] = 1*(A * s[i]*s[i]+B * s[i]+C)
        else:
            z[i] = 1*(-B ^ 2 / 4 / A + C)

    eqIn['z'] = z  # z coordinate of the path (projection of x,y on the raster)
    eqIn['indSplit'] = 50  # index of split point

    eqParams = com2AB.setEqParameters(smallAva=False)
    eqOut = com2AB.calcAB(eqIn, eqParameters)

    alpha = eqOut['alpha']
    alphaSD = eqOut['alphaSD']

    assert (alpha == alpha_ref) and (alphaSD == alphaSD_ref)
