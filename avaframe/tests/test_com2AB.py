"""Tests for module com2AB"""
import numpy as np
import math
import pytest
import avaframe.com2AB.com2AB as com2AB


def test_com2AB(capfd):
    '''Simple test for module com2AB'''

    eqIn = {}
    s =  np.linspace(0.0, 4000.0, num=401)
    eqIn['s'] = s  # curvilinear coordinate (of the x, y path)
    eqIn['x'] = []  # x coordinate of the path
    eqIn['y'] = []  # y coordinate of the path
    A = 0.000158
    B = -0.821131
    C = 1061.832865
    z = np.empty(np.shape(s))
    for i in range(401):
        if (s[i] < (-B/(2*A))):
            z[i] = 1*(A * s[i]*s[i]+B * s[i]+C)
        else:
            z[i] = 1*(-B*B / (4*A) + C)

    eqIn['z'] = z  # z coordinate of the path (projection of x,y on the raster)
    eqIn['indSplit'] = 2  # index of split point
    eqParams = com2AB.setEqParameters(smallAva=False)
    eqOut = com2AB.calcAB(eqIn, eqParams)

    alpha = eqOut['alpha']
    alphaSD = eqOut['alphaSD']
    print(alphaSD)
    assert (alpha == pytest.approx(24.700714632848285)) and (alphaSD[0] == pytest.approx(25.95071463)) and (alphaSD[1] == pytest.approx(23.45071463)) and (alphaSD[2] == pytest.approx(22.20071463))
