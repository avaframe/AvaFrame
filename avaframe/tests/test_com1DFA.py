"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
from avaframe.com1DFA import com1DFA
import avaframe.in2Trans.ascUtils as IOf
import pytest
import configparser



def test_setDEMOriginToZero():
    """ test if origin is set to zero """

    # setup required input
    tHeader = IOf.cASCheader()
    tHeader.xllcenter = 10.
    tHeader.yllcenter = 4.0
    dem = {'header': tHeader}

    # call function to be tested
    demTest = com1DFA.setDEMoriginToZero(dem)

    assert demTest['header'].xllcenter == 0.0
    assert demTest['header'].yllcenter == 0.0


def test_placeParticles():
    """ test placing of particles """

    # setup required inputs
    massCell = 10.
    indx = 0
    indy = 1
    csz = 5
    massPerPart = 2.
    initPartDistType = 'uniform'
    rng = np.random.default_rng(12345)

    # call funciton to be tested - uniform
    xpart, ypart, mPart, nPart =  com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest = np.asarray([-1.66666666, 0.0, 1.66666666, -1.66666666, 0., 1.66666666, -1.66666666,
                        0.0, 1.66666666])
    ypartTest = np.asarray([3.33333333, 3.33333333, 3.33333333, 5.0, 5., 5., 6.66666666, 6.66666666,
                        6.66666666])

    # call funciton to be tested - uniform
    massCell = 8.
    xpart2, ypart2, mPart2, nPart2 =  com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest2 = np.asarray([-1.25, 1.25, -1.25, 1.25])
    ypartTest2 = np.asarray([3.75, 3.75, 6.25, 6.25])

    # call funciton to be tested - random
    massCell = 11.5
    initPartDistType = 'random'
    xpart3, ypart3, mPart3, nPart3 =  com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest3 = np.asarray([-0.9162083, 1.48682729, 0.88127335, -0.54445225, -0.83593036, 0.49154377])
    ypartTest3 = np.asarray([3.43367093, 5.86378022, 7.20901433, 3.74122857, 7.24440576, 5.83618727])

    assert nPart == 9.0
    assert np.isclose(mPart, 1.111111)
    assert np.allclose(xpart, xpartTest)
    assert np.allclose(ypart, ypartTest)
    assert nPart2 == 4.0
    assert mPart2 == 2.
    assert np.allclose(xpart2, xpartTest2)
    assert np.allclose(ypart2, ypartTest2)
    assert nPart3 == 6.0
    assert np.isclose(mPart3, 1.9166666666666)
    assert np.allclose(xpart3, xpartTest3)
    assert np.allclose(ypart3, ypartTest3)
