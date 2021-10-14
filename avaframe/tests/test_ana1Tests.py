"""Tests for module ana1Tests"""
import numpy as np
import pytest

# Local imports
import avaframe.ana1Tests.analysisTools as anaTools


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
