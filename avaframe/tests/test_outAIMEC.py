''' Tests for module ana3AIMEC plots '''
import pandas as pd
import numpy as np
import pathlib
import configparser
import pytest

# Local imports
import avaframe.out3Plot.outAIMEC as oA
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils

def test_getIndicesVel():
    """ test if indices of a  vector exceeding a threshold are found """

    # setup required input
    a = np.zeros(10)
    a[2:8] = 10.
    a[2] = 2.
    a[5] = 1.2
    velocityThreshold = 1.

    # call function
    indStart, indStop = oA.getIndicesVel(a, velocityThreshold)

    assert indStart == 2
    assert indStop == 7

    a = np.zeros(10) * np.nan
    a[2:8] = 10.
    a[2] = 2.
    a[5] = 1.2
    velocityThreshold = 1.

    # call function
    indStart, indStop = oA.getIndicesVel(a, velocityThreshold)

    assert indStart == 2
    assert indStop == 7

    a = np.zeros(10)
    a[2:] = 10.
    a[2] = 0.9
    velocityThreshold = 1.

    # call function
    indStart, indStop = oA.getIndicesVel(a, velocityThreshold)

    assert indStart == 3
    assert indStop == 9

    a = np.zeros(10)
    a[2:] = 0.4
    a[2] = 0.9
    velocityThreshold = 1.

    with pytest.raises(AssertionError) as e:
        assert oA.getIndicesVel(a, velocityThreshold)
    assert 'No peak flow velocity max along thalweg found exceeding' in str(e.value)
