"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
from avaframe.com1DFA import com1DFA
import avaframe.in2Trans.ascUtils as IOf
import pytest



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
