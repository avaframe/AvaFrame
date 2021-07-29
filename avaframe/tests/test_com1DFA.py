"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
from avaframe.com1DFA import com1DFA
import pytest



def test_setDEMOriginToZero():
    """ test if origin is set to zero """

    # setup required input
    class testHeader:
        def __init__(self):
            self.xllcenter = 1.0
            self.yllcenter = 4.0

    tHeader = testHeader()
    dem = {'header': tHeader}

    demTest = com1DFA.setDEMoriginToZero(dem)

    assert demTest['header'].xllcenter == 0.0
    assert demTest['header'].yllcenter == 0.0
