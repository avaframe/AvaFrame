"""
    Pytest for module com1DFA
"""

import configparser
import copy
import logging
import pathlib
import pickle
import shutil

#  Load modules
import numpy as np
import pytest

from avaframe.com4FlowPy import flowClass

def test_add_os():
    cell = flowClass.Cell(1,1,
                          np.array([[10,10,10], [10,10,10], [10,10,10]]), 10,
                          1,0, None,
                          20, 8, 3e-4, 270,
                          startcell=True)
    cell.add_os(0.2)
    assert cell.flux == 1.2
