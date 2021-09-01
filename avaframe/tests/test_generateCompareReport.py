"""
    Pytest for generate Report

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.log2Report import generateCompareReport as gR

import pytest
import configparser
import pathlib


def test_makeLists():
    """ test creating lists for table creation """

    # setup required input
    simDict = {'type': 'list', 'maxP': '100.', 'minF': 0.0, 'testKey': 10., 'testIgnore': ''}
    benchDict = {'type': 'list', 'maxP': '200.', 'minF': 0.0, 'testIgnore2': 'A'}

    # call function to be tested
    parameterList, valuesSim, valuesBench = gR.makeLists(simDict, benchDict)

    assert parameterList == ['maxP', 'minF', 'testKey']
    assert valuesSim == ['100.', 0.0, 10.]
    assert valuesBench == ['200.', 0.0, 'non existent']
