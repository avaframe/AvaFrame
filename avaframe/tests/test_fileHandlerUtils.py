"""
    Pytest for fileHandlerUtils

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.in3Utils import fileHandlerUtils as fU
import pytest
import configparser


def test_readLogFile():
    """ Test if logDict is generated correctly """

    # Test function
    logName = os.path.join(os.getcwd(), 'avaframe', 'tests', 'data', 'testExpLog.txt')
    cfg = configparser.ConfigParser()
    cfg = {'varPar' : 'RelTh'}
    logDict = fU.readLogFile(logName, cfg)
    logDictMu = fU.readLogFile(logName)

    assert logDict['noSim'][4] == 5
    assert logDict['simName'][2] == 'release1HS2_null_dfa'
    assert logDict['RelTh'][2] == 4.0
    assert logDictMu['noSim'][4] == 5
    assert logDictMu['simName'][2] == 'release1HS2_null_dfa'
    assert logDictMu['Mu'][2] == 4.0
