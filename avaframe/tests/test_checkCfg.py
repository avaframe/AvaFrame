"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
import logging
import pytest
import configparser
import pathlib


from avaframe.com1DFA import checkCfg


def test_checkCfgConsistency():
    """ test check if sphOption and viscOption settings are consistent """

    # setup requuired input data
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'viscOption': 2, 'sphOption': 2}

    # call function to be tested
    testFlag = checkCfg.checkCfgConsistency(cfg)

    assert testFlag == True

    cfg['GENERAL'] = {'viscOption': 2, 'sphOption': 1}

    # call function to be tested
    with pytest.raises(AssertionError) as e:
        assert checkCfg.checkCfgConsistency(cfg)
    assert 'If viscOption is set to 2' in str(e.value)


def test_checkCellSizeKernelRadius():
    """ test check if sphOption and viscOption settings are consistent """

    # setup requuired input data
    cfg = {'GENERAL': {'sphKernelRadius': '5', 'meshCellSize': '10'}}

    # call function to be tested
    cfg = checkCfg.checkCellSizeKernelRadius(cfg)

    assert cfg['GENERAL']['sphKernelRadius'] == '5'

    # setup requuired input data
    cfg = {'GENERAL': {'sphKernelRadius': 'meshCellSize', 'meshCellSize': '10'}}

    # call function to be tested
    cfg = checkCfg.checkCellSizeKernelRadius(cfg)

    assert cfg['GENERAL']['sphKernelRadius'] == '10'
