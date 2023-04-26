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


def test_checkCfgFrictionModel():
    """ test if friction model parameters are set to nan if not used """

    cfg = {'GENERAL': {'frictModel': 'samosAT', 'musamosat': '0.155', 'tau0samosat': '0', 'Rs0samosat': '0.222',
        'kappasamosat': '0.43', 'Rsamosat': '0.05', 'Bsamosat': '4.13',
        'muvoellmy': '4000.', 'xsivoellmy': '4000.'}}

    cfg = checkCfg.checkCfgFrictionModel(cfg)

    assert cfg['GENERAL']['frictModel'] == 'samosAT'
    assert cfg['GENERAL']['musamosat'] == '0.155'
    assert cfg['GENERAL']['tau0samosat'] == '0'
    assert cfg['GENERAL']['kappasamosat'] == '0.43'
    assert cfg['GENERAL']['Rsamosat'] == '0.05'
    assert cfg['GENERAL']['Bsamosat'] == '4.13'
    assert np.isnan(cfg['GENERAL']['muvoellmy'])
    assert np.isnan(cfg['GENERAL']['xsivoellmy'])


    cfg = {'GENERAL': {'frictModel': 'samosAT', 'musamosat': '0.155', 'tau0samosat': '0', 'Rs0samosat': '0.222',
        'kappasamosat': '0.43', 'Rsamosat': '0.05', 'Bsamosat': '9.13',
        'muvoellmy': '4000.', 'xsivoellmy': 'nan'}}

    cfg = checkCfg.checkCfgFrictionModel(cfg)

    assert cfg['GENERAL']['frictModel'] == 'samosAT'
    assert cfg['GENERAL']['musamosat'] == '0.155'
    assert cfg['GENERAL']['tau0samosat'] == '0'
    assert cfg['GENERAL']['kappasamosat'] == '0.43'
    assert cfg['GENERAL']['Rsamosat'] == '0.05'
    assert cfg['GENERAL']['Bsamosat'] == '9.13'
    assert np.isnan(cfg['GENERAL']['muvoellmy'])
    assert np.isnan(cfg['GENERAL']['xsivoellmy'])

    cfg = {'GENERAL': {'frictModel': 'samosAT', 'musamosat': '0.155', 'tau0samosat': '0', 'Rs0samosat': '0.222',
        'kappasamosat': '0.43', 'Rsamosat': '0.05', 'Bsamosat': 'nan',
        'muvoellmy': '4000.', 'xsivoellmy': 'nan'}}

    # call function to be tested
    with pytest.raises(ValueError) as e:
        assert checkCfg.checkCfgFrictionModel(cfg)
    assert 'Friction model used samosAT, but Bsamosat is nan - not valid' in str(e.value)

    cfg = {'GENERAL': {'frictModel': 'samosAT', 'musamosat': '0.155', 'tau0samosat': '0', 'Rs0samosat': '0.222',
        'kappasamosat': '0.43', 'Rsamosat': '0.05', 'Bsamosat': 'test',
        'muvoellmy': '4000.', 'xsivoellmy': 'nan'}}

    # call function to be tested
    with pytest.raises(ValueError) as e:
        assert checkCfg.checkCfgFrictionModel(cfg)
    assert 'Friction model used samosAT, but test is not of valid' in str(e.value)
