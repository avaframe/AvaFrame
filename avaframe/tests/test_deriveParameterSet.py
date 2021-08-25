''' Tests for dfa2Aimec '''

import pathlib
import pytest
import configparser
import numpy as np

# Local imports
import avaframe.com1DFA.deriveParameterSet as dP


def test_getVariationDict():
    """ test creating a variation dictionary """

    # setup required input
    avaDir = 'test/avaTest'
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg['GENERAL'] = {'simTypeList' : 'null|ent', 'modelType' : 'dfa', 'resType' : 'ppr|pfd|pfv|particles|FD',
                      'tSteps' : '0:1', 'initPartDistType': 'random', 'initialiseParticlesFromFile': 'False',
                      'particleFile': '', 'seed': '12345', 'rho': '300|400', 'rhoEnt': '100', 'relTh': '1.',
                      'secRelArea': 'True', 'secondaryRelTh': '0.5', 'dt': '0.05', 'tEnd': '400'}
    modDict = {'GENERAL': {'simTypeList': ['null|ent', 'available'], 'resType': ['ppr|pfd|pfv', 'ppr|pfd|pfv|particles|FD'],
                'tSteps': ['0:1', '1'], 'rho': ['300|400', '200'], 'secRelArea': ['True', 'False']},
                'TEST': {'test': ['test1', '']}}

    # call function to be tested
    variations = dP.getVariationDict(avaDir, cfg, modDict)

    print('variations', variations)

    variationsTest = {'simTypeList': ['null', 'ent'], 'rho': np.asarray([300, 400])}

    assert len(variations.keys()) == len(variationsTest.keys())
    assert variations['simTypeList'][0] == 'null'
    assert variations['simTypeList'][1] == 'ent'
    assert np.array_equal(variations['rho'], np.asarray([300, 400]))


def test_validateVarDict():
    """ test if variation dict has only valid parameters """

    # setup required input
    variationDict = {'simTypeList': ['null', 'ent'], 'rho': np.asarray([300, 400]), 'relThAA': np.asarray([1.0, 2.0]),
                     'secRelAreA': ['False', 'True'], 'rhoEnt': '200.'}
    standardCfg = configparser.ConfigParser()
    standardCfg.optionxform = str
    standardCfg['GENERAL'] = {'simTypeList' : 'available', 'modelType' : 'dfa', 'resType' : 'ppr|pfd|pfv',
                  'tSteps' : '0:1', 'initPartDistType': 'random', 'initialiseParticlesFromFile': 'False',
                  'particleFile': '', 'seed': '12345', 'rho': '300|400', 'rhoEnt': '100', 'relTh': '1.',
                  'secRelArea': 'True', 'secondaryRelTh': '0.5', 'dt': '0.05', 'tEnd': '400'}

    # call function to be tested
    variationDictTest = dP.validateVarDict(variationDict, standardCfg)

    print('variationDictTest', variationDictTest)

    assert len(variationDictTest.keys()) == 3
    assert variationDictTest['simTypeList'][0] == 'null'
    assert variationDictTest['simTypeList'][1] == 'ent'
    assert variationDictTest['rhoEnt'] == ['200.']
    assert np.array_equal(variationDictTest['rho'], np.asarray([300, 400]))
    assert 'relThAA' not in variationDictTest.keys()
