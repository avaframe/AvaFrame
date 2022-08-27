''' Tests for dfa2Aimec '''

import pathlib
import pytest
import configparser
import numpy as np
import pandas as pd

# Local imports
import avaframe.com1DFA.deriveParameterSet as dP


def test_getVariationDict():
    """ test creating a variation dictionary """

    # setup required input
    avaDir = 'test/avaTest'
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg['GENERAL'] = {'simTypeList' : 'null|ent', 'modelType' : 'dfa', 'resType' : 'ppr|pft|pfv|particles|FT',
                      'tSteps' : '0:1', 'initPartDistType': 'random', 'initialiseParticlesFromFile': 'False',
                      'particleFile': '', 'seed': '12345', 'rho': '300|400', 'rhoEnt': '100', 'relTh': '1.',
                      'secRelArea': 'True', 'secondaryRelTh': '0.5', 'dt': '0.05', 'tEnd': '400'}
    modDict = {'GENERAL': {'simTypeList': ['null|ent', 'available'], 'resType': ['ppr|pft|pfv', 'ppr|pft|pfv|particles|FT'],
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


    cfg['GENERAL']['relThFromShp'] = 'True'
    cfg['GENERAL']['relThFromFile'] = 'False'
    cfg['GENERAL']['relTh'] = ''
    cfg['GENERAL']['relThPercentVariation'] = '40$5'

    modDict['GENERAL']['relThPercentVariation'] = '40$5'

    # call function to be tested
    variations2 = dP.getVariationDict(avaDir, cfg, modDict)

    print('variations2', variations2)

    variationsTest2 = {'simTypeList': ['null', 'ent'], 'rho': np.asarray([300, 400]),
        'relThPercentVariation': np.linspace(0.6, 1.4, 5)}

    assert len(variations2.keys()) == len(variationsTest2.keys())
    assert variations2['simTypeList'][0] == 'null'
    assert variations2['simTypeList'][1] == 'ent'
    assert np.array_equal(variations2['relThPercentVariation'], variationsTest2['relThPercentVariation'])
    assert np.array_equal(variations2['rho'], np.asarray([300, 400]))


def test_validateVarDict():
    """ test if variation dict has only valid parameters """

    # setup required input
    variationDict = {'simTypeList': ['null', 'ent'], 'rho': np.asarray([300, 400]), 'relThAA': np.asarray([1.0, 2.0]),
                     'secRelAreA': ['False', 'True'], 'rhoEnt': '200.'}
    standardCfg = configparser.ConfigParser()
    standardCfg.optionxform = str
    standardCfg['GENERAL'] = {'simTypeList': 'available', 'modelType': 'dfa', 'resType': 'ppr|pft|pfv',
                  'tSteps': '0:1', 'initPartDistType': 'random', 'initialiseParticlesFromFile': 'False',
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

    variationDict = {'simTypeList': ['null', 'ent'], 'rho': np.asarray([300, 400]), 'relThAA': np.asarray([1.0, 2.0]),
                     'secRelAreA': ['False', 'True'], 'rhoEnt': '100:200:2&400'}
    # call function to be tested
    variationDictTest = dP.validateVarDict(variationDict, standardCfg)

    print('variationDictTest', variationDictTest)

    assert len(variationDictTest.keys()) == 3
    assert variationDictTest['simTypeList'][0] == 'null'
    assert variationDictTest['simTypeList'][1] == 'ent'
    assert variationDictTest['rhoEnt'][0] == 100.
    assert variationDictTest['rhoEnt'][1] == 200.
    assert variationDictTest['rhoEnt'][2] == 400.
    assert len(variationDictTest['rhoEnt']) == 3
    assert np.array_equal(variationDictTest['rho'], np.asarray([300, 400]))
    assert 'relThAA' not in variationDictTest.keys()


def test_checkResType():
    """ test checking if desired result type is in availble result types """

    # setup required input
    fullCfg = configparser.ConfigParser()
    fullCfg['GENERAL'] = {'resType': 'pft|ppr|pfv|particles|test1|test2'}
    section = 'GENERAL'
    key = 'resType'
    value = 'pft|ppr|pfv|particles|test1|test2'

    # call function to be tested
    fullCfg = dP.checkResType(fullCfg, section, key, value)

    print('fullCfg', fullCfg)

    assert fullCfg['GENERAL']['resType'] == 'pft|ppr|pfv|particles'


def test_getThicknessValue():
    """ test setting of thickness values """

    inputSimFiles = {'release1HS': {'thickness': ['1.2', '1.4']}}
    inputSimFiles['release1HS']['id'] = ['0', '1']
    inputSimFiles['release1HS']['ci95'] = ['None', 'None']

    thType = 'relTh'
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'relThFromFile': 'False', 'relTh': '', 'relThFromShp': 'True',
        'relThPercentVariation': '40$3', 'relThDistVariation': ''}
    cfg['INPUT'] = {'releaseScenario': 'release1HS'}

    cfg = dP.getThicknessValue(cfg, inputSimFiles, 'release1HS', thType)

    assert cfg['INPUT']['relThId'] == '0|1'
    assert cfg['INPUT']['relThThickness'] == '1.2|1.4'
    assert cfg['GENERAL']['relThPercentVariation'] == '40$3'

    inputSimFiles = {'release1HS': {'thickness': ['1.2', 'None']}}
    inputSimFiles['release1HS']['id'] = ['0', '1']
    inputSimFiles['release1HS']['ci95'] = ['None', 'None']

    thType = 'relTh'
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'relThFromFile': 'False', 'relTh': '', 'relThFromShp': 'True',
        'relThPercentVariation': '40$3', 'relThDistVariation': ''}
    cfg['INPUT'] = {'releaseScenario': 'release1HS'}

    with pytest.raises(AssertionError) as e:
        assert dP.getThicknessValue(cfg, inputSimFiles, 'release1HS', thType)
    assert str(e.value) == "Not all features in shape file have a thickness value - check shape file attributes: %s" % ('release1HS')

    inputSimFiles = {'release1HS': {'thickness': ['1.2', 'None']}}
    inputSimFiles['release1HS']['id'] = ['0', '1']
    inputSimFiles['release1HS']['ci95'] = ['None', 'None']

    thType = 'relTh'
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'relThFromFile': 'False', 'relTh': '40$2', 'relThFromShp': 'False',
        'relThPercentVariation': '', 'relThDistVariation': ''}
    cfg['INPUT'] = {'releaseScenario': 'release1HS'}

    with pytest.raises(AssertionError) as e:
        assert dP.getThicknessValue(cfg, inputSimFiles, 'release1HS', thType)
    assert str(e.value) == "Format of relTh value in ini file is not correct - for variation from ini use refValue$percent$numberOfSteps"

    inputSimFiles = {'release1HS': {'thickness': ['1.2', 'None']}}
    inputSimFiles['release1HS']['id'] = ['0', '1']
    inputSimFiles['release1HS']['ci95'] = ['None', 'None']

    thType = 'relTh'
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'relThFromFile': 'False', 'relTh': '1.0', 'relThFromShp': 'False',
        'relThPercentVariation': '', 'relThDistVariation': ''}
    cfg['INPUT'] = {'releaseScenario': 'release1HS'}

    cfg = dP.getThicknessValue(cfg, inputSimFiles, 'release1HS', thType)

    assert cfg.has_option('INPUT', 'relThId') is False
    assert cfg.has_option('INPUT', 'relThThickness') is False
    assert cfg['GENERAL']['relThPercentVariation'] == ''
    assert cfg['GENERAL']['relTh'] == '1.0'


def test_checkThicknessSettings():
    """ test checking thickness settings function """

    # setup required inputs
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'entThFromShp': 'True', 'entTh': '', 'entThPercentVariation': '',
                      'entThRangeVariation': ''}

    thName = 'entTh'

    thicknessSettingsCorrect = dP.checkThicknessSettings(cfg, thName)

    assert thicknessSettingsCorrect

    cfg['GENERAL']['entTh'] = '0.3'

    with pytest.raises(AssertionError) as e:
        assert dP.checkThicknessSettings(cfg, thName)
    assert str(e.value) == "If %s is set to True - it is not allowed to set a value for %s" % ('entThFromShp', 'entTh')

    cfg['GENERAL']['entThFromShp'] = 'False'
    cfg['GENERAL']['entTh'] = ''
    cfg['GENERAL']['entThFromFile'] = 'False'


    with pytest.raises(AssertionError) as e:
        assert dP.checkThicknessSettings(cfg, thName)
    assert str(e.value) == "If %s is set to False - it is required to set a value for %s" % ('entThFromShp', 'entTh')


    cfg['GENERAL']['entThFromShp'] = ''


    with pytest.raises(AssertionError) as e:
        assert dP.checkThicknessSettings(cfg, thName)
    assert str(e.value) == "Check %s - needs to be True or False" % ('entThFromShp')

    cfg['GENERAL']['relThFromShp'] = 'False'
    cfg['GENERAL']['relThFromFile'] = 'True'
    cfg['GENERAL']['relTh'] = '1.0'

    with pytest.raises(AssertionError) as e:
        assert dP.checkThicknessSettings(cfg, 'relTh')
    assert str(e.value) == ("If %s is set to True - it is not allowed to set %s to True or provide a value in %s" %
        ('relThFromFile', 'relThFromShp', 'relTh'))


def test_appendShpThickness():
    """ test appending thickness values """

    # setup required inputs
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'secRelArea': 'False', 'simTypeActual': 'null', 'relThFromShp': 'True', 'relTh': '',
        'relThFromFile': 'False', 'relThPercentVariation': '', 'relThRangeVariation': '',
        'relThDistVariation': ''}
    cfg['INPUT'] = {'relThThickness': '1.2|1.4', 'relThId': '0|1', 'releaseScenario': 'release1HS'}


    # call function to be tested
    cfg = dP.appendShpThickness(cfg)

    assert cfg['GENERAL']['relTh0'] == '1.2'
    assert cfg['GENERAL']['relTh1'] == '1.4'

    # call function to be tested with different inputs again
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'secRelArea': 'False', 'simTypeActual': 'null', 'relThFromShp': 'True', 'relTh': '',
        'relThFromFile': 'False', 'relThPercentVariation': '40$3', 'relThRangeVariation': '',
        'relThDistVariation': ''}
    cfg['INPUT'] = {'relThThickness': '1.2|1.4', 'relThId': '0|1', 'releaseScenario': 'release1HS'}
    cfg = dP.appendShpThickness(cfg)

    assert cfg.has_option('GENERAL', 'relTh0') is False
    assert cfg.has_option('GENERAL', 'relTh1') is False


def test_setThicknessValueFromVariation():
    """ test setting thickness from variation """

    # setup required inputs
    key = 'relThPercentVariation'
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'secRelArea': 'False', 'relThFromShp': 'True', 'relTh': '',
        'relThFromFile': 'False', 'relThPercentVariation': '40$3', 'relThDistVariation': ''}
    cfg['INPUT'] = {'relThThickness': '1.2|1.4', 'relThId': '0|1', 'releaseScenario': 'release1HS',
        'relThCi95': 'None|None'}

    data = {'relThPercentVariation': [0.6, 1.0, 1.4], 'relTh': ['', '', ''], 'ind': [0, 1, 2]}
    simDF = pd.DataFrame.from_dict(data)

    for row in simDF.itertuples():
        if row._asdict()['ind'] == 0:
            cfg = dP.setThicknessValueFromVariation(key, cfg, 'null', row)

            assert cfg['GENERAL']['relThPercentVariation'] == '-40.0$1'
            assert cfg['GENERAL']['relTh0'] == '0.72'
            assert cfg['GENERAL']['relTh1'] == '0.84'
        elif row._asdict()['ind'] == 1:
            cfg = dP.setThicknessValueFromVariation(key, cfg, 'null', row)

            assert cfg['GENERAL']['relThPercentVariation'] == ''
            assert cfg['GENERAL']['relTh0'] == '1.2'
            assert cfg['GENERAL']['relTh1'] == '1.4'

        elif row._asdict()['ind'] == 2:
            cfg = dP.setThicknessValueFromVariation(key, cfg, 'null', row)

            assert cfg['GENERAL']['relThPercentVariation'] == '+39.99999999999999$1'
            assert cfg['GENERAL']['relTh0'] == '1.68'
            assert cfg['GENERAL']['relTh1'] == '1.9599999999999997'


    cfg['GENERAL'] = {'secRelArea': 'False', 'relThFromShp': 'False', 'relTh': '1.',
        'relThFromFile': 'False', 'relThPercentVariation': '40$3', 'relThDistVariation': ''}
    cfg['INPUT'] = {'releaseScenario': 'release1HS'}

    data = {'relThPercentVariation': [0.6, 1.0, 1.4], 'relTh': ['1.', '1.', '1.'], 'ind': [0, 1, 2]}
    simDF = pd.DataFrame.from_dict(data)

    for row in simDF.itertuples():
        if row._asdict()['ind'] == 0:
            cfg['GENERAL']['relTh'] = '1.'
            cfg = dP.setThicknessValueFromVariation(key, cfg, 'null', row)

            assert cfg['GENERAL']['relThPercentVariation'] == ''
            assert cfg['GENERAL']['relTh'] == '0.6'

        elif row._asdict()['ind'] == 1:
            cfg['GENERAL']['relTh'] = '1.'
            cfg = dP.setThicknessValueFromVariation(key, cfg, 'null', row)

            assert cfg['GENERAL']['relThPercentVariation'] == ''
            assert cfg['GENERAL']['relTh'] == '1.0'

        elif row._asdict()['ind'] == 2:
            cfg['GENERAL']['relTh'] = '1.'
            cfg = dP.setThicknessValueFromVariation(key, cfg, 'null', row)

            assert cfg['GENERAL']['relThPercentVariation'] == ''
            assert cfg['GENERAL']['relTh'] == '1.4'


def test_splitVariationToArraySteps():
    """ testing splitting variation to array steps """

    # setup required inputs
    value = '40$5'
    key = 'relThPercentVariation'

    fullCfg = configparser.ConfigParser()
    fullCfg['GENERAL'] = {'addStandardConfig': 'True'}

    itemsArray = dP.splitVariationToArraySteps(value, key, fullCfg)

    assert np.array_equal(itemsArray, np.linspace(0.6, 1.4, 5))

    value = '40$4'
    itemsTest = np.append(np.linspace(0.6, 1.4, 4), np.array([1.]))
    itemsArray = dP.splitVariationToArraySteps(value, key, fullCfg)
    print('itemsTest', itemsTest)
    print('itemsArray', itemsArray)

    assert np.array_equal(itemsArray, itemsTest)

    value = '2$5'
    key = 'relThRangeVariation'
    fullCfg['GENERAL'] = {'addStandardConfig': 'False'}
    itemsArray = dP.splitVariationToArraySteps(value, key, fullCfg)

    assert np.array_equal(itemsArray, np.linspace(-2, 2, 5))

    value = '+2$5'
    key = 'relThRangeVariation'
    fullCfg['GENERAL'] = {'addStandardConfig': 'False'}
    itemsArray = dP.splitVariationToArraySteps(value, key, fullCfg)

    assert np.array_equal(itemsArray, np.linspace(0, 2, 5))

    value = '-2$5'
    key = 'relThRangeVariation'
    fullCfg['GENERAL'] = {'addStandardConfig': 'False'}
    itemsArray = dP.splitVariationToArraySteps(value, key, fullCfg)

    assert np.array_equal(itemsArray, np.linspace(-2, 0, 5))

    value = 'normaldistribution$3$0.1$95$ci95$10000'
    key = 'relThDistVariation'
    fullCfg['GENERAL'] = {'addStandardConfig': 'False'}
    itemsArray = dP.splitVariationToArraySteps(value, key, fullCfg)

    assert np.array_equal(itemsArray, ['0$normaldistribution$3$0.1$95$ci95$10000',
        '1$normaldistribution$3$0.1$95$ci95$10000', '2$normaldistribution$3$0.1$95$ci95$10000'])

    value = '4$normaldistribution$3$0.1$95$ci95$10000'
    key = 'relThDistVariation'
    fullCfg['GENERAL'] = {'addStandardConfig': 'False'}
    itemsArray = dP.splitVariationToArraySteps(value, key, fullCfg)

    assert np.array_equal(itemsArray, ['4$normaldistribution$3$0.1$95$ci95$10000'])
