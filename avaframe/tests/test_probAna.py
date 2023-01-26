"""
    Pytest for module ana4Stats

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.ana4Stats import probAna as pA
import pytest
import configparser
import shutil
import pathlib


def test_probAna(tmp_path):
    """ test probAna function to compute mask for parameter exceeding threshold """

    # set input directory
    avaName = 'avaParabola'
    avaTestDir = 'avaParabolaStatsTest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir
    avaDirtmp = pathlib.Path(tmp_path, avaName)
    inputDir = pathlib.Path(tmp_path, avaName)
    inputDir1 = avaDir
    shutil.copytree(inputDir1, inputDir)

    # set configurations
    testCfg = os.path.join(inputDir, '%sProbAna_com1DFACfg.ini' % avaName)
    cfgMain = cfgUtils.getModuleConfig(com1DFA, testCfg)

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'peakLim': 1.0, 'peakVar': 'ppr'}
    cfg['FILTER'] = {}

    # provide optional filter criteria for simulations
    parametersDict = fU.getFilterDict(cfg, 'FILTER')

    print('PARAMETERSDIC', parametersDict)

    # call function to test
    pA.probAnalysis(avaDirtmp, cfg, com1DFA, parametersDict=parametersDict, inputDir='')
    probTest = np.loadtxt(os.path.join(avaDirtmp, 'Outputs', 'ana4Stats', 'avaParabola_probMap1.0.asc'), skiprows=6)

    # Load reference solution
    probSol = np.loadtxt(os.path.join(inputDir1, 'avaParabola_probMap1.0.txt'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(probTest, probSol, atol=1.e-6)

    # Test
    assert (testRes is True)

    # call function to test
    testInputDir = avaDir / 'Outputs' / 'com1DFA'
    avaDirtmp2 = pathlib.Path(tmp_path, 'avaTest')
    avaDirtmp2.mkdir()
    pA.probAnalysis(avaDirtmp2, cfg, com1DFA, parametersDict='', inputDir=testInputDir)
    probTest2 = np.loadtxt(os.path.join(avaDirtmp2, 'Outputs', 'ana4Stats', 'avaTest_probMap1.0.asc'), skiprows=6)

    # Compare result to reference solution
    testRes2 = np.allclose(probTest2, probSol, atol=1.e-6)

    # Test
    assert (testRes2 is True)


def test_createComModConfig(tmp_path):
    """ test creatig a config file """

    # set input directory
    avaName = 'avaParabola'
    avaDir = pathlib.Path(tmp_path, avaName)
    dirPath = pathlib.Path(__file__).parents[0]
    cfgFile = dirPath / 'data' / 'testCom1DFA' / 'probA_com1DFACfg.ini'

    cfgProb = configparser.ConfigParser()
    cfgProb['PROBRUN'] = {'varParList': 'mu|relTh', 'variationType': 'percent',
                          'variationValue': '60|50', 'numberOfSteps': '2|3',
                          'addStandardConfig': 'True',
                          'defaultSetup': 'True', 'samplingStrategy': '2',
                          'defaultComModuleCfg': 'True'}
    cfgProb['sampling_override'] = {'defaultConfig': 'True'}

    # call function to be tested
    cfgFiles, outDir = pA.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod='')

    print('cfgFiles', cfgFiles)

    # load cfg from file
    cfgMu = configparser.ConfigParser()
    cfgRelTh = configparser.ConfigParser()

    cfgMu.read(cfgFiles[0])
    cfgRelTh.read(cfgFiles[1])

    print(cfgMu['GENERAL']['mu'], cfgMu['GENERAL']['relTh'], cfgRelTh['GENERAL']['mu'], cfgRelTh['GENERAL']['relTh'])

    assert cfgMu['GENERAL']['mu'] == '0.155$60$2&0.155'
    assert cfgMu['GENERAL']['relTh'] == ''
    assert cfgRelTh['GENERAL']['mu'] == '0.155'
    assert cfgRelTh['GENERAL']['relTh'] == ''
    assert cfgRelTh['GENERAL']['relThFromShp'] == 'True'
    assert cfgRelTh['GENERAL']['relThPercentVariation'] == '50$3'
    assert cfgRelTh['GENERAL']['relThFromShp'] == 'True'
    assert cfgMu['GENERAL']['relThFromShp'] == 'True'
    assert cfgRelTh['GENERAL']['addStandardConfig'] == 'True'


    cfgProb['PROBRUN']['addStandardConfig'] = 'False'
    cfgProb['PROBRUN']['defaultComModuleCfg'] = 'False'
    # call function to be tested
    cfgFiles, outDir = pA.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod=cfgFile)

    # load cfg from file
    cfgMu = configparser.ConfigParser()
    cfgRelTh = configparser.ConfigParser()

    cfgMu.read(cfgFiles[0])
    cfgRelTh.read(cfgFiles[1])

    print(cfgMu['GENERAL']['mu'], cfgMu['GENERAL']['relTh'], cfgRelTh['GENERAL']['mu'], cfgRelTh['GENERAL']['relTh'])

    assert cfgMu['GENERAL']['mu'] == '0.155$60$2'
    assert cfgMu['GENERAL']['relTh'] == '2.'
    assert cfgRelTh['GENERAL']['mu'] == '0.155'
    assert cfgRelTh['GENERAL']['relTh'] == '2.'
    assert cfgRelTh['GENERAL']['relThFromShp'] == 'False'
    assert cfgRelTh['GENERAL']['relThPercentVariation'] == '50$3'
    assert cfgRelTh['GENERAL']['relThFromShp'] == 'False'
    assert cfgMu['GENERAL']['relThFromShp'] == 'False'
    assert cfgRelTh['GENERAL']['addStandardConfig'] == 'False'


    cfgProb = configparser.ConfigParser()
    cfgProb.optionxform = str
    cfgProb['PROBRUN'] = {'varParList': 'mu|relTh', 'variationType': 'percent',
                          'variationValue': '60|50', 'numberOfSteps': '2|3',
                          'addStandardConfig': 'True',
                          'defaultSetup': 'True', 'samplingStrategy': '1',
                          'defaultComModuleCfg': 'False',
                          'varParType': 'float|float', 'nSample': '40', 'sampleSeed': '12345',
                          'sampleMethod': 'latin'}


    # call function to be tested
    cfgFiles, outDir = pA.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod=cfgFile)

    print('cfgFiles', cfgFiles)

    cfgTest = configparser.ConfigParser()
    cfgTest.read(cfgFiles[0])
    print('cfgTest', cfgTest['GENERAL']['relThFromShp'], cfgTest['GENERAL']['relTh'],
        cfgTest['GENERAL']['relThPercentVariation'], cfgTest['GENERAL']['mu'])

    assert cfgTest['GENERAL']['relThFromShp'] == 'False'
    assert cfgTest['GENERAL']['relTh'] == '1.6341620830145123'
    assert len(cfgFiles) == 40

    cfgProb = configparser.ConfigParser()
    cfgProb.optionxform = str
    cfgProb['PROBRUN'] = {'varParList': 'mu|relTh', 'variationType': 'percent',
                          'variationValue': '60|50', 'numberOfSteps': '2|3',
                          'addStandardConfig': 'True',
                          'defaultSetup': 'True', 'samplingStrategy': '1', 'defaultComModuleCfg': 'True',
                          'varParType': 'float|float', 'nSample': '40', 'sampleSeed': '12345',
                          'sampleMethod': 'latin'}


    # call function to be tested
    cfgFiles2, outDir = pA.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod='')

    print('cfgFiles', cfgFiles2)

    cfgTest2 = configparser.ConfigParser()
    cfgTest2.read(cfgFiles2[0])
    print('cfgTest', cfgTest['GENERAL']['relThFromShp'], cfgTest['GENERAL']['relTh'],
        cfgTest['GENERAL']['relThPercentVariation'], cfgTest['GENERAL']['mu'])

    assert cfgTest2['GENERAL']['relThFromShp'] == 'True'
    assert cfgTest2['GENERAL']['relTh'] == ''
    assert cfgTest2['GENERAL']['relThPercentVariation'] == '-18.291895849274383$1'
    assert len(cfgFiles) == 40

    cfgProb = configparser.ConfigParser()
    cfgProb.optionxform = str
    cfgProb['PROBRUN'] = {'varParList': 'mu|relTh', 'variationType': 'range',
                          'variationValue': '0.2|1.2', 'numberOfSteps': '2|3',
                          'addStandardConfig': 'False',
                          'samplingStrategy': '2',
                          'defaultComModuleCfg': 'True'}

    # call function to be tested
    cfgFiles3, outDir = pA.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod='')

    print('cfgFiles', cfgFiles3)

    cfgTest3 = configparser.ConfigParser()
    cfgTest3.read(cfgFiles3[1])

    assert cfgTest3['GENERAL']['relThFromShp'] == 'True'
    assert cfgTest3['GENERAL']['relTh'] == ''
    assert cfgTest3['GENERAL']['relThRangeVariation'] == '1.2$3'
    assert len(cfgFiles3) == 2


def test_updateCfgRange():
    """ test updating cfg values """

    # setup inputs
    cfg = configparser.ConfigParser()
    cfg['PROBRUN'] = {'varParList': 'mu|relTh', 'variationType': 'percent','addStandardConfig': 'True',
                      'variationValue': '60|50', 'numberOfSteps': '2|2',
                      'defaultSetup': 'True'}

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'mu'

    varDict = pA.makeDictFromVars(cfg['PROBRUN'])
    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    assert cfgNew['GENERAL']['mu'] == '0.155$60$2&0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == ''


    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'relTh'

    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    assert cfgNew['GENERAL']['mu'] == '0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == '50$2'
    assert cfgNew['GENERAL']['addStandardConfig'] == 'True'


    cfg['PROBRUN']['addStandardConfig'] = 'False'

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'mu'
    varDict = {}
    varDict = pA.makeDictFromVars(cfg['PROBRUN'])

    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    assert cfgNew['GENERAL']['mu'] == '0.155$60$2'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == ''


    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'relTh'

    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    assert cfgNew['GENERAL']['mu'] == '0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == '50$2'
    assert cfgNew['GENERAL']['addStandardConfig'] == 'False'

    # setup inputs
    cfg = configparser.ConfigParser()
    cfg['PROBRUN'] = {'varParList': 'mu|relTh', 'variationType': 'range','addStandardConfig': 'False',
                      'variationValue': '0.1|0.5', 'numberOfSteps': '5|12',
                      'defaultSetup': 'True'}

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'mu'

    varDict = pA.makeDictFromVars(cfg['PROBRUN'])
    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    assert np.isclose(float(cfgNew['GENERAL']['mu'].split(':')[0]), 0.055, rtol=1.e-5)
    assert np.isclose(float(cfgNew['GENERAL']['mu'].split(':')[1]), 0.255, rtol=1.e-5)
    assert cfgNew['GENERAL']['mu'].split(':')[2] == '5'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == ''


    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'relTh'
    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    assert cfgNew['GENERAL']['mu'] == '0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThRangeVariation'] == '0.5$12'


    # setup inputs
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg['PROBRUN'] = {'varParList': 'mu|relTh', 'variationType': 'normaldistribution',
                      'addStandardConfig': 'False',
                      'variationValue': '0.1|0.3', 'numberOfSteps': '3|12',
                      'defaultSetup': 'True'}
    cfg['computeFromDistribution_override'] = {'buildType': 'ci95', 'minMaxInterval': '95',
        'defaultConfig': 'True'}

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'mu'

    varDict = pA.makeDictFromVars(cfg['PROBRUN'])
    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    assert np.isclose(float(cfgNew['GENERAL']['mu'].split('|')[0]), 0.055, rtol=1.e-3)
    assert np.isclose(float(cfgNew['GENERAL']['mu'].split('|')[1]), 0.155, rtol=1.e-3)
    assert np.isclose(float(cfgNew['GENERAL']['mu'].split('|')[2]), 0.255, rtol=1.e-3)
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThDistVariation'] == ''

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'relTh'

    varDict = pA.makeDictFromVars(cfg['PROBRUN'])
    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName, varDict[varName])

    print('value', cfgNew['GENERAL']['relThPercentVariation'])

    assert cfgNew['GENERAL']['mu'] == '0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThDistVariation'] == 'normaldistribution$12$0.3$95$ci95$10000'
