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
    cfgProb['PROBRUN'] = {'varParList': 'mu|relTh', 'percentVariation': 'True',
                          'muVariation': '60', 'muSteps': '2', 'addStandardConfig': 'True',
                          'relThVariation': '50', 'relThSteps': '3', 'defaultSetup': 'True'}

    # call function to be tested
    cfgFiles = pA.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod=cfgFile)

    # load cfg from file
    cfgMu = configparser.ConfigParser()
    cfgRelTh = configparser.ConfigParser()

    cfgMu.read(cfgFiles['mu']['cfgFile'])
    cfgRelTh.read(cfgFiles['relTh']['cfgFile'])

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
    # call function to be tested
    cfgFiles = pA.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod=cfgFile)

    # load cfg from file
    cfgMu = configparser.ConfigParser()
    cfgRelTh = configparser.ConfigParser()

    cfgMu.read(cfgFiles['mu']['cfgFile'])
    cfgRelTh.read(cfgFiles['relTh']['cfgFile'])

    print(cfgMu['GENERAL']['mu'], cfgMu['GENERAL']['relTh'], cfgRelTh['GENERAL']['mu'], cfgRelTh['GENERAL']['relTh'])

    assert cfgMu['GENERAL']['mu'] == '0.155$60$2'
    assert cfgMu['GENERAL']['relTh'] == ''
    assert cfgRelTh['GENERAL']['mu'] == '0.155'
    assert cfgRelTh['GENERAL']['relTh'] == ''
    assert cfgRelTh['GENERAL']['relThFromShp'] == 'True'
    assert cfgRelTh['GENERAL']['relThPercentVariation'] == '50$3'
    assert cfgRelTh['GENERAL']['relThFromShp'] == 'True'
    assert cfgMu['GENERAL']['relThFromShp'] == 'True'
    assert cfgRelTh['GENERAL']['addStandardConfig'] == 'False'


def test_updateCfgRange():
    """ test updating cfg values """

    # setup inputs
    cfg = configparser.ConfigParser()
    cfg['PROBRUN'] = {'varParList': 'mu|relTh', 'percentVariation': 'True','addStandardConfig': 'True',
                      'muVariation': '60', 'muSteps': '2', 'relThVariation': '50',
                      'relThSteps': '2', 'defaultSetup': 'True'}

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'mu'

    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName)

    assert cfgNew['GENERAL']['mu'] == '0.155$60$2&0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == ''


    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'relTh'

    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName)

    assert cfgNew['GENERAL']['mu'] == '0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == '50$2'
    assert cfgNew['GENERAL']['addStandardConfig'] == 'True'


    cfg['PROBRUN']['addStandardConfig'] = 'False'

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'mu'

    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName)

    assert cfgNew['GENERAL']['mu'] == '0.155$60$2'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == ''


    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'relTh'

    # call function
    cfgNew = pA.updateCfgRange(com1DFACfg, cfg, varName)

    assert cfgNew['GENERAL']['mu'] == '0.155'
    assert cfgNew['GENERAL']['relTh'] == ''
    assert cfgNew['GENERAL']['relThFromShp'] == 'True'
    assert cfgNew['GENERAL']['relThPercentVariation'] == '50$2'
    assert cfgNew['GENERAL']['addStandardConfig'] == 'False'
