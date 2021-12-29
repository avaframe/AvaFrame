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
    assert (testRes == True)

    # call function to test
    testInputDir = avaDir / 'Outputs' / 'com1DFA'
    avaDirtmp2 = pathlib.Path(tmp_path, 'avaTest')
    avaDirtmp2.mkdir()
    pA.probAnalysis(avaDirtmp2, cfg, com1DFA, parametersDict='', inputDir=testInputDir)
    probTest2 = np.loadtxt(os.path.join(avaDirtmp2, 'Outputs', 'ana4Stats', 'avaTest_probMap1.0.asc'), skiprows=6)

    # Compare result to reference solution
    testRes2 = np.allclose(probTest2, probSol, atol=1.e-6)

    # Test
    assert (testRes2 == True)


def test_createComModConfig(tmp_path):
    """ test creatig a config file """

    # set input directory
    avaName = 'avaParabola'
    avaDir = pathlib.Path(tmp_path, avaName)

    cfgProb = configparser.ConfigParser()
    cfgProb['PROBRUN'] = {'varParList': 'mu|relTh', 'muRange': '0.055', 'muSteps': '2',
        'relThRange': '0.5', 'relThSteps': '3', 'defaultSetup': 'True'}

    # call function to be tested
    cfgFiles = pA.createComModConfig(cfgProb, avaDir, com1DFA)

    # load cfg from file
    cfgMu = configparser.ConfigParser()
    cfgRelTh = configparser.ConfigParser()

    cfgMu.read(cfgFiles['mu']['cfgFile'])
    cfgRelTh.read(cfgFiles['relTh']['cfgFile'])

    print(cfgMu['GENERAL']['mu'], cfgMu['GENERAL']['relTh'], cfgRelTh['GENERAL']['mu'], cfgRelTh['GENERAL']['relTh'])

    assert cfgMu['GENERAL']['mu'] == '0.1:0.21:2&0.155'
    assert cfgMu['GENERAL']['relTh'] == '1.'
    assert cfgRelTh['GENERAL']['mu'] == '0.15500'
    assert cfgRelTh['GENERAL']['relTh'] == '0.5:1.5:3'
    assert cfgRelTh['GENERAL']['useRelThFromIni'] == 'True'
    assert cfgMu['GENERAL']['useRelThFromIni'] == 'True'


def test_updateCfgRange():
    """ test updating cfg values """

    # setup inputs
    cfg = configparser.ConfigParser()
    cfg['PROBRUN'] = {'varParList': 'mu|relTh', 'muRange': '0.055', 'muSteps': '2',
        'relThRange': '0.5', 'relThSteps': '3', 'defaultSetup': 'True'}

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'mu'

    # call function
    cfgNew, refIn = pA.updateCfgRange(com1DFACfg, cfg, varName)

    assert refIn == True
    assert cfgNew['GENERAL']['mu'] == '0.1:0.21:2&0.155'
    assert cfgNew['GENERAL']['relTh'] == '1.'

    com1DFACfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    varName = 'relTh'

    # call function
    cfgNew, refIn = pA.updateCfgRange(com1DFACfg, cfg, varName)

    assert refIn == False
    assert cfgNew['GENERAL']['mu'] == '0.15500'
    assert cfgNew['GENERAL']['relTh'] == '0.5:1.5:3'
