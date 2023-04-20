"""
    Pytest for module ana4Stats

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
import pathlib
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.ana4Stats import getStats
import pytest
import configparser
import shutil


def test_getStats(tmp_path):
    """ test get statistics of data """

    # define test directory
    avaName = 'avaStats'
    avaDirtmp = pathlib.Path(tmp_path, avaName)
    avaDirPeakFiles = avaDirtmp / 'Outputs' / 'com1DFA' / 'peakFiles'
    avaDirConfigFiles = avaDirtmp / 'Outputs' / 'com1DFA' / 'configurationFiles'
    os.makedirs(avaDirtmp)
    os.makedirs(avaDirPeakFiles)
    os.makedirs(avaDirConfigFiles)

    # set path to required input files
    dirPath = pathlib.Path(__file__).parents[0]
    testDataDir = pathlib.Path(dirPath, 'data')
    data1 = testDataDir / 'testGetStats_1000000_ppr.asc'
    data2 = testDataDir / 'testGetStats_2000000_ppr.asc'
    cfg1 = testDataDir / 'testGetStats_1000000.ini'
    cfg2 = testDataDir / 'testGetStats_2000000.ini'

    # define destination paths
    fNames = [str(avaDirPeakFiles / 'release1_null_dfa_1000000_ppr.asc'),
              str(avaDirPeakFiles / 'release1_null_dfa_2000000_ppr.asc')]

    configFiles = [str(avaDirConfigFiles / 'release1_null_dfa_1000000.ini'),
                   str(avaDirConfigFiles / 'release1_null_dfa_2000000.ini')]

    # copy files as for extractMaxValues avaDir folder structure is needed
    shutil.copy(data1, fNames[0])
    shutil.copy(data2, fNames[1])
    shutil.copy(cfg1, configFiles[0])
    shutil.copy(cfg2, configFiles[1])

    # parameter dictionary
    varPar = 'relTh'
    inputDir = avaDirPeakFiles
    peakValues = getStats.extractMaxValues(inputDir, avaDirtmp, varPar, restrictType='ppr', nameScenario='', parametersDict='')

    # call function to be tested test 2
    parametersDict = {'relTh': 1.0}
    peakValues2 = getStats.extractMaxValues(inputDir, avaDirtmp, varPar, restrictType='', nameScenario='releaseScenario', parametersDict=parametersDict)

    assert peakValues['release1_null_dfa_1000000']['varPar'] == 1.0
    assert peakValues['release1_null_dfa_2000000']['varPar'] == 2.0
    assert peakValues['release1_null_dfa_1000000']['ppr']['max'] == 4.0
    assert peakValues['release1_null_dfa_1000000']['ppr']['mean'] == 2.0
    assert peakValues['release1_null_dfa_1000000']['ppr']['min'] == 1.0
    assert peakValues['release1_null_dfa_2000000']['ppr']['max'] == 10.0
    assert peakValues['release1_null_dfa_2000000']['ppr']['min'] == 1.0
    assert peakValues['release1_null_dfa_2000000']['ppr']['mean'] == 3.0
    assert peakValues2['release1_null_dfa_1000000']['scenario'] == 'release1'
