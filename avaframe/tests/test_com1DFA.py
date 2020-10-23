"""
    Pytest for module com1DFA

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.out3SimpPlot import outPlotAllPeak as oP
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import initializeProject as initProj
import pytest
import configparser
import shutil


def test_initialiseRun(tmp_path):
    """ test check for input data """

    # get input data
    dirPath = os.path.dirname(__file__)
    avaName = 'avaHockeySmoothChannel'
    avaDir  = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)

    # flags for entrainment and resistance
    flagEnt = True
    flagRes = True

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['PARAMETERVAR'] = {'flagVarPar': 'True', 'varPar': 'RelTh'}
    cfgPar = cfg['PARAMETERVAR']

    # call function to be tested
    dem, rels, ent, res, flagEntRes = com1DFA.initialiseRun(avaDir, flagEnt, flagRes, cfgPar)

    # open ExpLog
    # Load simulation report
    testDir = os.path.join(avaDir, 'Work', 'com1DFA')
    logFile = open(os.path.join(testDir, 'ExpLog.txt'), 'r')
    lines = logFile.readlines()
    lineVals = []
    for line in lines:
        lineVals.append(line.split(','))

    # Test
    assert dem == os.path.join(avaDir, 'Inputs', 'myDEM_HS2_Topo.asc')
    assert rels == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS2.shp')]
    assert res == ''
    assert ent == os.path.join(avaDir, 'Inputs', 'ENT', 'entrainment1HS2.shp')
    assert flagEntRes == True
    assert lineVals[0][2] == 'RelTh\n'


def test_execCom1Exe(tmp_path):
    """ test call to com1DFA executable without performing simulation but check log generation """

    # Initialise inputs
    com1Exe = 'testPath'
    dirPath = os.path.dirname(__file__)
    avaName = 'avaHockey'
    avaDir  = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)
    cintFile = os.path.join(dirPath, 'data', 'runBasicST.cint')
    logName = 'release1HS2_entres_dfa_0.155'
    fullOut = False

    # create required Output folder
    os.makedirs(os.path.join(avaDir, 'Outputs', 'com1DFA'))

    # Call function
    com1DFA.execCom1Exe(com1Exe, cintFile, avaDir, fullOut, logName)

    # reference log File name
    logFile = 'startrelease1HS2_entres_dfa_0.155.log'

    # check if log file has been created with correct name
    flagFile = False
    if os.path.isfile(os.path.join(avaDir, 'Outputs', 'com1DFA', logFile)):
        flagFile = True

    # Test
    assert flagFile == True


@pytest.mark.skip(reason="com1DFA exe is missing, no way of testing this")
def test_com1DFAMain(tmp_path):
    """ test call to com1DFA module """

    # get path to executable
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
    com1Exe = cfgCom1DFA['GENERAL']['com1Exe']

    # get input data
    dirPath = os.path.dirname(__file__)
    avaName = 'avaHockeySmoothChannel'
    avaDir  = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    testCfg = os.path.join(dirPath, 'data', '%s_com1DFACfg.ini' % avaName)
    shutil.copytree(avaData, avaInputs)

    # get configuration
    cfg = cfgUtils.getModuleConfig(com1DFA, testCfg)
    cfg['GENERAL']['com1Exe'] = com1Exe

    # Run Standalone DFA
    reportDictList = com1DFA.com1DFAMain(cfg, avaDir)
    reportD = reportDictList[0]

    # Test module
    assert reportD['simName']['name'] == 'release1HS2_null_dfa_0.155'
    assert reportD['Simulation Parameters']['Entrainment Area'] == ''
    assert reportD['Simulation Parameters']['Mu'] == '0.155'
    assert reportD['Simulation Parameters']['Parameter value'] == ''
