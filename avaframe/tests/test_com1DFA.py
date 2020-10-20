"""
    Pytest for module com1DFA

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.com1DFA import com1DFA
import pytest
import configparser


def test_initialiseRun():
    """ test check for input data """

    avaDir  = os.path.join('avaframe', 'tests', 'data', 'avaTest')
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
    testDir = os.path.join(os.getcwd(), avaDir, 'Work', 'com1DFA')
    logFile = open(os.path.join(testDir, 'ExpLog.txt'), 'r')
    lines = logFile.readlines()
    lineVals = []
    for line in lines:
        lineVals.append(line.split(','))

    # Test
    assert dem == os.path.join(avaDir, 'Inputs', 'myDEM_IP_Topo.asc')
    assert rels == [os.path.join(avaDir, 'Inputs', 'REL', 'release1IP.shp')]
    assert res == ''
    assert ent == os.path.join(avaDir, 'Inputs', 'ENT', 'entrainement1IP.shp')
    assert flagEntRes == True
    assert lineVals[0][2] == 'RelTh\n'


def test_execSamos():
    """ test call to samosAT executable """

    # Initialise inputs
    samosAT = 'testPath'
    cintFile = os.path.join('avaframe', 'tests', 'data', 'runBasicST.cint')
    avaDir  = os.path.join('avaframe', 'tests', 'data', 'avaTest')
    simName = 'release1IP'
    fullOut = False

    # Call function
    com1DFA.execSamos(samosAT, cintFile, avaDir, fullOut, simName)

    # reference log File name
    logName = 'startrelease1IP.log'

    # check if log file has been created with correct name
    flagFile = False
    if os.path.isfile(os.path.join(avaDir, 'Outputs', 'com1DFA', logName)):
        flagFile = True

    # Test
    assert flagFile == True
