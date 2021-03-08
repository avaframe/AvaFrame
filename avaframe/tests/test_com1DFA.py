"""
    Pytest for module com1DFA

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.ana1Tests import testUtilities as tU
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import initializeProject as initProj
from benchmarks import simParameters
import pytest
import configparser
import shutil


def test_execCom1Exe(tmp_path):
    """ test call to com1DFA executable without performing simulation but check log generation """

    # Initialise inputs
    com1Exe = 'testPath'
    dirPath = os.path.dirname(__file__)
    avaName = 'avaParabola'
    avaDir  = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)
    cintFile = os.path.join(dirPath, 'data', 'runBasicST.cint')
    logName = 'release1HS_entres_dfa_0.155'
    fullOut = False

    # create required Output folder
    os.makedirs(os.path.join(avaDir, 'Outputs', 'com1DFA'))

    # Call function
    com1DFA.execCom1Exe(com1Exe, cintFile, avaDir, fullOut, logName)

    # reference log File name
    logFile = 'startrelease1HS_entres_dfa_0.155.log'

    # check if log file has been created with correct name
    flagFile = False
    if os.path.isfile(os.path.join(avaDir, 'Outputs', 'com1DFA', logFile)):
        flagFile = True

    # Test
    assert flagFile == True


@pytest.mark.skip(reason="com1DFA exe is missing, no way of testing this")
@pytest.mark.parametrize("avaName", ["avaBowl", "avaFlatPlane", "avaHelix", "avaHelixChannel",
                          "avaParabola", "avaHockeyChannel", "avaHockeySmall", "avaInclinedPlane"])
def test_com1DFAMain(tmp_path, avaName):
    """ test call to com1DFA module """

    # get path to executable
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
    com1Exe = cfgCom1DFA['GENERAL']['com1Exe']

    # get input data
    dirPath = os.path.dirname(__file__)
    avaDir  = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    # load all benchmark info as dictionaries from description files
    testDictList = tU.readAllBenchmarkDesDicts(info=False)
    type = 'TAGS'
    valuesList = ['standardTest']
    testList = tU.filterBenchmarks(testDictList, type, valuesList)
    type = 'AVANAME'
    valuesList = avaName
    testList = tU.filterBenchmarks(testList, type, valuesList)
    avaTestName = testList[0]['NAME']
    print('avaTEstNAme', avaTestName)

    testCfg = os.path.join(dirPath, '..', '..', 'benchmarks', avaTestName, '%s_com1DFACfg.ini' % avaName)
    shutil.copytree(avaData, avaInputs)

    # get configuration
    cfg = cfgUtils.getModuleConfig(com1DFA, testCfg)
    cfg['GENERAL']['com1Exe'] = com1Exe

    # Run Standalone DFA
    reportDictList = com1DFA.com1DFAMain(cfg, avaDir)
    reportD = reportDictList[0]
    pprCom1DFA = np.loadtxt(os.path.join(avaDir, 'Outputs', 'com1DFA', 'peakFiles',
                                         reportD['simName']['name'] + '_ppr.asc'), skiprows=6)

    # Fetch simulation info from benchmark results
    benchDictList = simParameters.fetchBenchParameters(avaDir)
    benchDict = benchDictList[0]
    simBench = benchDict['simName'].replace('dfa', 'ref')
    pprBench = np.loadtxt(os.path.join(dirPath, '..', '..', 'benchmarks', avaTestName,
                                       simBench + '_ppr.asc'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(pprCom1DFA, pprBench, atol=1.e-12)

    # Test module
    assert reportD['simName']['name'] == benchDict['simName']
    assert reportD['Simulation Parameters']['Release Area'] == benchDict['Simulation Parameters']['Release Area']
    assert reportD['Simulation Parameters']['Entrainment Area'] == benchDict['Simulation Parameters']['Entrainment Area']
    assert reportD['Simulation Parameters']['Resistance Area'] == benchDict['Simulation Parameters']['Resistance Area']
    assert reportD['Simulation Parameters']['Mu'] == benchDict['Simulation Parameters']['Mu']
    assert reportD['Simulation Parameters']['Parameter value'] == ''
    assert reportD['Simulation Parameters']['Release thickness [m]'] == benchDict['Simulation Parameters']['Release thickness [m]']
    assert testRes == True
