"""
    Pytest for fileHandlerUtils

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
import pytest
import shutil
import configparser


def test_readLogFile():
    """ Test if logDict is generated correctly """

    # Test function
    dirPath = os.path.dirname(__file__)
    logName = os.path.join(dirPath, 'data', 'testExpLog.txt')
    cfg = configparser.ConfigParser()
    cfg = {'varPar' : 'RelTh'}
    logDict = fU.readLogFile(logName, cfg)
    logDictMu = fU.readLogFile(logName)

    assert logDict['noSim'][4] == 5
    assert logDict['simName'][2] == 'release1HS2_null_dfa'
    assert logDict['RelTh'][2] == 4.0
    assert logDictMu['noSim'][4] == 5
    assert logDictMu['simName'][2] == 'release1HS2_null_dfa'
    assert logDictMu['Mu'][2] == 4.0


def test_makeSimDict():
    """ Test if simulation dictionary is generated correctly """

    # Test function
    dirPath = os.path.dirname(__file__)
    inputDir = os.path.join(dirPath, 'data', 'testSim')
    cfg = configparser.ConfigParser()
    cfg = {'varPar' : 'test'}
    data = fU.makeSimDict(inputDir, cfg['varPar'])

    assert data['names'][0] == 'releaseTest1_entres_dfa_0.888_ppr'
    assert data['releaseArea'][0] == 'releaseTest1'
    assert data['simType'][0] == 'entres'
    assert data['resType'][0] == 'ppr'
    assert data['cellSize'][0] == 5.0
    assert data['test'][0] == '0.888'


def test_exportcom1DFAOutput(tmp_path):
    """ Test if export of result files works """

    # Create input directoy structure
    dirPath = os.path.dirname(__file__)
    avaName = 'avaParabola'
    avaNameTest = 'avaParabolaPyest'
    avaDir  = os.path.join(tmp_path, avaName)
    outDir = os.path.join(avaDir, 'Work', 'com1DFA', 'FullOutput_RelTh_1.25000', 'release1PF_entres_dfa', 'raster')
    os.makedirs(avaDir)
    os.makedirs(outDir)

    # copy inut data from benchmarks folder to tmp_path and rename correctly
    resType = ['ppr', 'pfd', 'pfv']
    for m in resType:
        if m == 'pfv':
            avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                               'release1PF_entres_dfa_1.25000_pfv.asc')
            input = os.path.join(avaDir, 'Work', 'com1DFA', 'FullOutput_RelTh_1.25000',
                                     'release1PF_entres_dfa', 'raster', 'release1PF_entres_dfa_pv.asc')
        else:
            avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                                'release1PF_entres_dfa_1.25000_%s.asc' % m)
            input = os.path.join(avaDir, 'Work', 'com1DFA', 'FullOutput_RelTh_1.25000',
                                 'release1PF_entres_dfa', 'raster', 'release1PF_entres_dfa_%s.asc' % m)
        shutil.copy(avaData, input)
    avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                           'ExpLog.txt')
    input = os.path.join(avaDir, 'Work', 'com1DFA', 'ExpLog.txt')
    shutil.copy(avaData, input)
    avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                           'test.html')
    input = os.path.join(avaDir, 'Work', 'com1DFA', 'FullOutput_RelTh_1.25000',
                             'release1PF_entres_dfa.html')
    shutil.copy(avaData, input)

    # Set cfg
    cfg = configparser.ConfigParser()
    cfg = {'varPar' : 'RelTh'}

    # Call function to test
    fU.exportcom1DFAOutput(avaDir, cfg)
    # load exported file
    pprTest = np.loadtxt(os.path.join(avaDir, 'Outputs', 'com1DFA', 'peakFiles',
                                         'release1PF_entres_dfa_1.25000_ppr.asc'), skiprows=6)

    # load initial file
    pprBench = np.loadtxt(os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                                       'release1PF_entres_dfa_1.25000_ppr.asc'), skiprows=6)
    # Compare result to reference solution
    testRes = np.allclose(pprTest, pprBench, atol=1.e-12)

    assert testRes == True


def test_splitIniValueToArraySteps():
    """ Test if splitting into an array works fine  """

    cfgValues = '1.0|2.5|3.8'
    cfgValuesList = np.asarray([1.0, 2.5, 3.8])

    cfgValues2 = '0:10:5'
    cfgValuesList2 = np.asarray([0., 2.5, 5., 7.5, 10.])

    # call function to be tested
    items = fU.splitIniValueToArraySteps(cfgValues)
    items2 = fU.splitIniValueToArraySteps(cfgValues2)


    assert len(items) == len(cfgValuesList)
    assert items[0] == cfgValuesList[0]
    assert items[1] == cfgValuesList[1]
    assert items[2] == cfgValuesList[2]
    assert len(items2) == len(cfgValuesList2)
    assert items2[0] == cfgValuesList2[0]
    assert items2[1] == cfgValuesList2[1]
    assert items2[2] == cfgValuesList2[2]
    assert items2[3] == cfgValuesList2[3]
    assert items2[4] == cfgValuesList2[4]


def test_splitTimeValueToArrayInterval():
    """ Test if splitting into an array works fine  """

    cfgValues = '1.0|2.5|3.8'
    cfgValuesList = np.asarray([1.0, 2.5, 3.8])

    cfgValues1 = '0.|2.5|3.8'
    cfgValuesList1 = np.asarray([2.5, 3.8])

    cfgValues2 = '0:5'
    cfgValuesList2 = np.asarray([5., 10., 15.])

    cfgValues3 = ''
    cfgValuesList3 = np.asarray([40.])

    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'tEnd': '20'}
    cfgGen = cfg['GENERAL']

    # call function to be tested
    cfgGen['tSteps'] = cfgValues
    items = fU.splitTimeValueToArrayInterval(cfgGen)
    cfgGen['tSteps'] = cfgValues1
    items1 = fU.splitTimeValueToArrayInterval(cfgGen)
    cfgGen['tSteps'] = cfgValues2
    items2 = fU.splitTimeValueToArrayInterval(cfgGen)
    cfgGen['tSteps'] = cfgValues3
    items3 = fU.splitTimeValueToArrayInterval(cfgGen)

    assert len(items) == len(cfgValuesList)
    assert items[0] == cfgValuesList[0]
    assert items[1] == cfgValuesList[1]
    assert items[2] == cfgValuesList[2]
    assert len(items1) == len(cfgValuesList1)
    assert items1[0] == cfgValuesList1[0]
    assert items1[1] == cfgValuesList1[1]
    assert len(items2) == len(cfgValuesList2)
    assert items2[0] == cfgValuesList2[0]
    assert items2[1] == cfgValuesList2[1]
    assert items2[2] == cfgValuesList2[2]
    assert len(items3) == len(cfgValuesList3)
    assert items3[0] == cfgValuesList3[0]
