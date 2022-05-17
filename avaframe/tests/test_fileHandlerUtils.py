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
import pathlib
import configparser


def test_makeADir(tmp_path):
    """ test make directory """

    # create temporary directory
    avaName = 'testAva'
    avaDir = os.path.join(tmp_path, avaName)
    fU.makeADir(avaDir)
    avaDir2 = os.path.join(tmp_path, avaName, 'test')

    dirTrue = os.path.isdir(avaDir)
    dirFalse = os.path.isdir(avaDir2)

    assert dirTrue
    assert dirFalse is False


def test_readLogFile():
    """ Test if logDict is generated correctly """

    # Test function
    dirPath = os.path.dirname(__file__)
    logName = os.path.join(dirPath, 'data', 'testExpLog.txt')
    cfg = configparser.ConfigParser()
    cfg = {'varPar': 'RelTh'}
    logDict = fU.readLogFile(logName, cfg)
    logDictMu = fU.readLogFile(logName)

    assert logDict['noSim'][4] == 5
    assert logDict['simName'][2] == 'release1HS2_null_dfa'
    assert logDict['RelTh'][2] == 4.0
    assert logDictMu['noSim'][4] == 5
    assert logDictMu['simName'][2] == 'release1HS2_null_dfa'
    assert logDictMu['Mu'][2] == 4.0


def test_extractLogInfo():
    """ test extracting info from logFile """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    logName = dirPath / 'data' / 'logTest.log'

    # call function to be tested
    logDict = fU.extractLogInfo(logName)

    print('logDict', logDict)
    # define test results
    time = np.asarray([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 301.1])
    mass = np.asarray([1.99393e+07, 1.99393e+07, 1.99393e+07, 1.99393e+07, 1.99393e+07, 1.99393e+07,
                       2.02876e+07])
    entrMass = np.asarray([0., 0., 0., 0., 0., 0., 0.])
    stopTime = 78.4059
    stopCrit = 'kinetic energy 1.00 of peak KE'

    assert logDict['indRun'] == [0, 7]
    assert np.array_equal(logDict['time'], time)
    assert np.array_equal(logDict['mass'], mass)
    assert np.array_equal(logDict['entrMass'], entrMass)
    assert logDict['stopTime'] == stopTime
    assert logDict['stopCrit'] == stopCrit


def test_checkIfFileExists():
    """ test if a file exists and if not throw error """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    testPath = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathData = testPath / 'Outputs' / 'com1DFAOri' / 'peakFiles' / 'release1HS_ent_dfa_67dc2dc10a_pft.asc'

    # call function to be tested
    with pytest.raises(FileNotFoundError) as e:
        assert fU.checkIfFileExists(pathData, fileType='')
    assert str(e.value) == ('No  file found called: %s' % str(pathData))

    # call function to be tested
    with pytest.raises(FileNotFoundError) as e:
        assert fU.checkIfFileExists(pathData, fileType='log info')
    assert str(e.value) == ('No log info file found called: %s' % str(pathData))

    # call function to be tested
    pathData2 = 'test/dataTest'
    print('pathDatastr', pathData2)
    with pytest.raises(FileNotFoundError) as e:
        assert fU.checkIfFileExists(pathData, fileType='log info')
    assert str(e.value) == ('No log info file found called: %s' % str(pathData))


def test_makeSimDF():
    """ Test if simulation dataFrame is generated correctly """

    # Test function
    dirPath = os.path.dirname(__file__)
    inputDir = os.path.join(dirPath, 'data', 'testSim')
    cfg = configparser.ConfigParser()
    cfg = {'varPar': 'test'}
    dataDF = fU.makeSimDF(inputDir, simID=cfg['varPar'])

    assert dataDF['names'][0] == 'releaseTest1_0.888_entres_dfa_ppr'
    assert dataDF['releaseArea'][0] == 'releaseTest1'
    assert dataDF['simType'][0] == 'entres'
    assert dataDF['resType'][0] == 'ppr'
    assert dataDF['cellSize'][0] == 5.0
    assert dataDF['test'][0] == '0.888'

    inputDir = os.path.join(dirPath, 'data', 'testSim1')
    dataDF = fU.makeSimDF(inputDir, simID=cfg['varPar'])
    assert dataDF['names'][0] == 'releaseTest1_test_AF_0.888_entres_dfa_ppr'
    assert dataDF['releaseArea'][0] == 'releaseTest1_test'
    assert dataDF['simType'][0] == 'entres'
    assert dataDF['resType'][0] == 'ppr'
    assert dataDF['cellSize'][0] == 5.0
    assert dataDF['test'][0] == '0.888'


def test_makeSimFromResDF():
    """ Test if simulation dataFrame is generated correctly """

    # Test function
    dirPath = os.path.dirname(__file__)
    inputDir = os.path.join(dirPath, 'data', 'testSim')
    dataDF, resTypeList = fU.makeSimFromResDF(inputDir, 'comModule', inputDir=inputDir)
    print(resTypeList)
    print(dataDF.columns)
    print(dataDF.index)
    assert dataDF['simName'][0] == 'releaseTest1_0.888_entres_dfa'
    assert dataDF['releaseArea'][0] == 'releaseTest1'
    assert dataDF['simType'][0] == 'entres'
    assert dataDF['cellSize'][0] == 5.0

    inputDir = os.path.join(dirPath, 'data', 'testSim1')
    dataDF, resTypeList = fU.makeSimFromResDF(inputDir, 'comModule', inputDir=inputDir)
    assert dataDF['simName'][0] == 'releaseTest1_test_AF_0.888_entres_dfa'
    assert dataDF['releaseArea'][0] == 'releaseTest1_test'
    assert dataDF['simType'][0] == 'entres'
    assert dataDF['cellSize'][0] == 5.0


def test_exportcom1DFAOrigOutput(tmp_path):
    """ Test if export of result files works """

    # Create input directoy structure
    dirPath = os.path.dirname(__file__)
    avaName = 'avaParabola'
    avaNameTest = 'avaParabolaPyest'
    avaDir = os.path.join(tmp_path, avaName)
    outDir = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'FullOutput_RelTh_1.25000', 'release1PF_entres_dfa', 'raster')
    os.makedirs(avaDir)
    os.makedirs(outDir)

    # copy inut data from benchmarks folder to tmp_path and rename correctly
    resType = ['ppr', 'pft', 'pfv']
    for m in resType:
        if m == 'pfv':
            avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                               'release1PF_entres_dfa_1.25000_pfv.asc')
            input = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'FullOutput_RelTh_1.25000',
                                'release1PF_entres_dfa', 'raster', 'release1PF_entres_dfa_pv.asc')
        elif m == 'pft':
            avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                               'release1PF_entres_dfa_1.25000_pft.asc')
            input = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'FullOutput_RelTh_1.25000',
                                'release1PF_entres_dfa', 'raster', 'release1PF_entres_dfa_pfd.asc')
        else:
            avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                                'release1PF_entres_dfa_1.25000_%s.asc' % m)
            input = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'FullOutput_RelTh_1.25000',
                                 'release1PF_entres_dfa', 'raster', 'release1PF_entres_dfa_%s.asc' % m)
        shutil.copy(avaData, input)
    avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                           'ExpLog.txt')
    input = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'ExpLog.txt')
    shutil.copy(avaData, input)
    avaData = os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                           'test.html')
    input = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'FullOutput_RelTh_1.25000',
                        'release1PF_entres_dfa.html')
    shutil.copy(avaData, input)

    # Set cfg
    cfg = configparser.ConfigParser()
    cfg = {'varPar': 'RelTh'}

    # Call function to test
    fU.exportcom1DFAOrigOutput(avaDir, cfg)
    # load exported file
    pprTest = np.loadtxt(os.path.join(avaDir, 'Outputs', 'com1DFAOrig', 'peakFiles',
                         'release1PF_entres_dfa_1.25000_ppr.asc'), skiprows=6)

    # load initial file
    pprBench = np.loadtxt(os.path.join(dirPath, '..', '..', 'benchmarks', avaNameTest,
                                       'release1PF_entres_dfa_1.25000_ppr.asc'), skiprows=6)
    # Compare result to reference solution
    testRes = np.allclose(pprTest, pprBench, atol=1.e-12)

    assert testRes is True


def test_splitIniValueToArraySteps():
    """ Test if splitting into an array works fine  """

    cfgValues = '1.0|2.5|3.8'
    cfgValuesList = np.asarray([1.0, 2.5, 3.8])

    cfgValues2 = '0:10:5'
    cfgValuesList2 = np.asarray([0., 2.5, 5., 7.5, 10.])

    # call function to be tested
    items = fU.splitIniValueToArraySteps(cfgValues)
    items2 = fU.splitIniValueToArraySteps(cfgValues2)
    items3 = fU.splitIniValueToArraySteps(cfgValues, returnList=True)

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
    assert len(items3) == len(cfgValuesList)
    assert items3[0] == '1.0'
    assert items3[1] == '2.5'
    assert items3[2] == '3.8'
    assert isinstance(items3, list)


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

    cfgValues4 = '0:22'
    cfgValuesList4 = np.asarray([20.])

    cfgValues5 = '0'
    cfgValuesList5 = np.asarray([40.])

    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'tEnd': '20'}
    cfgGen = cfg['GENERAL']

    # call function to be tested
    cfgGen['tSteps'] = cfgValues
    items = fU.splitTimeValueToArrayInterval(cfgGen['tSteps'], cfgGen.getfloat('tEnd'))
    cfgGen['tSteps'] = cfgValues1
    items1 = fU.splitTimeValueToArrayInterval(cfgGen['tSteps'], cfgGen.getfloat('tEnd'))
    cfgGen['tSteps'] = cfgValues2
    items2 = fU.splitTimeValueToArrayInterval(cfgGen['tSteps'], cfgGen.getfloat('tEnd'))
    cfgGen['tSteps'] = cfgValues3
    items3 = fU.splitTimeValueToArrayInterval(cfgGen['tSteps'], cfgGen.getfloat('tEnd'))
    cfgGen['tSteps'] = cfgValues4
    items4 = fU.splitTimeValueToArrayInterval(cfgGen['tSteps'], cfgGen.getfloat('tEnd'))
    cfgGen['tSteps'] = cfgValues5
    items5 = fU.splitTimeValueToArrayInterval(cfgGen['tSteps'], cfgGen.getfloat('tEnd'))

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
    assert len(items4) == len(cfgValuesList4)
    assert items4[0] == cfgValuesList4[0]
    assert len(items5) == len(cfgValuesList5)
    assert items5[0] == cfgValuesList5[0]


def test_getFilterDict():
    """ test generation of filter dictionary """

    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg['GENERAL'] = {'tEnd': '20'}
    cfg['FILTER'] = {'relTh': '1:2:3', 'entH': 200, 'simType': '', 'secRelArea': 'True'}

    parametersDict = fU.getFilterDict(cfg, 'FILTER')

    noKey = 'simType' in parametersDict

    print('parametersDict', parametersDict)

    assert np.allclose(parametersDict['relTh'], np.asarray([1, 1.5, 2]), atol=1e-10)
    assert noKey is False
    assert parametersDict['entH'] == [200.]
    assert parametersDict['secRelArea'] == ['True']

    parametersDict = fU.getFilterDict(cfg, 'TESTS')

    assert parametersDict == {}


def test_fetchFlowFields():
    """ test fetching fields in a folder """

    # setup required input
    # get input data
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannelPytest'
    flowFieldsDir = dirPath / '..' / '..' / 'benchmarks' / avaName / 'Outputs' / 'com1DFA' / 'peakFiles'
    suffix = 'ppr'

    # call function to be tested
    flowFields = fU.fetchFlowFields(flowFieldsDir, suffix=suffix)
    flowFields = sorted(flowFields)
    print('flowFields', flowFields, sorted(flowFields))

    assert flowFields[0].stem == 'release1HS_d10bdc1e81_ent_dfa_ppr'
    assert flowFields[1].stem == 'release2HS_e2145362b7_ent_dfa_ppr'
    assert len(flowFields) == 2

    # call function to be tested
    flowFields = fU.fetchFlowFields(flowFieldsDir)
    flowFields = sorted(flowFields)
    print('flowFields', flowFields, sorted(flowFields))

    assert flowFields[0].stem == 'release1HS_d10bdc1e81_ent_dfa_pft'
    assert len(flowFields) == 6
