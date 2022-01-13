''' Tests for dfa2Aimec '''

import pathlib
import pytest
import configparser
import numpy as np
import pandas as pd

# Local imports
import avaframe.ana3AIMEC.dfa2Aimec as dfa2Aimec


def test_extractCom1DFAMBInfo():
    """ test extracting mass balance from log file """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    inputDF = pd.DataFrame()

    # call function to be tested and check for correct error if file does not exist
    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.extractCom1DFAMBInfo(avaDir, inputDF, simNameInput='testName')
    assert 'starttestName.log' in str(e.value)

    # call function to be tested
    inputDF = dfa2Aimec.extractCom1DFAMBInfo(avaDir, inputDF, simNameInput='')

    print('inputDF', inputDF)

    # read created mass file
    massFile = avaDir / 'Outputs' / 'com1DFAOrig' / 'mass_release1HS_ent_dfa_0.15500.txt'
    massTime = np.loadtxt(massFile, delimiter=',', skiprows=1)
    print('massTime', massTime[1:10, :])
    print('ent mass', massTime[np.where(massTime[:, 0] == 30.2)], np.where(massTime[:, 0] == 30.2))

    assert str(inputDF['massBal'][0]) == str(massFile)
    assert np.where(massTime[:, 2] > 0.0)[0][0] == 301.0

    # call function to be tested
    simNameInput = 'release1HS_ent_dfa_0.15500'
    inputDF = dfa2Aimec.extractCom1DFAMBInfo(avaDir, inputDF, simNameInput=simNameInput)

    # read created mass file
    massFile2 = avaDir / 'Outputs' / 'com1DFAOrig' / 'mass_release1HS_ent_dfa_0.15500.txt'
    massTime2 = np.loadtxt(massFile2, delimiter=',', skiprows=1)
    print('massTime', massTime2[1:10, :])
    print('ent mass', massTime2[np.where(massTime2[:, 0] == 30.2)], np.where(massTime2[:, 0] == 30.2))

    assert str(inputDF['massBal'][0]) == str(massFile2)
    assert np.where(massTime2[:, 2] > 0.0)[0][0] == 301.0


def test_mainDfa2Aimec(tmp_path):

    # Initialise inputs
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    testPath = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathData = testPath / 'Outputs' / 'com1DFA' / 'peakFiles'
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'varParList': 'releaseScenario', 'ascendingOrder': 'True'}
    cfg['FLAGS'] = {'flagMass': 'True'}
    inputDF, resTypeList = dfa2Aimec.mainDfa2Aimec(testPath, 'com1DFA', cfg)
    print('path', dirPath)
    # get path dictionary for test
    pathDTest = {}
    pathDTest['ppr'] = [pathData / 'release1HS_ent_dfa_d10bdc1e81_ppr.asc',
                        pathData / 'release2HS_ent_dfa_e2145362b7_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_ent_dfa_d10bdc1e81_pfd.asc',
                        pathData / 'release2HS_ent_dfa_e2145362b7_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_ent_dfa_d10bdc1e81_pfv.asc',
                        pathData / 'release2HS_ent_dfa_e2145362b7_pfv.asc']
    pathDTest['massBal'] = [testPath / 'Outputs' / 'com1DFA' / 'mass_release1HS_ent_dfa_d10bdc1e81.txt',
                            testPath / 'Outputs' / 'com1DFA' / 'mass_release2HS_ent_dfa_e2145362b7.txt']
    print(inputDF['ppr'].to_list())
    print(pathDTest['ppr'])
    diff = set(inputDF['ppr'].to_list()) ^ set(pathDTest['ppr'])
    assert not diff
    diff = set(inputDF['pfd'].to_list()) ^ set(pathDTest['pfd'])
    assert not diff
    diff = set(inputDF['pfv'].to_list()) ^ set(pathDTest['pfv'])
    assert not diff
    diff = set(inputDF['massBal'].to_list()) ^ set(pathDTest['massBal'])
    assert not diff

    inputDF, resTypeList = dfa2Aimec.mainDfa2Aimec(testPath, 'com1DFAOrig', cfg)

    pathData = testPath / 'Outputs' / 'com1DFAOrig' / 'peakFiles'
    pathDTest['ppr'] = [pathData / 'release1HS_ent_dfa_0.15500_ppr.asc',
                        pathData / 'release2HS_ent_dfa_0.15500_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_ent_dfa_0.15500_pfd.asc',
                        pathData / 'release2HS_ent_dfa_0.15500_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_ent_dfa_0.15500_pfv.asc',
                        pathData / 'release2HS_ent_dfa_0.15500_pfv.asc']
    pathDTest['massBal'] = [testPath / 'Outputs' / 'com1DFAOrig' / 'mass_release1HS_ent_dfa_0.15500.txt',
                            testPath / 'Outputs' / 'com1DFAOrig' / 'mass_release2HS_ent_dfa_0.15500.txt']

    diff = set(inputDF['ppr'].to_list()) ^ set(pathDTest['ppr'])
    assert not diff
    diff = set(inputDF['pfd'].to_list()) ^ set(pathDTest['pfd'])
    assert not diff
    diff = set(inputDF['pfv'].to_list()) ^ set(pathDTest['pfv'])
    assert not diff
    diff = set(inputDF['massBal'].to_list()) ^ set(pathDTest['massBal'])
    assert not diff


def test_dfaComp2Aimec(tmp_path):

    # Initialise inputs
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    testPath = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathData = testPath / 'Outputs' / 'com1DFAOrig' / 'peakFiles'
    pathData2 = testPath / 'Outputs' / 'com1DFA' / 'peakFiles'
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'comModules': 'com1DFAOrig|com1DFA'}
    cfg['FLAGS'] = {'flagMass': 'True'}
    inputDF, pathDict = dfa2Aimec.dfaBench2Aimec(testPath, cfg, 'release1HS_ent', 'release1HS_ent')

    # get path dictionary for test
    massNameRef = 'mass_release1HS_ent_dfa_0.15500.txt'
    massNameSim = 'mass_release1HS_ent_dfa_d10bdc1e81.txt'
    pathDTest = {}
    pathDTest['ppr'] = [pathData / 'release1HS_ent_dfa_0.15500_ppr.asc',
                        pathData2 / 'release1HS_ent_dfa_d10bdc1e81_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_ent_dfa_0.15500_pfd.asc',
                        pathData2 / 'release1HS_ent_dfa_d10bdc1e81_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_ent_dfa_0.15500_pfv.asc',
                        pathData2 / 'release1HS_ent_dfa_d10bdc1e81_pfv.asc']
    pathDTest['massBal'] = [testPath / 'Outputs' / 'com1DFAOrig' / massNameRef,
                            testPath / 'Outputs' / 'com1DFA' / massNameSim]
    print(inputDF['ppr'].to_list())
    print(pathDTest['ppr'])
    diff = set(inputDF['ppr'].to_list()) ^ set(pathDTest['ppr'])
    assert not diff
    diff = set(inputDF['pfd'].to_list()) ^ set(pathDTest['pfd'])
    assert not diff
    diff = set(inputDF['pfv'].to_list()) ^ set(pathDTest['pfv'])
    assert not diff
    diff = set(inputDF['massBal'].to_list()) ^ set(pathDTest['massBal'])
    assert not diff

    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.dfaBench2Aimec(testPath, cfg, 'release3HS_ent', 'release3HS_ent')
    assert 'Did not find the reference simulation : release3HS_ent' in str(e.value)


def test_getRefMB():
    """ test get reference mass balance info """

    # setup required input
    avaTestName = 'avaHockeyChannelPytest'
    pathDict = {'ppr': 'test', 'massBal': []}
    simName = 'release1HS_ent_dfa_67dc2dc10a'

    # call function to be tested
    pathDict = dfa2Aimec.getRefMB(avaTestName, pathDict, simName)

    print('pathDict', pathDict)

    assert 'mass_release1HS_ent_dfa_67dc2dc10a' in str(pathDict['massBal'][0])


def test_getMBInfo():
    """ test get mass balance info """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathDict = {'ppr': 'test', 'massBal': []}
    comMod = 'com1DFA'
    simName = ''
    inputsDF = pd.DataFrame()
    # call function to be tested and check for correct error if file does not exist
    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.getMBInfo(avaDir, inputsDF, comMod, simName='testName')
    assert 'mass_testName.txt' in str(e.value)

    with pytest.raises(FileNotFoundError) as e:
        comModTest = 'com1DFATest'
        assert dfa2Aimec.getMBInfo(avaDir, inputsDF, comModTest, simName='')
    assert 'avaHockeyChannelPytest/Outputs/com1DFATest' in str(e.value)

    # call fucntion to be tested
    inputsDF = pd.DataFrame()
    inputsDF = dfa2Aimec.getMBInfo(avaDir, inputsDF, comMod, simName=simName)

    print('pathDict', pathDict)
    assert 'mass_release1HS_ent_dfa_d10bdc1e81' in str(inputsDF['massBal'][0])
    assert 'mass_release2HS_ent_dfa_e2145362b7' in str(inputsDF['massBal'][1])

    # call fucntion to be tested
    simName = 'release1HS_ent_dfa_d10bdc1e81'
    inputsDF = pd.DataFrame()
    inputsDF = dfa2Aimec.getMBInfo(avaDir, inputsDF, comMod, simName=simName)

    print('pathDict2', inputsDF)
    assert 'mass_release1HS_ent_dfa_d10bdc1e81' in str(inputsDF['massBal'][0])
    assert 'mass_release2HS_ent_dfa_e2145362b7' not in str(inputsDF['massBal'][0])



def test_dfaBench2Aimec():
    """ test export data used for aimec """

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'comModules': 'com1DFAOrig|com1DFA'}
    cfg['FLAGS'] = {'flagMass': 'True'}
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    simNameRef = 'release1HS_ent_dfa_0.15500'
    simNameComp = 'release1HS_ent_dfa_d10bdc1e81'

    # call function to be tested
    inputsDF, pathDict = dfa2Aimec.dfaBench2Aimec(avaDir, cfg, simNameRef, simNameComp)

    print('inputsDF', inputsDF)

    assert 'mass_release1HS_ent_dfa_d10bdc1e81' in str(inputsDF['massBal'][1])
    assert 'release1HS_ent_dfa_d10bdc1e81_pfd' in str(inputsDF['pfd'][1])
    assert 'release1HS_ent_dfa_d10bdc1e81_ppr' in str(inputsDF['ppr'][1])
    assert 'release1HS_ent_dfa_d10bdc1e81_pfv' in str(inputsDF['pfv'][1])
    assert 'mass_release1HS_ent_dfa_0.15500' in str(inputsDF['massBal'][0])
    assert 'release1HS_ent_dfa_0.15500_pfd' in str(inputsDF['pfd'][0])
    assert 'release1HS_ent_dfa_0.15500_ppr' in str(inputsDF['ppr'][0])
    assert 'release1HS_ent_dfa_0.15500_pfv' in str(inputsDF['pfv'][0])


def test_getCompDirs():
    """ test get comparison dirs """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'comModules': 'com1DFAOrig|com1DFA'}
    pathDict = {'ppr': []}

    # call function to be tested
    inputDirRef, inputDirComp, pathDict = dfa2Aimec.getCompDirs(avaDir, cfg['AIMECSETUP'])

    print('inputDirRef', inputDirRef)
    print('inputDirComp', inputDirComp)
    print('pathDict', pathDict)
    print('refModule', pathDict['compType'][1])

    assert str(inputDirRef) == str(avaDir / 'Outputs' / 'com1DFAOrig' / 'peakFiles')
    assert str(inputDirComp) == str(avaDir / 'Outputs' / 'com1DFA' / 'peakFiles')
    assert pathDict['compType'][1] == 'com1DFAOrig'
    assert pathDict['compType'] == ['comModules', 'com1DFAOrig', 'com1DFA']
    assert pathDict['contCmap'] is True

    # call function to be tested
    cfg['AIMECSETUP'] = {'comModules': 'benchmarkReference|com1DFA'}
    pathDict = {'ppr': []}
    cfg['AIMECSETUP']['testName'] = 'avaHockeyChannelEntTest'
    inputDirRef, inputDirComp, pathDict = dfa2Aimec.getCompDirs(avaDir, cfg['AIMECSETUP'])

    testPath = dirPath / '..' / '..' / 'benchmarks' / 'avaHockeyChannelEntTest'

    assert str(inputDirRef) == '../benchmarks/avaHockeyChannelEntTest'
    assert str(inputDirComp) == str(avaDir / 'Outputs' / 'com1DFA' / 'peakFiles')
    assert pathDict['compType'][1] == 'benchmarkReference'
    assert pathDict['compType'] == ['comModules', 'benchmarkReference', 'com1DFA']
    assert pathDict['contCmap'] is True
