''' Tests for dfa2Aimec '''

import pathlib
import pytest
import configparser
import numpy as np

# Local imports
import avaframe.ana3AIMEC.dfa2Aimec as dfa2Aimec


def test_extractCom1DFAMBInfo():
    """ test extracting mass balance from log file """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathDict = {'massBal': []}

    # call function to be tested and check for correct error if file does not exist
    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.extractCom1DFAMBInfo(avaDir, pathDict, simNameInput='testName')
    assert 'starttestName.log' in str(e.value)

    # call function to be tested
    pathDict = dfa2Aimec.extractCom1DFAMBInfo(avaDir, pathDict, simNameInput='')

    print('pathDict', pathDict)

    # read created mass file
    massFile = avaDir / 'Outputs' / 'com1DFAOrig' / 'mass_release1HS_ent_dfa_0.15500.txt'
    massTime = np.loadtxt(massFile, delimiter=',', skiprows=1)
    print('massTime', massTime[1:10, :])
    print('ent mass', massTime[np.where(massTime[:, 0] == 30.2)], np.where(massTime[:, 0] == 30.2))

    assert str(pathDict['massBal'][0]) == str(massFile)
    assert np.where(massTime[:, 2] > 0.0)[0][0] == 301.0

    # call function to be tested
    simNameInput = 'release1HS_ent_dfa_0.15500'
    pathDict2 = dfa2Aimec.extractCom1DFAMBInfo(avaDir, pathDict, simNameInput=simNameInput)

    # read created mass file
    massFile2 = avaDir / 'Outputs' / 'com1DFAOrig' / 'mass_release1HS_ent_dfa_0.15500.txt'
    massTime2 = np.loadtxt(massFile2, delimiter=',', skiprows=1)
    print('massTime', massTime2[1:10, :])
    print('ent mass', massTime2[np.where(massTime2[:, 0] == 30.2)], np.where(massTime2[:, 0] == 30.2))

    assert str(pathDict2['massBal'][0]) == str(massFile2)
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
    pathDict = dfa2Aimec.mainDfa2Aimec(testPath, 'com1DFA', cfg)
    print('path', dirPath)
    # get path dictionary for test
    pathDTest = {}
    pathDTest['ppr'] = [pathData / 'release1HS_ent_dfa_67dc2dc10a_ppr.asc',
                        pathData / 'release2HS_ent_dfa_872f0101a4_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_ent_dfa_67dc2dc10a_pfd.asc',
                        pathData / 'release2HS_ent_dfa_872f0101a4_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_ent_dfa_67dc2dc10a_pfv.asc',
                        pathData / 'release2HS_ent_dfa_872f0101a4_pfv.asc']
    pathDTest['massBal'] = [testPath / 'Outputs' / 'com1DFA' / 'mass_release1HS_ent_dfa_67dc2dc10a.txt',
                            testPath / 'Outputs' / 'com1DFA' / 'mass_release2HS_ent_dfa_872f0101a4.txt']

    assert pathDict['ppr'] == pathDTest['ppr']
    assert pathDict['pfd'] == pathDTest['pfd']
    assert pathDict['pfv'] == pathDTest['pfv']
    for massName1, massName2 in zip(pathDict['massBal'], pathDTest['massBal']):
        assert str(massName1) == str(massName2)

    pathDict2 = dfa2Aimec.mainDfa2Aimec(testPath, 'com1DFAOrig', cfg)

    pathData = testPath / 'Outputs' / 'com1DFAOrig' / 'peakFiles'
    pathDTest['ppr'] = [pathData / 'release1HS_ent_dfa_0.15500_ppr.asc',
                        pathData / 'release2HS_ent_dfa_0.15500_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_ent_dfa_0.15500_pfd.asc',
                        pathData / 'release2HS_ent_dfa_0.15500_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_ent_dfa_0.15500_pfv.asc',
                        pathData / 'release2HS_ent_dfa_0.15500_pfv.asc']
    pathDTest['massBal'] = [testPath / 'Outputs' / 'com1DFAOrig' / 'mass_release1HS_ent_dfa_0.15500.txt',
                            testPath / 'Outputs' / 'com1DFAOrig' / 'mass_release2HS_ent_dfa_0.15500.txt']

    assert pathDict2['ppr'] == pathDTest['ppr']
    assert pathDict2['pfd'] == pathDTest['pfd']
    assert pathDict2['pfv'] == pathDTest['pfv']
    for massName1, massName2 in zip(pathDict2['massBal'], pathDTest['massBal']):
        assert str(massName1) == str(massName2)


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
    pathDict = dfa2Aimec.dfaComp2Aimec(testPath, cfg, 'release1HS', 'ent')

    # get path dictionary for test
    massNameRef = 'mass_release1HS_ent_dfa_0.15500.txt'
    massNameSim = 'mass_release1HS_ent_dfa_67dc2dc10a.txt'
    pathDTest = {}
    pathDTest['ppr'] = [pathData / 'release1HS_ent_dfa_0.15500_ppr.asc',
                        pathData2 / 'release1HS_ent_dfa_67dc2dc10a_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_ent_dfa_0.15500_pfd.asc',
                        pathData2 / 'release1HS_ent_dfa_67dc2dc10a_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_ent_dfa_0.15500_pfv.asc',
                        pathData2 / 'release1HS_ent_dfa_67dc2dc10a_pfv.asc']
    pathDTest['massBal'] = [testPath / 'Outputs' / 'com1DFAOrig' / massNameRef,
                            testPath / 'Outputs' / 'com1DFA' / massNameSim]

    assert pathDict['ppr'] == pathDTest['ppr']
    assert pathDict['pfd'] == pathDTest['pfd']
    assert pathDict['pfv'] == pathDTest['pfv']
    for massName1, massName2 in zip(pathDict['massBal'], pathDTest['massBal']):
        assert str(massName1) == str(massName2)

    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.dfaComp2Aimec(testPath, cfg, 'release3HS', 'ent')
    assert 'No matching simulations found for reference and comparison simulation' in str(e.value)


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

    # call function to be tested and check for correct error if file does not exist
    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.getMBInfo(avaDir, pathDict, comMod, simName='testName')
    assert 'mass_testName.txt' in str(e.value)

    with pytest.raises(FileNotFoundError) as e:
        comModTest = 'com1DFATest'
        assert dfa2Aimec.getMBInfo(avaDir, pathDict, comModTest, simName='')
    assert 'avaHockeyChannelPytest/Outputs/com1DFATest' in str(e.value)

    # call fucntion to be tested
    pathDict = dfa2Aimec.getMBInfo(avaDir, pathDict, comMod, simName=simName)

    print('pathDict', pathDict)
    assert 'mass_release1HS_ent_dfa_67dc2dc10a' in str(pathDict['massBal'][0])
    assert 'mass_release2HS_ent_dfa_872f0101a4' in str(pathDict['massBal'][1])

    # call fucntion to be tested
    simName = 'release1HS_ent_dfa_67dc2dc10a'
    pathDict2 = dfa2Aimec.getMBInfo(avaDir, pathDict, comMod, simName=simName)

    print('pathDict2', pathDict2)
    assert 'mass_release1HS_ent_dfa_67dc2dc10a' in str(pathDict['massBal'][0])
    assert 'mass_release2HS_ent_dfa_872f0101a4' not in str(pathDict['massBal'][0])


def test_getPathsFromSimName():
    """ test get paths from simName """

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'comModules': 'com1DFAOrig|com1DFA'}
    cfg['FLAGS'] = {'flagMass': 'True'}
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    inputDirRef = avaDir / 'Outputs' / 'com1DFAOrig' / 'peakFiles'
    simNameRef = 'release1HS_ent_dfa_0.15500'
    inputDirComp = avaDir / 'Outputs' / 'com1DFA' / 'peakFiles'
    simNameComp = 'release1HS_ent_dfa_67dc2dc10a'
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'massBal': []}

    # call function to be tested
    pathDict = dfa2Aimec.getPathsFromSimName(pathDict, avaDir, cfg, inputDirRef, simNameRef, inputDirComp, simNameComp)

    print('pathDict', pathDict)

    assert 'com1DFA/mass_release1HS_ent_dfa_67dc2dc10a' in str(pathDict['massBal'][1])
    assert 'com1DFA/peakFiles/release1HS_ent_dfa_67dc2dc10a_pfd' in str(pathDict['pfd'][1])
    assert 'com1DFA/peakFiles/release1HS_ent_dfa_67dc2dc10a_ppr' in str(pathDict['ppr'][1])
    assert 'com1DFA/peakFiles/release1HS_ent_dfa_67dc2dc10a_pfv' in str(pathDict['pfv'][1])
    assert 'com1DFAOrig/mass_release1HS_ent_dfa_0.15500' in str(pathDict['massBal'][0])
    assert 'com1DFAOrig/peakFiles/release1HS_ent_dfa_0.15500_pfd' in str(pathDict['pfd'][0])
    assert 'release1HS_ent_dfa_0.15500_ppr' in str(pathDict['ppr'][0])
    assert 'release1HS_ent_dfa_0.15500_pfv' in str(pathDict['pfv'][0])

    # call function to be tested
    cfg['AIMECSETUP']['comModules'] = 'benchmarkReference|com1DFA'
    cfg['AIMECSETUP']['testName'] = 'avaHockeyChannelEntTest'
    cfg['FLAGS'] = {'flagMass': 'False'}
    inputDirRef = dirPath / '..' / '..' / 'benchmarks' / 'avaHockeyChannelEntTest'
    simNameRef = 'release1HS_ent_ref_0.15500'
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'massBal': []}
    pathDict2 = dfa2Aimec.getPathsFromSimName(pathDict, avaDir, cfg, inputDirRef, simNameRef, inputDirComp, simNameComp)

    print('pathDict2', pathDict2)

    assert pathDict2['massBal'] == []
    assert 'release1HS_ent_dfa_67dc2dc10a_pfd' in str(pathDict2['pfd'][1])
    assert 'release1HS_ent_dfa_67dc2dc10a_ppr' in str(pathDict2['ppr'][1])
    assert 'release1HS_ent_dfa_67dc2dc10a_pfv' in str(pathDict2['pfv'][1])
    assert 'avaHockeyChannelEntTest/release1HS_ent_ref_0.15500_pfd' in str(pathDict2['pfd'][0])
    assert 'avaHockeyChannelEntTest/release1HS_ent_ref_0.15500_ppr' in str(pathDict2['ppr'][0])
    assert 'avaHockeyChannelEntTest/release1HS_ent_ref_0.15500_pfv' in str(pathDict2['pfv'][0])

    # call function to be tested
    cfg['AIMECSETUP']['comModules'] = 'com1DFAOrig|com1DFA'
    cfg['AIMECSETUP']['testName'] = 'avaHockeyChannelEntTest'
    cfg['FLAGS'] = {'flagMass': 'False'}
    inputDirRef = dirPath / '..' / '..' / 'benchmarks' / 'avaHockeyChannelEntTest'
    simNameRef = 'release1HS_null_ref_0.15500'
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'massBal': []}

    refFileTest = inputDirRef / simNameRef
    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.getPathsFromSimName(pathDict, avaDir, cfg, inputDirRef, simNameRef, inputDirComp, simNameComp)
    assert ('No reference simulation file found called: %s' % str(refFileTest))  in str(e.value)

    simNameRef = 'release1HS_ent_ref_0.15500'
    simNameComp = 'release1HS_null_dfa_67dc2dc10a'
    compFileTest = inputDirComp / simNameComp
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'massBal': []}
    with pytest.raises(FileNotFoundError) as e:
        assert dfa2Aimec.getPathsFromSimName(pathDict, avaDir, cfg, inputDirRef, simNameRef, inputDirComp, simNameComp)
    assert ('No comparison simulation file found called: %s' % str(compFileTest))  in str(e.value)


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
    simNameComp = 'release1HS_ent_dfa_67dc2dc10a'

    # call function to be tested
    pathDict = dfa2Aimec.dfaBench2Aimec(avaDir, cfg, simNameRef, simNameComp)

    print('pathDict', pathDict)

    assert 'mass_release1HS_ent_dfa_67dc2dc10a' in str(pathDict['massBal'][1])
    assert 'release1HS_ent_dfa_67dc2dc10a_pfd' in str(pathDict['pfd'][1])
    assert 'release1HS_ent_dfa_67dc2dc10a_ppr' in str(pathDict['ppr'][1])
    assert 'release1HS_ent_dfa_67dc2dc10a_pfv' in str(pathDict['pfv'][1])
    assert 'mass_release1HS_ent_dfa_0.15500' in str(pathDict['massBal'][0])
    assert 'release1HS_ent_dfa_0.15500_pfd' in str(pathDict['pfd'][0])
    assert 'release1HS_ent_dfa_0.15500_ppr' in str(pathDict['ppr'][0])
    assert 'release1HS_ent_dfa_0.15500_pfv' in str(pathDict['pfv'][0])


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
    inputDirRef, inputDirComp, pathDict, refModule = dfa2Aimec.getCompDirs(avaDir, cfg['AIMECSETUP'], pathDict)

    print('inputDirRef', inputDirRef)
    print('inputDirComp', inputDirComp)
    print('pathDict', pathDict)
    print('refModule', refModule)

    assert str(inputDirRef) == str(avaDir / 'Outputs' / 'com1DFAOrig' / 'peakFiles')
    assert str(inputDirComp) == str(avaDir / 'Outputs' / 'com1DFA' / 'peakFiles')
    assert refModule == 'com1DFAOrig'
    assert pathDict['compType'] == ['comModules', 'com1DFAOrig', 'com1DFA']
    assert pathDict['referenceFile'] == 0
    assert pathDict['contCmap'] is True

    # call function to be tested
    cfg['AIMECSETUP'] = {'comModules': 'benchmarkReference|com1DFA'}
    pathDict = {'ppr': []}
    cfg['AIMECSETUP']['testName'] = 'avaHockeyChannelEntTest'
    inputDirRef, inputDirComp, pathDict, refModule = dfa2Aimec.getCompDirs(avaDir, cfg['AIMECSETUP'], pathDict)

    testPath = dirPath / '..' / '..' / 'benchmarks' / 'avaHockeyChannelEntTest'

    assert str(inputDirRef) == '../benchmarks/avaHockeyChannelEntTest'
    assert str(inputDirComp) == str(avaDir / 'Outputs' / 'com1DFA' / 'peakFiles')
    assert refModule == 'benchmarkReference'
    assert pathDict['compType'] == ['comModules', 'benchmarkReference', 'com1DFA']
    assert pathDict['referenceFile'] == 0
    assert pathDict['contCmap'] is True


def test_indiDfa2Aimec():
    """ test only save paths to data for one result type for aimec analysis """

    # setup required path
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'varParList': 'releaseScenario|relTh', 'ascendingOrder': 'True'}
    suffix = 'ppr'
    cfg['AIMECSETUP']['anaMod'] = 'com1DFA'

    # call function to be tested
    pathDict = dfa2Aimec.indiDfa2Aimec(avaDir, suffix, cfg['AIMECSETUP'])
    print('pathDict', pathDict)

    assert 'pfv' not in pathDict.keys()
    assert 'release1HS_ent_dfa_67dc2dc10a_ppr' in str(pathDict['ppr'][0])
    assert 'pfd' not in pathDict.keys()
    assert pathDict['compType'] == ['singleModule', 'com1DFA']
    assert pathDict['colorParameter'] == ['release1HS', 'release2HS']
