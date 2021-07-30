''' Tests for dfa2Aimec '''

import pathlib
import configparser

# Local imports
import avaframe.ana3AIMEC.dfa2Aimec as dfa2Aimec


def test_mainDfa2Aimec(tmp_path):

    # Initialise inputs
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    testPath = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathData = testPath / 'Outputs' / 'com1DFA' / 'peakFiles'
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'varParList': 'releaseScenario', 'ascendingOrder': 'True'}
    pathDict = dfa2Aimec.mainDfa2Aimec(testPath, 'com1DFA', cfg['AIMECSETUP'])
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
    pathDTest = {}
    pathDTest['ppr'] = [pathData / 'release1HS_ent_dfa_0.15500_ppr.asc', pathData2 / 'release1HS_ent_dfa_67dc2dc10a_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_ent_dfa_0.15500_pfd.asc', pathData2 / 'release1HS_ent_dfa_67dc2dc10a_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_ent_dfa_0.15500_pfv.asc', pathData2 / 'release1HS_ent_dfa_67dc2dc10a_pfv.asc']
    pathDTest['massBal'] = [testPath / 'Outputs' / 'com1DFAOrig' / 'mass_release1HS_ent_dfa_0.15500.txt', testPath / 'Outputs' / 'com1DFA' / 'mass_release1HS_ent_dfa_67dc2dc10a.txt']

    assert pathDict['ppr'] == pathDTest['ppr']
    assert pathDict['pfd'] == pathDTest['pfd']
    assert pathDict['pfv'] == pathDTest['pfv']
    for massName1, massName2 in zip(pathDict['massBal'], pathDTest['massBal']):
        assert str(massName1) == str(massName2)
