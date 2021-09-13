''' Tests for module ana3AIMEC aimecTools '''
import numpy as np
import pathlib
import configparser

# Local imports
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


def test_fetchReferenceSimNo(tmp_path):
    """ test fetchReferenceSimNo"""

    # setup required input
    testPath = pathlib.Path(tmp_path, 'testDir')
    test1PPR = testPath / 'testSim_no1_ppr.asc'
    test1PFD = testPath / 'testSim_no1_pfd.asc'
    test1PFV = testPath / 'testSim_no1_pfv.asc'
    test2PPR = testPath / 'testSim_no2_ppr.asc'
    test2PFD = testPath / 'testSim_no2_pfd.asc'
    test2PFV = testPath / 'testSim_no2_pfv.asc'
    pathDict = {'ppr': [test1PPR, test2PPR],
                'pfv': [test1PFV, test2PFV],
                'pfd': [test1PFD, test2PFD]}
    pathDict['colorParameter'] = []

    cfgSetup = configparser.ConfigParser()
    cfgSetup['GENERAL'] = {'resType': 'pfv', 'referenceSimName': 'testSim_no2', 'referenceSimValue': '5'}

    pathDict = aT.fetchReferenceSimNo(pathDict, cfgSetup['GENERAL'])
    print('pathDict test1', pathDict)
    assert pathDict['referenceFile'] == 1
    assert pathDict[cfgSetup['GENERAL']['resType']][pathDict['referenceFile']] == test2PFV

    cfgSetup['GENERAL']['referenceSimName'] = ''
    pathDict = {'ppr': [test1PPR, test2PPR],
                'pfv': [test1PFV, test2PFV],
                'pfd': [test1PFD, test2PFD]}
    pathDict['colorParameter'] = []
    pathDict = aT.fetchReferenceSimNo(pathDict, cfgSetup['GENERAL'])
    print('pathDict test2', pathDict)
    assert pathDict['referenceFile'] == 0
    assert pathDict[cfgSetup['GENERAL']['resType']][pathDict['referenceFile']] == test1PFV

    cfgSetup['GENERAL']['referenceSimValue'] = '0.5'
    cfgSetup['GENERAL']['varParList'] = 'relTh'
    pathDict = {'ppr': [test1PPR, test2PPR],
                'pfv': [test1PFV, test2PFV],
                'pfd': [test1PFD, test2PFD]}
    pathDict['colorParameter'] = [0.5, 1.0]
    pathDict = aT.fetchReferenceSimNo(pathDict, cfgSetup['GENERAL'])
    print('pathDict test3', pathDict)
    assert pathDict['referenceFile'] == 0
    assert pathDict[cfgSetup['GENERAL']['resType']][pathDict['referenceFile']] == test1PFV

    cfgSetup['GENERAL'] = {'resType': 'pfv', 'referenceSimName': 'testSim_no2', 'referenceSimValue': '1.0'}
    cfgSetup['GENERAL']['varParList'] = 'relTh'
    pathDict = {'ppr': [test1PPR, test2PPR],
                'pfv': [test1PFV, test2PFV],
                'pfd': [test1PFD, test2PFD]}
    pathDict['colorParameter'] = [0.5, 1.0]
    pathDict = aT.fetchReferenceSimNo(pathDict, cfgSetup['GENERAL'])
    print('pathDict test4', pathDict)
    assert pathDict['referenceFile'] == 1
    assert pathDict[cfgSetup['GENERAL']['resType']][pathDict['referenceFile']] == test2PFV

    cfgSetup['GENERAL'] = {'resType': 'pfv', 'referenceSimName': 'testSim_no4', 'referenceSimValue': '1.0'}
    cfgSetup['GENERAL']['varParList'] = 'relTh'
    pathDict = {'ppr': [test1PPR, test2PPR],
                'pfv': [test1PFV, test2PFV],
                'pfd': [test1PFD, test2PFD]}
    pathDict['colorParameter'] = []
    pathDict = aT.fetchReferenceSimNo(pathDict, cfgSetup['GENERAL'])
    print('pathDict test5', pathDict)
    assert pathDict['referenceFile'] == 0
    assert pathDict[cfgSetup['GENERAL']['resType']][pathDict['referenceFile']] == test1PFV

    cfgSetup['GENERAL'] = {'resType': 'pfd', 'referenceSimName': '', 'referenceSimValue': 'Coulomb'}
    cfgSetup['GENERAL']['varParList'] = 'frictModel'
    pathDict = {'ppr': [test1PPR, test2PPR],
                'pfv': [test1PFV, test2PFV],
                'pfd': [test1PFD, test2PFD]}
    pathDict['colorParameter'] = ['samosAT', 'Coulomb']
    pathDict = aT.fetchReferenceSimNo(pathDict, cfgSetup['GENERAL'])
    print('pathDict test6', pathDict)
    assert pathDict['referenceFile'] == 1
    assert pathDict[cfgSetup['GENERAL']['resType']][pathDict['referenceFile']] == test2PFD

    cfgSetup['GENERAL'] = {'resType': 'pfv', 'referenceSimName': 'testSim_no4', 'referenceSimValue': '1.0'}
    cfgSetup['GENERAL']['varParList'] = ''
    pathDict = {'ppr': [test1PPR, test2PPR],
                'pfv': [test1PFV, test2PFV],
                'pfd': [test1PFD, test2PFD]}
    pathDict['colorParameter'] = []
    pathDict = aT.fetchReferenceSimNo(pathDict, cfgSetup['GENERAL'])
    print('pathDict test7', pathDict)
    assert pathDict['referenceFile'] == 0
    assert pathDict[cfgSetup['GENERAL']['resType']][pathDict['referenceFile']] == test1PFV
