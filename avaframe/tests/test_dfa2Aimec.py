''' Tests for dfa2Aimec '''


import numpy as np
import os
import glob
import configparser

# Local imports
import avaframe.ana3AIMEC.dfa2Aimec as dfa2Aimec
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils


def test_mainDfa2Aimec(tmp_path):

    # Initialise inputs
    dirPath = os.path.dirname(__file__)
    avaTestName = 'avaHockeyChannelPytest'
    testPath = os.path.join(dirPath, '..', '..', 'benchmarks', avaTestName)
    pathData = os.path.join(testPath, 'Outputs', 'com1DFA', 'peakFiles')
    pathDict = dfa2Aimec.mainDfa2Aimec(testPath)

    # get path dictionary for test
    pathDTest = {}
    pathDTest['ppr'] = [os.path.join(pathData, 'release1HS_entres_dfa_0.15500_ppr.asc'),
                        os.path.join(pathData, 'release2HS_entres_dfa_0.15500_ppr.asc')]
    pathDTest['pfd'] = [os.path.join(pathData, 'release1HS_entres_dfa_0.15500_pfd.asc'),
                        os.path.join(pathData, 'release2HS_entres_dfa_0.15500_pfd.asc')]
    pathDTest['pfv'] = [os.path.join(pathData, 'release1HS_entres_dfa_0.15500_pfv.asc'),
                        os.path.join(pathData, 'release2HS_entres_dfa_0.15500_pfv.asc')]
    pathDTest['mb'] = [os.path.join(testPath, 'Outputs', 'com1DFA', 'mass_release1HS_entres_dfa_0.15500.txt'),
                       os.path.join(testPath, 'Outputs', 'com1DFA', 'mass_release2HS_entres_dfa_0.15500.txt')]


    assert pathDict['ppr'] == pathDTest['ppr']
    assert pathDict['pfd'] == pathDTest['pfd']
    assert pathDict['pfv'] == pathDTest['pfv']
    assert pathDict['mb'] == pathDTest['mb']



def test_dfaComp2Aimec(tmp_path):

    # Initialise inputs
    dirPath = os.path.dirname(__file__)
    avaTestName = 'avaHockeyChannelPytest'
    testPath = os.path.join(dirPath, '..', '..', 'benchmarks', avaTestName)
    pathData = os.path.join(testPath, 'Outputs', 'com1DFA', 'peakFiles')
    pathData2 = os.path.join(testPath, 'Outputs', 'com1DFAPy', 'peakFiles')
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'comModules': 'com1DFA|com1DFAPy'}
    cfg['FLAGS'] = {'flagMass': 'True'}
    pathDict = dfa2Aimec.dfaComp2Aimec(testPath, cfg, 'release1HS', 'entres')

    # get path dictionary for test
    pathDTest = {}
    pathDTest['ppr'] = [os.path.join(pathData, 'release1HS_entres_dfa_0.15500_ppr.asc'), os.path.join(pathData2, 'release1HS_entres_dfa_0.15500_ppr.asc')]
    pathDTest['pfd'] = [os.path.join(pathData, 'release1HS_entres_dfa_0.15500_pfd.asc'), os.path.join(pathData2, 'release1HS_entres_dfa_0.15500_pfd.asc')]
    pathDTest['pfv'] = [os.path.join(pathData, 'release1HS_entres_dfa_0.15500_pfv.asc'), os.path.join(pathData2, 'release1HS_entres_dfa_0.15500_pfv.asc')]
    pathDTest['mb'] = [os.path.join(testPath, 'Outputs', 'com1DFA', 'mass_release1HS_entres_dfa_0.15500.txt'), os.path.join(testPath, 'Outputs', 'com1DFAPy', 'mass_release1HS_entres_dfa_0.15500.txt')]


    assert pathDict['ppr'] == pathDTest['ppr']
    assert pathDict['pfd'] == pathDTest['pfd']
    assert pathDict['pfv'] == pathDTest['pfv']
    assert pathDict['mb'] == pathDTest['mb']
