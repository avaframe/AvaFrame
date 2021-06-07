''' Tests for dfa2Aimec '''


import numpy as np
import os
import glob
import configparser
import pathlib

# Local imports
import avaframe.ana3AIMEC.dfa2Aimec as dfa2Aimec
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils


def test_mainDfa2Aimec(tmp_path):

    # Initialise inputs
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    testPath = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathData = testPath / 'Outputs' / 'com1DFA' / 'peakFiles'
    pathDict = dfa2Aimec.mainDfa2Aimec(testPath)

    # get path dictionary for test
    pathDTest = {}
    pathDTest['ppr'] = [pathData / 'release1HS_entres_dfa_0.15500_ppr.asc',
                        pathData / 'release2HS_entres_dfa_0.15500_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_entres_dfa_0.15500_pfd.asc',
                        pathData / 'release2HS_entres_dfa_0.15500_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_entres_dfa_0.15500_pfv.asc',
                        pathData / 'release2HS_entres_dfa_0.15500_pfv.asc']
    pathDTest['massBal'] = [os.path.join(testPath, 'Outputs', 'com1DFA', 'mass_release1HS_entres_dfa_0.15500.txt'),
                       os.path.join(testPath, 'Outputs', 'com1DFA', 'mass_release2HS_entres_dfa_0.15500.txt')]


    assert pathDict['ppr'] == pathDTest['ppr']
    assert pathDict['pfd'] == pathDTest['pfd']
    assert pathDict['pfv'] == pathDTest['pfv']
    assert pathDict['massBal'] == pathDTest['massBal']


def test_dfaComp2Aimec(tmp_path):

    # Initialise inputs
    dirPath = pathlib.Path(__file__).parents[0]
    avaTestName = 'avaHockeyChannelPytest'
    testPath = dirPath / '..' / '..' / 'benchmarks' / avaTestName
    pathData = testPath / 'Outputs' / 'com1DFA' / 'peakFiles'
    pathData2 = testPath / 'Outputs' / 'com1DFAPy' / 'peakFiles'
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'comModules': 'com1DFA|com1DFAPy'}
    cfg['FLAGS'] = {'flagMass': 'True'}
    pathDict = dfa2Aimec.dfaComp2Aimec(testPath, cfg, 'release1HS', 'entres')

    # get path dictionary for test
    pathDTest = {}
    pathDTest['ppr'] = [pathData / 'release1HS_entres_dfa_0.15500_ppr.asc', pathData2 / 'release1HS_entres_dfa_0.15500_ppr.asc']
    pathDTest['pfd'] = [pathData / 'release1HS_entres_dfa_0.15500_pfd.asc', pathData2 / 'release1HS_entres_dfa_0.15500_pfd.asc']
    pathDTest['pfv'] = [pathData / 'release1HS_entres_dfa_0.15500_pfv.asc', pathData2 / 'release1HS_entres_dfa_0.15500_pfv.asc']
    pathDTest['massBal'] = [os.path.join(testPath, 'Outputs', 'com1DFA', 'mass_release1HS_entres_dfa_0.15500.txt'), os.path.join(testPath, 'Outputs', 'com1DFAPy', 'mass_release1HS_entres_dfa_0.15500.txt')]


    assert pathDict['ppr'] == pathDTest['ppr']
    assert pathDict['pfd'] == pathDTest['pfd']
    assert pathDict['pfv'] == pathDTest['pfv']
    assert pathDict['massBal'] == pathDTest['massBal']
