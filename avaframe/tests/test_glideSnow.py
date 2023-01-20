"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
import logging
import pytest
import configparser
import pathlib
import shutil


from avaframe.com1DFA import checkCfg
from avaframe.com5GlideSnow import com5GlideSnow
from avaframe.in3Utils import cfgUtils


def test_glideSnow(tmp_path):
    """ test calling of glideSnow """

    testDir = pathlib.Path(__file__).parents[0]
    inputDir = testDir / 'data' / 'testCom1DFA'
    avaDir = pathlib.Path(tmp_path, 'testCom1DFA')
    shutil.copytree(inputDir, avaDir)

    # setup requuired input data
    cfgMain = configparser.ConfigParser()
    cfgMain['FLAGS'] = {'showPlot': 'False', 'savePlot': 'True', 'ReportDir': 'True', 'reportOneFile': 'True',
        'debugPlot': 'False'}
    cfgMain['MAIN'] = {'avalancheDir': str(avaDir)}

    glideSnowCfg = cfgUtils.getModuleConfig(com5GlideSnow)
    glideSnowCfg['com1DFA_override']['dt'] = '1'
    glideSnowCfg['com1DFA_override']['meshCellSize'] = '5'
    glideSnowCfg['com1DFA_override']['tEnd'] = '5'
    glideSnowCfg['com1DFA_override']['entThFromShp'] = 'False'
    glideSnowCfg['com1DFA_override']['entTh'] = '0.3'

    # call function to be tested
    _, _, _, simDF = com5GlideSnow.runGlideSnow(cfgMain, glideSnowCfg)

    print('simDF', simDF.to_string())

    assert simDF['tEnd'].iloc[0] == 5
    assert simDF['relThFromShp'].iloc[0] == 'True'
    assert simDF['glideSnow'].iloc[0] == 1
