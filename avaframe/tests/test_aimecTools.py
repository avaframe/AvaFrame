''' Tests for module ana3AIMEC aimecTools '''
import pandas as pd
import pathlib
import configparser
import pytest

# Local imports
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils


def test_fetchReferenceSimNo(tmp_path):
    """ test fetchReferenceSimNo"""

    # setup required input
    avaDir = pathlib.Path(tmp_path, 'testDir')
    testPath = avaDir / 'Outputs' / 'comModule' / 'peakFiles'
    fU.makeADir(testPath)
    test1PPR = testPath / 'testSim_no1_ppr.asc'
    test1PFT = testPath / 'testSim_no1_pft.asc'
    test1PFV = testPath / 'testSim_no1_pfv.asc'
    test2PPR = testPath / 'testSim_no2_ppr.asc'
    test2PFT = testPath / 'testSim_no2_pft.asc'
    test2PFV = testPath / 'testSim_no2_pfv.asc'
    d = {'simName': ['testSim_no1', 'testSim_no2'], 'ppr': [test1PPR, test2PPR],
         'pft': [test1PFT, test2PFT], 'pfv': [test1PFV, test2PFV]}
    inputsDF = pd.DataFrame(data=d, index=['testSim_no1', 'testSim_no2'])
    cfgSetup = configparser.ConfigParser()
    cfgSetup['GENERAL'] = {'resType': 'pfv', 'referenceSimName': 'testSim_no2', 'referenceSimValue': '',
                           'varParList': ''}
    refSimHash, refSimName, inputsDF, colorParameter = aT.fetchReferenceSimNo(avaDir, inputsDF, 'comModule',
                                                                              cfgSetup['GENERAL'])
    assert refSimName == 'testSim_no2'
    assert colorParameter is False

    cfgSetup['GENERAL']['referenceSimName'] = ''
    refSimHash, refSimName, inputsDF, colorParameter = aT.fetchReferenceSimNo(avaDir, inputsDF, 'comModule',
                                                                              cfgSetup['GENERAL'])
    assert refSimName == 'testSim_no1'
    assert colorParameter is False
    assert inputsDF.loc[refSimHash, cfgSetup['GENERAL']['resType']] == test1PFV


def test_computeCellSizeSL(tmp_path):
    """ test fetchReferenceSimNo"""
    cfg = configparser.ConfigParser()
    cfg['AIMECSETUP'] = {'cellSizeSL': ''}
    cfgSetup = cfg['AIMECSETUP']
    demHeader = {'cellsize': 1}

    # read the cell size from the header
    cellSizeSL = aT.computeCellSizeSL(cfgSetup, demHeader['cellsize'])
    assert cellSizeSL == 1

    # read the cell size from the cfg
    cfgSetup['cellSizeSL'] = '3'
    cellSizeSL = aT.computeCellSizeSL(cfgSetup, demHeader['cellsize'])
    assert cellSizeSL == 3

    # read the cell size from the cfg
    cfgSetup['cellSizeSL'] = '3.1'
    cellSizeSL = aT.computeCellSizeSL(cfgSetup, demHeader['cellsize'])
    assert cellSizeSL == 3.1

    # check error if no number provided but a character
    cfgSetup['cellSizeSL'] = 'c'
    message = ('cellSizeSL is read from the configuration file but should be a number, you provided: c')
    with pytest.raises(ValueError) as e:
        assert aT.computeCellSizeSL(cfgSetup, demHeader['cellsize'])
    assert str(e.value) == message
