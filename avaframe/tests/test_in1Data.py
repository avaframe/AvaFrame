"""
    Pytest for module in1Data

    This file is part of Avaframe.

 """

#  Load modules
import os
import pathlib
from avaframe.in1Data import getInput
import configparser
import shutil


def test_getInputData(tmp_path):
    """ test check for input data """

    # get input data
    dirPath = os.path.dirname(__file__)
    avaName = 'avaHockeyChannel'
    avaDir = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'flagEnt': 'True', 'flagRes': 'True', 'flagDev': 'False', 'releaseScenario': ''}
    cfgGen = cfg['GENERAL']

    # call function to be tested
    dem, rels, ent, res, entResInfo = getInput.getInputData(avaDir, cfgGen)
    # second option
    cfg['GENERAL']['releaseScenario'] = 'release1HS'
    dem2, rels2, ent2, res2, entResInfo2 = getInput.getInputData(avaDir, cfgGen)
    # Test
    assert str(dem) == str(pathlib.Path(avaDir, 'Inputs', 'DEM_HS_Topo.asc'))
    assert rels == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp'), os.path.join(avaDir, 'Inputs', 'REL', 'release2HS.shp'), os.path.join(avaDir, 'Inputs', 'REL', 'release3HS.shp')]
    assert rels2 == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp')]
    assert res == ''
    assert str(ent) == str(os.path.join(avaDir, 'Inputs', 'ENT', 'entrainment1HS.shp'))
    assert entResInfo['flagEnt'] == "Yes"
    assert entResInfo['flagRes'] == "No"


def test_getInputDataCom1DFA(tmp_path):
    """ test check for input data """

    # get input data
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannel'
    avaDir = pathlib.Path(tmp_path, avaName)
    avaInputs = pathlib.Path(avaDir, 'Inputs')
    avaData = pathlib.Path(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'flagEnt': 'True', 'flagRes': 'True', 'flagDev': 'False', 'releaseScenario': ''}
    cfgGen = cfg['GENERAL']

    # call function to be tested
    inputSimFiles = getInput.getInputDataCom1DFA(avaDir, cfgGen)
    # second option
    cfg['GENERAL']['releaseScenario'] = 'release1HS'
    inputSimFiles2 = getInput.getInputDataCom1DFA(avaDir, cfgGen)
    # Test
    print(inputSimFiles['demFile'])
    print(pathlib.Path(avaDir, 'Inputs', 'DEM_HS_Topo.asc'))
    print(inputSimFiles['relFiles'])
    print([pathlib.Path(avaDir, 'Inputs', 'REL', 'release1HS.shp'), pathlib.Path(avaDir, 'Inputs', 'REL', 'release2HS.shp'), pathlib.Path(avaDir, 'Inputs', 'REL', 'release3HS.shp')])
    assert inputSimFiles['demFile'] == pathlib.Path(avaDir, 'Inputs', 'DEM_HS_Topo.asc')
    assert inputSimFiles['relFiles'] == [pathlib.Path(avaDir, 'Inputs', 'REL', 'release1HS.shp'), pathlib.Path(avaDir, 'Inputs', 'REL', 'release2HS.shp'), pathlib.Path(avaDir, 'Inputs', 'REL', 'release3HS.shp')]
    assert inputSimFiles2['relFiles'] == [pathlib.Path(avaDir, 'Inputs', 'REL', 'release1HS.shp')]
    assert inputSimFiles['resFile'] == None
    assert inputSimFiles['entFile'] == pathlib.Path(avaDir, 'Inputs', 'ENT', 'entrainment1HS.shp')
    assert inputSimFiles['entResInfo']['flagEnt'] == "Yes"
    assert inputSimFiles['entResInfo']['flagRes'] == "No"
