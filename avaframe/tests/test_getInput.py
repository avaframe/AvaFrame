"""
    Pytest for module in1Data

    This file is part of Avaframe.

 """

#  Load modules
import os
import pathlib
from avaframe.in1Data import getInput
import configparser
import pytest
import shutil
import numpy as np
from scipy.interpolate import interp1d


def test_readDEM():
    """ test reading DEM """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannel'
    avaDir = dirPath / '..' / 'data' / avaName

    # call function to be tested
    dem = getInput.readDEM(avaDir)

    print('demHeader', dem['header'])

    assert dem['header']['ncols'] == 1001
    assert dem['header']['nrows'] == 401
    assert dem['header']['xllcenter'] == 1000
    assert dem['header']['yllcenter'] == -5000
    assert dem['header']['cellsize'] == 5
    assert dem['rasterData'].shape == (401, 1001)


def test_getDEMPath(tmp_path):
    """ test get path of DEM """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannel'
    avaDir = dirPath / '..' / 'data' / avaName
    avaDirInputs =  avaDir / 'Inputs'
    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / 'Inputs'
    shutil.copytree(avaDirInputs, avaTestDirInputs)

    # call function to be tested
    demPath = getInput.getDEMPath(avaTestDir)

    print('dem path', demPath)

    assert 'DEM_HS_Topo.asc' in str(demPath)

    # call function to be tested
    inputFile = avaDirInputs / 'DEM_HS_Topo.asc'
    testFile = avaTestDirInputs / 'DEM_HS_Topo2.asc'
    shutil.copyfile(inputFile, testFile)

    with pytest.raises(AssertionError) as e:
        assert getInput.getDEMPath(avaTestDir)
    assert str(e.value) == "There should be exactly one topography .asc file in %s/Inputs/" % (avaTestDir)


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

    # third option
    cfg['GENERAL']['releaseScenario'] = 'release1HS.shp'
    dem3, rels3, ent3, res3, entResInfo3 = getInput.getInputData(avaDir, cfg['GENERAL'])
    assert str(ent3) == str(os.path.join(avaDir, 'Inputs', 'ENT', 'entrainment1HS.shp'))
    assert entResInfo3['flagEnt'] == "Yes"
    assert entResInfo3['flagRes'] == "No"
    assert rels3 == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp')]

    # fourth option
    dem5, rels5, ent5, res5, entResInfo5 = getInput.getInputData(avaDir, cfgGen, flagDev=True)
    assert str(ent5) == str(os.path.join(avaDir, 'Inputs', 'ENT', 'entrainment1HS.shp'))
    assert entResInfo5['flagEnt'] == "Yes"
    assert entResInfo5['flagRes'] == "No"
    assert rels5 == []

    # call function to be tested
    cfg['GENERAL']['releaseScenario'] = 'release4HS'
    releaseF = os.path.join(avaDir, 'Inputs', 'REL', 'release4HS.shp')
    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getInputData(avaDir, cfg['GENERAL'])
    assert str(e.value) == ("No release scenario called: %s" % releaseF)

    # fifth option
    cfg['GENERAL']['flagDev'] = 'False'
    cfg['GENERAL']['releaseScenario'] = 'release1BL.shp'
    avaName = 'avaBowl'
    avaDir = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)
    dem6, rels6, ent6, res6, entResInfo6 = getInput.getInputData(avaDir, cfg['GENERAL'])
    assert ent6 == ''
    assert res6 == ''
    assert entResInfo6['flagEnt'] == "No"
    assert entResInfo6['flagRes'] == "No"


def test_getInputDataCom1DFA(tmp_path):
    """ test check for input data """

    # get input data
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannel'
    avaDir = pathlib.Path(tmp_path, avaName)
    avaInputs = avaDir / 'Inputs'
    avaData = dirPath / '..'/ 'data'/ avaName / 'Inputs'
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
    print(avaDir / 'Inputs' / 'DEM_HS_Topo.asc')
    print(inputSimFiles['relFiles'])
    print([avaDir / 'Inputs' / 'REL' / 'release1HS.shp', avaDir / 'Inputs' / 'REL' / 'release2HS.shp', avaDir / 'Inputs' / 'REL' / 'release3HS.shp'])
    assert inputSimFiles['demFile'] == avaDir / 'Inputs' / 'DEM_HS_Topo.asc'
    assert inputSimFiles['relFiles'] == [avaDir / 'Inputs' / 'REL' / 'release1HS.shp',
                                         avaDir / 'Inputs' / 'REL' / 'release2HS.shp',
                                         avaDir / 'Inputs' / 'REL' / 'release3HS.shp']
    assert inputSimFiles2['relFiles'] == [avaDir / 'Inputs' / 'REL' / 'release1HS.shp']
    assert inputSimFiles['resFile'] == None
    assert inputSimFiles['entFile'] == avaDir / 'Inputs' / 'ENT' / 'entrainment1HS.shp'
    assert inputSimFiles['entResInfo']['flagEnt'] == "Yes"
    assert inputSimFiles['entResInfo']['flagRes'] == "No"

    # call function to be tested
    cfg['GENERAL']['releaseScenario'] = 'release1HS.shp'
    inputSimFiles3 = getInput.getInputDataCom1DFA(avaDir, cfg['GENERAL'])

    assert inputSimFiles2['relFiles'] == [avaDir / 'Inputs' / 'REL' / 'release1HS.shp']

    # call function to be tested
    cfg['GENERAL']['flagDev'] = 'True'
    inputSimFiles4 = getInput.getInputDataCom1DFA(avaDir, cfg['GENERAL'])

    assert inputSimFiles4['relFiles'] == []

    # call function to be tested
    cfg['GENERAL']['releaseScenario'] = 'release4HS'
    cfg['GENERAL']['flagDev'] = 'False'
    releaseF = avaDir / 'Inputs' / 'REL' / 'release4HS.shp'
    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getInputDataCom1DFA(avaDir, cfg['GENERAL'])
    assert str(e.value) == ("No release scenario called: %s" % releaseF)


def test_getAndCheckInputFiles(tmp_path):
    """ test fetching input files and checking if exist """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannel'
    avaDir = dirPath / '..' / 'data' / avaName
    avaDirInputs =  avaDir / 'Inputs'
    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / 'Inputs'
    shutil.copytree(avaDirInputs, avaTestDirInputs)
    folder = 'ENT'
    inputType = 'entrainment'

    # call function to be tested
    outFile, available = getInput.getAndCheckInputFiles(avaTestDirInputs, folder, inputType)

    print('outfile', outFile)
    print('available', available)

    assert available == 'Yes'
    assert 'Inputs/ENT/entrainment1HS.shp' in str(outFile)

    # call function to be tested
    inputFile = avaDirInputs / 'ENT' / 'entrainment1HS.shp'
    testFile = avaTestDirInputs / 'ENT' / 'entrainment1HS2.shp'
    shutil.copyfile(inputFile, testFile)
    with pytest.raises(AssertionError) as e:
        assert getInput.getAndCheckInputFiles(avaTestDirInputs, folder, inputType)
    assert str(e.value) == ("More than one %s .shp file in %s/%s/ not allowed" %
                           (inputType, avaTestDirInputs, folder))
