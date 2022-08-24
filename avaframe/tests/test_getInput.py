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
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils


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

    # call function to be tested
    avaDirTest2 = pathlib.Path(tmp_path, 'avaTest')
    avaDirTest2Inputs = avaDirTest2 / 'Inputs'
    fU.makeADir(avaDirTest2Inputs)
    inputFile = avaDirInputs / 'DEM_HS_Topo.asc'
    testFile = avaDirTest2Inputs / 'DEM_HS_Topo2.txt'
    shutil.copyfile(inputFile, testFile)

    with pytest.raises(AssertionError) as e:
        assert getInput.getDEMPath(avaDirTest2)
    assert str(e.value) == "DEM file format not correct in %s/Inputs/ - only .asc is allowed but %s is provided" % (avaDirTest2, testFile.name)

    # call function to be tested
    avaDirTest3 = pathlib.Path(tmp_path, 'avaTest2')


    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getDEMPath(avaDirTest3)
    assert str(e.value) == "No topography .asc file in %s/Inputs/" % (avaDirTest3)


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
    cfg['INPUT'] = {'releaseScenario': ''}

    # call function to be tested
    dem, rels, ent, res, entResInfo = getInput.getInputData(avaDir, cfg['INPUT'])
    # second option
    cfg['INPUT']['releaseScenario'] = 'release1HS'
    dem2, rels2, ent2, res2, entResInfo2 = getInput.getInputData(avaDir, cfg['INPUT'])
    # Test
    assert str(dem) == str(pathlib.Path(avaDir, 'Inputs', 'DEM_HS_Topo.asc'))
    assert rels == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp'), os.path.join(avaDir, 'Inputs', 'REL', 'release2HS.shp'), os.path.join(avaDir, 'Inputs', 'REL', 'release3HS.shp')]
    assert rels2 == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp')]
    assert res == ''
    assert str(ent) == str(os.path.join(avaDir, 'Inputs', 'ENT', 'entrainment1HS.shp'))
    assert entResInfo['flagEnt'] == "Yes"
    assert entResInfo['flagRes'] == "No"

    # third option
    cfg['INPUT']['releaseScenario'] = 'release1HS.shp'
    dem3, rels3, ent3, res3, entResInfo3 = getInput.getInputData(avaDir, cfg['INPUT'])
    assert str(ent3) == str(os.path.join(avaDir, 'Inputs', 'ENT', 'entrainment1HS.shp'))
    assert entResInfo3['flagEnt'] == "Yes"
    assert entResInfo3['flagRes'] == "No"
    assert rels3 == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp')]


    # call function to be tested
    cfg['INPUT']['releaseScenario'] = 'release4HS'
    releaseF = os.path.join(avaDir, 'Inputs', 'REL', 'release4HS.shp')
    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getInputData(avaDir, cfg['INPUT'])
    assert str(e.value) == ("No release scenario called: %s" % releaseF)

    # fifth option
    cfg['INPUT']['releaseScenario'] = 'release1BL.shp'
    avaName = 'avaBowl'
    avaDir = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)
    dem6, rels6, ent6, res6, entResInfo6 = getInput.getInputData(avaDir, cfg['INPUT'])
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
    cfg['GENERAL'] = {'relThFromFile': 'False'}
    cfg['INPUT'] = {'releaseScenario': ''}

    # call function to be tested
    inputSimFiles = getInput.getInputDataCom1DFA(avaDir, cfg)
    # second option
    cfg['INPUT']['releaseScenario'] = 'release1HS'
    inputSimFiles2 = getInput.getInputDataCom1DFA(avaDir, cfg)
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
    assert inputSimFiles['relThFile'] == ''

    # call function to be tested
    cfg['INPUT']['releaseScenario'] = 'release1HS.shp'
    inputSimFiles3 = getInput.getInputDataCom1DFA(avaDir, cfg)

    assert inputSimFiles2['relFiles'] == [avaDir / 'Inputs' / 'REL' / 'release1HS.shp']

    # call function to be tested
    cfg['INPUT']['releaseScenario'] = 'release4HS'
    releaseF = avaDir / 'Inputs' / 'REL' / 'release4HS.shp'
    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getInputDataCom1DFA(avaDir, cfg)
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
    outFile, available = getInput.getAndCheckInputFiles(avaTestDirInputs, folder, inputType,
        fileExt='shp')

    print('outfile', outFile)
    print('available', available)

    assert available == 'Yes'
    assert 'Inputs/ENT/entrainment1HS.shp' in str(outFile)

    # call function to be tested
    inputFile = avaDirInputs / 'ENT' / 'entrainment1HS.shp'
    testFile = avaTestDirInputs / 'ENT' / 'entrainment1HS2.shp'
    shutil.copyfile(inputFile, testFile)
    with pytest.raises(AssertionError) as e:
        assert getInput.getAndCheckInputFiles(avaTestDirInputs, folder, inputType, fileExt='shp')
    assert str(e.value) == ("More than one %s .shp file in %s/%s/ not allowed" %
                           (inputType, avaTestDirInputs, folder))


def test_getThickness(tmp_path):
    """ test fetching thickness from shapefiles attributes """

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = 'avaHockeyChannel'
    avaDir = dirPath / '..' / 'data' / avaName
    avaDirInputs =  avaDir / 'Inputs'
    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / 'Inputs'
    shutil.copytree(avaDirInputs, avaTestDirInputs)
    outDir = avaTestDir / 'Outputs' / 'com1DFA'
    fU.makeADir(outDir)

    cfg = configparser.ConfigParser()
    cfg = cfgUtils.getModuleConfig(com1DFA, toPrint=False)

    cfg['GENERAL']['relThFromFile'] = 'False'
    cfg['GENERAL']['simTypeList'] = 'null|ent'
    cfg['GENERAL']['secRelAra'] = 'False'
    cfg['INPUT'] = {'releaseScenario': ''}

    demFile = avaTestDirInputs / 'DEM_HS_Topo.asc'
    relFile1 = avaTestDirInputs / 'REL' / 'release1HS.shp'
    relFile2 = avaTestDirInputs / 'REL' / 'release2HS.shp'
    entFile = avaTestDirInputs / 'ENT' / 'entrainment1HS.shp'
    inputSimFiles = {'demFile': demFile, 'relFiles': [relFile1, relFile2], 'entFile': entFile,
        'secondaryReleaseFile': None, 'entResInfo': {'flagRes': 'No', 'flagEnt': 'Yes',
        'flagSecondaryRelease': 'No'}}

    inputSimFiles, cfgFilesRels = getInput.getThickness(inputSimFiles, avaTestDir, com1DFA, cfg)

    print('inputSimFiles', inputSimFiles)
    print('cfgFilesRels', sorted(cfgFilesRels))

    assert inputSimFiles['release1HS']['thickness'] == ['1.0']
    assert inputSimFiles['release2HS']['thickness'] == ['1.0', '1.0']
    assert inputSimFiles['release1HS']['id'] == ['0']
    assert inputSimFiles['release2HS']['id'] == ['0', '1']
    assert inputSimFiles['entrainment1HS']['thickness'] == ['0.3']
    assert inputSimFiles['entrainment1HS']['id'] == ['0']
    assert cfgFilesRels[0].name == 'release1HS_com1DFACfg.ini'
    assert cfgFilesRels[1].name == 'release2HS_com1DFACfg.ini'
    assert len(cfgFilesRels) == 2

    cfgTest1 = configparser.ConfigParser()
    cfgTest1.read(cfgFilesRels[1])

    assert cfgTest1['GENERAL']['relTh'] == ''
    assert cfgTest1['GENERAL'].getboolean('relThFromShp') == True
    assert cfgTest1['GENERAL'].getboolean('relThFromFile') == False
    assert cfgTest1['GENERAL']['entTh'] == ''
    assert cfgTest1['GENERAL'].getboolean('entThFromShp') == True
    assert cfgTest1['INPUT']['releaseScenario'] == 'release2HS'
    assert cfgTest1['INPUT']['relThId'] == '0|1'
    assert cfgTest1['INPUT']['relThThickness'] == '1.0|1.0'
    assert cfgTest1['INPUT']['entrainmentScenario'] == 'entrainment1HS'
    assert cfgTest1['INPUT']['entThId'] == '0'
    assert cfgTest1['INPUT']['entThThickness'] == '0.3'


def test_selectReleaseScenario(tmp_path):
    """ testing selecting a release area scenario according to configuration settings """

    # setup the required inputs
    testPath = pathlib.Path(tmp_path, 'avaTest', 'Inputs', 'REL')
    rel1 = testPath / 'rel1.shp'
    rel2 = testPath / 'rel2.shp'

    inputSimFiles = {'relFiles': [rel1, rel2]}
    cfg = configparser.ConfigParser()
    cfg['INPUT'] = {'releaseScenario': 'rel1'}

    # call function to be tested
    inputSimFiles = getInput.selectReleaseScenario(inputSimFiles, cfg['INPUT'])

    assert inputSimFiles['relFiles'][0].name == 'rel1.shp'
    assert len(inputSimFiles['relFiles']) == 1

    cfg = configparser.ConfigParser()
    cfg['INPUT'] = {'releaseScenario': 'rel2'}
    inputSimFiles = {'relFiles': [rel1, rel2]}

    # call function to be tested
    inputSimFiles = getInput.selectReleaseScenario(inputSimFiles, cfg['INPUT'])


    assert inputSimFiles['relFiles'][0].name == 'rel2.shp'
    assert len(inputSimFiles['relFiles']) == 1


    cfg = configparser.ConfigParser()
    cfg['INPUT'] = {'releaseScenario': 'rel2'}
    inputSimFiles = {'relFiles': [rel1]}


    with pytest.raises(FileNotFoundError) as e:
        assert getInput.selectReleaseScenario(inputSimFiles, cfg['INPUT'])
    assert str(e.value) == ("Release area scenario %s not found - check input data" % ('rel2'))
