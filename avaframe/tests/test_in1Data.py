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
import numpy as np
from avaframe.in1Data import computeFromDistribution as cD


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
    assert inputSimFiles['relFiles'] == [avaDir / 'Inputs' / 'REL' / 'release1HS.shp', avaDir / 'Inputs' / 'REL' / 'release2HS.shp', avaDir / 'Inputs' / 'REL' / 'release3HS.shp']
    assert inputSimFiles2['relFiles'] == [avaDir / 'Inputs' / 'REL' / 'release1HS.shp']
    assert inputSimFiles['resFile'] == None
    assert inputSimFiles['entFile'] == avaDir / 'Inputs' / 'ENT' / 'entrainment1HS.shp'
    assert inputSimFiles['entResInfo']['flagEnt'] == "Yes"
    assert inputSimFiles['entResInfo']['flagRes'] == "No"


def test_computeParameters():
    """ test computing parameters """

    # setup required input
    a = 2
    b = 4
    c = 3

    # call function to be tested
    alpha, beta, mu = cD.computeParameters(a, b, c)

    # test
    muTest = 3.5
    alphaTest = 9
    betaTest = -3

    print('alpa, beta, mu', alpha, beta, mu)

    assert alpha == alphaTest
    assert beta == betaTest
    assert mu == muTest


def test_extractUniform():
    """ test extracting a uniform distribution """

    # setup required input
    a = 10
    c = 105
    sampleSize = 20
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'sampleSize': sampleSize, 'flagMinMax': 'True', 'support': 10000}
    steps = 10000

    # compute the support of the distribution
    x = np.linspace(a, c, steps)
    # call function to be tested
    CDF, CDFInt, sampleVect = cD.extractUniform(a, c, x, cfg['GENERAL'])

    print('CDF', CDF)
    print('sample', sampleVect)
    print('CDF', CDF)

    assert np.array_equal(sampleVect, np.linspace(10,105,20))
    assert len(sampleVect) == sampleSize
    assert len(CDF) == 10000
    assert CDFInt(57.5) == 0.5

    # call function to be tested
    cfg['GENERAL']['flagMinMax'] = 'False'
    CDF, CDFInt, sampleVect = cD.extractUniform(a, c, x, cfg['GENERAL'])

    print('CDF', CDF)
    print('sample', sampleVect)
    print('CDF', CDF)

    assert np.allclose(sampleVect, np.linspace(10,105,22)[1:-1], atol=1.e-6)
    assert len(sampleVect) == sampleSize
    assert len(CDF) == 10000
    assert CDFInt(57.5) == 0.5
