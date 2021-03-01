"""
    Pytest for module in1Data

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.in1Data import getInput
import pytest
import configparser
import shutil


def test_getInputData(tmp_path):
    """ test check for input data """

    # get input data
    dirPath = os.path.dirname(__file__)
    avaName = 'avaHockeyChannel'
    avaDir  = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, 'Inputs')
    avaData = os.path.join(dirPath, '..', 'data', avaName, 'Inputs')
    shutil.copytree(avaData, avaInputs)

    # flags for entrainment and resistance
    flagEnt = True
    flagRes = True

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'flagEnt': 'True', 'flagRes': 'True', 'flagDev': 'False', 'releaseScenario': ''}
    cfgGen = cfg['GENERAL']

    # call function to be tested
    dem, rels, ent, res, flagEntRes = getInput.getInputData(avaDir, cfgGen)
    # second option
    cfg['GENERAL']['releaseScenario'] = 'release1HS'
    dem2, rels2, ent2, res2, flagEntRes2 = getInput.getInputData(avaDir, cfgGen)
    # Test
    assert dem == os.path.join(avaDir, 'Inputs', 'DEM_HS_Topo.asc')
    assert rels == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp'), os.path.join(avaDir, 'Inputs', 'REL', 'release2HS.shp'), os.path.join(avaDir, 'Inputs', 'REL', 'release3HS.shp')]
    assert rels2 == [os.path.join(avaDir, 'Inputs', 'REL', 'release1HS.shp')]
    assert res == ''
    assert ent == os.path.join(avaDir, 'Inputs', 'ENT', 'entrainment1HS.shp')
    assert flagEntRes == True
