"""
    Pytest for module ana4Stats

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.ana4Stats import probAna as pA
import pytest
import configparser
import shutil


def test_probAna(tmp_path):
    """ test probAna function to compute mask for parameter exceeding threshold """

    # set input directory
    avaName = 'avaHockey'
    avaTestDir = 'avaHockeyStatsTest'
    dirPath = os.path.dirname(__file__)
    avaDir = os.path.join(dirPath, '..', '..', 'benchmarks', avaTestDir)
    avaDir2 = os.path.join(dirPath, '..', '..', 'benchmarks', avaName)
    inputDir = os.path.join(avaDir, 'ana4Stats')
    outDir = tmp_path

    # set configurations
    testCfg = os.path.join(inputDir, '%sProbAna_com1DFACfg.ini' % avaName)
    cfgMain = cfgUtils.getModuleConfig(com1DFA, testCfg)

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'peakLim': 10.0, 'peakVar': 'ppr'}

    # call function to test
    pA.probAnalysis(avaDir2, cfg, cfgMain, inputDir, outDir)
    probTest = np.loadtxt(os.path.join(tmp_path, 'avaHockey_probMap10.0.asc'), skiprows=6)

    # Load reference solution
    probSol = np.loadtxt(os.path.join(inputDir, 'avaHockey_probMap10.0.txt'), skiprows=6)

    # Compare result to reference solution
    testRes = np.allclose(probTest, probSol, atol=1.e-6)

    # Test
    assert (testRes == True)
