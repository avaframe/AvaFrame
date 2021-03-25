"""
    Pytest for module outQuickPlot

 """

#  Load modules
import numpy as np
import os
from avaframe.out3Plot import outQuickPlot as oP
from avaframe.in3Utils import cfgUtils
import pytest
import configparser
import shutil


def test_outQuickPlot(tmp_path):
    """ test  computation of mean, max and min of difference of two datasets """

    # Initialise inputs
    avaName = 'avaHockeyChannel'
    avaTestDir = 'avaHockeyChannelEntTest'
    dirPath = os.path.dirname(__file__)
    avaDir = os.path.join(dirPath, '..', '..', 'benchmarks', avaTestDir)
    outDir = tmp_path

    data1File = os.path.join(avaDir, 'release1HS_entres_ref_0.15500_pfd.asc')
    data2File = os.path.join(avaDir, 'release2HS_entres_ref_0.15500_pfd.asc')
    data1 = np.loadtxt(data1File, skiprows=6)
    data2 = np.loadtxt(data2File, skiprows=6)
    cellSize = 5.
    cfg = configparser.ConfigParser()
    cfg['FLAGS'] = {'showPlot': False}

    dataDict = {'data1': data1, 'data2': data2, 'name1': 'release1HS_entres_ref_0.15500_pfd',
                'name2': 'release2HS_entres_ref_0.15500_pfd.asc', 'compareType': 'compToRef',
                'simName': 'release1HS_entres_ref_0.15500', 'suffix': 'pfd', 'cellSize': cellSize, 'unit': 'm'}

    # Initialise plotList
    rel = 'release1HS'
    plotDict = {'relArea' : rel, 'plots': [], 'difference': [], 'stats': []}
    plotDictNew = oP.generatePlot(dataDict, avaName, outDir, cfg, plotDict)


    assert (plotDict['difference'][0] == 7.24068)
    assert (plotDict['difference'][1] == -0.24090283440839708)
    assert (plotDict['difference'][2] == -3.26828)
