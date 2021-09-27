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
import pathlib


def test_generatePlot(tmp_path):

    # Initialise inputs
    avaName = 'avaHockeyChannel'
    avaTestDir = 'avaPlotPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir
    outDir = tmp_path

    data1File = avaDir / 'release1HS_entres_ref_0.15500_pfd.asc'
    data2File = avaDir / 'release2HS_entres_ref_0.15500_pfd.asc'
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
    plotDict = {'relArea' : rel, 'plots': [], 'difference': [], 'differenceZoom': [], 'stats': []}
    plotDictNew = oP.generatePlot(dataDict, avaName, outDir, cfg, plotDict)


    assert (plotDict['difference'][0] == 4.0)
    assert (plotDict['difference'][1] == -0.5)
    assert (plotDict['difference'][2] == -8.0)
    assert (plotDict['differenceZoom'][0] == 0)
    assert (plotDict['differenceZoom'][1] == -0.023255813953488372)
    assert (plotDict['differenceZoom'][2] == -1.0)


def test_quickPlotSimple(tmp_path):
    """ test generating a comparison plot of two raster datasets """

    # setup required input
    testDir = pathlib.Path(__file__).parents[0]
    inputDir = testDir / '..' / '..' / 'benchmarks' / 'avaPlotPytest'
    avaDir = pathlib.Path(tmp_path, 'avaPlots')
    shutil.copytree(inputDir, avaDir)

    cfg = configparser.ConfigParser()
    cfg['FLAGS'] = {'showPlot': 'False'}

    plotDict = oP.quickPlotSimple(avaDir, avaDir, cfg)

    assert (plotDict['difference'][0] == 4.0)
    assert (plotDict['difference'][1] == -0.5)
    assert (plotDict['difference'][2] == -8.0)
