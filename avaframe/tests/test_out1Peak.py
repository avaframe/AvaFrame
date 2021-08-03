"""
    Pytest for module out1Peak

 """

#  Load modules
import numpy as np
from avaframe.out1Peak import outPlotAllPeak as oP
import pytest
import configparser
import pathlib
import shutil


def test_plotAllPeakFields(tmp_path):

    # Initialise inputs
    avaName = 'avaAlr'
    avaTestDir = 'avaAlrNullTest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks'/ avaTestDir
    peakFile1 = avaDir / 'relAlr_null_ref_0.15500_pfd.asc'
    peakFile2 = avaDir / 'relAlr_null_ref_0.15500_pfv.asc'
    demFile = dirPath / '..' / 'data' / 'avaAlr' / 'Inputs' / 'avaAlr.asc'
    avaDirTmp = pathlib.Path(tmp_path, avaTestDir)
    resultDir = avaDirTmp / 'Outputs' / 'com1DFA' / 'peakFiles'
    peakFileResult1 = resultDir / 'relAlr_null_ref_0.15500_pfd.asc'
    peakFileResult2 = resultDir / 'relAlr_null_ref_0.15500_pfv.asc'
    resultDir.mkdir(parents=True)
    shutil.copy(peakFile1, peakFileResult1)
    shutil.copy(peakFile2, peakFileResult2)


    # initialise DEM
    demData1 = np.loadtxt(demFile, skiprows=6)
    demData = {}
    demData['header'] = {'noDataValue': -9999}
    demData['rasterData'] = demData1

    # initialise configparser
    cfg = configparser.ConfigParser()
    cfg['FLAGS'] = {'showPlot': False}
    modName = 'com1DFA'

    # call function to be tested
    plotDict = oP.plotAllPeakFields(avaDirTmp, cfg['FLAGS'], modName, demData=demData)
    plotPath = avaDirTmp / 'Outputs' / 'out1Peak' / 'relAlr_null_ref_0.15500_pfd.png'

    assert 'relAlr_null_ref_0.15500' in plotDict
    assert plotDict['relAlr_null_ref_0.15500']['pfd'] == str(plotPath)
