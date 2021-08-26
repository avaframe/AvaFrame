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
    inputDir = avaDirTmp / 'Inputs'
    resultDir.mkdir(parents=True)
    inputDir.mkdir()
    demInputFile = inputDir / 'avaAlr.asc'
    shutil.copy(peakFile1, peakFileResult1)
    shutil.copy(peakFile2, peakFileResult2)
    shutil.copy(demFile, demInputFile)

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

    # call function to be tested
    plotDict2 = oP.plotAllPeakFields(avaDirTmp, cfg['FLAGS'], modName, demData='')
    plotPath = avaDirTmp / 'Outputs' / 'out1Peak' / 'relAlr_null_ref_0.15500_pfd.png'

    assert 'relAlr_null_ref_0.15500' in plotDict2
    assert plotDict2['relAlr_null_ref_0.15500']['pfd'] == str(plotPath)


def test_plotAllFields(tmp_path):

    # Initialise inputs
    avaName = 'avaAlr'
    avaTestDir = 'avaAlrNullTest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks'/ avaTestDir
    inputDir = avaDir
    outDir = pathlib.Path(tmp_path, 'avaTest')
    outDir.mkdir()

    # call function to be tested
    oP.plotAllFields(avaDir, inputDir, outDir, unit='', constrainData=True)
    plotPath = outDir / 'relAlr_null_ref_0.15500_ppr.png'

    assert plotPath.is_file()

    # call function to be tested
    outDir2 = pathlib.Path(tmp_path, 'avaTest2')
    outDir2.mkdir()
    oP.plotAllFields(avaDir, inputDir, outDir2, unit='', constrainData=False)
    plotPath2 = outDir2 / 'relAlr_null_ref_0.15500_ppr.png'

    assert plotPath2.is_file()
