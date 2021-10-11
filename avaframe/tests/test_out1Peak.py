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
    avaTestDir1 = 'avaAlrNullTest'
    avaTestDir2 = 'avaAlrNullTest2'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks'/ avaTestDir1
    peakFile1 = avaDir / 'relAlr_null_ref_0.15500_pfd.asc'
    peakFile2 = avaDir / 'relAlr_null_ref_0.15500_pfv.asc'
    demFile = dirPath / '..' / 'data' / 'avaAlr' / 'Inputs' / 'avaAlr.asc'
    avaDirTmp1 = pathlib.Path(tmp_path, avaTestDir1)
    resultDir1 = avaDirTmp1 / 'Outputs' / 'com1DFA' / 'peakFiles'
    peakFileResult1 = resultDir1 / 'relAlr_null_ref_0.15500_pfd.asc'
    peakFileResult2 = resultDir1 / 'relAlr_null_ref_0.15500_pfv.asc'
    inputDir1 = avaDirTmp1 / 'Inputs'
    resultDir1.mkdir(parents=True)
    inputDir1.mkdir()
    demInputFile1 = inputDir1 / 'avaAlr.asc'
    shutil.copy(peakFile1, peakFileResult1)
    shutil.copy(peakFile2, peakFileResult2)
    shutil.copy(demFile, demInputFile1)

    avaDirTmp2 = pathlib.Path(tmp_path, avaTestDir2)
    resultDir2 = avaDirTmp2 / 'Outputs' / 'com1DFA' / 'peakFiles'
    peakFileResult1 = resultDir2 / 'relAlr_null_ref_0.15500_pfd.asc'
    peakFileResult2 = resultDir2 / 'relAlr_null_ref_0.15500_pfv.asc'
    inputDir2 = avaDirTmp2 / 'Inputs'
    resultDir2.mkdir(parents=True)
    inputDir2.mkdir()
    demInputFile2 = inputDir2 / 'avaAlr.asc'
    shutil.copy(peakFile1, peakFileResult1)
    shutil.copy(peakFile2, peakFileResult2)
    shutil.copy(demFile, demInputFile2)

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
    plotDict = oP.plotAllPeakFields(avaDirTmp1, cfg['FLAGS'], modName, demData=demData)
    plotPath = avaDirTmp1 / 'Outputs' / 'out1Peak' / 'relAlr_null_ref_0.15500_pfd.png'

    assert 'relAlr_null_ref_0.15500' in plotDict
    assert plotDict['relAlr_null_ref_0.15500']['pfd'] == plotPath

    # call function to be tested
    plotDict2 = oP.plotAllPeakFields(avaDirTmp2, cfg['FLAGS'], modName, demData='')
    plotPath = avaDirTmp2 / 'Outputs' / 'out1Peak' / 'relAlr_null_ref_0.15500_pfd.png'
    print(plotDict2)
    assert 'relAlr_null_ref_0.15500' in plotDict2
    assert plotDict2['relAlr_null_ref_0.15500']['pfd'] == plotPath


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
