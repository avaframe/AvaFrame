"""
    Pytest for module ana4Stats

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.ana4Stats import getStats
import pytest
import configparser
import shutil


def test_getStats(tmp_path):
    """ test get statistics of data """

    # define input directory
    avaName = 'avaStats'
    avaDirtmp = os.path.join(tmp_path, avaName)
    avaDirPeakFiles = os.path.join(tmp_path, avaName, 'Outputs', 'com1DFA', 'peakFiles')
    avaDirConfigFiles = os.path.join(tmp_path, avaName, 'Outputs', 'com1DFA', 'configurationFiles')
    os.makedirs(avaDirtmp)
    os.makedirs(avaDirPeakFiles)
    os.makedirs(avaDirConfigFiles)

    data1 = np.asarray([[0., 1., 0.], [0., 2., 2.], [0., 4., 1.]])
    data2 = np.asarray([[0., 10., 0.], [2., 2., 1.], [1., 4., 1.]])

    ncols = data1.shape[1]
    nrows = data1.shape[0]
    noDataValue = -9999
    xllcenter = 0.0
    yllcenter = 0.0
    cellSize = 1.

    resultArrayList = [data1, data2]
    fNames = [os.path.join(avaDirPeakFiles, 'release1_null_dfa_1000000_ppr.asc'),
              os.path.join(avaDirPeakFiles, 'release1_null_dfa_2000000_ppr.asc')]

    # write configuration files
    cfg1 = configparser.ConfigParser()
    cfg1.optionxform = str
    cfg1['GENERAL'] = {'relTh': '1.0', 'releaseScenario': 'release1'}
    cfg2 = configparser.ConfigParser()
    cfg2.optionxform = str
    cfg2['GENERAL'] = {'relTh': '2.0', 'releaseScenario': 'release1'}
    configurationSettings = [cfg1, cfg2]

    configFiles =  [os.path.join(avaDirConfigFiles, 'release1_null_dfa_1000000.ini'),
                    os.path.join(avaDirConfigFiles, 'release1_null_dfa_2000000.ini')]

    for k in range(2):
        # Open outfile
        with open(fNames[k], 'w') as outFile:

            # write the header and array values to file
            outFile.write("ncols %d\n" % ncols)
            outFile.write("nrows %d\n" % nrows)
            outFile.write("xllcenter %.2f\n" % xllcenter)
            outFile.write("yllcenter %.2f\n" % yllcenter)
            outFile.write("cellsize %.2f\n" % cellSize)
            outFile.write("nodata_value %.2f\n" % noDataValue)

            resultArray = resultArrayList[k]
            M = resultArray.shape[0]
            for m in range(M):
                line = np.array([resultArray[m,:]])
                np.savetxt(outFile, line, fmt='%.16g')

            outFile.close()

        cfg = configurationSettings[k]
        cfg.optionxform = str
        with open(os.path.join(configFiles[k]), 'w') as conf:
            cfg.write(conf)

    # parameter dictionary
    varPar = 'relTh'
    inputDir = avaDirPeakFiles
    peakValues = getStats.extractMaxValues(inputDir, avaDirtmp, varPar, restrictType='ppr', nameScenario='', parametersDict='')

    assert peakValues['release1_null_dfa_1000000']['varPar'] == 1.0
    assert peakValues['release1_null_dfa_2000000']['varPar'] == 2.0
    assert peakValues['release1_null_dfa_1000000']['ppr']['max']== 4.0
    assert peakValues['release1_null_dfa_1000000']['ppr']['mean']== 2.0
    assert peakValues['release1_null_dfa_1000000']['ppr']['min']== 1.0
    assert peakValues['release1_null_dfa_2000000']['ppr']['max'] == 10.0
    assert peakValues['release1_null_dfa_2000000']['ppr']['min'] == 1.0
    assert peakValues['release1_null_dfa_2000000']['ppr']['mean']== 3.0
