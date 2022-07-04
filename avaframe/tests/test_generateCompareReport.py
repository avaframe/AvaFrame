"""
    Pytest for generate Report

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.log2Report import generateCompareReport as gR

import pytest
import configparser
import pathlib


def test_makeLists():
    """ test creating lists for table creation """

    # setup required input
    simDict = {'type': 'list', 'maxP': '100.', 'minF': 0.0, 'testKey': 10., 'testIgnore': ''}
    benchDict = {'type': 'list', 'maxP': '200.', 'minF': 0.0, 'testIgnore2': 'A'}

    # call function to be tested
    parameterList, valuesSim, valuesBench = gR.makeLists(simDict, benchDict)

    assert parameterList == ['maxP', 'minF', 'testKey']
    assert valuesSim == ['100.', 0.0, 10.]
    assert valuesBench == ['200.', 0.0, 'non existent']


def test_copyQuickPlots(tmp_path):
    """ test copy plots """

    # setup required input
    outDir = pathlib.Path(tmp_path)
    dirPath = pathlib.Path(__file__).parents[0]
    pathPlot = dirPath / 'data' / 'release1HX_ent_dfa_3f8d2d9327_pft.png'
    plotListRep = {'pft': pathPlot}
    avaName = 'avaTest'
    testName = 'avaTestName'

    # call function to be tested
    plotPaths = gR.copyQuickPlots(avaName, testName, outDir, plotListRep, rel='')
    testFile = outDir / 'avaTestName__pft.png'

    assert str(plotPaths['pft']) == str(testFile)
    assert len(plotPaths) == 1


def test_copyAimecPlots(tmp_path):
    """ test copying aimec plots to report location """

    # setup required path
    outDir = pathlib.Path(tmp_path)
    dirPath = pathlib.Path(__file__).parents[0]
    pathPlot = dirPath / 'data' / 'release1HX_ent_dfa_3f8d2d9327_pft.png'
    plotFiles = [{'aimec plot 1': {'sim1': pathPlot}}, {'aimec plot 2': {'sim2': pathPlot}}]
    testName = 'avaTest'
    plotPaths = {}

    # call function to be tested
    plotPaths = gR.copyAimecPlots(plotFiles, testName, outDir, plotPaths)
    testFile = outDir / 'avaTest__release1HX_ent_dfa_3f8d2d9327_pft.png'
    print(plotPaths)
    assert str(plotPaths['aimec plot 1']) == str(testFile)
    assert len(plotPaths) == 2


def test_writeCompareReport(tmp_path):
    """ test writing report file """

    # setup required input
    reportFile = pathlib.Path(tmp_path, 'reportTest.md')
    avaName = 'avaTest'
    plot1 = reportFile / 'testplot.png'
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'diffLim': '0.01', 'perDiff': '0.01'}
    reportD = {'avaName': {'type': 'avaName', 'name': 'data/avaKotST'},
        'simName': {'type': 'simName', 'name': 'relKot_null_dfa_a50bf70f08'},
        'time': {'type': 'time', 'time': '02/09/2021 13:20:17'},
        'Simulation Parameters': {'type': 'list',
                                  'Program version': 'development',
                                  'Mu': '0.15500',
                                  'Density [kgm-3]': '200'},
        'Simulation Results': {'ppr': plot1,
                               'Aimec comparison of mean and max values along path': plot1},
        'Aimec analysis': {'type': 'list', 'runout [m]': 16.,
                           'max peak pressure [kPa]': 361.,
                           'max peak flow thickness [m]': 5.4,
                           'max peak flow velocity [ms-1]': 42.50},
        'Simulation Difference': {'ppr': [243., -2.2, -166.],
                                  'pft': [0.7, -0.01, -0.9],
                                  'pfv': [34., -0.2, -28.]},
        'Simulation Stats': {'ppr': [364., 0.0], 'pft': [5., 0.0], 'pfv': [42., 0.0]}}
    benchD = {'avaName': {'type': 'avaName', 'name': 'data/avaKotST'},
        'simName': {'type': 'simName', 'name': 'relKot_null_dfa_0.15500'},
        'time': {'type': 'time', 'time': '02/09/2021 13:19:50'},
        'Simulation Parameters': {'type': 'list',
                                  'Mu': '0.15500',
                                  'Stop criterion': 'end time reached: 112.40 s',
                                  'Avalanche run time [s]': 112.4},

        'Release Area': {'type': 'columns', 'Release area scenario': 'relKot',
                         'Release features': ['KoT'],
                         'Release thickness [m]': [1.0]},
        'Simulation Results': {'pft': plot1,
                               'type': 'image'},
        'Aimec analysis': {'type': 'list', 'runout [m]': 16.,
                           'max peak pressure [kPa]': 358.,
                           'max peak flow thickness [m]': 5.49,
                           'max peak flow velocity [ms-1]': 42.},
        'Test Info': {'type': 'text',
                      'Test Info': 'Compare com1DFAOrig (Reference) to com1DFA (Simulation) results.'}}

    # call function to be tested
    gR.writeCompareReport(reportFile, reportD, benchD, avaName, cfg)

    # Load simulation report
    reportFile = open(os.path.join(tmp_path, 'reportTest.md'), 'r')
    lines = reportFile.readlines()
    lineVals = []
    for line in lines:
        lineVals.append(line)

    reportFile.close()

    assert lineVals[0] == '## *avaTest* \n'
    assert lineVals[1] == '### Simulation name: *relKot_null_dfa_a50bf70f08* \n'
    assert '#### <span style="color:red"> Reference simulation name is different' in lineVals[2]
    assert 'relKot_null_dfa_0.15500  </span> \n' in lineVals[2]
    assert lineVals[4] == '#### Test Info \n'
    assert lineVals[6] == 'Compare com1DFAOrig (Reference) to com1DFA (Simulation) results. \n'
    assert lineVals[8] == '#### Simulation Parameters \n'
    assert lineVals[10] == '| Parameter | Reference | Simulation | \n'
    assert lineVals[12] == '| Program version | non existent | <span style="color:red"> development </span> | \n'
    assert lineVals[17] == '#### Aimec Analysis \n'
    assert lineVals[20] == '| --------- | --------- | ---------- | \n'
    assert lineVals[21] == '| runout [m] | 16.0 | 16.0 | \n'
    assert lineVals[24] == '| max peak flow velocity [ms-1] | 42.0 | <span style="color:red"> 42.5 </span> | \n'
    assert lineVals[27] == '#### Comparison Plots \n'
    assert lineVals[29] == '##### Figure:   ppr \n'
    assert ' <span style="color:red"> Warning absolute difference' in lineVals[31]
    assert 'exceeds the tolerance of 1 percent of ppr-max value' in lineVals[31]
    assert ' </span> \n' in lineVals[31]
    assert ' <span style="color:red"> Difference is: Max = 243.00,' in lineVals[32]
    assert 'Mean = -2.20 and Min = -166.00 </span> \n' in lineVals[32]
    assert lineVals[34] == '![ppr](testplot.png) \n'
    assert lineVals[37] == '##### Figure:   Aimec comparison of mean and max values along path \n'
    assert lineVals[39] == '![Aimec comparison of mean and max values along path](testplot.png) \n'
