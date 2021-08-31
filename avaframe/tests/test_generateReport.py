"""
    Pytest for generate Report

    This file is part of Avaframe.

 """

#  Load modules
import numpy as np
import os
from avaframe.log2Report import generateReport as gR

import pytest
import configparser


def test_writeReport(tmp_path):
    """ Test generation markdown report file """

    # Initialise input dictionary
    testDict = {'title info': {'name': 'This is my report title', 'type': 'title',
                'testNo': 'this should not appear in the report if working fine!'},
        'avaName': {'type': 'avaName', 'name': 'data/avaTest'},
        'time': {'type': 'time', 'time': '1/01/2021 12:46:56'},
        'simName': {'type': 'simName', 'name': 'this is my simulation name'},
        'Testing Parameters': {'release area': 'release1HS2', 'Mu': '0.155', 'type': 'list'},
        'Release area': {'release area scenario': 'my test release area', 'type': 'columns'},
        'Special Input': {'Additional snow-covered area': 'Area 1', 'type': 'columns'},
        'Additional Info': {'type': 'text', 'Report Info': 'This simulations were performed with a developement version.'},
        'Fancy plots': {'type': 'image', 'Peak Pressure Field of my test' : 'release1HS2_entres_dfa_0.750_pfd.png'}}

    # make a list of input dictionaries
    reportDictList = [testDict]

    # initialise cfg object and test directory
    cfg = configparser.ConfigParser()
    cfg['FLAGS'] = {'reportOneFile': True}
    plotDict = {'this is my simulation name': {'pfd': 'testPath/pfd.png'}}
    # Call function to be tested
    gR.writeReport(tmp_path, reportDictList, cfg['FLAGS'], plotDict)

    # Load simulation report
    reportFile = open(os.path.join(tmp_path, 'fullSimulationReport.md'), 'r')
    lines = reportFile.readlines()
    lineVals = []
    for line in lines:
        lineVals.append(line)

    # Test
    assert lineVals[0] == '# This is my report title \n'
    assert lineVals[2] == '### Avalanche directory: *data/avaTest* \n'
    assert lineVals[4] == '### Date: 1/01/2021 12:46:56 \n'
    assert lineVals[11] == '| ---------- | ------ | \n'
    assert lineVals[12] == '| release area | release1HS2 | \n'
    assert lineVals[23] == '| Additional snow-covered area | \n'
    assert lineVals[24] == '| ----------| \n'
    assert lineVals[-3] == '![pfd](testPath/pfd.png) \n'
    assert lineVals[37] == '![Peak Pressure Field of my test](release1HS2_entres_dfa_0.750_pfd.png) \n'
