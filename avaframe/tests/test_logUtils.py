"""Tests for module logUtilis"""
import avaframe.in3Utils.logUtils as logUtils
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import os
import pathlib


def test_initiateLogger(tmp_path):
    '''Simple test for module initiateLogger'''
    logName = 'testCFG'
    logUtils.initiateLogger(tmp_path, logName)
    # get all files in tmp_path
    pathToLogs = pathlib.Path(tmp_path)
    logFiles = list(pathToLogs.glob('testCFG*.log'))
    assert(len(logFiles) == 1)


def test_writeCfg2Log(tmp_path):
    '''Simple test for module writeCfg2Log'''
    dirname = os.path.dirname(__file__)
    avalancheDir = dirname
    logName = 'testCFG'
    log, logFileName = logUtils.initiateLogger(tmp_path, logName)
    cfg = cfgUtils.getModuleConfig(test_logUtils)
    logFileNameRef = os.path.join(avalancheDir, 'data', 'testCFGRef.tog')

    f = open(logFileName).readlines()
    for i in range(4):
        firstLine = f.pop(0)

    fref = open(logFileNameRef).readlines()
    assert f == fref
