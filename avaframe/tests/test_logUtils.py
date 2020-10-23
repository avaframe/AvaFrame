"""Tests for module logUtilis"""
import avaframe.in3Utils.logUtils as logUtils
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import os


def test_initiateLogger(tmp_path):
    '''Simple test for module initiateLogger'''
    logName = 'testCFG'
    logUtils.initiateLogger(tmp_path, logName)
    logFileName = os.path.join(tmp_path, 'testCFG.log')
    assert os.path.isfile(logFileName)


def test_writeCfg2Log(tmp_path):
    '''Simple test for module writeCfg2Log'''
    dirname = os.path.dirname(__file__)
    avalancheDir = dirname
    logName = 'testCFG'
    logUtils.initiateLogger(tmp_path, logName)
    cfg = cfgUtils.getModuleConfig(test_logUtils)

    logFileName = os.path.join(tmp_path, 'testCFG.log')
    logFileNameRef = os.path.join(avalancheDir, 'data', 'testCFGRef.tog')
    f = open(logFileName).readlines()
    for i in range(3):
        firstLine = f.pop(0)

    fref = open(logFileNameRef).readlines()
    assert f == fref
