"""Tests for module logUtilis"""
import avaframe.in3Utils.logUtils as logUtils
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import os
import pathlib
import logging


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
    log = logUtils.initiateLogger(tmp_path, logName)
    logFileName = log.handlers[-1].baseFilename

    cfg = cfgUtils.getModuleConfig(test_logUtils)

    logFileNameRef = os.path.join(avalancheDir, 'data', 'testCFGRef.tog')

    f = open(logFileName).readlines()
    for i in range(4):
        firstLine = f.pop(0)
    

    fref = open(logFileNameRef).readlines()

    for line,refLine in zip(f, fref):
        line = line[line.find('cfgUtils'):]
        assert(line == refLine)

def test_silentLogger(tmp_path):
    '''Simple test for silentLogger context manager functionality'''
    logName = 'testSilent'
    log = logUtils.initiateLogger(tmp_path, logName)
    originalLevel = logging.INFO
    log.setLevel(originalLevel)
    
    with logUtils.silentLogger():
        assert log.level == logging.ERROR
        log.info("This should not appear in log")
        
    assert log.level == originalLevel
    
    logFile = log.handlers[-1].baseFilename
    with open(logFile, 'r') as f:
        logContent = f.read()
    
    assert "This should not appear in log" not in logContent