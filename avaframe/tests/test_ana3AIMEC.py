"""Tests for module ana3AIMEC"""
import numpy as np
import pytest
import os
import configparser

# Local imports
import avaframe.ana3AIMEC.ana3AIMEC as ana3AIMEC
import avaframe.in3Utils.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils

def test_analyzeArea(capfd):
    '''Simple test for module analyzeArea'''
    # get input data
    dirname = os.path.dirname(__file__)
    dataRef = os.path.join(dirname, 'data', 'refTestAimecTopo.asc')
    dataSim = os.path.join(dirname, 'data', 'simTestAimecTopo.asc')
    dataMass = os.path.join(dirname, 'data', '000000.txt')
    dataMass1 = os.path.join(dirname, 'data', '000001.txt')
    cfgPath = {}
    cfgPath['pressurefileList'] = [dataRef, dataSim]
    cfgPath['massfileList'] = [dataMass, dataMass1]
    pathResult = os.path.join(dirname, 'data')
    cfgPath['pathResult'] = pathResult
    cfgPath['dirName'] = 'testAIMEC'

    cfg = cfgUtils.getModuleConfig(test_logUtils)
    cfgFlags = cfg['FLAGS']

    pLim = 0.5

    avalData = np.array(([None] * 2))
    data = IOf.readRaster(dataRef)
    avalData[0] = np.transpose(data['rasterData'])
    data = IOf.readRaster(dataSim)
    avalData[1] = np.transpose(data['rasterData'])

    newRasters = {}
    newRasters['newRasterPressure'] = avalData
    newRasters['newRasterDepth'] = avalData
    newRasters['newRasterSpeed'] = avalData
    newRasters['newRasterDEM'] = np.transpose(data['rasterData'])
    rasterTransfo = {}
    rasterTransfo['s'] = np.linspace(0, 499, 500)
    rasterTransfo['l'] = np.linspace(0, 99, 100)
    rasterTransfo['x'] = rasterTransfo['s']
    rasterTransfo['y'] = 50*np.ones(np.shape(rasterTransfo['s']))
    rasterTransfo['rasterArea'] = np.ones((500,100))
    rasterTransfo['indRunoutPoint'] = 400
    rasterTransfo['runoutAngle'] = 10
    resAnalysis = {}
    runout = np.array([100])
    resAnalysis['runout'] = runout

    # testing analyzeFields function
    resAnalysis = ana3AIMEC.analyzeFields(rasterTransfo, pLim, newRasters, cfgPath)

    assert (resAnalysis['runout'][0][0] == 49) and (resAnalysis['runout'][1][1] == 419) and (
    resAnalysis['runout'][2][0] == 50) and (resAnalysis['MMPP'][1] == 1)

    # testing analyzeArea function
    resAnalysis = ana3AIMEC.analyzeArea(rasterTransfo, resAnalysis, pLim, newRasters, cfgPath, cfgFlags)

    assert (resAnalysis['TP'][1] == 800) and (resAnalysis['FN'][1] == 1700) and (
        resAnalysis['FP'][1] == 200) and (resAnalysis['TN'][1] == 7300)
