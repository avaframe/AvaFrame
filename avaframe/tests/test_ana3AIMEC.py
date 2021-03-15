''' Tests for module ana3AIMEC '''
import numpy as np
import os
import glob

# Local imports
import avaframe.ana3AIMEC.ana3AIMEC as ana3AIMEC
import avaframe.in2Trans.ascUtils as IOf
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
    cfgPath['projectName'] = 'NameOfAvalanche'
    cfgPath['ppr'] = [dataRef, dataSim]
    cfgPath['mb'] = [dataMass, dataMass1]
    pathResult = os.path.join(dirname, 'data')
    cfgPath['pathResult'] = pathResult
    cfgPath['dirName'] = 'testAIMEC'
    cfgPath['referenceFile'] = 0

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
    rasterTransfo['rasterArea'] = np.ones((500, 100))
    rasterTransfo['indStartOfRunout'] = 400
    rasterTransfo['startOfRunoutAngle'] = 10
    resAnalysis = {}
    runout = np.array([100])
    resAnalysis['runout'] = runout

    # testing analyzeFields function
    resAnalysis = ana3AIMEC.postProcessAIMEC(rasterTransfo, pLim, newRasters, cfgPath, cfgFlags)

    assert (resAnalysis['runout'][0][0] == 449) and (
            resAnalysis['runout'][1][1] == 419) and (
            resAnalysis['runout'][2][0] == 50) and (
            resAnalysis['MMPPR'][1] == 1)

    assert (resAnalysis['TP'][1] == 800) and (
            resAnalysis['FN'][1] == 1700) and (
            resAnalysis['FP'][1] == 200) and (resAnalysis['TN'][1] == 7300)


def test_makeDomainTransfo(capfd):
    '''Simple test for module makeDomainTransfo'''
    # Extract input file locations
    cfgPath = {}
    dir = os.path.dirname(__file__)
    dirname = os.path.join(dir, 'data', 'testAna3Aimec')
    pathData = os.path.join(dirname, 'data')

    profileLayer = glob.glob(os.path.join(dirname, 'LINES', '*aimec*.shp'))
    cfgPath['profileLayer'] = ''.join(profileLayer)

    splitPointLayer = glob.glob(os.path.join(dirname, 'POINTS', '*.shp'))
    cfgPath['splitPointSource'] = ''.join(splitPointLayer)

    demSource = glob.glob(os.path.join(dirname, '*.asc'))
    cfgPath['demSource'] = ''.join(demSource)

    cfgPath['ppr'] = [os.path.join(pathData, 'testAimec_0.asc'), os.path.join(pathData, 'testAimec_1.asc'),
                      os.path.join(pathData, 'testAimec_2.asc'), os.path.join(pathData, 'testAimec_3.asc'),
                      os.path.join(pathData, 'testAimec_4.asc')]
    cfgPath['pfd'] = [os.path.join(pathData, 'testAimec_0.asc'), os.path.join(pathData, 'testAimec_1.asc'),
                      os.path.join(pathData, 'testAimec_2.asc'), os.path.join(pathData, 'testAimec_3.asc'),
                      os.path.join(pathData, 'testAimec_4.asc')]
    cfgPath['pfv'] = [os.path.join(pathData, 'testAimec_0.asc'), os.path.join(pathData, 'testAimec_1.asc'),
                      os.path.join(pathData, 'testAimec_2.asc'), os.path.join(pathData, 'testAimec_3.asc'),
                      os.path.join(pathData, 'testAimec_4.asc')]

    cfgPath['mb'] = [os.path.join(dirname, '000001.txt')]*5

    pathResult = os.path.join(dirname, 'results')
    cfgPath['pathResult'] = pathResult

    cfgPath['projectName'] = 'testAna3Aimec'
    pathName = os.path.basename(profileLayer[0])
    cfgPath['pathName'] = pathName
    cfgPath['dirName'] = 'com1DFA'
    cfgPath['referenceFile'] = 0

    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    cfgFlags['savePlot'] = 'False'
    cfgSetup['startOfRunoutAngle'] = '0'
    cfgSetup['domainWidth'] = '160'
    cfgSetup['pressureLimit'] = '0.9'

    rasterTransfo = ana3AIMEC.makeDomainTransfo(cfgPath, cfgSetup)

    assert rasterTransfo['gridx'][-1, 0] == 60
    assert rasterTransfo['gridx'][-1, -1] == 220
    assert rasterTransfo['gridy'][0, 0] == 180
    assert rasterTransfo['gridy'][0, -1] == 20
    assert rasterTransfo['gridy'][-1, -1] == 258

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    interpMethod = cfgSetup['interpMethod']
    newRasters['newRasterPressure'] = ana3AIMEC.assignData(cfgPath['ppr'], rasterTransfo, interpMethod)
    newRasters['newRasterDepth'] = newRasters['newRasterPressure']
    newRasters['newRasterSpeed'] = newRasters['newRasterPressure']
    newRasterDEM = ana3AIMEC.assignData([cfgPath['demSource']], rasterTransfo,
                                        interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    pressureLimit = float(cfgSetup['pressureLimit'])
    resAnalysis = ana3AIMEC.postProcessAIMEC(rasterTransfo, pressureLimit, newRasters, cfgPath, cfgFlags)

    for i in range(5):
        rasterSource = cfgPath['ppr'][i]
        sourceData = IOf.readRaster(rasterSource)
        rasterdata = sourceData['rasterData']
        error = (resAnalysis['TP'][i]+resAnalysis['FP'][i]-np.nansum(rasterdata))/(
                np.nansum(rasterdata)*100)
        assert error < 0.4
        assert np.abs(resAnalysis['runout'][0, i] - (240 + 10*(i+1))) < 5
