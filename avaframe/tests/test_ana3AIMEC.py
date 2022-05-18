''' Tests for module ana3AIMEC '''
import numpy as np
import pandas as pd
import pathlib

# Local imports
import avaframe.ana3AIMEC.ana3AIMEC as ana3AIMEC
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


def test_getAimecInputs(capfd):
    """ test readAIMECinputs and fetchReferenceSimNo"""
    dirname = pathlib.Path(__file__).parents[0]

    avalancheDir = pathlib.Path(dirname, '..', 'data', 'avaParabola')
    pathDict = {}
    pathDict = aT.readAIMECinputs(avalancheDir, pathDict, dirName='com1DFA')
    print(pathDict)
    assert 'avaframe/tests/../data/avaParabola/Inputs/LINES/path_aimec.shp' in str(pathDict['profileLayer'])
    assert 'avaframe/tests/../data/avaParabola/Inputs/POINTS/splitPoint.shp' in str(pathDict['splitPointSource'])
    assert 'avaframe/tests/../data/avaParabola/Inputs/DEM_PF_Topo.asc' in str(pathDict['demSource'])
    assert 'avaframe/tests/../data/avaParabola/Outputs/ana3AIMEC/com1DFA' in str(pathDict['pathResult'])
    assert pathDict['projectName'] == 'avaParabola'
    assert pathDict['pathName'] == 'path_aimec'


def test_analyzeArea(capfd):
    '''Simple test for module analyzeArea'''
    # get input data
    dirname = pathlib.Path(__file__).parents[0]
    dataRef = dirname / 'data' / 'refTestAimecTopo.asc'
    dataSim = dirname / 'data' / 'simTestAimecTopo.asc'
    dataMass = dirname / 'data' / '000000.txt'
    dataMass1 = dirname / 'data' / '000001.txt'
    pathDict = {}
    pathDict['projectName'] = 'NameOfAvalanche'
    d = {}
    d['simName'] = ['refTestAimecTopo', 'simTestAimecTopo']
    d['ppr'] = [dataRef, dataSim]
    d['pft'] = [dataRef, dataSim]
    d['pfv'] = [dataRef, dataSim]
    d['massBal'] = [dataMass, dataMass1]
    inputsDF = pd.DataFrame(data=d, index=['refTestAimecTopo', 'simTestAimecTopo'])
    pathResult = dirname / 'data'
    pathDict['pathResult'] = pathResult
    pathDict['dirName'] = 'testAIMEC'
    pathDict['refSimulation'] = 'refTestAimecTopo'
    pathDict['compType'] = ['singleModule', 'com1DFA']
    pathDict['contCmap'] = True
    pathDict['resTypeList'] = ['ppr', 'pft', 'pfv']

    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    cfgSetup['resType'] = 'ppr'
    cfgSetup['thresholdValue'] = '0.9'
    cfgSetup['contourLevels'] = '0.1|0.5|1'
    cfgSetup['domainWidth'] = '600'
    cfgSetup['startOfRunoutAreaAngle'] = '10'

    avalData = np.array(([None] * 2))
    data = IOf.readRaster(dataRef)
    avalData[0] = np.transpose(data['rasterData'])
    data = IOf.readRaster(dataSim)
    avalData[1] = np.transpose(data['rasterData'])

    newRasters = {}
    newRasters['newRasterDEM'] = np.transpose(data['rasterData'])
    rasterTransfo = {}
    rasterTransfo['s'] = np.linspace(0, 499, 500)
    rasterTransfo['l'] = np.linspace(0, 99, 100)
    gridy, gridx = np.meshgrid(rasterTransfo['l'], rasterTransfo['s'])
    rasterTransfo['x'] = rasterTransfo['s']
    rasterTransfo['y'] = 50*np.ones(np.shape(rasterTransfo['s']))
    rasterTransfo['gridx'] = gridx
    rasterTransfo['gridy'] = gridy
    rasterTransfo['rasterArea'] = np.ones((500, 100))
    rasterTransfo['indStartOfRunout'] = 400
    rasterTransfo['startOfRunoutAreaAngle'] = 10

    # Analyze data
    # postprocess reference
    inputsDFrow = inputsDF.loc[inputsDF['simName'] == pathDict['refSimulation']].squeeze()
    timeMass = None
    resAnalysisDF = inputsDF[['simName']].copy()
    resAnalysisDF, newRasters, timeMass = ana3AIMEC.postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDFrow,
                                                                     newRasters, timeMass, pathDict['refSimulation'],
                                                                     resAnalysisDF)

    # postprocess other simulations
    for index, inputsDFrow in inputsDF.iterrows():
        simName = inputsDFrow['simName']
        if simName != pathDict['refSimulation']:
            resAnalysisDF, newRasters, timeMass = ana3AIMEC.postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDFrow,
                                                                             newRasters, timeMass, simName,
                                                                             resAnalysisDF)
    print(resAnalysisDF['sRunout'])
    print(resAnalysisDF['xRunout'])
    print(resAnalysisDF['yRunout'])
    assert (resAnalysisDF['sRunout'][0] == 449) and (
            resAnalysisDF['xRunout'][1] == 419) and (
            resAnalysisDF['yRunout'][0] == 31) and (
            resAnalysisDF['maxpprCrossMax'][1] == 1)
    print(resAnalysisDF['TP'])
    print(resAnalysisDF['FN'])
    print(resAnalysisDF['FP'])
    print(resAnalysisDF['TN'])
    print(resAnalysisDF['TN']+resAnalysisDF['FP']+resAnalysisDF['FN']+resAnalysisDF['TP'])
    assert (resAnalysisDF['TP'][1] == 780) and (
            resAnalysisDF['FN'][1] == 1670) and (
            resAnalysisDF['FP'][1] == 200) and (resAnalysisDF['TN'][1] == 7350)


def test_makeDomainTransfo(capfd):
    '''Simple test for module makeDomainTransfo'''
    # Extract input file locations
    pathDict = {}
    dir = pathlib.Path(__file__).parents[0]
    dirname = dir / 'data' / 'testAna3Aimec'
    pathData = dirname / 'data'

    refDir = dirname / 'LINES'
    profileLayer = list(refDir.glob('*aimec*.shp'))
    pathDict['profileLayer'] = profileLayer[0]

    refDir = dirname / 'POINTS'
    splitPointLayer = list(refDir.glob('*.shp'))
    pathDict['splitPointSource'] = splitPointLayer[0]

    demSource = list(dirname.glob('*.asc'))
    pathDict['demSource'] = demSource[0]
    demSource = pathDict['demSource']
    dem = IOf.readRaster(demSource)

    d = {}
    d['simName'] = ['testAimec_0', 'testAimec_1', 'testAimec_2', 'testAimec_3', 'testAimec_4']
    d['ppr'] = [pathData / 'testAimec_0.asc', pathData / 'testAimec_1.asc',
                pathData / 'testAimec_2.asc', pathData / 'testAimec_3.asc',
                pathData / 'testAimec_4.asc']
    d['pft'] = [pathData / 'testAimec_0.asc', pathData / 'testAimec_1.asc',
                pathData / 'testAimec_2.asc', pathData / 'testAimec_3.asc',
                pathData / 'testAimec_4.asc']
    d['pfv'] = [pathData / 'testAimec_0.asc', pathData / 'testAimec_1.asc',
                pathData / 'testAimec_2.asc', pathData / 'testAimec_3.asc',
                pathData / 'testAimec_4.asc']
    d['massBal'] = [dirname / '000001.txt']*5
    inputsDF = pd.DataFrame(data=d, index=['testAimec_0', 'testAimec_1', 'testAimec_2', 'testAimec_3', 'testAimec_4'])
    pathDict['contCmap'] = True

    pathResult = dirname / 'results'
    pathDict['pathResult'] = pathResult

    pathDict['projectName'] = 'testAna3Aimec'
    pathName = profileLayer[0].stem
    pathDict['pathName'] = pathName
    pathDict['dirName'] = 'com1DFA'
    pathDict['refSimulation'] = 'testAimec_0'
    pathDict['compType'] = ['singleModule', 'com1DFA']
    pathDict['resTypeList'] = ['ppr', 'pft', 'pfv']

    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    cfgSetup['startOfRunoutAreaAngle'] = '0'
    cfgSetup['domainWidth'] = '160'
    cfgSetup['resType'] = 'ppr'
    cfgSetup['thresholdValue'] = '0.9'
    cfgSetup['contourLevels'] = '0.1|0.5|1'

    refSimulationName = pathDict['refSimulation']
    refResultSource = inputsDF[inputsDF['simName'] == refSimulationName][cfgSetup['resType']].to_list()[0]
    refRaster = IOf.readRaster(refResultSource)
    refHeader = refRaster['header']
    rasterTransfo = aT.makeDomainTransfo(pathDict, dem, refHeader['cellsize'], cfgSetup)

    assert rasterTransfo['gridx'][-1, 0] == 60
    assert rasterTransfo['gridx'][-1, -1] == 220
    assert rasterTransfo['gridy'][0, 0] == 180
    assert rasterTransfo['gridy'][0, -1] == 20
    assert rasterTransfo['gridy'][-1, -1] == 258

    # transform pressure_data, thickness_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    interpMethod = cfgSetup['interpMethod']
    newRasterDEM = aT.transform(dem, pathDict['demSource'], rasterTransfo, interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM

    # Analyze data
    # postprocess reference
    inputsDFrow = inputsDF.loc[inputsDF['simName'] == pathDict['refSimulation']].squeeze()
    timeMass = None
    resAnalysisDF = inputsDF[['simName']].copy()
    resAnalysisDF, newRasters, timeMass = ana3AIMEC.postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDFrow,
                                                                     newRasters, timeMass, pathDict['refSimulation'],
                                                                     resAnalysisDF)

    # postprocess other simulations
    for index, inputsDFrow in inputsDF.iterrows():
        simName = inputsDFrow['simName']
        if simName != pathDict['refSimulation']:
            resAnalysisDF, newRasters, timeMass = ana3AIMEC.postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDFrow,
                                                                             newRasters, timeMass, simName,
                                                                             resAnalysisDF)

    for i in range(5):
        rasterSource = inputsDF['ppr'][i]
        sourceData = IOf.readRaster(rasterSource)
        rasterdata = sourceData['rasterData']
        error = (resAnalysisDF['TP'][i]+resAnalysisDF['FP'][i]-np.nansum(rasterdata))/(
                np.nansum(rasterdata)*100)
        assert error < 0.4
        assert np.abs(resAnalysisDF['sRunout'][i] - (240 + 10*(i+1))) < 5

    # resAnalysis = ana3AIMEC.postProcessAIMECIndi(rasterTransfo, newRasters, cfgSetup, pathDict)
    # for i in range(5):
    #     rasterSource = pathDict['ppr'][i]
    #     sourceData = IOf.readRaster(rasterSource)
    #     rasterdata = sourceData['rasterData']
    #     error = (resAnalysis['TP'][i]+resAnalysis['FP'][i]-np.nansum(rasterdata))/(
    #             np.nansum(rasterdata)*100)
    #     assert error < 0.4
    #     assert np.abs(resAnalysis['runout'][0, i] - (240 + 10*(i+1))) < 5
