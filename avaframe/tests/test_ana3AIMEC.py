''' Tests for module ana3AIMEC '''
import numpy as np
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
    cfgPath = {}
    cfgPath = aT.readAIMECinputs(avalancheDir, cfgPath, dirName='com1DFA')
    print(cfgPath)
    assert 'avaframe/tests/../data/avaParabola/Inputs/LINES/path_aimec.shp' in str(cfgPath['profileLayer'])
    assert 'avaframe/tests/../data/avaParabola/Inputs/POINTS/splitPoint.shp' in str(cfgPath['splitPointSource'])
    assert 'avaframe/tests/../data/avaParabola/Inputs/DEM_PF_Topo.asc' in str(cfgPath['demSource'])
    assert 'avaframe/tests/../data/avaParabola/Outputs/ana3AIMEC/com1DFA' in str(cfgPath['pathResult'])
    assert cfgPath['projectName'] == 'avaParabola'
    assert cfgPath['pathName'] == 'path_aimec'


def test_analyzeArea(capfd):
    '''Simple test for module analyzeArea'''
    # get input data
    dirname = pathlib.Path(__file__).parents[0]
    dataRef = pathlib.Path(dirname, 'data', 'refTestAimecTopo.asc')
    dataSim = pathlib.Path(dirname, 'data', 'simTestAimecTopo.asc')
    dataMass = pathlib.Path(dirname, 'data', '000000.txt')
    dataMass1 = pathlib.Path(dirname, 'data', '000001.txt')
    cfgPath = {}
    cfgPath['projectName'] = 'NameOfAvalanche'
    cfgPath['ppr'] = [dataRef, dataSim]
    cfgPath['massBal'] = [dataMass, dataMass1]
    pathResult = pathlib.Path(dirname, 'data')
    cfgPath['pathResult'] = pathResult
    cfgPath['dirName'] = 'testAIMEC'
    cfgPath['referenceFile'] = 0
    cfgPath['compType'] = ['singleModule', 'com1DFA']
    cfgPath['numSim'] = 2
    cfgPath['contCmap'] = True

    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    cfgFlags['savePlot'] = 'True'
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
    newRasters['newRasterPPR'] = avalData
    newRasters['newRasterPFD'] = avalData
    newRasters['newRasterPFV'] = avalData
    newRasters['newRasterDEM'] = np.transpose(data['rasterData'])
    rasterTransfo = {}
    rasterTransfo['s'] = np.linspace(0, 499, 500)
    rasterTransfo['l'] = np.linspace(0, 99, 100)
    rasterTransfo['x'] = rasterTransfo['s']
    rasterTransfo['y'] = 50*np.ones(np.shape(rasterTransfo['s']))
    rasterTransfo['rasterArea'] = np.ones((500, 100))
    rasterTransfo['indStartOfRunout'] = 400
    rasterTransfo['startOfRunoutAreaAngle'] = 10

    # testing analyzeFields function
    resAnalysis = ana3AIMEC.postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, cfgPath, cfgFlags)

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
    dir = pathlib.Path(__file__).parents[0]
    dirname = pathlib.Path(dir, 'data', 'testAna3Aimec')
    pathData = pathlib.Path(dirname, 'data')

    refDir = pathlib.Path(dirname, 'LINES')
    profileLayer = list(refDir.glob('*aimec*.shp'))
    cfgPath['profileLayer'] = profileLayer[0]

    refDir = pathlib.Path(dirname, 'POINTS')
    splitPointLayer = list(refDir.glob('*.shp'))
    cfgPath['splitPointSource'] = splitPointLayer[0]

    refDir = pathlib.Path(dirname)
    demSource = list(refDir.glob('*.asc'))
    cfgPath['demSource'] = demSource[0]

    cfgPath['ppr'] = [pathlib.Path(pathData, 'testAimec_0.asc'), pathlib.Path(pathData, 'testAimec_1.asc'),
                      pathlib.Path(pathData, 'testAimec_2.asc'), pathlib.Path(pathData, 'testAimec_3.asc'),
                      pathlib.Path(pathData, 'testAimec_4.asc')]
    cfgPath['pfd'] = [pathlib.Path(pathData, 'testAimec_0.asc'), pathlib.Path(pathData, 'testAimec_1.asc'),
                      pathlib.Path(pathData, 'testAimec_2.asc'), pathlib.Path(pathData, 'testAimec_3.asc'),
                      pathlib.Path(pathData, 'testAimec_4.asc')]
    cfgPath['pfv'] = [pathlib.Path(pathData, 'testAimec_0.asc'), pathlib.Path(pathData, 'testAimec_1.asc'),
                      pathlib.Path(pathData, 'testAimec_2.asc'), pathlib.Path(pathData, 'testAimec_3.asc'),
                      pathlib.Path(pathData, 'testAimec_4.asc')]

    cfgPath['massBal'] = [pathlib.Path(dirname, '000001.txt')]*5

    cfgPath['contCmap'] = True

    pathResult = pathlib.Path(dirname, 'results')
    cfgPath['pathResult'] = pathResult

    cfgPath['projectName'] = 'testAna3Aimec'
    pathName = profileLayer[0].stem
    cfgPath['pathName'] = pathName
    cfgPath['dirName'] = 'com1DFA'
    cfgPath['referenceFile'] = 0
    cfgPath['compType'] = ['singleModule', 'com1DFA']

    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    cfgFlags['savePlot'] = 'False'
    cfgSetup['startOfRunoutAreaAngle'] = '0'
    cfgSetup['domainWidth'] = '160'
    cfgSetup['resType'] = 'ppr'
    cfgSetup['thresholdValue'] = '0.9'
    cfgSetup['contourLevels'] = '0.1|0.5|1'
    cfgPath['numSim'] = 5

    rasterTransfo = aT.makeDomainTransfo(cfgPath, cfgSetup)

    assert rasterTransfo['gridx'][-1, 0] == 60
    assert rasterTransfo['gridx'][-1, -1] == 220
    assert rasterTransfo['gridy'][0, 0] == 180
    assert rasterTransfo['gridy'][0, -1] == 20
    assert rasterTransfo['gridy'][-1, -1] == 258

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    interpMethod = cfgSetup['interpMethod']
    newRasters['newRasterPPR'] = aT.assignData(cfgPath['ppr'], rasterTransfo, interpMethod)
    newRasters['newRasterPFD'] = newRasters['newRasterPPR']
    newRasters['newRasterPFV'] = newRasters['newRasterPPR']
    newRasterDEM = aT.assignData([cfgPath['demSource']], rasterTransfo,
                                        interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    resAnalysis = ana3AIMEC.postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, cfgPath, cfgFlags)

    for i in range(5):
        rasterSource = cfgPath['ppr'][i]
        sourceData = IOf.readRaster(rasterSource)
        rasterdata = sourceData['rasterData']
        error = (resAnalysis['TP'][i]+resAnalysis['FP'][i]-np.nansum(rasterdata))/(
                np.nansum(rasterdata)*100)
        assert error < 0.4
        assert np.abs(resAnalysis['runout'][0, i] - (240 + 10*(i+1))) < 5

    resAnalysis = ana3AIMEC.postProcessAIMECIndi(rasterTransfo, newRasters, cfgSetup, cfgPath, cfgFlags)
    for i in range(5):
        rasterSource = cfgPath['ppr'][i]
        sourceData = IOf.readRaster(rasterSource)
        rasterdata = sourceData['rasterData']
        error = (resAnalysis['TP'][i]+resAnalysis['FP'][i]-np.nansum(rasterdata))/(
                np.nansum(rasterdata)*100)
        assert error < 0.4
        assert np.abs(resAnalysis['runout'][0, i] - (240 + 10*(i+1))) < 5
