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
    pathDict['ppr'] = [dataRef, dataSim]
    pathDict['massBal'] = [dataMass, dataMass1]
    pathResult = dirname / 'data'
    pathDict['pathResult'] = pathResult
    pathDict['dirName'] = 'testAIMEC'
    pathDict['referenceFile'] = 0
    pathDict['compType'] = ['singleModule', 'com1DFA']
    pathDict['numSim'] = 2
    pathDict['contCmap'] = True

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
    resAnalysis = ana3AIMEC.postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags)

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

    pathDict['ppr'] = [pathData / 'testAimec_0.asc', pathData / 'testAimec_1.asc',
                      pathData / 'testAimec_2.asc', pathData / 'testAimec_3.asc',
                      pathData / 'testAimec_4.asc']
    pathDict['pfd'] = [pathData / 'testAimec_0.asc', pathData / 'testAimec_1.asc',
                      pathData / 'testAimec_2.asc', pathData / 'testAimec_3.asc',
                      pathData / 'testAimec_4.asc']
    pathDict['pfv'] = [pathData / 'testAimec_0.asc', pathData / 'testAimec_1.asc',
                      pathData / 'testAimec_2.asc', pathData / 'testAimec_3.asc',
                      pathData / 'testAimec_4.asc']

    pathDict['massBal'] = [dirname / '000001.txt']*5

    pathDict['contCmap'] = True

    pathResult = dirname / 'results'
    pathDict['pathResult'] = pathResult

    pathDict['projectName'] = 'testAna3Aimec'
    pathName = profileLayer[0].stem
    pathDict['pathName'] = pathName
    pathDict['dirName'] = 'com1DFA'
    pathDict['referenceFile'] = 0
    pathDict['compType'] = ['singleModule', 'com1DFA']

    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    cfgFlags['savePlot'] = 'False'
    cfgSetup['startOfRunoutAreaAngle'] = '0'
    cfgSetup['domainWidth'] = '160'
    cfgSetup['resType'] = 'ppr'
    cfgSetup['thresholdValue'] = '0.9'
    cfgSetup['contourLevels'] = '0.1|0.5|1'
    pathDict['numSim'] = 5

    rasterTransfo = aT.makeDomainTransfo(pathDict, cfgSetup)

    assert rasterTransfo['gridx'][-1, 0] == 60
    assert rasterTransfo['gridx'][-1, -1] == 220
    assert rasterTransfo['gridy'][0, 0] == 180
    assert rasterTransfo['gridy'][0, -1] == 20
    assert rasterTransfo['gridy'][-1, -1] == 258

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    interpMethod = cfgSetup['interpMethod']
    newRasters['newRasterPPR'] = aT.assignData(pathDict['ppr'], rasterTransfo, interpMethod)
    newRasters['newRasterPFD'] = newRasters['newRasterPPR']
    newRasters['newRasterPFV'] = newRasters['newRasterPPR']
    newRasterDEM = aT.assignData([pathDict['demSource']], rasterTransfo,
                                        interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    resAnalysis = ana3AIMEC.postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags)

    for i in range(5):
        rasterSource = pathDict['ppr'][i]
        sourceData = IOf.readRaster(rasterSource)
        rasterdata = sourceData['rasterData']
        error = (resAnalysis['TP'][i]+resAnalysis['FP'][i]-np.nansum(rasterdata))/(
                np.nansum(rasterdata)*100)
        assert error < 0.4
        assert np.abs(resAnalysis['runout'][0, i] - (240 + 10*(i+1))) < 5

    resAnalysis = ana3AIMEC.postProcessAIMECIndi(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags)
    for i in range(5):
        rasterSource = pathDict['ppr'][i]
        sourceData = IOf.readRaster(rasterSource)
        rasterdata = sourceData['rasterData']
        error = (resAnalysis['TP'][i]+resAnalysis['FP'][i]-np.nansum(rasterdata))/(
                np.nansum(rasterdata)*100)
        assert error < 0.4
        assert np.abs(resAnalysis['runout'][0, i] - (240 + 10*(i+1))) < 5
