"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
from avaframe.com1DFA import com1DFA
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
import pytest
import configparser
import pathlib
import copy
import pickle
import pandas as pd


def test_prepareInputData():
    """ test preparing input data """

    # setup requuired input data
    inputSimFiles = {'entResInfo': {'flagEnt': 'Yes', 'flagRes': 'No', 'flagSecondaryRelease': 'No'}}
    dirName = pathlib.Path(__file__).parents[0]
    avaDir = dirName / '..' / 'data' / 'avaAlr'
    relFile = avaDir / 'Inputs' / 'REL' / 'relAlr.shp'
    inputSimFiles['releaseScenario'] = relFile
    inputSimFiles['demFile'] = avaDir / 'Inputs' / 'avaAlr.asc'
    inputSimFiles['entFile'] = avaDir / 'Inputs' / 'ENT' / 'entAlr.shp'

    # call function to be tested
    demOri, inputSimLines = com1DFA.prepareInputData(inputSimFiles)

    assert demOri['header']['ncols'] == 417
    assert demOri['header']['nrows'] == 915
    assert inputSimLines['releaseLine']['d0'] == ['None']
    assert inputSimLines['releaseLine']['Start'] == np.asarray([0.])
    assert inputSimLines['releaseLine']['Length'] == np.asarray([33.])
    assert inputSimLines['releaseLine']['Name'] == ['AlR']
    assert inputSimLines['entLine']['d0'] == ['None']
    assert inputSimLines['entLine']['Start'] == np.asarray([0.])
    assert inputSimLines['entLine']['Length'] == np.asarray([48.])
    assert inputSimLines['entLine']['Name'] == ['entAlr']
    assert inputSimLines['resLine'] is None
    assert inputSimLines['entrainmentArea'] == 'entAlr.shp'

    # call function to be tested
    inputSimFiles = {'entResInfo': {'flagEnt': 'No', 'flagRes': 'Yes', 'flagSecondaryRelease': 'No'}}
    dirName = pathlib.Path(__file__).parents[0]
    avaDir = dirName / '..' / 'data' / 'avaParabola'
    relFile = avaDir / 'Inputs' / 'REL' / 'release1PF.shp'
    inputSimFiles['releaseScenario'] = relFile
    inputSimFiles['demFile'] = avaDir / 'Inputs' / 'DEM_PF_Topo.asc'
    inputSimFiles['resFile'] = avaDir / 'Inputs' / 'RES' / 'resistance1PF.shp'
    demOri, inputSimLines = com1DFA.prepareInputData(inputSimFiles)

    print('inputSimLines', inputSimLines)

    assert inputSimLines['entLine'] is None
    assert inputSimLines['resLine']['Start'] == np.asarray([0.])
    assert inputSimLines['resLine']['Length'] == np.asarray([5.])
    assert inputSimLines['resLine']['Name'] == ['']


def test_prepareRelase(tmp_path):
    """ test preparing release areas """

    # setup required inputs
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'secRelArea': 'True', 'relTh': '1.32', 'secondaryRelTh': '2.5'}
    inputSimLines = {}
    inputSimLines['entResInfo'] = {'flagSecondaryRelease': 'Yes'}
    inputSimLines['releaseLine'] = {'d0': ['None', 'None']}
    inputSimLines['secondaryReleaseLine'] = {'d0': ['1.789']}
    rel = pathlib.Path(tmp_path, 'release1PF_test.shp')

    # call function to be tested
    relName, inputSimLines, badName = com1DFA.prepareRelase(cfg, rel, inputSimLines)

    assert relName == 'release1PF_test'
    assert inputSimLines['entResInfo']['flagSecondaryRelease'] == 'Yes'
    assert inputSimLines['releaseLine']['d0'] == [1.32, 1.32]
    assert inputSimLines['secondaryReleaseLine']['d0'] == [1.789]
    assert badName is True

    # setup required inputs
    inputSimLines = {}
    inputSimLines['entResInfo'] = {'flagSecondaryRelease': 'Yes'}
    inputSimLines['releaseLine'] = {'d0': ['1.78', '4.328']}
    inputSimLines['secondaryReleaseLine'] = {'d0': ['None']}
    rel = pathlib.Path(tmp_path, 'release1PF_test.shp')

    # call function to be tested
    relName2, inputSimLines2, badName2 = com1DFA.prepareRelase(cfg, rel, inputSimLines)

    assert relName2 == 'release1PF_test'
    assert inputSimLines2['entResInfo']['flagSecondaryRelease'] == 'Yes'
    assert inputSimLines2['releaseLine']['d0'] == [1.78, 4.328]
    assert inputSimLines2['secondaryReleaseLine']['d0'] == [2.5]
    assert badName2 is True

    # setup required inputs
    inputSimLines = {}
    inputSimLines['entResInfo'] = {'flagSecondaryRelease': 'No'}
    inputSimLines['releaseLine'] = {'d0': ['1.78', '4.328']}
    rel = pathlib.Path(tmp_path, 'release1PF_test.shp')

    # call function to be tested
    with pytest.raises(FileNotFoundError) as e:
        assert com1DFA.prepareRelase(cfg, rel, inputSimLines)
    assert str(e.value) == "No secondary release file found"

    # call function to be tested
    cfg['GENERAL']['secRelArea'] = 'False'
    relName3, inputSimLines3, badName3 = com1DFA.prepareRelase(cfg, rel, inputSimLines)

    assert relName3 == 'release1PF_test'
    assert inputSimLines3['entResInfo']['flagSecondaryRelease'] == 'No'
    assert inputSimLines3['releaseLine']['d0'] == [1.78, 4.328]
    assert inputSimLines3['secondaryReleaseLine'] == None
    assert badName3 is True



def test_createReportDict():
    """ test creating a report dictionary """

    # setup required input
    avaDir = 'data/avaTest'
    logName = 'testName'
    relName = 'relTest'
    inputSimLines = {'entrainmentArea': 'entTest', 'resistanceArea': 'resTest', 'releaseLine':
                    {'Name': 'relTestFeature', 'd0': '1.45'}}
    reportAreaInfo = {'entrainment': 'Yes', 'resistance': 'Yes', 'Release area info':
                      {'Projected Area [m2]': 'm2'}}
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'mu': '0.15500', 'rho': '200.', 'frictModel': 'samosAT', 'hEnt': '0.3',
                      'rhoEnt': '100.0'}

    # call function to be tested
    reportST = com1DFA.createReportDict(avaDir, logName, relName, inputSimLines, cfg['GENERAL'], reportAreaInfo)

    assert 'Simulation Parameters' in reportST
    assert 'Program version' in reportST['Simulation Parameters']
    assert reportST['avaName']['name'] == avaDir
    assert reportST['simName']['name'] == logName
    assert reportST['Simulation Parameters']['Release Area Scenario'] == relName
    assert reportST['Simulation Parameters']['Entrainment'] == 'Yes'
    assert reportST['Simulation Parameters']['Resistance'] == 'Yes'
    assert reportST['Simulation Parameters']['Mu'] == '0.15500'
    assert reportST['Simulation Parameters']['Density [kgm-3]'] == '200.'
    assert reportST['Simulation Parameters']['Friction model'] == 'samosAT'
    assert reportST['Release Area']['Release area scenario'] == relName
    assert reportST['Release Area']['Release Area'] == 'relTestFeature'
    assert reportST['Release Area']['Release thickness [m]'] == '1.45'
    assert reportST['Entrainment area']['Entrainment area scenario'] == 'entTest'
    assert 'Projected Area [m2]' in reportST['Release Area']


def test_reportAddTimeMassInfo():
    """ test adding mass and time info to report dict """

    # setup required input
    reportDict = {'Simulation Parameters': {'testItem': 1.0}}
    tcpuDFA = 400.
    infoDict = {'initial mass': 400000.2345, 'final mass': 400000.8345, 'entrained mass': 0.8,
                'entrained volume': 0.2, 'stopInfo': {'Stop criterion': '0.1 percent of PKE'}}

    # call function to be tested
    reportDict = com1DFA.reportAddTimeMassInfo(reportDict, tcpuDFA, infoDict)

    assert reportDict['Simulation Parameters']['testItem'] == 1.0
    assert reportDict['Simulation Parameters']['Initial mass [kg]'] == '400000.23'
    assert reportDict['Simulation Parameters']['Final mass [kg]'] == '400000.83'
    assert reportDict['Simulation Parameters']['Entrained mass [kg]'] == '0.80'
    assert reportDict['Simulation Parameters']['Entrained volume [m3]'] == '0.20'
    assert reportDict['Simulation Parameters']['Stop criterion'] == '0.1 percent of PKE'


def test_prepareArea():
    """ test converting a polygon from a shape file to a raster """

    # setup required input
    releaseLine = {'Name': ['testRel', 'test2'], 'Start': np.asarray([0., 5]), 'Length': np.asarray([5, 5]),
                   'x': np.asarray([0, 10., 10.0, 0., 0., 20., 26., 26., 20., 20.]),
                   'y': np.asarray([0., 0., 10.0, 10., 0.0, 21., 21., 27., 27., 21.])}
    demHeader = {}
    demHeader['xllcenter'] = 0.0
    demHeader['yllcenter'] = 0.0
    demHeader['cellsize'] = 5.0
    demHeader['noDataValue'] = -9999
    demHeader['nrows'] = 7
    demHeader['ncols'] = 7
    dem = {'header': demHeader}
    dem['rasterData'] = np.ones((demHeader['nrows'], demHeader['ncols']))
    radius = 0.01
    thList = [1.234, 7.8]
    combine = True
    checkOverlap = True

    # call function to be tested
    # test 1
    line = com1DFA.prepareArea(releaseLine, dem, radius, thList='', combine=True, checkOverlap=True)

    # test 2
    releaseLine2 = {'Name': ['testRel', 'test2'], 'Start': np.asarray([0., 5]), 'Length': np.asarray([5, 5]),
                   'x': np.asarray([0, 10., 10.0, 0., 0., 20., 26., 26., 20., 20.]),
                   'y': np.asarray([0., 0., 10.0, 10., 0.0, 21., 21., 27., 27., 21.])}
    line2 = com1DFA.prepareArea(releaseLine2, dem, 0.6, thList=thList, combine=True, checkOverlap=True)

    # test 3
    releaseLine3 = {'Name': ['testRel', 'test2'], 'Start': np.asarray([0., 5]), 'Length': np.asarray([5, 5]),
                   'x': np.asarray([0, 10., 10.0, 0., 0., 5, 15., 15., 5., 5]),
                   'y': np.asarray([0., 0., 10.0, 10., 0.0, 5, 5, 15., 15., 5.])}

    with pytest.raises(AssertionError) as e:
        assert com1DFA.prepareArea(releaseLine3, dem, 0.6, thList=thList, combine=True, checkOverlap=True)
    assert str(e.value) == "Features are overlaping - this is not allowed"

    line5 = com1DFA.prepareArea(releaseLine3, dem, 0.6, thList=thList, combine=True, checkOverlap=False)

    print('line5', line5)



    # test 4
    releaseLine4 = {'Name': ['testRel', 'test2'], 'Start': np.asarray([0., 5]), 'Length': np.asarray([5, 5]),
                    'x': np.asarray([0, 10., 10.0, 0., 0., 20., 26., 26., 20., 20.]),
                    'y': np.asarray([0., 0., 10.0, 10., 0.0, 21., 21., 27., 27., 21.])}
    line4 = com1DFA.prepareArea(releaseLine4, dem, 0.6, thList=thList, combine=False, checkOverlap=True)

    # test results
    testRaster = np.zeros((demHeader['nrows'], demHeader['ncols']))
    testRaster[0:3, 0:3] = 1.0
    testRaster[5, 4:6] = 1.0
    testRaster2 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    testRaster2[0:3, 0:3] = 1.234
    testRaster2[4, 4:6] = 7.8
    testRaster2[5, 4:6] = 7.8
    testRaster4 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    testRaster5 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    testRaster4[0:3, 0:3] = 1.234
    testRaster5[4, 4:6] = 7.8
    testRaster5[5, 4:6] = 7.8
    testList = [testRaster4, testRaster5]
    testRaster6 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    testRaster6[0:3, 0] = 1.234
    testRaster6[0, 0:3] = 1.234
    testRaster6[1:3, 1:3] = (1.234 + 7.8) / 2.
    testRaster6[3, 1:4] = 7.8
    testRaster6[1:4, 3] = 7.8

    print('line Raster', line['rasterData'])
    print('testRaster', testRaster)
    print('line Raster', line2['rasterData'])
    print('testRaster', testRaster2)
    print('test 4 line Raster', line4['rasterData'])
    print('testRaster', testList)
    print('test raster 6', testRaster6)

    assert np.array_equal(line['rasterData'], testRaster)
    assert np.array_equal(line2['rasterData'], testRaster2)
    assert np.array_equal(line4['rasterData'], testList)
    assert np.array_equal(line5['rasterData'], testRaster6)

def test_checkParticlesInRelease():
    """ test if particles are within release polygon and removed if not """

    # setup required input
    releaseLine = {'Name': ['testRel'], 'Start': np.asarray([0.]), 'Length': np.asarray([5]),
                   'x': np.asarray([0, 10., 10.0, 0., 0.]), 'y': np.asarray([0., 0., 10.0, 10., 0.0])}
    demHeader = {}
    demHeader['xllcenter'] = 0.0
    demHeader['yllcenter'] = 0.0
    demHeader['cellsize'] = 1.0
    demHeader['noDataValue'] = -9999
    demHeader['nrows'] = 5
    demHeader['ncols'] = 5
    releaseLine['header'] = demHeader
    particles = {'x': np.asarray([2.4, 9.7, 10.02, 11.5]), 'y': np.asarray([2.4, 9.7, 10.2, 11.5]),
                 'Npart': 4, 'm': np.asarray([1.4, 1.7, 1.4, 1.8])}
    radius = np.sqrt(2)

    # call function to be tested
    particles = com1DFA.checkParticlesInRelease(particles, releaseLine, radius)
    # call function to be tested
    particles2 = {'x': np.asarray([2.4, 9.7, 9.997, 10.09, 0.0]), 'y': np.asarray([2.4, 9.7, 9.994, 10.8, 0.0]),
                 'Npart': 5, 'm': np.asarray([1.4, 1.7, 1.4, 1.8, 1.1])}
    particles2 = com1DFA.checkParticlesInRelease(particles2, releaseLine, 0.01)

    print('particles', particles, particles2)
    assert np.array_equal(particles['x'], np.asarray([2.4, 9.7, 10.02]))
    assert np.array_equal(particles['y'], np.asarray([2.4, 9.7, 10.2]))
    assert particles['mTot'] == (1.4+1.7+1.4)
    assert particles['Npart'] == 3
    assert np.array_equal(particles2['x'], np.asarray([2.4, 9.7, 9.997, 0.0]))
    assert np.array_equal(particles2['y'], np.asarray([2.4, 9.7, 9.994, 0.0]))
    assert particles2['mTot'] == (1.4+1.7+1.4+1.1)
    assert particles2['Npart'] == 4


def test_pointInPolygon():
    """ test if a point is inside polygon   """

    # setup required input
    demHeader = {}
    demHeader['xllcenter'] = 0.0
    demHeader['yllcenter'] = 0.0
    radius = np.sqrt(2)
    line = {'x': np.asarray([0, 10., 10., 0., 0.]), 'y': np.asarray([0., 0., 10., 10., 0.0])}
    points = {'x': np.asarray([2.4, 9.7, -0.04, 10.09, 0.0]), 'y': np.asarray([2.4, 9.7, -0.11, 10.8, 0.0])}

    # call function to be tested
    mask = com1DFA.pointInPolygon(demHeader, points, line, radius)
    # call function to be tested
    mask2 = com1DFA.pointInPolygon(demHeader, points, line, 0.01)

    # test mask
    testMask = np.asarray([True, True, True, False, True])
    testMask2 = np.asarray([True, True, False, False, True])

    assert np.array_equal(mask, testMask)
    assert np.array_equal(mask2, testMask2)

    # call function to be tested
    line = {'x': np.asarray([0, 10., 10., 0.]), 'y': np.asarray([0., 0., 10., 10.])}
    mask3 = com1DFA.pointInPolygon(demHeader, points, line, 0.01)
    testMask3 = np.asarray([True, True, False, False, True])

    assert np.array_equal(mask3, testMask3)


def test_initializeMassEnt():
    """ test initializing entrainment area """

    # setup required input
    nrows = 110
    ncols = 150
    demHeader = {}
    demHeader['xllcenter'] = 0.0
    demHeader['yllcenter'] = 0.0
    demHeader['cellsize'] = 1.0
    demHeader['noDataValue'] = -9999
    demHeader['nrows'] = nrows
    demHeader['ncols'] = ncols
    dem = {'header': demHeader}
    dem['rasterData'] = np.ones((nrows, ncols))

    simTypeActual = 'entres'
    dirName = pathlib.Path(__file__).parents[0]
    fileName = dirName / 'testEnt.shp'
    entLine = {'fileName': fileName, 'Name': ['testEnt'], 'Start': np.asarray([0.]), 'Length': np.asarray([5]),
               'x': np.asarray([0, 10., 10.0, 0., 0.]), 'y': np.asarray([0., 0., 10.0, 10., 0.0])}
    reportAreaInfo = {'entrainment': '', }
    thresholdPointInPoly = 0.001

    # call function to be tested
    entrMassRaster, reportAreaInfo = com1DFA.initializeMassEnt(dem, simTypeActual, entLine, reportAreaInfo, thresholdPointInPoly)
    testData = np.zeros((nrows, ncols))
    testData[0:11, 0:11] = 1.0

    print('data', testData)
    print('ent', entrMassRaster, entLine)

    assert np.array_equal(entrMassRaster, testData)
    assert np.sum(entrMassRaster) == 121
    assert entrMassRaster.shape[0] == nrows
    assert reportAreaInfo['entrainment'] == 'Yes'

    # call function to be tested
    simTypeActual = 'res'
    entrMassRaster, reportAreaInfo = com1DFA.initializeMassEnt(dem, simTypeActual, entLine, reportAreaInfo, thresholdPointInPoly)

    assert np.array_equal(entrMassRaster, np.zeros((nrows, ncols)))
    assert reportAreaInfo['entrainment'] == 'No'


def test_initializeResistance():
    """ test initializing resistance area """

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'dRes': '0.3', 'cw': '0.5', 'sres': '5'}

    nrows = 11
    ncols = 15
    demHeader = {}
    demHeader['nrows'] = nrows
    demHeader['ncols'] = ncols
    demHeader['xllcenter'] = 0.0
    demHeader['yllcenter'] = 0.0
    demHeader['cellsize'] = 1.0
    demHeader['noDataValue'] = -9999
    dem = {'header': demHeader}
    dem['rasterData'] = np.ones((nrows, ncols))

    simTypeActual = 'entres'
    resLine = {'fileName': 'resTest', 'Name': 'resFeature', 'Start': np.asarray([0]), 'Length': np.asarray([5]),
               'Name': ['resTestFeat'],
               'x': np.asarray([0, 10., 10.0, 0., 0.]), 'y': np.asarray([0., 0., 10.0, 10., 0.0])}
    reportAreaInfo = {'entrainment': 'Yes', 'resistance': 'No'}
    thresholdPointInPoly = 0.01

    # call function to be tested
    cResRaster, reportAreaInfo = com1DFA.initializeResistance(cfg['GENERAL'], dem, simTypeActual, resLine, reportAreaInfo, thresholdPointInPoly)
    testArray = np.zeros((nrows, ncols))
    testArray[0:11, 0:11] = 0.003

    print('cResRaster', cResRaster)
    print('reportAreaInfo', reportAreaInfo)

    assert np.array_equal(cResRaster, testArray)
    assert np.sum(cResRaster) == 0.363
    assert reportAreaInfo['resistance'] == 'Yes'


def test_setDEMOriginToZero():
    """ test if origin is set to zero """

    # setup required input
    tHeader = {}
    tHeader['xllcenter'] = 10.
    tHeader['yllcenter'] = 4.0
    dem = {'header': tHeader}

    # call function to be tested
    demTest = com1DFA.setDEMoriginToZero(dem)

    assert demTest['header']['xllcenter'] == 0.0
    assert demTest['header']['yllcenter'] == 0.0


def test_placeParticles():
    """ test placing of particles """

    # setup required inputs
    massCell = 10.
    indx = 0
    indy = 1
    csz = 5
    massPerPart = 2.
    initPartDistType = 'uniform'
    rng = np.random.default_rng(12345)

    # call funciton to be tested - uniform
    xpart, ypart, mPart, nPart = com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest = np.asarray([-1.66666666, 0.0, 1.66666666, -1.66666666, 0., 1.66666666, -1.66666666,
                        0.0, 1.66666666])
    ypartTest = np.asarray([3.33333333, 3.33333333, 3.33333333, 5.0, 5., 5., 6.66666666, 6.66666666,
                        6.66666666])

    # call funciton to be tested - uniform
    massCell = 8.
    xpart2, ypart2, mPart2, nPart2 = com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest2 = np.asarray([-1.25, 1.25, -1.25, 1.25])
    ypartTest2 = np.asarray([3.75, 3.75, 6.25, 6.25])

    # call funciton to be tested - random
    massCell = 11.5
    initPartDistType = 'random'
    xpart3, ypart3, mPart3, nPart3 = com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)
    xpartTest3 = np.asarray([-0.9162083, 1.48682729, 0.88127335, -0.54445225, -0.83593036, 0.49154377])
    ypartTest3 = np.asarray([3.43367093, 5.86378022, 7.20901433, 3.74122857, 7.24440576, 5.83618727])

    assert nPart == 9.0
    assert np.isclose(mPart, 1.111111)
    assert np.allclose(xpart, xpartTest)
    assert np.allclose(ypart, ypartTest)
    assert nPart2 == 4.0
    assert mPart2 == 2.
    assert np.allclose(xpart2, xpartTest2)
    assert np.allclose(ypart2, ypartTest2)
    assert nPart3 == 6.0
    assert np.isclose(mPart3, 1.9166666666666)
    assert np.allclose(xpart3, xpartTest3)
    assert np.allclose(ypart3, ypartTest3)

    # call funciton to be tested - random
    massCell = 8
    initPartDistType = 'semiRandom'
    csz = 4
    xpart4, ypart4, mPart4, nPart4 = com1DFA.placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType)

    print('xpart4', xpart4)
    print('ypart4', ypart4)

    assert nPart4 == 4.0
    assert mPart2 == 2.
    assert -2.0 < xpart4[0] < 0.0
    assert 2.0 < ypart4[0] < 4.0
    assert 0.0 < xpart4[1] < 2.0
    assert 2.0 < ypart4[1] < 4.0
    assert -2.0 < xpart4[2] < 0.0
    assert 4.0 < ypart4[2] < 6.0
    assert 0.0 < xpart4[3] < 2.0
    assert 4.0 < ypart4[3] < 6.0


def test_initializeMesh():
    """ test mesh initialization """

    # setup required input
    demHeader = {}
    demHeader['xllcenter'] = 101.23
    demHeader['yllcenter'] = 24.54
    demHeader['cellsize'] = 1.0
    demHeader['noDataValue'] = -9999
    demHeader['nrows'] = 5
    demHeader['ncols'] = 5

    # define plane with constant slope of 45°
    demData = np.asarray([[1., 2., 3., 4., np.nan],
                          [1., 2., 3., 4., np.nan],
                          [1., 2., 3., 4., np.nan],
                          [1., 2., 3., 4., np.nan],
                          [1., 2., 3., 4., np.nan]])

    demOri = {'header': demHeader, 'rasterData': demData}
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'sphKernelRadius': '0.5', 'meshCellSizeThreshold': '0.0001',
                      'meshCellSize': '1.'}
    num = 1

    # setup testResults
    demNewHeader = {}
    demNewHeader['xllcenter'] = 0.
    demNewHeader['yllcenter'] = 0.
    demNewHeader['cellsize'] = 1.0
    demNewHeader['noDataValue'] = -9999
    demNewHeader['nrows'] = 5
    demNewHeader['ncols'] = 5
    demTest = {'header': demNewHeader}
    demTest['originOri'] = {'xllcenter': 101.23, 'yllcenter': 24.54}
    demTest['outOfDEM'] = np.asarray([[False, False, False, False, True],
                                      [False, False, False, False, True],
                                      [False, False, False, False, True],
                                      [False, False, False, False, True],
                                      [False, False, False, False, True]])
    # normal vector of plane
    demTest['Nx'] = np.asarray([[0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.]]) - 25.
    demTest['Ny'] = np.asarray([[0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.]])
    demTest['Nz'] = np.asarray([[0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0.]]) + 25.
    # setup neighbour grid
    headerNeighbourGrid = {}
    headerNeighbourGrid['cellsize'] = 0.5
    headerNeighbourGrid['ncols'] = 10
    headerNeighbourGrid['nrows'] = 10
    headerNeighbourGrid['xllcenter'] = 0
    headerNeighbourGrid['yllcenter'] = 0
    demTest['headerNeighbourGrid'] = headerNeighbourGrid
    areaCell = 1 / np.cos(np.deg2rad(45))
    demTest['areaRaster'] = np.zeros((5, 5)) + areaCell
    demTest['rasterData'] = demData

    # call function to be tested
    demOri, dem = com1DFA.initializeMesh(cfg['GENERAL'], demOri, num)

    assert dem['header']['xllcenter'] == demTest['header']['xllcenter']
    assert dem['header']['yllcenter'] == demTest['header']['yllcenter']
    assert dem['header']['ncols'] == demTest['header']['ncols']
    assert dem['header']['nrows'] == demTest['header']['nrows']
    assert dem['header']['cellsize'] == demTest['header']['cellsize']
    assert dem['header']['yllcenter'] == demTest['header']['yllcenter']
    assert np.array_equal(dem['rasterData'][0:4, 0:4], demTest['rasterData'][0:4, 0:4])
    assert np.all(np.isnan(dem['rasterData'][0:5, 4]))
    assert abs(dem['Nx'][2, 2]) == abs(dem['Nz'][2, 2])
    assert np.isclose(dem['areaRaster'][2, 2], demTest['areaRaster'][2, 2])
    assert dem['headerNeighbourGrid']['xllcenter'] == demTest['headerNeighbourGrid']['xllcenter']
    assert dem['headerNeighbourGrid']['yllcenter'] == demTest['headerNeighbourGrid']['yllcenter']
    assert dem['headerNeighbourGrid']['ncols'] == demTest['headerNeighbourGrid']['ncols']
    assert dem['headerNeighbourGrid']['nrows'] == demTest['headerNeighbourGrid']['nrows']
    assert dem['headerNeighbourGrid']['cellsize'] == demTest['headerNeighbourGrid']['cellsize']
    assert dem['headerNeighbourGrid']['yllcenter'] == demTest['headerNeighbourGrid']['yllcenter']


def test_polygon2Raster():
    """ test if polygon is converted to raster properly """

    # setup required inputs
    demHeader = {}
    demHeader['cellsize'] = 1
    demHeader['ncols'] = 10
    demHeader['nrows'] = 10
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0

    Line = {'x': np.asarray([0, 1., 0.989, 0., 0.]), 'y': np.asarray([0., 0., 0.989, 1., 0.0])}
    radius = 0.0001
    th = 1.2

    # call function to be tested
    Mask = com1DFA.polygon2Raster(demHeader, Line, radius, th=th)

    # setup test output
    maskTest = np.zeros((10, 10))
    maskTest[0, 0:2] = 1.2
    maskTest[1, 0] = 1.2

    # call function to be tested
    Mask2 = com1DFA.polygon2Raster(demHeader, Line, 0.1, th=th)
    maskTest2 = maskTest.copy()
    maskTest2[1, 1] = 1.2

    assert np.array_equal(maskTest, Mask)
    assert np.array_equal(maskTest2, Mask2)

    # call function to be tested
    Line = {'x': np.asarray([0, 1., 0.989, 0.]), 'y': np.asarray([0., 0., 0.989, 1.])}
    Mask3 = com1DFA.polygon2Raster(demHeader, Line, 0.1, th=th)
    maskTest3 = maskTest.copy()
    maskTest3[1, 1] = 1.2

    assert np.array_equal(maskTest3, Mask3)


def test_getSimTypeList():
    """ test create list of simTypes """

    # setup required input
    simTypeList = ['ent', 'res', 'null', 'available', 'entres']
    inputSimFiles = {'entResInfo': {'flagEnt': 'Yes', 'flagRes': 'Yes'}}

    # call function to be tested
    simTypeList = com1DFA.getSimTypeList(simTypeList, inputSimFiles)

    # setup test result
    simTypeListTest = ['ent', 'null', 'res', 'entres']

    assert set(simTypeListTest).issubset(simTypeList)
    assert 'available' not in simTypeList

    # call function to be tested
    simTypeList = ['ent', 'null', 'available']
    inputSimFiles['entResInfo']['flagRes'] = 'No'
    simTypeList2 = com1DFA.getSimTypeList(simTypeList, inputSimFiles)

    # setup test result
    simTypeListTest2 = ['ent', 'null']

    assert set(simTypeListTest2).issubset(simTypeList2)
    assert 'available' not in simTypeList2
    assert 'entres' not in simTypeList2
    assert 'res' not in simTypeList2

    # call function to be tested
    simTypeList = ['res', 'null', 'available']
    inputSimFiles['entResInfo']['flagEnt'] = 'No'
    inputSimFiles['entResInfo']['flagRes'] = 'Yes'
    simTypeList3 = com1DFA.getSimTypeList(simTypeList, inputSimFiles)

    # setup test result
    simTypeListTest3 = ['res', 'null']

    assert set(simTypeListTest3).issubset(simTypeList3)
    assert 'available' not in simTypeList3
    assert 'entres' not in simTypeList3
    assert 'ent' not in simTypeList3

    # call function to be tested
    simTypeList = ['ent', 'null', 'available', 'entres', 'res']
    inputSimFiles['entResInfo']['flagEnt'] = 'Yes'
    inputSimFiles['entResInfo']['flagRes'] = 'No'
    with pytest.raises(FileNotFoundError) as e:
        assert com1DFA.getSimTypeList(simTypeList, inputSimFiles)
    assert str(e.value) == "No resistance file found"

    # call function to be tested
    inputSimFiles['entResInfo']['flagEnt'] = 'No'
    inputSimFiles['entResInfo']['flagRes'] = 'Yes'
    with pytest.raises(FileNotFoundError) as e:
        assert com1DFA.getSimTypeList(simTypeList, inputSimFiles)
    assert str(e.value) == "No entrainment file found"


def test_appendFieldsParticles():
    """ test if correct fields and particles list are created for export """

    # setup required input
    fieldsListIn = [{'ppr': np.zeros((3, 3)), 'pfv': np.zeros((3, 3))}]
    particlesListIn = [{'x': np.asarray([0., 4., 0.]), 'y': np.asarray([0., 4., 0.]), 'm': np.asarray([0., 4., 0.])}]
    particles = {'x': np.asarray([0., 5., 0.]), 'y': np.asarray([0., 5., 0.]), 'm': np.asarray([0., 4., 0.])}
    fields = {'ppr': np.ones((3,3)), 'pfd': np.ones((3, 3)), 'pfv': np.ones((3, 3)), 'FD': np.ones((3, 3))}
    resTypes = ['ppr', 'pfv', 'pfd', 'particles']

    # call function to be tested
    fieldsList, particlesList = com1DFA.appendFieldsParticles(fieldsListIn, particlesListIn, particles, fields, resTypes)
    print('fieldsList', fieldsList[1])

    assert np.array_equal(fieldsList[1]['ppr'], np.ones((3, 3)))
    assert np.array_equal(fieldsList[1]['pfv'], np.ones((3, 3)))
    assert np.array_equal(fieldsList[1]['pfd'], np.ones((3, 3)))
    assert resTypes[0:3] == list(fieldsList[1].keys())
    assert len(fieldsList) == 2
    assert np.array_equal(particlesList[1]['x'], particles['x'])
    assert np.array_equal(particlesList[1]['y'], particles['y'])
    assert np.array_equal(particlesList[1]['m'], particles['m'])
    assert ['x', 'y', 'm'] == list(particlesList[1].keys())
    assert fieldsList[1].get('FD') is None


def test_releaseSecRelArea():
    """ test if secondary release area is triggered """

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'rho': '200.', 'gravAcc': '9.81', 'massPerParticleDeterminationMethod': 'MPPDH',
                      'interpOption': '2', 'sphKernelRadius': '1', 'deltaTh': '0.25', 'seed': '12345',
                      'initPartDistType': 'uniform', 'thresholdPointInPoly': '0.001', 'avalancheDir': 'data/avaTest'}
    demHeader = {}
    demHeader['cellsize'] = 1
    demHeader['ncols'] = 12
    demHeader['nrows'] = 12
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0
    demRaster = np.ones((demHeader['nrows'], demHeader['ncols']))
    areaRaster = np.ones((demHeader['nrows'], demHeader['ncols']))
    dem = {'header': demHeader, 'rasterData': demRaster, 'areaRaster': areaRaster}
    dem['originOri'] = {'xllcenter': 1.0, 'yllcenter': 1.0}
    secRelRaster2 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    secRelRaster2[6:8, 7] = 1.0
    secRelRaster3 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    secRelRaster3[9, 9] = 0.5
    secRelRaster1 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    secRelRaster1[1, 1] = 0.5
    secondaryReleaseInfo = {'x': np.asarray([1.5, 2.5, 2.5, 1.5, 1.5, 7.4, 8.5, 8.5, 7.4, 7.4, 9.5, 10.5, 10.5, 9.5, 9.5]),
                            'y': np.asarray([1.5, 1.5, 2.5, 2.5, 1.5, 7.4, 7.4, 8.5, 8.5, 7.4, 9.5, 9.5, 10.5, 10.5, 9.5]),
                            'Start': np.asarray([0, 5, 10]), 'Length': np.asarray([5, 5, 5]),
                            'Name': ['secRel1', 'secRel2', 'secRel3'], 'd0': [0.5, 1.0, 0.5],
                            'rasterData': [secRelRaster1, secRelRaster2, secRelRaster3]}
    secondaryReleaseInfo['header'] = demHeader
    secondaryReleaseInfo['header']['xllcenter'] = dem['originOri']['xllcenter']
    secondaryReleaseInfo['header']['yllcenter'] = dem['originOri']['yllcenter']
    secondaryReleaseInfo2 = copy.deepcopy(secondaryReleaseInfo)
    particlesIn = {'secondaryReleaseInfo': secondaryReleaseInfo}
    particlesIn['x'] = np.asarray([6., 7.])
    particlesIn['y'] = np.asarray([6., 7.])
    particlesIn['m'] = np.asarray([1250., 1250.])
    particlesIn['t'] = 1.0
    particlesIn['Npart'] = 2.
    fieldsFD = np.zeros((demHeader['nrows'], demHeader['ncols']))
    fieldsFD[7:9, 7:9] = 1.0
    fields = {'FD': fieldsFD}

    # call function to be tested
    particles = com1DFA.releaseSecRelArea(cfg['GENERAL'], particlesIn, fields, dem)
    print('particles IN pytest 1', particles)

    # call function to be tested test 2
    particlesIn2 = {'secondaryReleaseInfo': secondaryReleaseInfo2}
    particlesIn2['x'] = np.asarray([6., 7., 9.1])
    particlesIn2['y'] = np.asarray([6., 7., 9.1])
    particlesIn2['m'] = np.asarray([1250., 1250., 1250.])
    particlesIn2['t'] = 1.0
    particlesIn2['Npart'] = 3
    fieldsFD2 = np.zeros((demHeader['nrows'], demHeader['ncols']))
    fieldsFD2[7:9, 7:9] = 1.0
    fieldsFD2[9, 9] = 0.4
    fields2 = {'FD': fieldsFD2}

    particles2 = com1DFA.releaseSecRelArea(cfg['GENERAL'], particlesIn2, fields2, dem)

    print('particles IN pytest socond', particles2)

    assert particles['Npart'] == 6
    assert np.array_equal(particles['x'], np.asarray([6., 7., 6.75, 7.25, 6.75, 7.25]))
    assert np.array_equal(particles['y'], np.asarray([6., 7., 6.75, 6.75, 7.25, 7.25]))
    assert np.array_equal(particles['m'], np.asarray([1250., 1250., 50., 50., 50., 50.]))
    assert particles['mTot'] == 2700.0
    assert particles2['Npart'] == 11
    assert np.array_equal(particles2['x'], np.asarray([6., 7., 9.1, 6.75, 7.25, 6.75, 7.25, 8.75, 9.25, 8.75, 9.25]))
    assert np.array_equal(particles2['y'], np.asarray([6., 7., 9.1, 6.75, 6.75, 7.25, 7.25, 8.75, 8.75, 9.25, 9.25]))
    assert np.array_equal(particles2['m'], np.asarray([1250., 1250., 1250., 50., 50., 50., 50., 25., 25., 25., 25.]))
    assert particles2['mTot'] == 4050.0


def test_initializeParticles():
    """ test initialising particles """

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'rho': '200.', 'gravAcc': '9.81', 'massPerParticleDeterminationMethod': 'MPPDH',
                      'interpOption': '2', 'sphKernelRadius': '1', 'deltaTh': '0.25', 'seed': '12345',
                      'initPartDistType': 'uniform', 'thresholdPointInPoly': '0.001', 'avalancheDir': 'data/avaTest'}
    demHeader = {}
    demHeader['cellsize'] = 1
    demHeader['ncols'] = 12
    demHeader['nrows'] = 12
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0
    demRaster = np.ones((demHeader['nrows'], demHeader['ncols']))
    areaRaster = np.ones((demHeader['nrows'], demHeader['ncols']))
    dem = {'header': demHeader, 'rasterData': demRaster, 'areaRaster': areaRaster}
    dem['originOri'] = {'xllcenter': 1.0, 'yllcenter': 1.0}

    relRaster = np.zeros((12, 12))
    relRaster[6:8, 6:8] = 1.0
    releaseLine = {'x': np.asarray([6.9, 8.5, 8.5, 6.9, 6.9]), 'y': np.asarray([6.9, 6.9, 8.5, 8.5, 6.9]),
                    'Start': np.asarray([0]), 'Length': np.asarray([5]), 'Name': [''], 'd0': [1.0],
                    'rasterData': relRaster}

    releaseLine['header'] = demHeader
    releaseLine['header']['xllcenter'] = dem['originOri']['xllcenter']
    releaseLine['header']['yllcenter'] = dem['originOri']['yllcenter']

    # call function to be tested
    particles = com1DFA.initializeParticles(cfg['GENERAL'], releaseLine, dem)

    cfg['GENERAL']['massPerParticleDeterminationMethod'] = 'MPPDIR'
    cfg['GENERAL'].update({'massPerPart': '60.'})
    particles2 = com1DFA.initializeParticles(cfg['GENERAL'], releaseLine, dem)

    print('particles', particles2)

    assert particles['Npart'] == 9
    assert np.array_equal(particles['x'], np.asarray([6.25, 6.75, 7.25, 6.25, 6.25, 6.75, 7.25, 6.75, 7.25]))
    assert np.array_equal(particles['y'], np.asarray([6.25, 6.25, 6.25, 6.75, 7.25, 6.75, 6.75, 7.25, 7.25]))
    assert np.array_equal(particles['m'], np.asarray([50., 50., 50., 50., 50., 50., 50., 50., 50.]))
    assert particles['mTot'] == 450.0
    assert np.sum(particles['ux']) == 0.0

    assert particles2['Npart'] == 9
    assert np.array_equal(particles2['x'], np.asarray([6.25, 6.75, 7.25, 6.25, 6.25, 6.75, 7.25, 6.75, 7.25]))
    assert np.array_equal(particles2['y'], np.asarray([6.25, 6.25, 6.25, 6.75, 7.25, 6.75, 6.75, 7.25, 7.25]))
    assert np.array_equal(particles2['m'], np.asarray([50., 50., 50., 50., 50., 50., 50., 50., 50.]))
    assert particles2['mTot'] == 450.0
    assert np.sum(particles2['ux']) == 0.0


def test_writeMBFile(tmp_path):
    """ test writing of mass balance info to file """

    # setup required input
    infoDict = {'timeStep': np.asarray([0, 1, 2, 3, 4])}
    infoDict['massEntrained'] = np.asarray([0, 0, 10, 20, 30])
    infoDict['massTotal'] = np.asarray([60., 60., 70., 90., 120.])
    avaName = 'data/avaTest'
    avaDir = pathlib.Path(tmp_path, avaName)
    logName = 'simTestName'

    # call function to be tested
    com1DFA.writeMBFile(infoDict, avaDir, logName)

    mbFilePath = avaDir / 'Outputs' / 'com1DFA' / 'mass_simTestName.txt'
    mbInfo = np.loadtxt(mbFilePath, delimiter=',', skiprows=1)

    print('mbInfo', mbInfo)

    assert np.array_equal(mbInfo[:, 0], infoDict['timeStep'])
    assert np.array_equal(mbInfo[:, 2], infoDict['massEntrained'])
    assert np.array_equal(mbInfo[:, 1], infoDict['massTotal'])
    assert mbInfo.shape[0] == 5
    assert mbInfo.shape[1] == 3


def test_savePartToPickle(tmp_path):
    """ test saving particles info to pickle """

    # setup required input
    particles1 = {'x': np.asarray([1., 2., 3.]), 'y': np.asarray([1., 4., 5.]),
                  'm': np.asarray([10., 11., 11.]), 't': 0.}
    particles2 = {'x': np.asarray([10., 20., 30.]), 'y': np.asarray([10., 40., 50.]),
                  'm': np.asarray([100., 110., 110.]), 't': 2.}
    dictList = [particles1, particles2]
    outDir = pathlib.Path(tmp_path, 'particles')
    outDir.mkdir()
    logName = 'simNameTest'

    # call function to be tested
    com1DFA.savePartToPickle(dictList, outDir, logName)

    # read pickle
    picklePath1 = outDir / 'particles_simNameTest_0000.0000.p'
    picklePath2 = outDir / 'particles_simNameTest_0002.0000.p'
    particlesRead1 = pickle.load(open(picklePath1, "rb"))
    particlesRead2 = pickle.load(open(picklePath2, "rb"))

    print('particklesRead1', particlesRead1)
    print('particklesRead2', particlesRead2)

    assert np.array_equal(particlesRead1['x'], particles1['x'])
    assert np.array_equal(particlesRead1['y'], particles1['y'])
    assert np.array_equal(particlesRead1['m'], particles1['m'])
    assert particlesRead1['t'] == 0.
    assert  np.array_equal(particlesRead2['x'], particles2['x'])
    assert np.array_equal(particlesRead2['y'], particles2['y'])
    assert np.array_equal(particlesRead2['m'], particles2['m'])
    assert particlesRead2['t'] == 2.

    # call function to be tested
    logName = 'simNameTest3'
    com1DFA.savePartToPickle(particles1, outDir, logName)

    # read pickle
    picklePath3 = outDir / 'particles_simNameTest3_0000.0000.p'
    particlesRead3 = pickle.load(open(picklePath3, "rb"))

    print('particklesRead3', particlesRead3)
    print('particklesRead2', particlesRead2)

    assert np.array_equal(particlesRead3['x'], particles1['x'])
    assert np.array_equal(particlesRead3['y'], particles1['y'])
    assert np.array_equal(particlesRead3['m'], particles1['m'])
    assert particlesRead3['t'] == 0.


def test_readPartFromPickle(tmp_path):
    """ test reading particle properties from pickle """

    # setup required inputs
    inDir = pathlib.Path(tmp_path, 'avaTest')
    inDir.mkdir()
    particlesTestDict = {'x': np.asarray([1., 2., 3.]), 'y': np.asarray([1., 4., 5.]),
                          'm': np.asarray([10., 11., 11.]), 't': 0.}
    pickle.dump(particlesTestDict, open(inDir / 'test.p', "wb"))
    testDir = inDir / 'Outputs' / 'com1DFA' / 'particles'
    testDir.mkdir(parents=True)
    pickle.dump(particlesTestDict, open(testDir / 'test.p', "wb"))

    # call function to be tested
    Particles, TimeStepInfo = com1DFA.readPartFromPickle(inDir, flagAvaDir=False)
    # call function to be tested
    Particles2, TimeStepInfo2 = com1DFA.readPartFromPickle(inDir, flagAvaDir=True)

    print('Particles', Particles)
    print('TimeStepInfo', TimeStepInfo)

    assert np.array_equal(Particles[0]['x'], particlesTestDict['x'])
    assert TimeStepInfo == [0.]
    assert np.array_equal(Particles2[0]['x'], particlesTestDict['x'])
    assert TimeStepInfo2 == [0.]


def test_savePartToCsv(tmp_path):
    """ test saving particle infos to csv file """

    # setup required input
    particleProperties = 'm|x|y|velocityMagnitude|test'
    particles1 = {'x': np.asarray([1., 2., 3.]), 'y': np.asarray([1., 4., 5.]),
                  'z': np.asarray([1., 4., 5.]),
                  'm': np.asarray([10., 11., 11.]), 't': 0., 'simName': 'simNameTest',
                  'xllcenter': 11., 'yllcenter': 12., 'ux': np.asarray([0., 0., 0.]),
                  'uy': np.asarray([0., 0., 0.]), 'uz': np.asarray([0., 0., 0.])}
    particles2 = {'x': np.asarray([10., 20., 30.]), 'y': np.asarray([10., 40., 50.]),
                  'z': np.asarray([1., 4., 5.]),
                  'm': np.asarray([100., 110., 110.]), 't': 2., 'simName': 'simNameTest',
                  'xllcenter': 4., 'yllcenter': 2., 'ux': np.asarray([4., 4., 4.]),
                  'uy': np.asarray([4., 4., 4.]), 'uz': np.asarray([4., 4., 4.])}
    dictList = [particles1, particles2]
    outDir = pathlib.Path(tmp_path, 'testDir')
    outDir.mkdir()

    # call function to be tested
    com1DFA.savePartToCsv(particleProperties, dictList, outDir)

    # read csv file
    partCsv1 = outDir / 'particlesCSV' / 'particlessimNameTest.csv.0'
    DF1 = pd.read_csv(partCsv1)
    partCsv2 = outDir / 'particlesCSV' / 'particlessimNameTest.csv.1'
    DF2 = pd.read_csv(partCsv2)
    velMag = np.sqrt(4**2 + 4**2 + 4**2)

    print('csv df1', DF1.to_string())
    print('csv df2', DF2.to_string())

    assert np.array_equal(DF1['X'], (np.asarray([12., 13., 14.])))
    assert np.array_equal(DF1['Y'], (np.asarray([13., 16., 17.])))
    assert np.array_equal(DF1['m'], np.asarray([10., 11., 11.]))
    assert np.array_equal(DF1['velocityMagnitude'], np.asarray([0., 0., 0.]))
    assert DF1['time'][0] == 0.0
    assert np.array_equal(DF2['X'], (np.asarray([14., 24., 34.])))
    assert np.array_equal(DF2['Y'], (np.asarray([12., 42., 52.])))
    assert np.array_equal(DF2['m'], np.asarray([100., 110., 110.]))
    assert np.array_equal(DF2['velocityMagnitude'], np.asarray([velMag, velMag, velMag]))
    assert DF2['time'][0] == 2.0


def test_exportFields(tmp_path):
    """ test exporting fields to ascii files """

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'resType': 'ppr|pfd|FD'}
    cfg['REPORT'] = {'plotFields': 'ppr|pfd|pfv'}
    Tsave = [0, 10, 15, 25, 40]
    demHeader = {}
    demHeader['cellsize'] = 1
    demHeader['ncols'] = 10
    demHeader['nrows'] = 10
    demHeader['xllcenter'] = 0
    demHeader['yllcenter'] = 0
    demHeader['noDataValue'] = -9999
    dem = {'header': demHeader}
    outDir = pathlib.Path(tmp_path, 'testDir')
    outDir.mkdir()
    logName = 'simNameTest'
    FD = np.zeros((5, 5))
    pfd = np.zeros((5, 5))
    pfv = np.zeros((5, 5))
    ppr = np.zeros((5, 5))
    fields1 = {'ppr': ppr+1, 'pfd': pfd+1, 'pfv': pfv+1, 'FD': FD+1}
    fields2 = {'ppr': ppr+2, 'pfd': pfd+2, 'pfv': pfv+2, 'FD': FD+2}
    fields3 = {'ppr': ppr+4, 'pfd': pfd+4, 'pfv': pfv+4, 'FD': FD+4}
    fields4 = {'ppr': ppr+5, 'pfd': pfd+5, 'pfv': pfv+5, 'FD': FD+5}
    fields5 = {'ppr': ppr+6, 'pfd': pfd+6, 'pfv': pfv+6, 'FD': FD+6}
    fieldsList = [fields1, fields2, fields3, fields4, fields5]

    # call function to be tested
    com1DFA.exportFields(cfg, Tsave, fieldsList, dem, outDir, logName)

    # read fields
    fieldDir = outDir / 'peakFiles'
    fieldDirTSteps = outDir / 'peakFiles' / 'timeSteps'
    fieldFiles = list(fieldDirTSteps.glob('*.asc'))
    fieldsListTest = []

    for f in fieldFiles:
        fieldsListTest.append(f.name)

    field1 = fieldDir / 'simNameTest_ppr.asc'
    field2 = fieldDirTSteps / 'simNameTest_pfd_t10.00.asc'
    fieldFinal = np.loadtxt(field1, skiprows=6)
    field10 = np.loadtxt(field2, skiprows=6)
    pprFinal = ppr+0.006
    pfdt10 = pfd+2

    print('field1', fieldFinal)
    print('pprFinal', pprFinal)
    print('fields', fieldsListTest)

    assert np.array_equal(fieldFinal, pprFinal)
    assert np.array_equal(field10, pfdt10)
    assert len(fieldsListTest) == 16

    # call function to be tested
    outDir2 = pathlib.Path(tmp_path, 'testDir2')
    outDir2.mkdir()
    cfg['GENERAL']['resType'] = ''
    cfg['REPORT'] = {'plotFields': 'ppr|pfd|pfv|particles'}
    com1DFA.exportFields(cfg, Tsave, fieldsList, dem, outDir2, logName)

    # read fields
    fieldDir = outDir2 / 'peakFiles'
    fieldDirTSteps = outDir2 / 'peakFiles' / 'timeSteps'
    fieldFiles = list(fieldDirTSteps.glob('*.asc'))
    fieldsListTest2 = []

    for f in fieldFiles:
        fieldsListTest2.append(f.name)

    assert len(fieldsListTest2) == 15


def test_initializeFields():
    """ test initializing fieldgetSimTypeLists """

    # setup required inputs
    demHeader = {'nrows': 11, 'ncols': 12, 'cellsize': 1}
    areaRaster = np.ones((11, 12))
    dem = {'header': demHeader, 'headerNeighbourGrid': demHeader, 'areaRaster': areaRaster}
    particles = {'x': np.asarray([1., 2., 3.]), 'y': np.asarray([1., 2., 3.]), 'Npart': 3,
                 'ux': np.asarray([0., 0., 0.]), 'uy': np.asarray([0., 0., 0.]),
                 'uz': np.asarray([0., 0., 0.]), 'm': np.asarray([10., 10., 10.])}
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'rho': '200.', 'interpOption': '2'}

    # call function to be tested
    particles, fields = com1DFA.initializeFields(cfg['GENERAL'], dem, particles)

    print('particles', particles)
    print('fields', fields)

    assert len(fields) == 9
    assert np.sum(fields['pfv']) == 0.0
    assert np.sum(fields['ppr']) == 0.0
    assert np.sum(fields['FV']) == 0.0
    assert np.sum(fields['P']) == 0.0
    assert np.sum(fields['Vx']) == 0.0
    assert np.sum(fields['Vy']) == 0.0
    assert np.sum(fields['Vz']) == 0.0
    assert np.sum(fields['pfd']) != 0.0
    assert np.sum(fields['FD']) != 0.0


def test_prepareVarSimDict():
    """ test prepare variation sim dictionary """

    # setup required input
    standardCfg = configparser.ConfigParser()
    standardCfg.optionxform = str
    standardCfg['GENERAL'] = {'simTypeList': 'entres|null', 'modelType': 'dfa', 'simTypeActual': 'entres'}
    relPath = pathlib.Path('test', 'relTest.shp')
    inputSimFiles = {'relFiles': [relPath], 'entResInfo': {'flagEnt': 'Yes', 'flagRes': 'Yes'}}
    variationDict = {'rho': np.asarray([200., 150.])}

    # call function to be tested
    simDict = com1DFA.prepareVarSimDict(standardCfg, inputSimFiles, variationDict)
    testCfg = configparser.ConfigParser()
    testCfg.optionxform = str
    testCfg['GENERAL'] = {'simTypeList': 'entres', 'modelType': 'dfa', 'rho': '200.0',
                          'simTypeActual': 'entres', 'releaseScenario': 'relTest'}
    simHash = cfgUtils.cfgHash(testCfg)
    simName1 = 'relTest_entres_dfa_' + simHash
    testDict = {'relTest_entres_dfa_f466369a03': {'simHash': 'f466369a03', 'releaseScenario': 'relTest',
                'simType': 'entres', 'relFile': relPath, 'cfgSim': testCfg}}

    print('simDict', simDict)
    print('simName1', simName1)

    for key in testDict[simName1]:
        assert simDict[simName1][key] == testDict[simName1][key]

    for section in testCfg.sections():
        for key in testCfg[section]:
            print('section', section, 'key', key)
            assert simDict[simName1]['cfgSim'][section][key] == testCfg[section][key]

    # call function to be tested
    relPath = pathlib.Path('test', 'relTest_extended.shp')
    inputSimFiles = {'relFiles': [relPath], 'entResInfo': {'flagEnt': 'Yes', 'flagRes': 'Yes'}}
    variationDict = {'rho': np.asarray([200., 150.]), 'simTypeList': ['entres', 'ent']}
    simDict2 = com1DFA.prepareVarSimDict(standardCfg, inputSimFiles, variationDict)

    testCfg2 = configparser.ConfigParser()
    testCfg2.optionxform = str
    testCfg2['GENERAL'] = {'simTypeList': 'entres', 'modelType': 'dfa', 'rho': '200.0',
                          'simTypeActual': 'entres', 'releaseScenario': 'relTest_extended'}
    simHash2 = cfgUtils.cfgHash(testCfg2)
    simName2 = 'relTest_extended_AF_entres_dfa_' + simHash2
    testDict2 = {simName2: {'simHash': simHash2, 'releaseScenario': 'relTest_extended',
                'simType': 'entres', 'relFile': relPath, 'cfgSim': testCfg2}}

    print('simDict', simDict2)
    print('simName1', simName2)

    for key in testDict2[simName2]:
        assert simDict2[simName2][key] == testDict2[simName2][key]

    for section in testCfg2.sections():
        for key in testCfg2[section]:
            print('section', section, 'key', key)
            assert simDict2[simName2]['cfgSim'][section][key] == testCfg2[section][key]


def test_initializeSimulation():
    """ test initializing a simulation """

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'methodMeshNormal': '1', 'thresholdPointInPoly': '0.001',
                      'sphKernelRadius': '1.', 'meshCellSizeThreshold': '0.0001',
                      'meshCellSize': '1.', 'simTypeActual': 'ent', 'rhoEnt': '100.', 'hEnt': '0.3',
                      'rho': '200.', 'gravAcc': '9.81', 'massPerParticleDeterminationMethod': 'MPPDH',
                      'interpOption': '2', 'sphKernelRadius': '1', 'deltaTh': '0.25', 'seed': '12345',
                      'initPartDistType': 'uniform', 'thresholdPointInPoly': '0.001', 'avalancheDir': 'data/avaTest',
                      'dRes': '0.3', 'cw': '0.5', 'sres': '5', 'initialiseParticlesFromFile': 'False'}
    # setup dem input
    demHeader = {}
    demHeader['xllcenter'] = 1.0
    demHeader['yllcenter'] = 2.0
    demHeader['cellsize'] = 1.0
    demHeader['noDataValue'] = -9999
    demHeader['nrows'] = 12
    demHeader['ncols'] = 12
    demData = np.ones((12, 12))
    demOri = {'header': demHeader, 'rasterData': demData}

    # setup release line, entrainment line
    releaseLine = {'x': np.asarray([6.9, 8.5, 8.5, 6.9, 6.9]), 'y': np.asarray([7.9, 7.9, 9.5, 9.5, 7.9]),
                    'Start': np.asarray([0]), 'Length': np.asarray([5]), 'Name': [''], 'd0': [1.0]}
    entLine = {'fileName': 'test/entTest.shp', 'Name': ['testEnt'], 'Start': np.asarray([0.]), 'Length': np.asarray([5]),
               'x': np.asarray([4, 5., 5.0, 4., 4.]), 'y': np.asarray([4., 4., 5.0, 5., 4.0])}
    inputSimLines = {'releaseLine': releaseLine, 'entResInfo': {'flagSecondaryRelease': 'No'}, 'entLine': entLine,
                     'resLine': ''}
    # set release thickness read from file or not
    relThField = ''
    logName = 'simLog'

    # call function to be tested
    particles, fields, dem, reportAreaInfo = com1DFA.initializeSimulation(cfg, demOri, inputSimLines, logName, relThField)

    print('particles', particles)
    print('fields', fields)
    print('dem', dem)
    print('reportAreaInfo', reportAreaInfo)

    assert np.array_equal(particles['y'], np.asarray([6.25, 6.25, 6.25, 6.75, 7.25, 6.75, 6.75, 7.25, 7.25]))
    assert np.sum(fields['pfv']) == 0.0
    assert np.sum(fields['pfd']) != 0.0
    assert dem['header']['xllcenter'] == 0.0
    assert dem['header']['yllcenter'] == 0.0
    assert dem['originOri']['xllcenter'] == 1.0
    assert dem['originOri']['yllcenter'] == 2.0
    assert particles['Npart'] == 9
    assert np.array_equal(particles['x'], np.asarray([6.25, 6.75, 7.25, 6.25, 6.25, 6.75, 7.25, 6.75, 7.25]))
    assert np.array_equal(particles['m'], np.asarray([50., 50., 50., 50., 50., 50., 50., 50., 50.]))
    assert particles['mTot'] == 450.
    assert np.sum(particles['ux']) == 0.0
    assert reportAreaInfo['Release area info']['Projected Area [m2]'] == '9.00'
    assert reportAreaInfo['entrainment'] == 'Yes'
    assert reportAreaInfo['resistance'] == 'No'

    # call function to be tested
    inputSimLines['entResInfo']['flagSecondaryRelease'] = 'Yes'
    inputSimLines['secondaryReleaseLine'] = {'x': np.asarray([1.5, 2.5, 2.5, 1.5, 1.5]),
                                             'y': np.asarray([2.5, 2.5, 3.5, 3.5, 2.5]),
                                             'Start': np.asarray([0]), 'Length': np.asarray([5]),
                                             'Name': ['secRel1'], 'd0': [0.5]}
    relThField = np.zeros((12, 12)) + 0.5
    particles2, fields2, dem2, reportAreaInfo2 = com1DFA.initializeSimulation(cfg, demOri, inputSimLines, logName, relThField)

    print('secRel', particles2['secondaryReleaseInfo'])
    print('particles', particles2)
    print('fields', fields2['pfd'])

    assert np.sum(fields2['pfv']) == 0.0
    assert np.sum(fields2['pfd']) != np.sum(fields['pfd'])
    assert dem2['header']['xllcenter'] == 0.0
    assert dem2['header']['yllcenter'] == 0.0
    assert dem2['originOri']['xllcenter'] == 1.0
    assert dem2['originOri']['yllcenter'] == 2.0
    assert particles2['Npart'] == 9
    assert particles2['mTot'] == 225.
    assert np.sum(particles['ux']) == 0.0
    assert reportAreaInfo['Release area info']['Projected Area [m2]'] == '9.00'
    assert reportAreaInfo['entrainment'] == 'Yes'
    assert reportAreaInfo['resistance'] == 'No'
    assert np.sum(particles2['secondaryReleaseInfo']['rasterData']) == 4.5
