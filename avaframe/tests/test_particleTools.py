"""Tests for module DFAtools"""
import numpy as np
import pickle
import configparser
import pathlib
import pandas as pd

# Local imports
import avaframe.com1DFA.particleTools as particleTools


def test_placeParticles():
    """ test placing of particles """

    # setup required inputs
    indx = 0
    indy = 1
    csz = 5
    aCell = csz * csz
    hCell = 10/25
    massPerPart = 2.
    thresholdMassSplit = 1.5
    initPartDistType = 'uniform'
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'rho': '1', 'thresholdMassSplit': '0.5', 'initPartDistType': 'uniform',
                      'massPerParticleDeterminationMethod': 'MPPDH'}
    nPPK = 10
    rng = np.random.default_rng(12345)
    # call funciton to be tested - uniform
    xpart, ypart, mPart, nPart, aPart = particleTools.placeParticles(hCell, aCell, indx, indy, csz, massPerPart, nPPK,
                                                                     rng, cfg['GENERAL'])
    xpartTest = np.asarray([-1.66666666, 0.0, 1.66666666, -1.66666666, 0., 1.66666666, -1.66666666,
                            0.0, 1.66666666])
    ypartTest = np.asarray([3.33333333, 3.33333333, 3.33333333, 5.0, 5., 5., 6.66666666, 6.66666666,
                            6.66666666])

    assert nPart == 9.0
    assert np.isclose(mPart, 1.111111)
    assert np.allclose(xpart, xpartTest)
    assert np.allclose(ypart, ypartTest)

    # call funciton to be tested - uniform
    hCell = 8/25
    xpart, ypart, mPart, nPart, aPart = particleTools.placeParticles(hCell, aCell, indx, indy, csz, massPerPart, nPPK,
                                                                     rng, cfg['GENERAL'])
    xpartTest = np.asarray([-1.25, 1.25, -1.25, 1.25])
    ypartTest = np.asarray([3.75, 3.75, 6.25, 6.25])

    assert nPart == 4.0
    assert mPart == 2.
    assert np.allclose(xpart, xpartTest)
    assert np.allclose(ypart, ypartTest)

    # call funciton to be tested - random
    hCell = 11.5/25
    cfg['GENERAL']['initPartDistType'] = 'random'
    xpart, ypart, mPart, nPart, aPart = particleTools.placeParticles(hCell, aCell, indx, indy, csz, massPerPart, nPPK,
                                                                     rng, cfg['GENERAL'])
    xpartTest = np.asarray(
        [-0.9162083, 1.48682729, 0.88127335, -0.54445225, -0.83593036, 0.49154377, -1.56632907])
    ypartTest = np.asarray(
        [5.86378022, 7.20901433, 3.74122857, 7.24440576, 5.83618727, 2.97948968, 4.70919833])

    print('xpart', xpart)
    print('ypart', ypart)
    assert nPart == 7.0
    assert np.isclose(mPart, 1.6428571428571428)
    assert np.allclose(xpart, xpartTest)
    assert np.allclose(ypart, ypartTest)

    # call funciton to be tested - random
    csz = 4
    aCell = csz * csz
    hCell = 8/16
    cfg['GENERAL']['initPartDistType'] = 'semiRandom'
    xpart, ypart, mPart, nPart, aPart = particleTools.placeParticles(hCell, aCell, indx, indy, csz, massPerPart, nPPK,
                                                                     rng, cfg['GENERAL'])

    print('xpart', xpart)
    print('ypart', ypart)

    assert nPart == 4.0
    assert -2.0 < xpart[0] < 0.0
    assert 2.0 < ypart[0] < 4.0
    assert 0.0 < xpart[1] < 2.0
    assert 2.0 < ypart[1] < 4.0
    assert -2.0 < xpart[2] < 0.0
    assert 4.0 < ypart[2] < 6.0
    assert 0.0 < xpart[3] < 2.0
    assert 4.0 < ypart[3] < 6.0


def test_removePart(capfd):
    particles = {}
    particles['nPart'] = 10
    particles['ID'] = np.arange(particles['nPart'])
    particles['parentID'] = np.arange(particles['nPart'])
    particles['nID'] = 10
    particles['m'] = np.linspace(0, 9, 10)
    particles['x'] = np.linspace(0, 9, 10)
    particles['ux'] = np.linspace(0, 9, 10)
    particles['mTot'] = np.sum(particles['m'])
    mask = np.array([True, True, False, True, True,
                     True, False, False, True, True])
    nRemove = 3
    particles = particleTools.removePart(particles, mask, nRemove)

    res = np.array([0, 1, 3, 4, 5, 8, 9])
    atol = 1e-10
    assert particles['nPart'] == 7
    assert np.allclose(particles['m'], res, atol=atol)
    assert np.allclose(particles['x'], res, atol=atol)
    assert np.allclose(particles['ux'], res, atol=atol)
    assert particles['mTot'] == np.sum(res)
    assert particles['nID'] == 10
    assert np.allclose(particles['parentID'], res, atol=atol)
    assert np.allclose(particles['ID'], res, atol=atol)


def test_splitPartMass(capfd):
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'rho': '1', 'thresholdMassSplit': '1.5', 'distSplitPart': '0.5'}
    particles = {}
    particles['nPart'] = 10
    particles['massPerPart'] = 1
    particles['ID'] = np.arange(particles['nPart'])
    particles['parentID'] = np.arange(particles['nPart'])
    particles['m'] = np.array([1, 2, 1.4, 3.6, 1, 1.6, 5, 1, 1, 1])
    particles['h'] = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    particles['x'] = np.linspace(0, 9, 10)
    particles['y'] = np.zeros((10))
    particles['z'] = np.zeros((10))
    particles['ux'] = np.linspace(0, 9, 10)
    particles['uy'] = np.zeros((10))
    particles['uz'] = np.zeros((10))
    particles['mTot'] = np.sum(particles['m'])
    particles['nID'] = 10
    particles = particleTools.splitPartMass(particles, cfg['GENERAL'])
    print(particles)
    massNew = np.array([1, 1, 1.4, 1.8, 1, 0.8, 2.5, 1, 1, 1,
                        1, 1.8, 0.8, 2.5])
    xNew = np.array([0., 0.60105772, 2., 2.46476277, 4., 4.64317518, 5.36921687, 7., 8., 9.,
                     1.39894228, 3.53523723, 5.35682482, 6.63078313])
    res = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 3, 5, 6])
    print(particles['m'])
    print(particles['x'])
    print(massNew)
    print(particles['nID'])
    print(particles['parentID'])
    print(particles['ID'])
    atol = 1e-10
    assert particles['nPart'] == 14
    assert np.allclose(particles['m'], massNew, atol=atol)
    assert np.allclose(particles['x'], xNew, atol=atol)
    assert np.allclose(particles['ux'], res, atol=atol)
    assert particles['mTot'] == np.sum(massNew)
    assert particles['nID'] == 14
    assert np.allclose(particles['parentID'], res, atol=atol)
    assert np.allclose(particles['ID'], np.arange(14), atol=atol)


def test_mergeParticleDict(capfd):

    particles1 = {}
    particles1['nPart'] = 5
    particles1['m'] = np.linspace(0, 4, 5)
    particles1['x'] = np.linspace(0, 4, 5)
    particles1['ux'] = np.linspace(0, 4, 5)
    particles1['mTot'] = np.sum(particles1['m'])
    particles1['ID'] = np.arange(particles1['nPart'])
    particles1['parentID'] = np.arange(particles1['nPart'])
    particles1['nID'] = particles1['nPart']

    particles2 = {}
    particles2['nPart'] = 4
    particles2['m'] = np.linspace(5, 8, 4)
    particles2['x'] = np.linspace(5, 8, 4)
    particles2['ux'] = np.linspace(5, 8, 4)
    particles2['mTot'] = np.sum(particles2['m'])
    particles2['ID'] = np.arange(particles2['nPart'])
    particles2['parentID'] = np.arange(particles2['nPart'])
    particles2['nID'] = particles2['nPart']

    particles = particleTools.mergeParticleDict(particles1, particles2)
    res = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
    atol = 1e-10
    assert particles['nPart'] == 9
    assert np.allclose(particles['m'], res, atol=atol)
    assert np.allclose(particles['x'], res, atol=atol)
    assert np.allclose(particles['ux'], res, atol=atol)
    assert particles['mTot'] == np.sum(res)
    assert particles['nID'] == 9
    assert np.allclose(particles['parentID'], res, atol=atol)
    assert np.allclose(particles['ID'], np.arange(9), atol=atol)


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
    Particles, TimeStepInfo = particleTools.readPartFromPickle(
        inDir, flagAvaDir=False)
    # call function to be tested
    Particles2, TimeStepInfo2 = particleTools.readPartFromPickle(
        inDir, flagAvaDir=True)

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
    particleTools.savePartToCsv(particleProperties, dictList, outDir)

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
    assert np.array_equal(DF2['velocityMagnitude'],
                          np.asarray([velMag, velMag, velMag]))
    assert DF2['time'][0] == 2.0
