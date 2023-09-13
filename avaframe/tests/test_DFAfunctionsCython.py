""
import numpy as np
import math
import configparser
import pytest
import matplotlib.tri as tri

# Local imports
import avaframe.com1DFA.DFAfunctionsCython as DFAfunC
import avaframe.com1DFA.DFAtools as DFAtls


def test_getNeighborsC(capfd):
    """ Test the grid search/particle location method"""
    header = {}
    header['ncols'] = 5
    header['nrows'] = 6
    header['cellsize'] = 1
    dem = {}
    dem['header'] = header
    dem['headerNeighbourGrid'] = header
    particles = {}
    particles['nPart'] = 18
    particles['x'] = np.array(
        [1.6, 0.4, 1, 2, 1, 2, 0, 1, 0, 2, 0, 2, 1, 2, 3, 3, 4, 0])
    particles['y'] = np.array(
        [2.6, 1.4, 0, 1, 3, 3, 2, 1, 0, 0, 3, 2, 2, 1, 1, 4, 5, 5])
    particles['z'] = np.array(
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    particles['m'] = particles['z']
    atol = 1e-10
    indPCell = np.array([0.,  # always an extra zero at the begining
                         1,  # found 1 particle
                         2,  # found 1 particle
                         3,  # found 1 particle
                         3,  # nothing happens
                         3,  # nothing happens
                         4,  # found 1 particle
                         5,  # found 1 particle
                         7,   # found 2 particles
                         8,  # found 1 particle
                         8,  # nothing happens
                         9,  # found 1 particle
                         10,  # found 1 particle
                         11,   # found 1 particles
                         11, 11,  # nothing happens
                         12,  # found 1 particle
                         13,  # found 1 particle
                         15,  # found 2 particle
                         15, 15, 15, 15, 15,  # nothing happens
                         16,  # found 1 particle
                         16,  # nothing happens
                         17,  # found 1 particle
                         17, 17, 17,  # nothing happens
                         18])  # found 1 particle
    pInC = np.array([8,
                     2,
                     9,
                     1,
                     7,
                     13, 3,
                     14,
                     6,
                     12,
                     11,
                     10,
                     4,
                     5, 0,
                     15,
                     17,
                     16])

    particles = DFAfunC.getNeighborsC(particles, dem)
    print(particles['inCellDEM'])
    print(particles['indPartInCell'])
    print(particles['partInCell'])
    assert np.allclose(particles['indPartInCell'], indPCell, atol=atol)
    assert np.allclose(particles['partInCell'], pInC, atol=atol)


def test_computeEntMassAndForce(capfd):
    """ Test the computeEntMassAndForce function"""
    dt = 0.1
    entrMassCell = 0
    areaPart = 1
    uMag = 10
    tau = 1000
    entEroEnergy = 5000
    rhoEnt = 100
    dm, areaEntrPart = DFAfunC.computeEntMassAndForce(
        dt, entrMassCell, areaPart, uMag, tau, entEroEnergy, rhoEnt)
    print(dm, areaEntrPart)
    assert dm == 0
    assert areaEntrPart == 1

    entrMassCell = 200
    entEroEnergy = 0
    dm, areaEntrPart = DFAfunC.computeEntMassAndForce(
        dt, entrMassCell, areaPart, uMag, tau, entEroEnergy, rhoEnt)
    print(dm, areaEntrPart)
    assert dm == 200
    assert areaEntrPart == 2

    entEroEnergy = 5000
    dm, areaEntrPart = DFAfunC.computeEntMassAndForce(
        dt, entrMassCell, areaPart, uMag, tau, entEroEnergy, rhoEnt)
    print(dm, areaEntrPart)
    assert dm == 0.2
    assert areaEntrPart == 1


def test_computeResForce(capfd):
    """ Test the computeResForce function"""
    hRes = 2
    h = 1
    areaPart = 1
    rho = 200
    cResCell = 1
    uMag = 10
    explicitFriction = 0
    cResPart = DFAfunC.computeResForce(
        hRes, h, areaPart, rho, cResCell, uMag, explicitFriction)
    print(cResPart)
    assert cResPart == -2000

    h = 3
    cResPart = DFAfunC.computeResForce(
        hRes, h, areaPart, rho, cResCell, uMag, explicitFriction)
    print(cResPart)
    assert cResPart == -4000

    explicitFriction = 1
    cResPart = DFAfunC.computeResForce(
        hRes, h, areaPart, rho, cResCell, uMag, explicitFriction)
    print(cResPart)
    assert cResPart == -40000


def test_account4FrictionForce(capfd):
    """ Test the account4FrictionForce function"""
    uxNew = 10
    uyNew = 0
    uzNew = 0
    m = 10
    dt = 0.1
    forceFrict = 100
    uMag = 10
    explicitFriction = 0
    uxNew, uyNew, uzNew, dtStop = DFAfunC.account4FrictionForce(
        uxNew, uyNew, uzNew, m, dt, forceFrict, uMag, explicitFriction)
    print(uxNew, uyNew, uzNew, dtStop)
    assert dtStop == 0.1
    assert uxNew == 5

    uxNew = 10
    uyNew = 0
    uzNew = 0
    explicitFriction = 1
    m = 0.5
    uxNew, uyNew, uzNew, dtStop = DFAfunC.account4FrictionForce(
        uxNew, uyNew, uzNew, m, dt, forceFrict, uMag, explicitFriction)
    print(uxNew, uyNew, uzNew, dtStop)
    print(dt*forceFrict/m)
    assert dtStop == 0.05
    assert uxNew == 0

    uxNew = 10
    uyNew = 0
    uzNew = 0
    m = 10
    uxNew, uyNew, uzNew, dtStop = DFAfunC.account4FrictionForce(
        uxNew, uyNew, uzNew, m, dt, forceFrict, uMag, explicitFriction)
    print(uxNew, uyNew, uzNew, dtStop)
    print(dt*forceFrict/m)
    assert dtStop == 0.1
    assert uxNew == 9.0


def test_updatePositionC():
    """ test updating the position of particles """

    # TODO: make test also if velocity in z not zero!! - so to test also when reprojecting onto surface

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'stopCrit': '0.01', 'stopCritIni': '0.1', 'stopCritIniSmall': '1.001', 'stopCritType': 'kinEnergy',
                      'uFlowingThreshold': '0.1', 'gravAcc': '9.81', 'velMagMin': '1.e-6',  'rho': '100.',
                      'interpOption': '2',   'explicitFriction': '0', 'centeredPosition': '1',
                      'reprojMethodPosition': '2', 'reprojectionIterations': '5', 'thresholdProjection': '0.001',
                      'dissDam': '1', 'snowSlide': '1', 'wetSnow': '1'}

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 'trajectoryLengthXY': np.asarray([0., 0., 0.]),
                  'trajectoryLengthXYCor': np.asarray([0., 0., 0.]), 'trajectoryLengthXYZ': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'uAcc': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 0.0,
                  'peakForceSPH': 0.0, 'forceSPHIni': 0.0, 'nPart': 3,
                  'peakMassFlowing': 0.0, 'iterate': True, 'totalEnthalpy': np.asarray([0., 0., 0.]),
                  'velocityMag': np.asarray([1., 1., 1.])}
    particles['potentialEne'] = np.sum(9.81 * particles['z'] * particles['m'])

    demHeader = {}
    demHeader['xllcenter'] = 0.0
    demHeader['yllcenter'] = 0.0
    demHeader['cellsize'] = 5.0
    demHeader['nodata_value'] = -9999
    demHeader['nrows'] = 10
    demHeader['ncols'] = 10
    dem = {'header': demHeader}
    dem['rasterData'] = np.ones((demHeader['nrows'], demHeader['ncols']))
    dem['outOfDEM'] = np.zeros((demHeader['nrows'], demHeader['ncols'])).flatten()
    dem['Nx'] = np.zeros((demHeader['nrows'], demHeader['ncols']))
    dem['Ny'] = np.zeros((demHeader['nrows'], demHeader['ncols']))
    dem['Nz'] = np.ones((demHeader['nrows'], demHeader['ncols']))

    force = {'forceZ': np.asarray([0., 0., 0.]), 'forceFrict': np.asarray([10., 10., 10.]),
             'dM': np.asarray([0., 0., 0.]),
             'forceX': np.asarray([50., 50., 50.]), 'forceY': np.asarray([50., 50., 50.]),
             'forceSPHX': np.asarray([50., 50., 50.]),
             'forceSPHY': np.asarray([50., 50., 50.]), 'forceSPHZ': np.asarray([0., 0., 0.])}
    fields = {'FT': np.zeros((2, 2))}
    # crete a dummy dict (needed so that cython runs)
    wallLineDict = {'dam': 0, 'nIterDam': 1, 'cellsCrossed': np.zeros((dem['header']['ncols']*dem['header']['nrows'])).astype(int)}
    for key in ['x', 'y', 'z', 'xCrown', 'yCrown', 'zCrown', 'xTangent', 'yTangent', 'zTangent']:
        wallLineDict[key] = np.ones((1))*1.0
    for key in ['nPoints', 'height', 'slope', 'restitutionCoefficient']:
        wallLineDict[key] = 0
    dem['damLine'] = wallLineDict
    typeStop = 0

    # kinetic energy new
    kinEneNew = 0.0
    potEneNew = 0.0
    for k in range(3):
        kinEneNew = kinEneNew + particles['m'][k] * np.sqrt(5.5**2 + 5.5**2 + 0**2)**2 * 0.5
        potEneNew = potEneNew + particles['m'][k] * 9.81 + 0.0

    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, fields, typeStop=typeStop)
    uAcc = (np.sqrt(5.5**2 + 5.5**2 +0.) - np.sqrt(1.**2 + 1.**2 + 0.)) / 1.
    velocityMag = np.asarray([np.sqrt((5.5**2) + (5.5**2)), np.sqrt((5.5**2) + (5.5**2)), np.sqrt((5.5**2) + (5.5**2))])

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([3.25, 4.25, 5.25]))
    assert np.array_equal(particles['y'], np.asarray([5.25, 6.25, 7.25]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.allclose(particles['uAcc'], np.asarray([uAcc, uAcc, uAcc]), atol=1.e-4)
    assert np.allclose(particles['velocityMag'], velocityMag, atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == True

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 'trajectoryLengthXY': np.asarray([0., 0., 0.]),
                  'trajectoryLengthXYCor': np.asarray([0., 0., 0.]), 'trajectoryLengthXYZ': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'uAcc': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 100000.0,
                  'peakForceSPH': 0.0, 'forceSPHIni': 0.0, 'nPart': 3,
                  'velocityMag': np.asarray([1., 1., 1.]),
                  'peakMassFlowing': 0.0, 'iterate': True, 'totalEnthalpy': np.asarray([0., 0., 0.])}
    particles['potentialEne'] = np.sum(9.81 * particles['z'] * particles['m'])

    # call function to be tested
    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, fields, typeStop=typeStop)

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([3.25, 4.25, 5.25]))
    assert np.array_equal(particles['y'], np.asarray([5.25, 6.25, 7.25]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.allclose(particles['uAcc'], np.asarray([uAcc, uAcc, uAcc]), atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == False

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 'trajectoryLengthXY': np.asarray([0., 0., 0.]),
                  'trajectoryLengthXYCor': np.asarray([0., 0., 0.]), 'trajectoryLengthXYZ': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'uAcc': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 10000.0,
                  'peakForceSPH': 100000.0, 'forceSPHIni': 1.e5, 'nPart': 3,
                  'velocityMag': np.asarray([1., 1., 1.]),
                  'peakMassFlowing': 0.0, 'iterate': True,  'totalEnthalpy': np.asarray([0., 0., 0.])}
    particles['potentialEne'] = np.sum(9.81 * particles['z'] * particles['m'])
    typeStop = 1

    sphForceNew = 0.0
    kinEneNew = 0.0
    potEneNew = 0.0
    for k in range(3):
        sphForceNew = sphForceNew + particles['m'][k] * np.sqrt(50.**2 +50.**2 + 0.**2)**2.
        kinEneNew = kinEneNew + particles['m'][k] * np.sqrt(11.**2 +11.**2 + 0**2)**2 * 0.5
        potEneNew = potEneNew + particles['m'][k] * 9.81 + 0.0

    # call function to be tested
    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, fields, typeStop=typeStop)
    print('sph', particles['peakForceSPH'], sphForceNew)

    uAcc = (np.sqrt(11.**2 + 11.0**2 +0.) - np.sqrt(1.**2 + 1.**2 + 0.)) / 1.
    velocityMag = np.asarray([np.sqrt((11.**2) + (11.**2)), np.sqrt((11.**2) + (11.**2)), np.sqrt((11.**2) + (11.**2))])

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([6., 7., 8.]))
    assert np.array_equal(particles['y'], np.asarray([8., 9., 10.]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.allclose(particles['velocityMag'], velocityMag, atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert np.allclose(particles['uAcc'], np.asarray([uAcc, uAcc, uAcc]), atol=1.e-4)
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == False

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 'trajectoryLengthXY': np.asarray([0., 0., 0.]),
                  'trajectoryLengthXYCor': np.asarray([0., 0., 0.]), 'trajectoryLengthXYZ': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'uAcc': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 10000.0,
                  'peakForceSPH': 1000.0, 'forceSPHIni': 1.e5, 'nPart': 3,
                  'velocityMag': np.asarray([1., 1., 1.]),
                  'peakMassFlowing': 0.0, 'iterate': True, 'totalEnthalpy': np.asarray([0., 0., 0.])}
    particles['potentialEne'] = np.sum(9.81 * particles['z'] * particles['m'])
    typeStop = 1

    sphForceNew = 0.0
    kinEneNew = 0.0
    potEneNew = 0.0
    for k in range(3):
        sphForceNew = sphForceNew + particles['m'][k] * np.sqrt(50.**2 +50.**2 + 0.**2)**2.
        kinEneNew = kinEneNew + particles['m'][k] * np.sqrt(11.**2 +11.**2 + 0**2)**2 * 0.5
        potEneNew = potEneNew + particles['m'][k] * 9.81 + 0.0

    # call function to be tested
    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, fields, typeStop=typeStop)
    print('sph', particles['peakForceSPH'], sphForceNew)

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([6., 7., 8.]))
    assert np.array_equal(particles['y'], np.asarray([8., 9., 10.]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.allclose(particles['uAcc'], np.asarray([uAcc, uAcc, uAcc]), atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == True


def test_computeTrajectoryAngle():
    # first compute travel angle for each particle
    # get parent Id in order to  get the first z position
    parentID = np.array([0, 1, 2, 0])
    nPart = 4
    # get z0
    zPartArray0 = np.array([10.0, 9.0, 8.0])
    s = np.array([10.0, 10.0, 0.0, 10.0])
    z = np.array([0.0, 0.0, 0.0, 1.0])
    particles = {'nPart': nPart, 'parentID': parentID, 'trajectoryLengthXY': s, 'z': z}
    particles = DFAfunC.computeTrajectoryAngleC(particles, zPartArray0)
    print(particles['trajectoryAngle'])
    gamma = particles['trajectoryAngle']
    assert gamma[2] == 0
    assert gamma[0] == 45
    assert gamma[1] == pytest.approx(41.9872125, rel=1e-6)
    assert gamma[3] == pytest.approx(41.9872125, rel=1e-6)


def test_initializeBondsC():
    nPart = 3
    x = np.array([0., 1., 0.])
    y = np.array([0., 0., 1.])
    z = np.array([0., 0., 0.])
    # original triangulation
    triangles = tri.Triangulation(x, y)
    particles = {'nPart': nPart, 'x': x, 'y': y, 'z': z}
    particles = DFAfunC.initializeBondsC(particles, triangles)
    print(triangles.triangles)
    print(triangles.edges)
    print(particles['bondStart'])
    print(particles['bondDist'])
    print(particles['bondPart'])
    bondStart = particles['bondStart']
    bondDist = particles['bondDist']
    bondPart = particles['bondPart']
    assert np.array_equal(bondStart, np.asarray([0, 2, 4, 6]))
    for k in range(nPart):
        # loop on all bonded particles
        neighbors = list()
        for ib in range(bondStart[k], bondStart[k + 1]):
            l = bondPart[ib]
            neighbors.append(l)

        neighbors.sort()
        if k == 0:
            assert neighbors == [1, 2]
        if k == 1:
            assert neighbors == [0, 2]
        if k == 2:
            assert neighbors == [0, 1]
    bondDist.sort()
    assert np.array_equal(bondDist, np.asarray([1, 1, 1, 1, math.sqrt(2), math.sqrt(2)]))


def test_removeBondsC():
    nPart = 3
    x = np.array([0., 1., 0.])
    y = np.array([0., 0., 1.])
    z = np.array([0., 0., 0.])
    # original triangulation
    triangles = tri.Triangulation(x, y)
    particles = {'nPart': nPart, 'x': x, 'y': y, 'z': z}
    particles = DFAfunC.initializeBondsC(particles, triangles)
    print(particles['bondStart'])
    print(particles['bondDist'])
    print(particles['bondPart'])

    # now remove one particle (here particle 1):
    keepParticle = np.array([1., 0., 1.])
    nRemove = 1
    nBondRemove = DFAfunC.countRemovedBonds(particles, keepParticle, nRemove)
    print(nBondRemove)
    assert nBondRemove == 4
    particles = DFAfunC.removedBonds(particles, keepParticle, nRemove, nBondRemove)
    print(particles['bondStart'])
    print(particles['bondDist'])
    print(particles['bondPart'])
    bondStart = particles['bondStart']
    bondDist = particles['bondDist']
    bondPart = particles['bondPart']
    assert np.array_equal(bondStart, np.asarray([0, 1, 1, 2]))
    for k in range(nPart):
        # loop on all bonded particles
        neighbors = list()
        for ib in range(bondStart[k], bondStart[k + 1]):
            l = bondPart[ib]
            neighbors.append(l)

        neighbors.sort()
        if k == 0:
            assert neighbors == [2]
        if k == 1:
            assert neighbors == []
        if k == 2:
            assert neighbors == [0]
    bondDist.sort()
    assert np.array_equal(bondDist, np.asarray([1, 1]))


def test_computeCohesionForceC():
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'cohesiveSurfaceTension': '50000', 'cohesionMaxStrain': '0.2', 'minDistCohesion' : '1.0e-3'}
    nPart = 3
    x = np.array([0., 1., 0.])
    y = np.array([0., 0., 1.])
    z = np.array([0., 0., 0.])
    ux = np.array([0., 0., 0.])
    uy = np.array([0., 0., 0.])
    uz = np.array([0., 0., 0.])
    m = np.array([1., 1., 1.])
    h = np.array([1., 1., 1.])
    force = {}
    force['forceSPHX'] = np.zeros(np.shape(x))
    force['forceSPHY'] = np.zeros(np.shape(x))
    force['forceSPHZ'] = np.zeros(np.shape(x))
    # original triangulation
    triangles = tri.Triangulation(x, y)
    particles = {'nPart': nPart, 'x': x, 'y': y, 'z': z, 'ux': ux, 'uy': uy, 'uz': uz, 'm': m, 'h': h, 'dt': 0.05}
    particles = DFAfunC.initializeBondsC(particles, triangles)
    # bond breaking
    particles['x'][1] = 1.21
    print('Here')
    force, particles = DFAfunC.computeCohesionForceC(cfg['GENERAL'], particles, force)
    print('Here')
    print(particles['bondStart'])
    print(particles['bondDist'])
    print(particles['bondPart'])
    print(force['forceSPHX'])
    print(force['forceSPHY'])
    print(force['forceSPHZ'])
    bondStart = particles['bondStart']
    bondDist = particles['bondDist']
    bondPart = particles['bondPart']
    for k in range(nPart):
        # loop on all bonded particles
        neighbors = list()
        for ib in range(bondStart[k], bondStart[k + 1]):
            l = bondPart[ib]
            neighbors.append(l)
            if k == 0:
                if l == 1:
                    assert bondDist[ib] == -1
            if k == 0:
                if l == 1:
                    assert bondDist[ib] == -1

        neighbors.sort()
        if k == 0:
            assert neighbors == [1, 2]
        if k == 1:
            assert neighbors == [0, 2]
        if k == 2:
            assert neighbors == [0, 1]
    bondDist.sort()
    assert np.array_equal(bondDist, np.asarray([-1, -1, 1, 1, math.sqrt(2), math.sqrt(2)]))
    assert force['forceSPHX'][0] == 0
    assert force['forceSPHY'][0] == 0
    assert force['forceSPHX'][1] < 0
    assert force['forceSPHY'][1] > 0
    assert force['forceSPHX'][2] > 0
    assert force['forceSPHY'][2] < 0
    assert force['forceSPHX'][1] == -force['forceSPHX'][2]
    assert force['forceSPHY'][1] == -force['forceSPHY'][2]

    # bound not breaking
    particles['x'][1] = 1
    particles = {'nPart': nPart, 'x': x, 'y': y, 'z': z, 'ux': ux, 'uy': uy, 'uz': uz, 'm': m, 'h': h, 'dt': 0.05}
    particles = DFAfunC.initializeBondsC(particles, triangles)
    particles['x'][1] = 1.09
    force, particles = DFAfunC.computeCohesionForceC(cfg['GENERAL'], particles, force)
    print(particles['bondStart'])
    print(particles['bondDist'])
    print(particles['bondPart'])
    print(force['forceSPHX'])
    print(force['forceSPHY'])
    print(force['forceSPHZ'])
    bondStart = particles['bondStart']
    bondDist = particles['bondDist']
    bondPart = particles['bondPart']
    for k in range(nPart):
        # loop on all bonded particles
        neighbors = list()
        for ib in range(bondStart[k], bondStart[k + 1]):
            l = bondPart[ib]
            neighbors.append(l)

        neighbors.sort()
        if k == 0:
            assert neighbors == [1, 2]
        if k == 1:
            assert neighbors == [0, 2]
        if k == 2:
            assert neighbors == [0, 1]
    bondDist.sort()
    assert np.array_equal(bondDist, np.asarray([1, 1, 1, 1, math.sqrt(2), math.sqrt(2)]))

    assert force['forceSPHX'][0] > 0
    assert force['forceSPHY'][0] == 0
    assert force['forceSPHX'][1] < 0
    assert force['forceSPHY'][1] > 0
    assert force['forceSPHX'][2] > 0
    assert force['forceSPHY'][2] < 0
    assert force['forceSPHY'][1] == -force['forceSPHY'][2]
