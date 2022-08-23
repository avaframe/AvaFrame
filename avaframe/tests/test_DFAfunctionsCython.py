"""Tests for module com1DFAtools"""
import numpy as np
import matplotlib.pyplot as plt
import configparser
import pytest

# Local imports
import avaframe.com1DFA.DFAfunctionsCython as DFAfunC
import avaframe.com1DFA.DFAtools as DFAtls


def test_normalizeC(capfd):
    '''test DFAfunctions tools
    norm, norm2, normalize, crossProd and scalProd'''
    x = np.array([1.])
    y = np.array([1.])
    z = np.array([1.])
    norme = DFAfunC.norm(x, y, z)
    norme2 = DFAfunC.norm2(x, y, z)
    xn, yn, zn = DFAfunC.normalize(x, y, z)
    print(xn, yn, zn)
    atol = 1e-10
    assert norme == np.sqrt(3.)
    assert norme2 == 3.
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1., rel=atol)
    assert xn == 1/np.sqrt(3.)
    assert yn == 1/np.sqrt(3.)
    assert zn == 1/np.sqrt(3.)

    x = np.array([0.])
    y = np.array([0.])
    z = np.array([1e-18])
    xn, yn, zn = DFAfunC.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1, rel=atol)

    x = np.array([0.])
    y = np.array([0.])
    z = np.array([0.])
    xn, yn, zn = DFAfunC.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(0, rel=atol)

    x = np.array([1.])
    y = np.array([2.])
    z = np.array([3.])
    xn, yn, zn = DFAfunC.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1., rel=atol)

    x = np.array([1.])
    y = np.array([0.])
    z = np.array([1.])
    xn, yn, zn = DFAfunC.normalize(x, y, z)
    assert np.sqrt(xn*xn + yn*yn + zn*zn) == pytest.approx(1, rel=atol)
    assert xn == pytest.approx(1/np.sqrt(2.), rel=atol)
    assert yn == pytest.approx(0, rel=atol)
    assert zn == pytest.approx(1/np.sqrt(2.), rel=atol)

    x = np.array([1.])
    y = np.array([2.])
    z = np.array([3.])
    x1 = np.array([4.])
    y1 = np.array([5.])
    z1 = np.array([6.])
    xn, yn, zn = DFAfunC.crossProd(x, y, z, x1, y1, z1)
    assert xn == -3
    assert yn == 6
    assert zn == -3

    x = np.array([1.])
    y = np.array([2.])
    z = np.array([3.])
    x1 = np.array([4.])
    y1 = np.array([5.])
    z1 = np.array([6.])
    scal = DFAfunC.scalProd(x, y, z, x1, y1, z1)
    assert scal == 32


def test_getWeightsC(capfd):
    '''getWeights'''
    ncols = 4
    nrows = 3
    X = np.array([0., 1, 2.5, 4., 4., 4.9])
    Y = np.array([0., 0., 2.5, 2.5, 3., 4.9])
    F00 = np.array([[1., 1., 0., 0., 0., 0.],
                    [1./4, 1./4, 1./4, 1./4, 1./4, 1./4],
                    [1., 0.8, 1./4, 0.09999999999999998, 0.07999999999999999, 0.0003999999999999963]])
    F10 = np.array([[0., 0., 0., 0., 0., 0.],
                    [1./4, 1./4, 1./4, 1./4, 1./4, 1./4],
                    [0., 0.2, 1./4, 0.4, 0.32, 0.01959999999999991]])
    F01 = np.array([[0., 0., 0., 0., 0., 0.],
                    [1./4, 1./4, 1./4, 1./4, 1./4, 1./4],
                    [0., 0., 1./4, 0.09999999999999998, 0.11999999999999997, 0.01959999999999991]])
    F11 = np.array([[0., 0., 1., 1., 1., 1.],
                    [1./4, 1./4, 1./4, 1./4, 1./4, 1./4],
                    [0., 0., 1./4, 0.4, 0.48, 0.9604000000000001]])
    csz = 5.
    atol = 1e-10
    # option 0: nearest, 1, unifom, 2 bilinear
    for interpOption in range(3):
        FF00 = F00[interpOption]
        FF10 = F10[interpOption]
        FF01 = F01[interpOption]
        FF11 = F11[interpOption]
        for x, y, ff00, ff10, ff01, ff11 in zip(X, Y, FF00, FF10, FF01, FF11):
            Lx0, Ly0, iCell, f00, f10, f01, f11 = DFAfunC.getCellAndWeights(x, y, ncols, nrows,
                                                                            csz, interpOption)
            print(Lx0, Ly0)
            print(f00, f10, f01, f11)
            assert Lx0 == 0.
            assert Ly0 == 0.
            assert iCell == 0.
            assert ff00 == pytest.approx(f00, rel=atol)
            assert ff10 == pytest.approx(f10, rel=atol)
            assert ff01 == pytest.approx(f01, rel=atol)
            assert ff11 == pytest.approx(f11, rel=atol)

    # same as before but shift 2 in x dir and 1 in y dir
    X = X + 2*csz
    Y = Y + csz
    for interpOption in range(3):
        FF00 = F00[interpOption]
        FF10 = F10[interpOption]
        FF01 = F01[interpOption]
        FF11 = F11[interpOption]
        for x, y, ff00, ff10, ff01, ff11 in zip(X, Y, FF00, FF10, FF01, FF11):
            Lx0, Ly0, iCell, f00, f10, f01, f11 = DFAfunC.getCellAndWeights(x, y, ncols, nrows,
                                                                            csz, interpOption)
            print(Lx0, Ly0, iCell)
            print(f00, f10, f01, f11)
            assert Lx0 == 2.
            assert Ly0 == 1.
            assert iCell == 6.
            assert ff00 == pytest.approx(f00, rel=atol)
            assert ff10 == pytest.approx(f10, rel=atol)
            assert ff01 == pytest.approx(f01, rel=atol)
            assert ff11 == pytest.approx(f11, rel=atol)

    ncols = 4
    nrows = 3
    csz = 1.
    X = np.array([0., 1, -0.5, 0, 0, 1, 3, 3, 3.1, 2])
    Y = np.array([0., -0.5, 1, 1, 2, 0, 0, 2, 1, 2.1])
    cellInd = np.array([0, -1, -1, 4, -1, 1, -1, -1, -1, -1])
    for x, y, res in zip(X, Y, cellInd):
        iCell = DFAfunC.getCells(x, y, ncols, nrows, csz)
        print(iCell)
        assert iCell == res


def test_getScalVect(capfd):
    '''getScalar and getVector'''
    ncols = 100
    nrows = 10
    csz = 5
    interpOption = 2
    X = np.linspace(0, csz*(ncols-1), ncols)
    Y = np.linspace(0, csz*(nrows-1), nrows)
    XX, YY = np.meshgrid(X, Y)
    ZZ = -2*XX + 1000 - YY + 500

    xpExpected = [0, 12.5, 12.5, 52, 95, 494]
    ypExpected = [0, 0, 12.5, 41, 37, 44]
    for xpExp, ypExp in zip(xpExpected, ypExpected):
        Lx0, Ly0, iCell, w0, w1, w2, w3 = DFAfunC.getCellAndWeights(xpExp, ypExp, ncols, nrows,
                                                                    csz, interpOption)
        zScalRes = DFAfunC.getScalar(Lx0, Ly0, w0, w1, w2, w3, ZZ)
        xRes, yRes, zRes = DFAfunC.getVector(Lx0, Ly0, w0, w1, w2, w3, ZZ,
                                             2*ZZ, 3*ZZ)
        res = -2*xpExp + 1000-ypExp + 500
        print(xRes, res)
        atol = 1e-10
        assert zScalRes == pytest.approx(res, rel=atol)
        assert xRes == pytest.approx(res, rel=atol)
        assert yRes == pytest.approx(2*res, rel=atol)
        assert zRes == pytest.approx(3*res, rel=atol)

    zRes = DFAfunC.projOnRaster(np.array(xpExpected), np.array(ypExpected), ZZ, csz, ncols, nrows,
                                interpOption)
    res = -2*np.array(xpExpected) + 1000-np.array(ypExpected) + 500
    atol = 1e-10
    assert np.allclose(zRes, res, atol=atol)


def test_reprojectionC(capfd):
    '''test reprojection'''
    ncols = 100
    nrows = 10
    csz = 5
    header = {}
    header['ncols'] = ncols
    header['nrows'] = nrows
    header['cellsize'] = csz
    dem = {}
    dem['header'] = header
    x0 = 500
    z0 = 1000
    X = np.linspace(0, csz*(ncols-1), ncols)
    Y = np.linspace(0, csz*(nrows-1), nrows)
    XX, YY = np.meshgrid(X, Y)
    ZZ = z0/(x0*x0) * (XX-x0)*(XX-x0)
    dem['rasterData'] = ZZ
    num = 1
    interpOption = 2
    threshold = 0.001

    # build normals corresponding to dem
    dem = DFAtls.getNormalMesh(dem, num)
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']

    # expected normal projection
    xpExpected = 242
    ypExpected = 25
    Lx0, Ly0, iCell, w0, w1, w2, w3 = DFAfunC.getCellAndWeights(xpExpected, ypExpected, ncols,
                                                                nrows, csz, interpOption)
    zpExpected = DFAfunC.getScalar(Lx0, Ly0, w0, w1, w2, w3, ZZ)
    nx, ny, nz = DFAfunC.getVector(Lx0, Ly0, w0, w1, w2, w3, Nx, Ny, Nz)
    nx, ny, nz = DFAfunC.normalize(nx, ny, nz)
    # previous position
    xPrev = 230
    yPrev = 25
    Lx0, Ly0, iCell, w0, w1, w2, w3 = DFAfunC.getCellAndWeights(xPrev, yPrev, ncols, nrows,
                                                                csz, interpOption)
    zPrev = DFAfunC.getScalar(Lx0, Ly0, w0, w1, w2, w3, ZZ)

    # make a point above the parabola at dist d:
    d = 15

    print(xpExpected, zpExpected)
    x1 = xpExpected + d*nx
    y1 = ypExpected + d*ny
    z1 = zpExpected + d*nz
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(X, ZZ[5, :], 'ok-')
    ax.plot(x1, z1, 'ro')
    ax.plot(xpExpected, zpExpected, 'go')

    for reprojectionIterations in range(10):
        xpn, ypn, iCell, Lx0, Ly0, w0, w1, w2, w3 = DFAfunC.samosProjectionIteratrive(x1, y1, z1,
                                                                                      ZZ, Nx, Ny, Nz, csz, ncols, nrows, interpOption, reprojectionIterations)
        zpn = DFAfunC.getScalar(Lx0, Ly0, w0, w1, w2, w3, ZZ)
        print(xpn, zpn)
        ax.plot(xpn, zpn, 'bo')
    ax.set_aspect('equal', 'box')
    # plt.show()
    plt.close(fig)

    atol = 1e-2

    assert xpn == pytest.approx(xpExpected, rel=atol)
    assert ypn == pytest.approx(ypExpected, rel=atol)
    assert zpn == pytest.approx(zpExpected, rel=atol)

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(X, ZZ[5, :], 'ok-')
    ax.plot(x1, z1, 'ro')
    ax.plot(xpExpected, zpExpected, 'go')

    for reprojectionIterations in range(10):
        xpn, ypn, iCell, Lx0, Ly0, w0, w1, w2, w3 = DFAfunC.normalProjectionIteratrive(x1, y1, z1,
                                                                                       ZZ, Nx, Ny, Nz, csz, ncols, nrows, interpOption, reprojectionIterations, threshold)
        zpn = DFAfunC.getScalar(Lx0, Ly0, w0, w1, w2, w3, ZZ)
        print(xpn, zpn)
        ax.plot(xpn, zpn, 'bo')

    ax.plot([xpExpected, x1], [zpExpected, z1], 'b')
    ax.set_aspect('equal', 'box')
    # plt.show()
    plt.close(fig)

    atol = 1e-2

    dist = DFAfunC.norm(xpn-xpExpected, ypn-ypExpected, zpn-zpExpected)
    assert dist <= csz*threshold
    assert xpn == pytest.approx(xpExpected, rel=atol)
    assert ypn == pytest.approx(ypExpected, rel=atol)
    assert zpn == pytest.approx(zpExpected, rel=atol)

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(X, ZZ[5, :], 'ok-')
    ax.plot(x1, z1, 'ro')
    ax.plot(xPrev, zPrev, 'mo')
    ax.plot(xpExpected, zpExpected, 'go')

    for reprojectionIterations in range(10):
        xpn, ypn, zpn, iCell, Lx0, Ly0, w0, w1, w2, w3 = DFAfunC.distConservProjectionIteratrive(
            xPrev, yPrev, zPrev, ZZ, Nx, Ny, Nz, x1, y1, z1, csz, ncols, nrows, interpOption,
            reprojectionIterations, threshold)
        print(xpn, zpn)
        ax.plot(xpn, zpn, 'bo')

    ax.plot([xpExpected, x1], [zpExpected, z1], 'b')
    dist0 = DFAfunC.norm(x1-xPrev, y1-yPrev, z1-zPrev)
    circ = plt.Circle([xPrev, zPrev], dist0, fill=True)
    ax.add_patch(circ)
    ax.set_aspect('equal', 'box')
    # plt.show()
    plt.close(fig)

    atol = 1e-2

    dist0 = DFAfunC.norm(x1-xPrev, y1-yPrev, z1-zPrev)
    dist = DFAfunC.norm(xpn-xPrev, ypn-yPrev, zpn-zPrev)
    assert abs(dist-dist0) <= threshold*(dist0 + csz)


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


def test_SamosATfric(capfd):
    """ Test the account4FrictionForce function"""
    tau0 = 0
    Rs0 = 0.222
    kappa = 0.43
    R = 0.05
    B = 4.13
    rho = 0.1
    mu = 0.155
    uMag = 10
    sigmaB = 10
    h = 1
    tau = DFAfunC.SamosATfric(rho, tau0, Rs0, mu, kappa, B, R, uMag, sigmaB, h)
    print(tau)
    assert tau == 1.9128193823277053
    tau0 = 2
    tau = DFAfunC.SamosATfric(rho, tau0, Rs0, mu, kappa, B, R, uMag, sigmaB, h)
    print(tau)
    assert tau == 3.9128193823277053


def test_updatePositionC():
    """ test updating the position of particles """

    # TODO: make test also if velocity in z not zero!! - so to test also when reprojecting onto surface

    # setup required input
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'stopCrit': '0.01', 'stopCritIni': '0.1', 'stopCritIniSmall': '1.001', 'stopCritType': 'kinEnergy',
                      'uFlowingThreshold': '0.1', 'gravAcc': '9.81', 'velMagMin': '1.e-6',  'rho': '100.',
                      'interpOption': '2',   'explicitFriction': '0', 'centeredPosition': '1',
                      'reprojMethodPosition': '2', 'reprojectionIterations': '5', 'thresholdProjection': '0.001' }

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 's': np.asarray([0., 0., 0.]),
                  'sCor': np.asarray([0., 0., 0.]), 'l': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 0.0,
                  'peakForceSPH': 0.0, 'forceSPHIni': 0.0, 'nPart': 3,
                  'peakMassFlowing': 0.0, 'iterate': True}
    particles['potentialEne'] = np.sum(9.81 * particles['z'] * particles['m'])

    demHeader = {}
    demHeader['xllcenter'] = 0.0
    demHeader['yllcenter'] = 0.0
    demHeader['cellsize'] = 5.0
    demHeader['noDataValue'] = -9999
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

    typeStop = 0

    # kinetic energy new
    kinEneNew = 0.0
    potEneNew = 0.0
    for k in range(3):
        kinEneNew = kinEneNew + particles['m'][k] * np.sqrt(5.5**2 +5.5**2 + 0**2)**2 * 0.5
        potEneNew = potEneNew + particles['m'][k] * 9.81 + 0.0

    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, typeStop=typeStop)

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([3.25, 4.25, 5.25]))
    assert np.array_equal(particles['y'], np.asarray([5.25, 6.25, 7.25]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == True

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 's': np.asarray([0., 0., 0.]),
                  'sCor': np.asarray([0., 0., 0.]), 'l': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 100000.0,
                  'peakForceSPH': 0.0, 'forceSPHIni': 0.0, 'nPart': 3,
                  'peakMassFlowing': 0.0, 'iterate': True}
    particles['potentialEne'] = np.sum(9.81 * particles['z'] * particles['m'])

    # call function to be tested
    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, typeStop=typeStop)

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([3.25, 4.25, 5.25]))
    assert np.array_equal(particles['y'], np.asarray([5.25, 6.25, 7.25]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([5.5, 5.5, 5.5]), atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == False

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 's': np.asarray([0., 0., 0.]),
                  'sCor': np.asarray([0., 0., 0.]), 'l': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 10000.0,
                  'peakForceSPH': 100000.0, 'forceSPHIni': 1.e5, 'nPart': 3,
                  'peakMassFlowing': 0.0, 'iterate': True}
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
    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, typeStop=typeStop)
    print('sph', particles['peakForceSPH'], sphForceNew)

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([6., 7., 8.]))
    assert np.array_equal(particles['y'], np.asarray([8., 9., 10.]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == False

    particles = {'dt': 1.0, 'm': np.asarray([10., 10., 10.]), 'idFixed': np.asarray([0., 0., 0.]), 's': np.asarray([0., 0., 0.]),
                  'sCor': np.asarray([0., 0., 0.]), 'l': np.asarray([0., 0., 0.]), 'x': np.asarray([0., 1., 2.]), 'y': np.asarray([2., 3., 4.]),
                  'z': np.asarray([1., 1., 1.]), 'ux': np.asarray([1., 1., 1.]), 'uy': np.asarray([1., 1., 1.]),
                  'uz': np.asarray([0., 0., 0.]), 'kineticEne': 0.0, 'peakKinEne': 10000.0,
                  'peakForceSPH': 1000.0, 'forceSPHIni': 1.e5, 'nPart': 3,
                  'peakMassFlowing': 0.0, 'iterate': True}
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
    particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, typeStop=typeStop)
    print('sph', particles['peakForceSPH'], sphForceNew)

    assert np.array_equal(particles['m'], np.asarray([10., 10., 10.]))
    assert np.array_equal(particles['x'], np.array([6., 7., 8.]))
    assert np.array_equal(particles['y'], np.asarray([8., 9., 10.]))
    assert np.array_equal(particles['z'], np.asarray([1., 1., 1.]))
    assert np.allclose(particles['ux'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.allclose(particles['uy'], np.asarray([11., 11., 11.]), atol=1.e-4)
    assert np.array_equal(particles['uz'], np.asarray([0., 0., 0.]))
    assert particles['massEntrained'] == 0.0
    assert particles['nPart'] == 3
    assert (kinEneNew- 1.e-4) < particles['kineticEne'] < (kinEneNew+1.e-4)
    assert (potEneNew-1.e-4) < particles['potentialEne'] < (potEneNew +1.e-4)
    assert particles['iterate'] == True


def test_computeTravelAngle():
    # first compute travel angle for each particle
    # get parent Id in order to  get the first z position
    parentID = np.array([0, 1, 2, 0])
    nPart = 5
    # get z0
    zPartArray0 = np.array([10.0, 9.0, 8.0])
    s = np.array([10.0, 10.0, 0.0, 10.0])
    z = np.array([0.0, 0.0, 0.0, 1.0])
    particles = {'nPart': nPart, 'parentID': parentID, 's': s, 'z': z}
    particles = DFAfunC.computeTravelAngleC(particles, zPartArray0)
    print(particles['travelAngle'])
    gamma = particles['travelAngle']
    assert gamma[2] == 0
    assert gamma[0] == 45
    assert gamma[1] == pytest.approx(41.9872125, rel=1e-6)
    assert gamma[3] == pytest.approx(41.9872125, rel=1e-6)
