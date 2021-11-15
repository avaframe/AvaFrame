import numpy as np
import math
import shapefile
import matplotlib.pyplot as plt

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls
from avaframe.in3Utils import geoTrans
import avaframe.out3Plot.plotUtils as pU


def getCom1DFAPath(particlesList, tSaveList, dem):
    # get DEM Path
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']

    relField = np.zeros((nrows, ncols))
    V2Path = np.empty((0, 1))
    EkinPath = np.empty((0, 1))
    # EpotPath = np.empty((0, 1))
    # ax = plt.subplot(111)
    part_p = np.empty((0, 5))
    part_m = np.empty((0, 5))
    part_kE = np.empty((0, 5))
    count = 0
    for t in tSaveList:
        particles = particlesList[count]
        m = particles['m']
        x = particles['x'] + xllc
        y = particles['y'] + yllc
        z = particles['z']
        ux = particles['ux']
        uy = particles['uy']
        uz = particles['uz']
        s = particles['s']
        sCor = particles['sCor']
        u = DFAtls.norm(ux, uy, uz)
        U2 = u*u
        nPart = particles['nPart']
        kineticEne = 0.5*m*u*u
        kineticEneSum = np.sum(kineticEne)
        # mass-averaged path
        if kineticEneSum <= 100:
            pond = np.ones(np.shape(kineticEne))
            pondSum = nPart
        else:
            pond = kineticEne
            pondSum = kineticEneSum

        pond = m
        pondSum = np.sum(m)
        v2coE = np.sum(pond*U2)/pondSum
        V2Path = np.append(V2Path, v2coE)

        # update energy
        EkinPath = np.append(EkinPath, kineticEneSum)
        # kinetic energy-averaged path
        if kineticEneSum > 0:
            temp_kE = np.array([[np.sum(kineticEne*x)/kineticEneSum, np.sum(kineticEne*y)/kineticEneSum,
                    np.sum(kineticEne*z)/kineticEneSum, np.sum(kineticEne*s)/kineticEneSum,
                    np.sum(kineticEne*sCor)/kineticEneSum]])
        else:
            temp_kE = np.array([[np.sum(m*x)/np.sum(m), np.sum(m*y)/np.sum(m),
                    np.sum(m*z)/np.sum(m), np.sum(m*s)/np.sum(m),
                    np.sum(m*sCor)/np.sum(m)]])
        # mass-averaged path
        temp_m = np.array([[np.sum(m*x)/np.sum(m), np.sum(m*y)/np.sum(m),
                np.sum(m*z)/np.sum(m), np.sum(m*s)/np.sum(m),
                    np.sum(m*sCor)/np.sum(m)]])
        # particle-averaged path
        temp_p = np.array([[np.sum(x)/nPart, np.sum(y)/nPart, np.sum(z)/nPart,
                np.sum(s)/nPart, np.sum(sCor)/nPart]])
        # add a line to part
        part_p = np.append(part_p, temp_p, axis=0)
        part_m = np.append(part_m, temp_m, axis=0)
        part_kE = np.append(part_kE, temp_kE, axis=0)
        count = count + 1

        avaProfilePart = {'x': part_p[:, 0], 'y': part_p[:, 1], 'z': part_p[:, 2], 's': part_p[:, 3], 'sCor': part_p[:, 4]}
        avaProfileMass = {'x': part_m[:, 0], 'y': part_m[:, 1], 'z': part_m[:, 2], 's': part_m[:, 3], 'sCor': part_m[:, 4]}
        avaProfileKE = {'x': part_kE[:, 0], 'y': part_kE[:, 1], 'z': part_kE[:, 2], 's': part_kE[:, 3], 'sCor': part_kE[:, 4]}
    return avaProfilePart, avaProfileMass, avaProfileKE, V2Path, EkinPath


def elongateCom1DFAPath(dem, particlesIni, avaProfilePart, avaProfileMass, avaProfileKE):

    for profile in [avaProfilePart, avaProfileMass, avaProfileKE]:
        profile, _ = geoTrans.prepareLine(dem, profile, distance=dem['header']['cellsize'], Point=None)
        profile = extendProfileTop(dem, particlesIni, profile)
        profile = extendProfileBottom(dem, profile)
        profile, _ = geoTrans.prepareLine(dem, profile, distance=dem['header']['cellsize'], Point=None)
        isNotNan = ~(np.isnan(profile['z']))
        profile['x'] = profile['x'][isNotNan]
        profile['y'] = profile['y'][isNotNan]
        profile['s'] = profile['s'][isNotNan]
        profile['z'] = profile['z'][isNotNan]
    return avaProfilePart, avaProfileMass, avaProfileKE


def extendProfileTop(dem, particlesIni, profile):
    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    zRaster = dem['rasterData']
    # get highest particle
    indHighest = np.argmax(particlesIni['z'])
    xHighest = particlesIni['x'][indHighest] + xllc
    yHighest = particlesIni['y'][indHighest] + yllc
    zHighest = particlesIni['z'][indHighest]

    xFirst = profile['x'][0]
    yFirst = profile['y'][0]
    zFirst = profile['z'][0]
    sFirst = profile['s'][0]
    # compute distance from first point:
    r = DFAtls.norm(profile['x']-xFirst, profile['y']-yFirst, profile['z']-zFirst)
    pointsOfInterestFirst = np.where(r < 10 * csz)[0]
    xInterest = profile['x'][pointsOfInterestFirst]
    yInterest = profile['y'][pointsOfInterestFirst]
    zInterest = profile['z'][pointsOfInterestFirst]

    vDirX = xInterest - xFirst
    vDirY = yInterest - yFirst
    vDirZ = zInterest - zFirst
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    vDirX = np.sum(vDirX)
    vDirY = np.sum(vDirY)
    vDirZ = np.sum(vDirZ)
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    gamma = (zHighest - zFirst) / vDirZ
    xExtTop = np.array([xFirst + gamma * vDirX])
    yExtTop = np.array([yFirst + gamma * vDirY])
    zExtTop = np.array([zFirst + gamma * vDirZ])
    sExtTop = np.array([sFirst - gamma**2 * np.sqrt(vDirX**2 + vDirY**2)])
    zExtTop, _ = geoTrans.projectOnGrid(xExtTop, yExtTop, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    profile['x'] = np.append(xExtTop, profile['x'])
    profile['y'] = np.append(yExtTop, profile['y'])
    profile['z'] = np.append(zExtTop, profile['z'])
    profile['s'] = np.append(sExtTop, profile['s'])
    # make sure the s coordinate starts at 0
    profile['s'] = profile['s'] - profile['s'][0]

    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Extend path')
    ax.plot(particlesIni['x'] + xllc, particlesIni['y'] + yllc, '.k')
    ax.plot(xHighest, yHighest, '.r')
    ax.plot(profile['x'], profile['y'], '.b')
    ax.plot(profile['x'][pointsOfInterestFirst], profile['y'][pointsOfInterestFirst], '.m')
    ax.plot(profile['x'][0], profile['y'][0], '.g')
    plt.legend()
    plt.show()
    return profile


def extendProfileBottom(dem, profile):
    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    zRaster = dem['rasterData']

    xLast = profile['x'][-1]
    yLast = profile['y'][-1]
    zLast = profile['z'][-1]
    sLast = profile['s'][-1]
    # compute distance from first point:
    r = DFAtls.norm(profile['x']-xLast, profile['y']-yLast, profile['z']-zLast)
    pointsOfInterestLast = np.where(r < 10 * csz)[0]
    xInterest = profile['x'][pointsOfInterestLast]
    yInterest = profile['y'][pointsOfInterestLast]
    zInterest = profile['z'][pointsOfInterestLast]

    vDirX = xLast - xInterest
    vDirY = yLast - yInterest
    vDirZ = zLast - zInterest
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    vDirX = np.sum(vDirX)
    vDirY = np.sum(vDirY)
    vDirZ = np.sum(vDirZ)
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    gamma = 0.01 * sLast / np.sqrt(vDirX**2 + vDirY**2)
    xExtBottom = np.array([xLast + gamma * vDirX])
    yExtBottom = np.array([yLast + gamma * vDirY])
    zExtBottom = np.array([zLast + gamma * vDirZ])
    sExtBottom = np.array([sLast + gamma**2 * np.sqrt(vDirX**2 + vDirY**2)])
    zExtBottom, _ = geoTrans.projectOnGrid(xExtBottom, yExtBottom, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    profile['x'] = np.append(profile['x'], xExtBottom)
    profile['y'] = np.append(profile['y'], yExtBottom)
    profile['z'] = np.append(profile['z'], zExtBottom)
    profile['s'] = np.append(profile['s'], sExtBottom)

    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Extend path')
    ax.plot(profile['x'], profile['y'], '.b')
    ax.plot(profile['x'][pointsOfInterestLast], profile['y'][pointsOfInterestLast], '.m')
    ax.plot(profile['x'][-1], profile['y'][-1], '.g')
    plt.legend()
    plt.show()
    return profile


def writeLine2SHPfile(avaProfile, lineName, fileName):
    nPnts = np.size(avaProfile['x'])
    lineArray = np.zeros((nPnts, 3))
    for ind, x, y, z in zip(range(nPnts), avaProfile['x'], avaProfile['y'], avaProfile['z']):
        lineArray[ind, :] = np.array([x, y, z])
    w = shapefile.Writer(str(fileName))
    w.field('name', 'C')
    w.line([lineArray])
    w.record(lineName)
    w.close()
