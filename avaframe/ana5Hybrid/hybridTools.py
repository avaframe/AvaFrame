import pathlib
import numpy as np
import shapefile
import matplotlib.pyplot as plt

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls
from avaframe.in3Utils import geoTrans
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outCom1DFA as outCom1DFA

debug = False


def getCom1DFAPath(particlesList, dem):
    """ compute part, mass and energy averaged path from particles

    Also returns the averaged velocity ant kinetik energy associated

    Parameters
    -----------
    particlesList: list
        list of particles dict
    dem: dict
        dem dict

    Returns
    --------
    avaProfilePart: dict
        particle averaged profile
    avaProfileMass: dict
        mass averaged profile
    avaProfileKE: dict
        kinetic energy averaged profile
    """
    # get DEM Path
    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']

    proList = ['x', 'y', 'z', 's', 'sCor']
    avaProfilePart = {'x': np.empty((0, 1)), 'y': np.empty((0, 1)), 'z': np.empty((0, 1)), 's': np.empty((0, 1)),
                      'sCor': np.empty((0, 1)), 'v2': np.empty((0, 1)), 'ekin': np.empty((0, 1))}
    avaProfileMass = {'x': np.empty((0, 1)), 'y': np.empty((0, 1)), 'z': np.empty((0, 1)), 's': np.empty((0, 1)),
                      'sCor': np.empty((0, 1)), 'v2': np.empty((0, 1)), 'ekin': np.empty((0, 1))}
    avaProfileKE = {'x': np.empty((0, 1)), 'y': np.empty((0, 1)), 'z': np.empty((0, 1)), 's': np.empty((0, 1)),
                    'sCor': np.empty((0, 1)), 'v2': np.empty((0, 1)), 'ekin': np.empty((0, 1))}
    # loop on each particle dictionary (ie each time step saved)
    for particles in particlesList:
        if particles['nPart'] > 0:
            m = particles['m']
            ux = particles['ux']
            uy = particles['uy']
            uz = particles['uz']
            u = DFAtls.norm(ux, uy, uz)
            U2 = u*u
            kineticEne = 0.5*m*u*u
            kineticEneSum = np.nansum(kineticEne)

            # kinetic energy-averaged path
            if kineticEneSum > 0:
                w = kineticEne
                for prop in proList:
                    avaProfileKE[prop] = np.append(avaProfileKE[prop], np.average(particles[prop], weights=w))
                avaProfileKE['v2'] = np.append(avaProfileKE['v2'], np.average(U2, weights=w))
                avaProfileKE['ekin'] = np.append(avaProfileKE['ekin'], kineticEneSum)
            else:
                w = m
                for prop in proList:
                    avaProfileKE[prop] = np.append(avaProfileKE[prop], np.average(particles[prop], weights=w))
                avaProfileKE['v2'] = np.append(avaProfileKE['v2'], np.average(U2, weights=w))
                avaProfileKE['ekin'] = np.append(avaProfileKE['ekin'], kineticEneSum)

            # mass-averaged path
            w = m
            for prop in proList:
                avaProfileMass[prop] = np.append(avaProfileMass[prop], np.average(particles[prop], weights=w))
            avaProfileMass['v2'] = np.append(avaProfileMass['v2'], np.average(U2, weights=w))
            avaProfileMass['ekin'] = np.append(avaProfileMass['ekin'], kineticEneSum)

            # particle-averaged path
            w = None
            for prop in proList:
                avaProfilePart[prop] = np.append(avaProfilePart[prop], np.average(particles[prop], weights=w))
            avaProfilePart['v2'] = np.append(avaProfilePart['v2'], np.average(U2, weights=w))
            avaProfilePart['ekin'] = np.append(avaProfilePart['ekin'], kineticEneSum)

    avaProfilePart['x'] = avaProfilePart['x'] + xllc
    avaProfilePart['y'] = avaProfilePart['y'] + yllc
    avaProfileMass['x'] = avaProfileMass['x'] + xllc
    avaProfileMass['y'] = avaProfileMass['y'] + yllc
    avaProfileKE['x'] = avaProfileKE['x'] + xllc
    avaProfileKE['y'] = avaProfileKE['y'] + yllc
    return avaProfilePart, avaProfileMass, avaProfileKE


def extendCom1DFAPath(dem, particlesIni, avaProfilePart, avaProfileMass, avaProfileKE):
    """ extend the DFA path at the top and bottom

    Parameters
    -----------
    dem: dict
        dem dict
    particlesIni: dict
        initial particles dict
    avaProfilePart: dict
        particle averaged profile
    avaProfileMass: dict
        mass averaged profile
    avaProfileKE: dict
        kinetic energy averaged profile

    Returns
    --------
    avaProfilePart: dict
        extended particle averaged profile
    avaProfileMass: dict
        extended mass averaged profile
    avaProfileKE: dict
        extended kinetic energy averaged profile
    """
    for profile in [avaProfilePart, avaProfileMass, avaProfileKE]:
        # project the profile on the dem
        profile, _ = geoTrans.projectOnRaster(dem, profile, interp='bilinear')
        indNotNan = np.where(~np.isnan(profile['z']))
        profile['x'] = profile['x'][indNotNan]
        profile['y'] = profile['y'][indNotNan]
        profile['z'] = profile['z'][indNotNan]
        # resample the profile
        profile, _ = geoTrans.prepareLine(dem, profile, distance=dem['header']['cellsize']*5, Point=None)
        profile = extendProfileTop(dem, particlesIni, profile)
        # resample the profile
        profile, _ = geoTrans.prepareLine(dem, profile, distance=dem['header']['cellsize']*5, Point=None)
        profile = extendProfileBottom(dem, profile)
        # resample the profile
        profile, _ = geoTrans.prepareLine(dem, profile, distance=dem['header']['cellsize']*5, Point=None)
        # remove points that lay outside of the dem
        isNotNan = ~(np.isnan(profile['z']))
        profile['x'] = profile['x'][isNotNan]
        profile['y'] = profile['y'][isNotNan]
        profile['s'] = profile['s'][isNotNan]
        profile['z'] = profile['z'][isNotNan]
    return avaProfilePart, avaProfileMass, avaProfileKE


def extendProfileTop(dem, particlesIni, profile):
    """ extend the DFA path at the top

    Find the direction in which to extend considerind the first point of the profile
    and a few following ones (distFromFirt <= 30 * csz). Extend in this diretion until
    the z of the highest particle in the release is reached.

    Parameters
    -----------
    dem: dict
        dem dict
    particlesIni: dict
        initial particles dict
    profile: dict
        profile to extend

    Returns
    --------
    profile: dict
        extended profile
    """
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
    # get first particle of the path
    xFirst = profile['x'][0]
    yFirst = profile['y'][0]
    zFirst = profile['z'][0]
    # compute distance from first point:
    r = DFAtls.norm(profile['x']-xFirst, profile['y']-yFirst, profile['z']-zFirst)
    # find the following first points
    pointsOfInterestFirst = np.where(r <= 30 * csz)[0]
    xInterest = profile['x'][pointsOfInterestFirst]
    yInterest = profile['y'][pointsOfInterestFirst]
    zInterest = profile['z'][pointsOfInterestFirst]
    # find the direction in which we need to extend the path
    vDirX = xInterest - xFirst
    vDirY = yInterest - yFirst
    vDirZ = zInterest - zFirst
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    vDirX = np.sum(vDirX)
    vDirY = np.sum(vDirY)
    vDirZ = np.sum(vDirZ)
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    # find the point in this direction that has the same z as the highest particle
    xExtTop, yExtTop, zExtTop = findPointOnDEM(dem, vDirX, vDirY, vDirZ, zHighest, xFirst, yFirst, zFirst)
    zExtTop, _ = geoTrans.projectOnGrid(xExtTop, yExtTop, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    # extend profile
    profile['x'] = np.append(xExtTop, profile['x'])
    profile['y'] = np.append(yExtTop, profile['y'])
    profile['z'] = np.append(zExtTop, profile['z'])

    if debug:
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        ax.set_title('Extend path towards the top')
        ax.plot(particlesIni['x'] + xllc, particlesIni['y'] + yllc, '.c', label='particles at t=0s')
        ax.plot(xHighest, yHighest, '.r', label='highest particle at t=0s')
        ax.plot(profile['x'][1:], profile['y'][1:], '.k', label='mass averaged path')
        ax.plot(xInterest[1:], yInterest[1:], '.m', markersize=10, label='points considered to find drection')
        ax.plot(xInterest[0], yInterest[0], '.b', markersize=10, label='top point of the mass averaged path')
        ax.plot(profile['x'][0], profile['y'][0], '.g', label='point at the same hight as the highest point in the release \n and in the extention direction')
        ax.plot(profile['x'][0:2], profile['y'][0:2], 'k--', label='extended path')
        plt.legend()
        plt.show()
    return profile


def findPointOnDEM(dem, vDirX, vDirY, vDirZ, zHighest, xFirst, yFirst, zFirst):
    """ find point on dem given a direction and a z value to reach

    Parameters
    -----------
    dem: dict
        dem dict
    vDirX, vDirY, vDirZ: floats
        x, y and z components of the direction in which to extend
    zHighest: float
        z value to reach
    xFirst, yFirst, zFirst: floats
        x, y and z coordinates of the starting point

    Returns
    --------
    xExtTop, yExtTop, zExtTop:floats
        x, y and z coordinates of the point found
    """
    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    zRaster = dem['rasterData']
    gamma = (zHighest - zFirst) / vDirZ * np.linspace(0.25, 2, 100)
    xArray = xFirst + gamma * vDirX
    yArray = yFirst + gamma * vDirY
    zArray, _ = geoTrans.projectOnGrid(xArray, yArray, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    idx = np.nanargmin(np.abs(zArray - np.array([zHighest])))
    xExtTop = np.array([xFirst + gamma[idx] * vDirX])
    yExtTop = np.array([yFirst + gamma[idx] * vDirY])
    zExtTop = np.array([zArray[idx]])
    return xExtTop, yExtTop, zExtTop


def extendProfileBottom(dem, profile):
    """ extend the DFA path at the bottom

    Find the direction in which to extend considerind the last point of the profile
    and a few previous ones but discardings the one that are too close ( 2* csz < distFromLast <= 30 * csz).
    Extend in this diretion for a distance 0.2 * length of the path.

    Parameters
    -----------
    dem: dict
        dem dict
    profile: dict
        profile to extend

    Returns
    --------
    profile: dict
        extended profile
    """
    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    zRaster = dem['rasterData']
    # get last point
    xLast = profile['x'][-1]
    yLast = profile['y'][-1]
    zLast = profile['z'][-1]
    sLast = profile['s'][-1]
    # compute distance from last point:
    r = DFAtls.norm(profile['x']-xLast, profile['y']-yLast, profile['z']-zLast)
    # find the previous points
    pointsOfInterestLast = np.where((r < 20 * csz) & (r > 2* csz))[0]
    xInterest = profile['x'][pointsOfInterestLast]
    yInterest = profile['y'][pointsOfInterestLast]
    zInterest = profile['z'][pointsOfInterestLast]
    # find the direction in which we need to extend the path
    vDirX = xLast - xInterest
    vDirY = yLast - yInterest
    vDirZ = zLast - zInterest
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    vDirX = np.sum(vDirX)
    vDirY = np.sum(vDirY)
    vDirZ = np.sum(vDirZ)
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), np.array([vDirZ]))
    # extend in this direction
    gamma = 0.2 * sLast / np.sqrt(vDirX**2 + vDirY**2)
    xExtBottom = np.array([xLast + gamma * vDirX])
    yExtBottom = np.array([yLast + gamma * vDirY])
    # project on cszDEM
    zExtBottom, _ = geoTrans.projectOnGrid(xExtBottom, yExtBottom, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    # extend profile
    profile['x'] = np.append(profile['x'], xExtBottom)
    profile['y'] = np.append(profile['y'], yExtBottom)
    profile['z'] = np.append(profile['z'], zExtBottom)

    if debug:
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        ax.set_title('Extend path towards the bottom')
        ax.plot(profile['x'][:-1], profile['y'][:-1], '.k', label='mass averaged path')
        ax.plot(xInterest, yInterest, '.m', markersize=10, label='points considered to find drection')
        ax.plot(xInterest[-1], yInterest[-1], '.b', markersize=10, label='bottom point of the mass averaged path')
        ax.plot(profile['x'][0], profile['y'][0], '.g', label='point in the extention direction at distance \n 0.2 x path length from the bottom point')
        ax.plot(profile['x'][-2:], profile['y'][-2:], 'k--', label='extended path')
        plt.legend()
        plt.show()
    return profile


def plotHybridRes(avalancheDir, resAB, resABNew, name, pathDict, simID, rasterTransfo, resAnalysis, dem, demOri,
                  particlesList, fieldsList, avaProfileMass, avaProfileMassNew):

    headerOri = demOri['header']
    xllcOri = headerOri['xllcenter']
    yllcOri = headerOri['yllcenter']
    # Plots
    indSim = pathDict['simID'].index(simID)
    # sAIMEC = resAnalysis['runout'][0][indSim]
    xAIMEC = resAnalysis['runout'][1][indSim]
    yAIMEC = resAnalysis['runout'][2][indSim]
    # zAIMEC = rasterTransfo['z'][indSim]

    alpha = resAB[name]['alpha']
    ids_alpha = resAB[name]['ids_alpha']
    # indSplit = resAB[name]['indSplit']
    xAB = resAB[name]['x'][ids_alpha]
    yAB = resAB[name]['y'][ids_alpha]
    # sAB = resAB[name]['s'][ids_alpha]
    # zAB = resAB[name]['z'][ids_alpha]
    # fAB = resAB[name]['f']
    # splitPoint = resAB[name]['splitPoint']

    alphaNew = resABNew[name]['alpha']
    ids_alphaNew = resABNew[name]['ids_alpha']
    xABNew = resABNew[name]['x'][ids_alphaNew]
    yABNew = resABNew[name]['y'][ids_alphaNew]

    V2Path = avaProfileMassNew['v2']
    EkinPath = avaProfileMassNew['ekin']

    g = 9.81
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    ax1 = outCom1DFA.addDem2Plot(ax1, dem, what='slope')
    ax1 = outCom1DFA.addResult2Plot(ax1, dem, fieldsList[-1]['pta'], 'pta')
    ax1.plot(xAIMEC - xllcOri, yAIMEC - yllcOri, 'X', markersize=8, color='0.7', label='Runout point com1DFA (AIMEC pfd=0m)')
    ax1.plot(xAB - xllcOri, yAB - yllcOri, 'x', markersize=8, label=r'com2AB $\alpha$ point')
    ax1.plot(xABNew - xllcOri, yABNew - yllcOri, '+', markersize=8, label=r'com2AB $\alpha$ point with new path')
    ax1.plot(avaProfileMass['x'] - xllcOri, avaProfileMass['y'] - yllcOri, 'b-', label='Initial center of mass\navalanche path extended')
    ax1.plot(avaProfileMassNew['x'] - xllcOri, avaProfileMassNew['y'] - yllcOri, 'k--', label='final center of mass\navalanche path extended')
    # ax1 = outCom1DFA.addParticles2Plot(particlesList[0], ax1, dem, whatS='h')
    ax1 = outCom1DFA.addParticles2Plot(particlesList[-1], ax1, dem, whatS='h')
    # set titel of output png
    title = ('myTitle1')
    ax1.legend(loc='lower left')
    plt.show()
    path = pathlib.Path(avalancheDir, 'Outputs')
    pU.saveAndOrPlot({'pathResult': path}, title, fig1)

    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    cmap = pU.cmapPlasma
    cmap.set_under(color='w')
    projectedZ = 'yes'
    if projectedZ == 'yes':
        avaProfileMassNew, _ = geoTrans.projectOnRaster(demOri, avaProfileMassNew, interp='bilinear')
    Zene = avaProfileMassNew['z'] + V2Path/(2*g)
    # Zene = z_m[0] + V2Path/(2*g)
    # Colorbar: kinietische Energie [J]
    scat = ax2.scatter(avaProfileMassNew['sCor'], Zene, marker='s', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Gesamtenergie(s_mod)')
    scat = ax2.scatter(avaProfileMassNew['s'], Zene, marker='o', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Gesamtenergie(s_real)')
    cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('kinetische Energie [J]')

    ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'], 'b-', label='Z_av(s_real)')
    ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'], 'b-.', label='Z_av(s_mod)')
    ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'], 'k-', label='Z_true(s_real)')
    ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'], 'k-.', label='Z_true(s_mod)')
    # ax2.plot(part_aimec[:, 3], z_aimec, 'k--', label='lin. extrapol. Lp_m_projZ')
    # ax2.plot(part_aimec[:, 3], part_aimec[:, 2], 'k-.', label='lin. extrapol. Lp_m')
    GK = avaProfileMassNew['sCor'][-1] * np.tan(alpha*np.pi/180)
    z_enda = avaProfileMassNew['z'][0] - GK
    s_geomL = [avaProfileMassNew['sCor'][0], avaProfileMassNew['sCor'][-1]]
    z_geomL = [avaProfileMassNew['z'][0], z_enda]
    ax2.plot(s_geomL, z_geomL, 'r-', linewidth=0.3, label='alpha line from Z_av func s_mod')
    GK = avaProfileMassNew['s'][-1] * np.tan(alpha*np.pi/180)
    z_enda = avaProfileMassNew['z'][0] - GK
    s_geomL = [avaProfileMassNew['s'][0], avaProfileMassNew['s'][-1]]
    z_geomL = [avaProfileMassNew['z'][0], z_enda]
    ax2.plot(s_geomL, z_geomL, 'g-', linewidth=0.3, label='alpha line from Z_true func s_mod')
    # ax2.plot(sAB, fAB, 'b-', linewidth=0.3, label='Alphalinie com2AB')
    # ax2.plot(sAB, zAB, 'x', markersize=8, label='A com2AB   = %1.2f' % zAB + 'm')
    # ax2.plot(sAIMEC, zAB, 'X', markersize=8, color='0.7', label='A ana3AIMEC   = %1.2f' % zAIMEC + 'm')

    ax2.set_xlabel('s [m]', fontsize=pU.fs)
    ax2.set_ylabel('z [m]', fontsize=pU.fs)
    ax2.legend(loc='lower left')

    # set titel of output png
    title = ('myTitle')

    plt.show()
    path = pathlib.Path(avalancheDir, 'Outputs')
    pU.saveAndOrPlot({'pathResult': path}, title, fig2)
