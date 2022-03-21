"""
    Tools for generating an avalanche simulation from a DFA simulation
"""

# Load modules
import math
import numpy as np
import logging

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.in3Utils.geoTrans as gT
import avaframe.com1DFA.particleTools as particleTools
import avaframe.com1DFA.DFAtools as DFAtls
from avaframe.com1DFA import com1DFA
import avaframe.out3Plot.outDebugPlots as debPlot
# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')


def generateAvragePath(avalancheDir, pathFromPart, simName, dem):
    """ extract path from fileds or particles

    Parameters
    -----------
    avalancheDir: pathlib
        avalanche directory pathlib path
    pathFromPart: boolean
        should the path be extracted from particles?
    simName: str
        simulation name
    dem: dict
        com1DFA simulation dictionary

    Returns
    --------
    avaProfileMass: dict
        mass averaged profile
    """
    if pathFromPart:
        particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName=simName, flagAvaDir=True,
                                                                       comModule='com1DFA')
        # postprocess to extract path and energy line
        avaProfileMass = getDFAPathFromPart(particlesList, dem)
    else:
        particlesList = ''
        # read field
        fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['FD', 'FV', 'FM'], simName=simName,
                                                               flagAvaDir=True, comModule='com1DFA')
        # postprocess to extract path and energy line
        avaProfileMass = getDFAPathFromField(fieldsList, fieldHeader, dem)

    return avaProfileMass


def getDFAPathFromPart(particlesList, dem):
    """ compute mass averaged path from particles

    Also returns the averaged velocity and kinetic energy associated

    Parameters
    -----------
    particlesList: list
        list of particles dict
    dem: dict
        dem dict

    Returns
    --------
    avaProfileMass: dict
        mass averaged profile
    """
    proList = ['x', 'y', 'z', 's', 'sCor']
    avaProfileMass = {'v2': np.empty((0, 1)), 'ekin': np.empty((0, 1))}
    for prop in proList:
        avaProfileMass[prop] = np.empty((0, 1))
        avaProfileMass[prop + 'std'] = np.empty((0, 1))

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

            # mass-averaged path
            avaProfileMass = appendAverageStd(proList, avaProfileMass, particles, U2, kineticEneSum, m)

    return avaProfileMass


def getDFAPathFromField(fieldsList, fieldHeader, dem):
    """ compute mass averaged path from fields

    Also returns the averaged velocity and kinetic energy associated
    The dem and fieldsList (FD, FV and FM) need to have identical dimentions and cell size.

    Parameters
    -----------
    fieldsList: list
        time sorted list of fields dict
    fieldHeader: dict
        field header dict
    dem: dict
        dem dict

    Returns
    --------
    avaProfileMass: dict
        mass averaged profile
    """
    # get DEM
    demRaster = dem['rasterData']
    # get fields header
    ncols = fieldHeader['ncols']
    nrows = fieldHeader['nrows']
    xllc = fieldHeader['xllcenter']
    yllc = fieldHeader['yllcenter']
    csz = fieldHeader['cellsize']
    X, Y = gT.makeCoordinateGrid(xllc, yllc, csz, ncols, nrows)

    # get normal vector of the grid mesh

    proList = ['x', 'y', 'z', 'm']
    avaProfileMass = {'v2': np.empty((0, 1)), 'ekin': np.empty((0, 1))}
    for prop in proList:
        avaProfileMass[prop] = np.empty((0, 1))
        avaProfileMass[prop + 'std'] = np.empty((0, 1))

    # loop on each particle dictionary (ie each time step saved)
    for field in fieldsList:
        nonZeroIndex = np.where(field['FD'] > 0)
        xArray = X[nonZeroIndex]
        yArray = Y[nonZeroIndex]
        zArray, _ = gT.projectOnGrid(xArray, yArray, demRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
        mArray = field['FM'][nonZeroIndex]
        uArray = field['FV'][nonZeroIndex]
        u2Array = uArray*uArray
        kineticEneArray = 0.5*mArray*u2Array
        kineticEne = np.nansum(kineticEneArray)
        particles = {'m': mArray, 'x': xArray, 'y': yArray, 'z': zArray}

        # mass-averaged path
        avaProfileMass = appendAverageStd(proList, avaProfileMass, particles, u2Array, kineticEne, mArray)

    avaProfileMass['x'] = avaProfileMass['x']
    avaProfileMass['y'] = avaProfileMass['y']

    # compute s
    avaProfileMass = gT.computeS(avaProfileMass)
    return avaProfileMass


def extendDFAPath(avalancheDir, cfg, dem, simName, avaProfile):
    """ extend the DFA path at the top and bottom
    Only returns the x, y, z path (if an s or sCor was provided..., it is removed)

    Parameters
    -----------
    dem: dict
        dem dict
    particlesIni: dict
        initial particles dict
    avaProfileExt: dict
        profile

    Returns
    --------
    avaProfile: dict
        extended profile
    """
    if cfg.getboolean('pathFromPart'):
        # read particles
        particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName=simName, flagAvaDir=True,
                                                                       comModule='com1DFA')
        particlesIni = particlesList[0]
    else:
        # read field
        fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['FD', 'FV', 'FM'], simName=simName,
                                                               flagAvaDir=True, comModule='com1DFA', timeStep=0)
        # get fields header
        ncols = fieldHeader['ncols']
        nrows = fieldHeader['nrows']
        xllc = fieldHeader['xllcenter']
        yllc = fieldHeader['yllcenter']
        csz = fieldHeader['cellsize']
        X, Y = gT.makeCoordinateGrid(xllc, yllc, csz, ncols, nrows)
        particlesIni = {'x': np.flatten(X), 'y':np.flatten(Y)}
        particlesIni = gT.projectOnRaster(dem, particlesIni)

    avaProfileExt = {key: avaProfile[key] for key in ['x', 'y', 'z']}
    # resample the profile
    resampleDistance = cfg.getfloat('nCellsResample') * dem['header']['cellsize']
    avaProfileExt, _ =  gT.prepareLine(dem, avaProfileExt, distance=resampleDistance, Point=None)
    avaProfileExt = extendProfileTop(cfg.getint('extTopOption'), particlesIni, avaProfileExt)
    avaProfileExt = extendProfileBottom(cfg, dem, avaProfileExt)
    avaProfileExt, _ =  gT.prepareLine(dem, avaProfileExt, distance=resampleDistance, Point=None)
    # remove points that lay outside of the dem# project the profile on the dem
    avaProfileExt, _ =  gT.projectOnRaster(dem, avaProfileExt, interp='bilinear')
    isNotNan = ~(np.isnan(avaProfileExt['z']))
    avaProfileExt['x'] = avaProfileExt['x'][isNotNan]
    avaProfileExt['y'] = avaProfileExt['y'][isNotNan]
    avaProfileExt['z'] = avaProfileExt['z'][isNotNan]
    return avaProfileExt


def extendProfileTop(extTopOption, particlesIni, profile):
    """ extend the DFA path at the top

    Find the direction in which to extend considering the first point of the profile
    and a few following ones (distFromFirt <= 30 * csz). Extend in this direction until
    the z of the highest particle in the release is reached.

    Parameters
    -----------
    extTopOption: boolean
        decid how to extend towards the top
    particlesIni: dict
        initial particles dict
    profile: dict
        profile to extend

    Returns
    --------
    profile: dict
        extended profile
    """
    if extTopOption == 0:
        # get highest particle
        indHighest = np.argmax(particlesIni['z'])
        xHighest = particlesIni['x'][indHighest]
        yHighest = particlesIni['y'][indHighest]
        zHighest = particlesIni['z'][indHighest]
        # get highest point
        xExtTop = xHighest
        yExtTop = yHighest
        zExtTop = zHighest
    elif extTopOption == 1:
        # get point with the most important runout gain
        # get first particle of the path
        xFirst = profile['x'][0]
        yFirst = profile['y'][0]
        zFirst = profile['z'][0]
        # get last particle of the path
        sLast = profile['s'][-1]
        zLast = profile['z'][-1]
        # compute runout angle for averaged path
        tanAngle = (zFirst-zLast)/sLast
        # compute ds
        dx = particlesIni['x'] - xFirst
        dy = particlesIni['y'] - yFirst
        ds = np.sqrt(dx**2 + dy**2)
        # compute dz
        dz = particlesIni['z'] - zFirst
        # remove the elevation needed to match the runout angle
        dz1 = dz - tanAngle * ds
        # get the particle with the highest potential
        indTop = np.argmax(dz1)
        xExtTop = particlesIni['x'][indTop]
        yExtTop = particlesIni['y'][indTop]
        zExtTop = particlesIni['z'][indTop]

    # extend profile
    profile['x'] = np.append(xExtTop, profile['x'])
    profile['y'] = np.append(yExtTop, profile['y'])
    profile['z'] = np.append(zExtTop, profile['z'])
    if debugPlot:
        debPlot.plotPathExtTop(profile, particlesIni, xFirst, yFirst, zFirst, dz1)
    return profile


def extendProfileBottom(cfg, dem, profile):
    """ extend the DFA path at the bottom

    Find the direction in which to extend considering the last point of the profile
    and a few previous ones but discarding the ones that are too close ( 2* csz < distFromLast <= 30 * csz).
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
    extendMinDistance = cfg.getfloat('nCellsMinExtend') * csz
    extendMaxDistance = cfg.getfloat('nCellsMaxExtend') * csz
    pointsOfInterestLast = np.where((r < extendMaxDistance) & (r > extendMinDistance))[0]
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
    gamma = cfg.getfloat('factBottomExt') * sLast / np.sqrt(vDirX**2 + vDirY**2)
    xExtBottom = np.array([xLast + gamma * vDirX])
    yExtBottom = np.array([yLast + gamma * vDirY])
    # project on cszDEM
    zExtBottom, _ =  gT.projectOnGrid(xExtBottom, yExtBottom, zRaster, csz=csz)
    # extend profile
    profile['x'] = np.append(profile['x'], xExtBottom)
    profile['y'] = np.append(profile['y'], yExtBottom)
    profile['z'] = np.append(profile['z'], zExtBottom)

    if debugPlot:
        debPlot.plotPathExtBot(profile, xInterest, yInterest, zInterest, xLast, yLast)
    return profile


def weightedAvgAndStd(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


def appendAverageStd(proList, avaProfile, particles, U2, kineticEneSum, weights):
    """ append averaged to path

    Parameters
    -----------
    proList: list
        list of properties to average and append
    avaProfile: dict
        path
    particles: dict
        particles dict
    U2: numpy array
        array of particles squared velocities (same size as particles)
    kineticEneSum: numpy array
        array of particles kinetic energy (same size as particles)
    weights: numpy array
        array of weights (same size as particles)

    Returns
    --------
    avaProfile: dict
        averaged profile
    """
    for prop in proList:
        avg, std = weightedAvgAndStd(particles[prop], weights)
        avaProfile[prop] = np.append(avaProfile[prop], avg)
        avaProfile[prop + 'std'] = np.append(avaProfile[prop + 'std'], std)
    avaProfile['v2'] = np.append(avaProfile['v2'], np.average(U2, weights=weights))
    avaProfile['ekin'] = np.append(avaProfile['ekin'], kineticEneSum)
    return avaProfile
