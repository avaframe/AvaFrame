"""
    Tools for generating an avalanche path from a DFA simulation
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


def generateAveragePath(avalancheDir, pathFromPart, simName, dem, addVelocityInfo=False, flagAvaDir=True,
                        comModule='com1DFA'):
    """ extract path from fileds or particles

    Parameters
    -----------
    avalancheDir: pathlib
        avalanche directory pathlib path
    pathFromPart: boolean
        compute path from particles if True, from fields (FT, FM, FV) if False
    simName: str
        simulation name
    dem: dict
        com1DFA simulation dictionary
    addVelocityInfo: boolean
        True to add (u2, ekin, totEKin) to result
        Will only work if the particles (ux, uy, uz) exist or if 'FV' exists
    flagAvaDir: bool
        if True avalancheDir corresponds to an avalanche directory and data is
        read from avaDir/Outputs/comModule/particles or avaDir/Outputs/comModule/peakFiles
        depending on if pathFromPart is True or False
    comModule: str
        module that computed the particles or fields

    Returns
    --------
    avaProfileMass: dict
        mass averaged profile (x, y, z, 's')
        if addVelocityInfo is True, kinetic energy and velocity information are added to
        the avaProfileMass dict (u2, ekin, totEKin)
    """
    if pathFromPart:
        particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName=simName, flagAvaDir=True,
                                                                       comModule='com1DFA')
        log.info('Using particles to generate avalanche path profile')
        # postprocess to extract path and energy line
        avaProfileMass = getDFAPathFromPart(particlesList, addVelocityInfo=addVelocityInfo)
    else:
        particlesList = ''
        # read field
        fieldName = ['FT', 'FM']
        if addVelocityInfo:
            fieldName.append['FV']
        fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, fieldName, simName=simName,
                                                               flagAvaDir=True, comModule='com1DFA')
        log.info('Using fields to generate avalanche path profile')
        # postprocess to extract path and energy line
        avaProfileMass = getDFAPathFromField(fieldsList, fieldHeader, dem)

    return avaProfileMass


def getDFAPathFromPart(particlesList, addVelocityInfo=False):
    """ compute mass averaged path from particles

    Also returns the averaged velocity and kinetic energy associated
    If addVelocityInfo is True, information about velocity and kinetic energy is computed

    Parameters
    -----------
    particlesList: list
        list of particles dict
    addVelocityInfo: boolean
        True to add (u2, ekin, totEKin) to result

    Returns
    --------
    avaProfileMass: dict
        mass averaged profile (x, y, z, 's', 'sCor')
        if addVelocityInfo is True, kinetic energy and velocity information are added to
        the avaProfileMass dict (u2, ekin, totEKin)
    """

    propList = ['x', 'y', 'z', 's', 'sCor']
    avaProfileMass = {}
    # do we have velocity info?
    if addVelocityInfo:
        propList.append('u2')
        propList.append('ekin')
        avaProfileMass['totEKin'] = np.empty((0, 1))
    # initialize other properties
    for prop in propList:
        avaProfileMass[prop] = np.empty((0, 1))
        avaProfileMass[prop + 'std'] = np.empty((0, 1))

    # loop on each particle dictionary (ie each time step saved)
    for particles in particlesList:
        if particles['nPart'] > 0:
            m = particles['m']
            if addVelocityInfo:
                ux = particles['ux']
                uy = particles['uy']
                uz = particles['uz']
                u = DFAtls.norm(ux, uy, uz)
                u2Array = u*u
                kineticEneArray = 0.5*m*u2Array
                particles['u2'] = u2Array
                particles['ekin'] = kineticEneArray

            # mass-averaged path
            avaProfileMass = appendAverageStd(propList, avaProfileMass, particles, m)

            if addVelocityInfo:
                avaProfileMass['totEKin'] = np.append(avaProfileMass['totEKin'], np.nansum(kineticEneArray))

    return avaProfileMass


def getDFAPathFromField(fieldsList, fieldHeader, dem):
    """ compute mass averaged path from fields

    Also returns the averaged velocity and kinetic energy associated
    The dem and fieldsList (FT, FV and FM) need to have identical dimentions and cell size.
    If FV is not provided, information about velocity and kinetic energy is not computed

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
        mass averaged profile (x, y, z, 's')
        if 'FV' in fieldsList, kinetic energy and velocity information are added to
        the avaProfileMass dict (u2, ekin, totEKin)
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

    propList = ['x', 'y', 'z']
    avaProfileMass = {}
    # do we have velocity info?
    addVelocityInfo = False
    if 'FV' in fieldsList[0]:
        propList.append('u2')
        propList.append('ekin')
        avaProfileMass['totEKin'] = np.empty((0, 1))
        addVelocityInfo = True
    # initialize other properties
    for prop in propList:
        avaProfileMass[prop] = np.empty((0, 1))
        avaProfileMass[prop + 'std'] = np.empty((0, 1))
    # loop on each field dictionary (ie each time step saved)
    for field in fieldsList:
        # find cells with snow
        nonZeroIndex = np.where(field['FT'] > 0)
        xArray = X[nonZeroIndex]
        yArray = Y[nonZeroIndex]
        zArray, _ = gT.projectOnGrid(xArray, yArray, demRaster, csz=csz, xllc=xllc, yllc=yllc)
        mArray = field['FM'][nonZeroIndex]
        particles = {'x': xArray, 'y': yArray, 'z': zArray}
        if addVelocityInfo:
            uArray = field['FV'][nonZeroIndex]
            u2Array = uArray*uArray
            kineticEneArray = 0.5*mArray*u2Array
            particles['u2'] = u2Array
            particles['ekin'] = kineticEneArray

        # mass-averaged path
        avaProfileMass = appendAverageStd(propList, avaProfileMass, particles, mArray)

        if addVelocityInfo:
            avaProfileMass['totEKin'] = np.append(avaProfileMass['totEKin'], np.nansum(kineticEneArray))

    avaProfileMass['x'] = avaProfileMass['x'] - xllc
    avaProfileMass['y'] = avaProfileMass['y'] - yllc

    # compute s
    avaProfileMass = gT.computeS(avaProfileMass)
    return avaProfileMass


def extendDFAPath(avalancheDir, cfg, dem, simName, avaProfile):
    """ extend the DFA path at the top (release) and bottom (runout area)
    cfg['pathFromPart'] decides if the path is extended from particles dict or fields
    Call extendDFAPathKernel for extending the path avaProfile with x, y, z, s information

    Parameters
    -----------
    avalancheDir: str or pathlib path
        path to avalanche directory
    cfg: configParser
        configuration object with:
        pathFromPart: boolean
            read information from particles file or fields file
        extTopOption: int
        how to extend towards the top?
            0 for highest point method
            a for largest runout method
        nCellsResample: int
            resampling length is given by nCellsResample*demCellSize
        nCellsMinExtend: int
            when extending towards the bottom, take points at more than nCellsMinExtend*demCellSize
            from last point to get the direction
        nCellsMaxExtend: int
            when extending towards the bottom, take points at less than nCellsMaxExtend*demCellSize
            from last point to get the direction
        factBottomExt: float
            extend the profile from factBottomExt*sMax
    dem: dict
        dem dict
    simName: str
        simHash or name
    avaProfile: dict
        profile to be extended
    Returns
    --------
    avaProfileExt: dict
        extended profile at top and bottom (x, y, z).
    """
    # read inputs from particles or fields
    if cfg.getboolean('pathFromPart'):
        # read particles
        particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName=simName, flagAvaDir=True,
                                                                       comModule='com1DFA')
        particlesIni = particlesList[0]
        log.info('Using particles to generate avalanche path profile')
    else:
        # read field
        fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['FT'], simName=simName,
                                                               flagAvaDir=True, comModule='com1DFA', timeStep=0)
        # get fields header
        ncols = fieldHeader['ncols']
        nrows = fieldHeader['nrows']
        csz = fieldHeader['cellsize']
        # we want the origin to be in (0, 0) as it is in the avaProfile that comes in
        X, Y = gT.makeCoordinateGrid(0, 0, csz, ncols, nrows)
        indNonZero = np.where(fieldsList[0]['FT'] > 0)
        # convert this data in a particles style (dict with x, y, z info)
        particlesIni = {'x': X[indNonZero], 'y': Y[indNonZero]}
        particlesIni, _ = gT.projectOnRaster(dem, particlesIni)
        log.info('Using fields to generate avalanche path profile')

    avaProfileExt = extendDFAPathKernel(cfg, avaProfile, dem, particlesIni)
    return avaProfileExt


def extendDFAPathKernel(cfg, avaProfile, dem, particlesIni):
    """ extend the DFA path at the top and bottom
    avaProfile with x, y, z, s information

    Parameters
    -----------
    cfg: configParser
        configuration object with:
        extTopOption: int
        how to extend towards the top?
            0 for heighst point method
            a for largest runout method
        nCellsResample: int
            resampling length is given by nCellsResample*demCellSize
        nCellsMinExtend: int
            when extending towards the bottom, take points at more than nCellsMinExtend*demCellSize
            from last point to get the direction
        nCellsMaxExtend: int
            when extending towards the bottom, take points at less than nCellsMaxExtend*demCellSize
            from last point to get the direction
        factBottomExt: float
            extend the profile from factBottomExt*sMax
    avaProfile: dict
        profile to be extended
    dem: dict
        dem dict
    particlesIni: dict
        initial particles dict

    Returns
    --------
    avaProfileExt: dict
        extended profile at top and bottom (x, y, z).
    """
    # resample the profile
    resampleDistance = cfg.getfloat('nCellsResample') * dem['header']['cellsize']
    avaProfile = extendProfileTop(cfg.getint('extTopOption'), particlesIni, avaProfile)
    avaProfile, _ =  gT.prepareLine(dem, avaProfile, distance=resampleDistance, Point=None)
    avaProfile = extendProfileBottom(cfg, dem, avaProfile)
    avaProfile, _ =  gT.prepareLine(dem, avaProfile, distance=resampleDistance, Point=None)
    # project the profile on the dem
    avaProfile, _ =  gT.projectOnRaster(dem, avaProfile, interp='bilinear')
    # remove points that lay outside of the dem
    isNotNan = ~(np.isnan(avaProfile['z']))
    avaProfileExt = {}
    avaProfileExt['x'] = avaProfile['x'][isNotNan]
    avaProfileExt['y'] = avaProfile['y'][isNotNan]
    avaProfileExt['z'] = avaProfile['z'][isNotNan]
    return avaProfileExt


def extendProfileTop(extTopOption, particlesIni, profile):
    """ extend the DFA path at the top (release)

    Either towards the highest point in particlesIni (extTopOption = 0)
    or the point leading to the longest runout (extTopOption = 1)

    Parameters
    -----------
    extTopOption: int
        decide how to extend towards the top
        if 0, extrapolate towards the highest point in the release
        if 1, extrapolate towards the point leading to the lonest runout
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
        xExtTop = particlesIni['x'][indHighest]
        yExtTop = particlesIni['y'][indHighest]
        zExtTop = particlesIni['z'][indHighest]
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
    """ extend the DFA path at the bottom (runout area)

    Find the direction in which to extend considering the last point of the profile
    and a few previous ones but discarding the ones that are too close
    (nCellsMinExtend* csz < distFromLast <= nCellsMaxExtend * csz).
    Extend in this diretion for a distance factBottomExt * length of the path.

    Parameters
    -----------
    cfg: configParser
        nCellsMinExtend: int
            when extending towards the bottom, take points at more than nCellsMinExtend*demCellSize
            from last point to get the direction
        nCellsMaxExtend: int
            when extending towards the bottom, take points at less than nCellsMaxExtend*demCellSize
            from last point to get the direction
        factBottomExt: float
            extend the profile from factBottomExt*sMax
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
    sLast = profile['s'][-1]
    # compute distance from last point:
    r = DFAtls.norm(profile['x']-xLast, profile['y']-yLast, 0)
    # find the previous points
    extendMinDistance = cfg.getfloat('nCellsMinExtend') * csz
    extendMaxDistance = cfg.getfloat('nCellsMaxExtend') * csz
    pointsOfInterestLast = np.where((r < extendMaxDistance) & (r > extendMinDistance))[0]
    xInterest = profile['x'][pointsOfInterestLast]
    yInterest = profile['y'][pointsOfInterestLast]
    # find the direction in which we need to extend the path
    vDirX = xLast - xInterest
    vDirY = yLast - yInterest
    vDirX, vDirY, vDirZ = DFAtls.normalize(np.array([vDirX]), np.array([vDirY]), 0*np.array([vDirY]))
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
        debPlot.plotPathExtBot(profile, xInterest, yInterest, 0*yInterest, xLast, yLast)
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


def appendAverageStd(propList, avaProfile, particles, weights):
    """ append averaged to path

    Parameters
    -----------
    propList: list
        list of properties to average and append
    avaProfile: dict
        path
    particles: dict
        particles dict
    weights: numpy array
        array of weights (same size as particles)

    Returns
    --------
    avaProfile: dict
        averaged profile
    """
    for prop in propList:
        avg, std = weightedAvgAndStd(particles[prop], weights)
        avaProfile[prop] = np.append(avaProfile[prop], avg)
        avaProfile[prop + 'std'] = np.append(avaProfile[prop + 'std'], std)
    return avaProfile
