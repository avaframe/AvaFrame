"""
    Tools specific to the com1DFA computational kernel
"""

# Load modules
import logging
import math
import pathlib
import numpy as np

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls
from avaframe.in3Utils import geoTrans
import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils
import avaframe.out3Plot.outDebugPlots as debPlot

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')


def getPartInitMethod(cfg, csz, relThForPart):
    """ Get particle initialization parameters

    Get the massPerPart and nPPK corresponding to the desired initialization method

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    csz: float
        cell size
    relThForPart: float
        relTh value

    Returns
    -------
    massPerPart : float
        mass per particle desired
    nPPK : int
        number of particles per kernel radius desired
    """
    rho = cfg.getfloat('rho')
    massPerParticleDeterminationMethod = cfg['massPerParticleDeterminationMethod']
    nPPK = 0
    # derive mass per particle to define number of particles per cell:
    if massPerParticleDeterminationMethod == 'MPPDIR':
        massPerPart = cfg.getfloat('massPerPart')
        log.debug('Number of particles defined by: mass per particle %s' % cfg['massPerPart'])
    elif massPerParticleDeterminationMethod == 'MPPDH':
        deltaTh = cfg.getfloat('deltaTh')
        ds = min(csz, cfg.getfloat('sphKernelRadius'))
        massPerPart = rho * ds * ds * deltaTh
        log.debug('Number of particles defined by: release thickness per particle: %s' % cfg['deltaTh'])
        log.debug('mass per particle is %.2f' % massPerPart)
    elif massPerParticleDeterminationMethod == 'MPPKR':
        sphKernelRadius = cfg.getfloat('sphKernelRadius')
        cszMin = min(csz, sphKernelRadius)
        nPPK0 = cfg.getfloat('nPPK0')
        sphKR0 = cfg.getfloat('sphKR0')
        aPPK = cfg.getfloat('aPPK')
        nPPK = round(nPPK0 * (cszMin/sphKR0)**aPPK)
        massPerPart = rho * math.pi * cszMin * cszMin * relThForPart / nPPK

    return massPerPart, nPPK


def extendCom1DFAPath(cfg, dem, particlesIni, avaProfile):
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
    avaProfileExt = {key: avaProfile[key] for key in ['x', 'y', 'z']}
    # resample the profile
    resampleDistance = cfg.getfloat('PATH', 'nCellsResample') * dem['header']['cellsize']
    avaProfileExt, _ = geoTrans.prepareLine(dem, avaProfileExt, distance=resampleDistance, Point=None)
    avaProfileExt = extendProfileTop(cfg, dem, particlesIni, avaProfileExt)
    avaProfileExt = extendProfileBottom(cfg, dem, avaProfileExt)
    avaProfileExt, _ = geoTrans.prepareLine(dem, avaProfileExt, distance=resampleDistance, Point=None)
    # remove points that lay outside of the dem# project the profile on the dem
    avaProfileExt, _ = geoTrans.projectOnRaster(dem, avaProfileExt, interp='bilinear')
    isNotNan = ~(np.isnan(avaProfileExt['z']))
    avaProfileExt['x'] = avaProfileExt['x'][isNotNan]
    avaProfileExt['y'] = avaProfileExt['y'][isNotNan]
    avaProfileExt['z'] = avaProfileExt['z'][isNotNan]
    return avaProfileExt


def extendProfileTop(cfg, dem, particlesIni, profile):
    """ extend the DFA path at the top

    Find the direction in which to extend considering the first point of the profile
    and a few following ones (distFromFirt <= 30 * csz). Extend in this direction until
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
    if cfg.getint('PATH', 'extTopOption') == 1:
        # compute distance from first point:
        r = DFAtls.norm(profile['x']-xFirst, profile['y']-yFirst, profile['z']-zFirst)
        # find the following first points
        extendMinDistance = cfg.getfloat('PATH', 'nCellsMinExtend') * csz
        extendMaxDistance = cfg.getfloat('PATH', 'nCellsMaxExtend') * csz
        pointsOfInterestFirst = np.where((r < extendMaxDistance) & (r > extendMinDistance))[0]
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
        xExtTop, yExtTop, zExtTop = geoTrans.findPointOnDEM(dem, vDirX, vDirY, vDirZ, zHighest, xFirst, yFirst, zFirst)
        zExtTop, _ = geoTrans.projectOnGrid(xExtTop, yExtTop, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    else:
        xExtTop = xHighest
        yExtTop = yHighest
        zExtTop = zHighest
    # extend profile
    profile['x'] = np.append(xExtTop, profile['x'])
    profile['y'] = np.append(yExtTop, profile['y'])
    profile['z'] = np.append(zExtTop, profile['z'])
    if debugPlot:
        debPlot.plotPathExtTop(cfg, profile, xHighest, yHighest, zHighest, xInterest, yInterest, zInterest,
                                     particlesIni, xllc, yllc, xFirst, yFirst, zFirst)
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
    extendMinDistance = cfg.getfloat('PATH', 'nCellsMinExtend') * csz
    extendMaxDistance = cfg.getfloat('PATH', 'nCellsMaxExtend') * csz
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
    gamma = 0.2 * sLast / np.sqrt(vDirX**2 + vDirY**2)
    xExtBottom = np.array([xLast + gamma * vDirX])
    yExtBottom = np.array([yLast + gamma * vDirY])
    # project on cszDEM
    zExtBottom, _ = geoTrans.projectOnGrid(xExtBottom, yExtBottom, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    # extend profile
    profile['x'] = np.append(profile['x'], xExtBottom)
    profile['y'] = np.append(profile['y'], yExtBottom)
    profile['z'] = np.append(profile['z'], zExtBottom)

    if debugPlot:
        debPlot.plotPathExtBot(cfg, profile, xInterest, yInterest, zInterest, xllc, yllc, xLast, yLast)
    return profile
