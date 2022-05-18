"""
  Create generic/idealised topographies
"""


# load modules
import logging
import numpy as np
from scipy.stats import norm
from scipy.interpolate import griddata
import pathlib

from avaframe.in3Utils import geoTrans

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getParabolaParams(cfg):
    """ Compute parameters for parabola """

    # input parameters
    C = float(cfg['TOPO']['C'])
    fLens = float(cfg['TOPO']['fLens'])
    meanAlpha = float(cfg['TOPO']['meanAlpha'])

    # If mean slope is given or distance to the start of the flat plane
    if meanAlpha != 0:
        fLen = C / np.tan(np.radians(meanAlpha))
        log.info('fLen computed from mean alpha: %.2f meters' % fLen)
    else:
        fLen = fLens
        log.info('flen directly set to: %.2f meters' % fLen)
    A = C / (fLen**2)
    B = (-C * 2.) / fLen

    return A, B, fLen


def getGridDefs(cfg):
    # determine number of rows and columns to define domain
    dx = float(cfg['TOPO']['dx'])
    xEnd = float(cfg['TOPO']['xEnd'])
    yEnd = float(cfg['TOPO']['yEnd'])

    return dx, xEnd, yEnd


def computeCoordGrid(dx, xEnd, yEnd):

    # Compute coordinate grid
    xv = np.arange(0, xEnd+dx, dx)
    yv = np.arange(-0.5 * yEnd, 0.5 * (yEnd+dx), dx)
    nRows = len(yv)
    nCols = len(xv)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nRows, nCols))

    return xv, yv, zv, x, y, nRows, nCols


def flatplane(cfg):
    """ Compute coordinates of flat plane topography """

    dx, xEnd, yEnd = getGridDefs(cfg)

    zElev = float(cfg['TOPO']['zElev'])

    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # Set elevation of surface
    zv = zv + zElev

    # Log info here
    log.info('Flatplane coordinates computed')

    return x, y, zv


def inclinedplane(cfg):
    """ Compute coordinates of inclined plane with given slope (meanAlpha)"""

    # input parameters
    dx, xEnd, yEnd = getGridDefs(cfg)

    z0 = float(cfg['TOPO']['z0'])
    meanAlpha = float(cfg['TOPO']['meanAlpha'])

    cFf = float(cfg['CHANNELS']['cff'])
    cRadius = float(cfg['CHANNELS']['cRadius'])

    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # Set surface elevation from slope and max. elevation
    zv = z0 - np.tan(np.radians(meanAlpha)) * x

    # If a channel shall be introduced
    if cfg['TOPO'].getboolean('channel'):
        # Compute cumulative distribution function and set horizontal extent of channel
        c0 = norm.cdf(xv, 0, cFf)
        cExtent = cRadius
        yv = np.reshape(yv, (nRows , 1))

        # if location within horizontal extent of channel,
        # make half sphere shaped channel with radius given by channel horizontal extent
        mask = np.zeros(np.shape(yv))
        mask[np.where(abs(yv) < cExtent)] = 1
        if cfg['TOPO'].getboolean('topoAdd'):
            zv = zv + cExtent*c0*(1. - np.sqrt(np.abs(1. - (np.square(yv) / (cExtent**2)))))*mask
            mask = np.zeros(np.shape(yv))
            mask[np.where(abs(yv) >= cExtent)] = 1
            zv = zv + cExtent*c0*mask
        else:
            zv = zv - cExtent*c0*np.sqrt(np.abs(1. - (np.square(yv) / (cExtent**2))))*mask

    # Log info here
    log.info('Inclined plane coordinates computed')

    return x, y, zv


def addDrop(cfg, x, y, zv):
    """ Add a drop to a given topography

    The drop is added in the x drection

    Parameters
    ----------
    cfg: configparser
        configuration for generateTopo
    x: 2D numpy array
        x coordinate of the raster
    y: 2D numpy array
        y coordinate of the raster
    zv: 2D numpy array
        z coordinate of the raster

    Returns
    -------
    zv: 2D numpy array
        z coordinate of the raster taking the drop into account
    """

    # input parameters
    dx, xEnd, yEnd = getGridDefs(cfg)

    xStartDrop = float(cfg['DROP']['xStartDrop'])
    dxDrop = float(cfg['DROP']['dxDrop'])
    alphaDrop = float(cfg['DROP']['alphaDrop'])

    # get zcoord
    # deduce drop height from the drop steepness and length in x direction
    dzDrop = dxDrop * np.tan(np.radians(alphaDrop))
    xEndDrop = xStartDrop + dxDrop

    nrows, ncols = np.shape(x)
    # get the z coordinate of the point at the beginning of the drop
    zIniDrop, _ = geoTrans.projectOnGrid(xStartDrop*np.ones((nrows)), y[:, 0],
                                         np.vstack((zv[0,:], zv)), csz=dx,
                                         xllc=x[0, 0], yllc=y[0, 0], interp='bilinear')
    zIniDrop = np.tile(zIniDrop, (ncols, 1)).transpose()
    # get the z coordinate of the point at the end of the drop
    zEndDrop, _ = geoTrans.projectOnGrid(xEndDrop*np.ones((nrows)), y[:, 0],
                                         np.vstack((zv[0,:], zv)), csz=dx,
                                         xllc=x[0, 0], yllc=y[0, 0], interp='bilinear')
    zEndDrop = np.tile(zEndDrop, (ncols, 1)).transpose()
    # Set surface elevation from slope and max. elevation
    zv = np.where(((x >= xStartDrop) & (x <= xEndDrop)),
                  zIniDrop - (x - xStartDrop)*np.tan(np.radians(alphaDrop)),
                  zv)
    zv = np.where(x > xEndDrop, zv - (dzDrop + zEndDrop - zIniDrop), zv)

    # Log info here
    log.info('Added drop to the topography')

    return zv


def hockey(cfg):
    """
        Compute coordinates of an inclined plane with a flat foreland  defined by
        total fall height z0, angle to flat foreland (meanAlpha) and a radius (rCirc) to
        smooth the transition from inclined plane to flat foreland
    """

    # input parameters
    rCirc = float(cfg['TOPO']['rCirc'])
    meanAlpha = float(cfg['TOPO']['meanAlpha'])
    z0 = float(cfg['TOPO']['z0'])

    cff = float(cfg['CHANNELS']['cff'])
    cRadius = float(cfg['CHANNELS']['cRadius'])
    cInit = float(cfg['CHANNELS']['cInit'])
    cMustart = float(cfg['CHANNELS']['cMustart'])
    cMuendFP = float(cfg['CHANNELS']['cMuendFP'])

    dx, xEnd, yEnd = getGridDefs(cfg)

    # Compute coordinate grid
    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # Compute distance to flat foreland for given meanAlpha
    x1 = z0 / np.tan(np.radians(meanAlpha))
    if x1 >= xEnd * 0.9:
        log.warning('Your domain (xEnd) is to small or the slope angle (meanAlpha) to'
                    'shallow to produce a signifcant (>10 percent of domain, in your case:'
                    ' %.2f m) flat foreland!' % (0.1 * (xEnd - dx)))

    # Compute circle parameters for smoothing the transition
    beta = (0.5 * (180. - (meanAlpha)))
    xc = rCirc / np.tan(np.radians(beta))
    yc = xc * np.cos(np.radians(meanAlpha))
    xCirc = x1 + xc
    # for plotting
    d1 = np.tan(np.radians(beta)) * x1

    # Set surface elevation
    zv = np.zeros((nRows, nCols))
    mask = np.zeros(np.shape(x))
    mask[np.where(x < (x1 - yc))] = 1
    zv = zv + (z0 - np.tan(np.radians(meanAlpha)) * x)*mask

    mask = np.zeros(np.shape(x))
    mask[np.where(((x1 - yc) <= x) & (x <= (x1 + xc)))] = 1
    # rCirc + np.sqrt(rCirc**2 - (xv[m] - xCirc)**2)
    zv = zv + (rCirc - np.sqrt(np.abs(rCirc**2 - (xCirc - x)**2)))*mask

    # If a channel shall be introduced
    if cfg['TOPO'].getboolean('channel'):
        # Compute cumulative distribution function - c1 for upper part (start)
        # of channel and c2 for lower part (end) of channel
        c1 = norm.cdf(xv, cMustart * (x1), cff)
        c2 = 1. - norm.cdf(xv, cMuendFP * (x1), cff)

        # combine both into one function separated at the the middle of
        #  the channel longprofile location
        mask = np.zeros(np.shape(xv))
        mask[np.where(xv < (x1 * (0.5 * (cMustart + cMuendFP))))] = 1
        c0 = c1 * mask

        mask = np.zeros(np.shape(xv))
        mask[np.where(xv >= (x1 * (0.5 * (cMustart + cMuendFP))))] = 1
        c0 = c0 + c2 * mask

        # Is the channel of constant width or narrowing
        if cfg['TOPO'].getboolean('narrowing'):
            cExtent = cInit * (1 - c0[:]) + (c0[:] * cRadius)
        else:
            cExtent = np.zeros(np.shape(xv)) + cRadius

        # Set surface elevation
        mask = np.zeros(np.shape(y))
        mask[np.where(abs(y) < cExtent)] = 1
        # Add surface elevation modification introduced by channel
        if cfg['TOPO'].getboolean('topoAdd'):
            zv = zv + cExtent*c0*(1. - np.sqrt(np.abs(1. - (np.square(y) / (cExtent**2)))))*mask
            # outside of the channel, add layer of channel thickness
            mask = np.zeros(np.shape(y))
            mask[np.where(abs(y) >= cExtent)] = 1
            zv = zv + cExtent*c0*mask
        else:
            zv = zv - cExtent*c0*np.sqrt(np.abs(1. - (np.square(y) / (cExtent**2))))*mask

    # Log info here
    log.info('Hockeystick coordinates computed')

    return x, y, zv


def parabola(cfg):
    """
        Compute coordinates of a parabolically-shaped slope with a flat foreland
        defined by total fall height C, angle (meanAlpha) or distance (fLen) to flat foreland
    """

    C = float(cfg['TOPO']['C'])
    cff = float(cfg['CHANNELS']['cff'])
    cRadius = float(cfg['CHANNELS']['cRadius'])
    cInit = float(cfg['CHANNELS']['cInit'])
    cMustart = float(cfg['CHANNELS']['cMustart'])
    cMuend = float(cfg['CHANNELS']['cMuend'])

    # Get grid definitons
    dx, xEnd, yEnd = getGridDefs(cfg)

    # Compute coordinate grid
    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # Get parabola Parameters
    [A, B, fLen] = getParabolaParams(cfg)

    # Set surface elevation
    zv = np.ones((nRows, nCols))
    # initialize superimposed channel
    superChannel = np.zeros(np.shape(xv))
    superDam = np.zeros(np.shape(xv))

    zv = zv * ((-B**2) / (4. * A) + C)
    mask = np.zeros(np.shape(xv))
    mask[np.where(xv < fLen)] = 1

    zv = zv + (A * xv**2 + B * xv + C)*mask

    # If a channel shall be introduced
    if cfg['TOPO'].getboolean('channel'):
        c1 = norm.cdf(xv, cMustart * fLen, cff)
        c2 = 1. - norm.cdf(xv, cMuend * fLen, cff)

        # combine both into one function separated at the the middle of
        #  the channel longprofile location
        mask = np.zeros(np.shape(xv))
        mask[np.where(xv < (fLen * (0.5 * (cMustart + cMuend))))] = 1
        c0 = c1 * mask

        mask = np.zeros(np.shape(xv))
        mask[np.where(xv >= (fLen * (0.5 * (cMustart + cMuend))))] = 1
        c0 = c0 + c2 * mask

        # Is the channel of constant width or narrowing
        if cfg['TOPO'].getboolean('narrowing'):
            cExtent = cInit * (1 - c0[:]) + (c0[:] * cRadius)
        else:
            cExtent = np.zeros(nCols) + cRadius


        # Add surface elevation modification introduced by channel
        mask = np.zeros(np.shape(y))
        mask[np.where(abs(y) < cExtent)] = 1
        if cfg['TOPO'].getboolean('topoAdd'):
            superChannel = superChannel + cExtent*c0*(1. - np.sqrt(np.abs(1. - (np.square(y) / (cExtent**2)))))*mask
            # outside of the channel, add layer of channel thickness
            mask = np.zeros(np.shape(y))
            mask[np.where(abs(y) >= cExtent)] = 1
            superChannel = superChannel + cExtent*c0*mask
        else:
            superChannel = superChannel - cExtent*c0*np.sqrt(np.abs(1. - (np.square(y) / (cExtent**2))))*mask


    if cfg['TOPO'].getboolean('dam'):
        #damPos = float(cfg['DAMS']['damPos'])
        damPos = cfg['DAMS'].getfloat('damPos')
        damHeight = cfg['DAMS'].getfloat('damHeight')
        damWidth = cfg['DAMS'].getfloat('damWidth')
        superDam = norm.pdf(xv, damPos * (-B/2/A) , damWidth)
        superDam = superDam / np.max(superDam) * damHeight
        if not cfg['TOPO'].getboolean('topoAdd'):
            superDam = superDam -  cExtent*c0

    # add channel and dam
    zv = zv + np.maximum( superDam , superChannel)

    # Log info here
    log.info('Parabola with flat outrun coordinates computed')

    return x, y, zv


def bowl(cfg):
    """ Compute coordinates of sphere with given radius (rBwol) """

    # input parameters
    rBwol = float(cfg['TOPO']['rBowl'])

    # Get grid definitions
    dx, xEnd, yEnd = getGridDefs(cfg)

    # Compute coordinate grid
    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # recompute xv yv and x, y as they are shifted
    xv = np.arange(-0.5 * xEnd, 0.5 * (xEnd+dx), dx)
    yv = np.arange(-0.5 * yEnd, 0.5 * (yEnd+dx), dx)
    x, y = np.meshgrid(xv, yv)

    # Set surface elevation
    zv = rBwol*np.ones((nRows, nCols))
    if cfg['TOPO'].getboolean('curvedSlope') == True:
        radius = np.sqrt(x**2)
    else:
        radius = np.sqrt(x**2 + y**2)
    mask = np.zeros(np.shape(x))
    mask[np.where(radius <= rBwol)] = 1
    zv = zv - (rBwol * np.sqrt(np.abs(1 - (radius / rBwol)**2)))*mask
    if cfg['TOPO'].getboolean('curvedSlope') == True:
        zv[x>=0] = 0.

    # Log info here
    log.info('Bowl coordinates computed')

    return x, y, zv


def helix(cfg):
    """ Compute coordinates of helix-shaped topography with given radius (rHelix) """

    # input parameters
    rHelix = float(cfg['TOPO']['rHelix'])
    C = float(cfg['TOPO']['C'])
    cff = float(cfg['CHANNELS']['cff'])
    cRadius = float(cfg['CHANNELS']['cRadius'])
    cInit = float(cfg['CHANNELS']['cInit'])
    cMustart = float(cfg['CHANNELS']['cMustart'])
    cMuend = float(cfg['CHANNELS']['cMuend'])

    # Get grid definitions
    dx, xEnd, yEnd = getGridDefs(cfg)

    # Compute coordinate grid
    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # recompute xv yv and x, y as they are shifted
    xv = np.arange(-0.5 * xEnd, 0.5 * (xEnd+dx), dx)
    yv = np.arange(-yEnd, 0+dx, dx)
    x, y = np.meshgrid(xv, yv)

    # Get parabola Parameters
    [A, B, fLen] = getParabolaParams(cfg)

    # Set surface elevation
    zv = np.ones((nRows, nCols))
    radius = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x) + np.pi
    zv = zv * ((-B**2) / (4. * A) + C)
    mask = np.zeros(np.shape(x))
    mask[np.where((theta * rHelix) < fLen)] = 1

    zv = zv + (A * (theta * rHelix)**2 + B * (theta * rHelix) + C)*mask

    # If channel is introduced to topography
    if cfg['TOPO'].getboolean('channel'):
        c0 = np.zeros(np.shape(x))
        mask = np.zeros(np.shape(x))
        mask[np.where((theta * rHelix) < (0.5 * (cMustart + cMuend) * fLen))] = 1
        c0 = c0 + norm.cdf(theta * rHelix, cMustart * fLen, cff)*mask
        mask = np.zeros(np.shape(x))
        mask[np.where((theta * rHelix) >= (0.5 * (cMustart + cMuend) * fLen))] = 1
        c0 = c0 + (1. - norm.cdf(theta * rHelix, cMuend * fLen, cff))*mask
        # c0 = np.ones(np.shape(zv))
        # If channel of constant width or becoming narrower in the middle
        if cfg['TOPO'].getboolean('narrowing'):
            cExtent = cInit * (1. - c0) + c0 * cRadius
        else:
            cExtent = cRadius

        if cfg['TOPO'].getboolean('topoAdd'):
            zv = zv + c0 * cExtent

        # Inner and outer boundaries of the channel
        boundIn = rHelix - cExtent
        boundExt = rHelix + cExtent

        # Set channel
        mask = np.zeros(np.shape(x))
        mask[np.where((radius >= rHelix) & (radius < boundExt))] = 1
        radius1 = radius - rHelix
        zv = zv - cExtent*c0*np.sqrt(np.abs(1. - (np.square(radius1) / np.square(cExtent))))*mask

        mask = np.zeros(np.shape(x))
        mask[np.where((radius < rHelix) & (radius > boundIn))] = 1
        radius2 = rHelix - radius
        zv = zv - cExtent*c0*np.sqrt(np.abs(1. - (np.square(radius1) / np.square(cExtent))))*mask


    # set last row at Center to fall height
    indCols = int(0.5*nCols)
    zv[-1, 0:indCols] = C

    # Log info here
    log.info('Helix coordinates computed')

    return x, y, zv


def pyramid(cfg):
    """ Generate a pyramid topography - in this case rectangular domain """

    # get parameters from ini
    meanAlpha = float(cfg['TOPO']['meanAlpha'])
    z0 = float(cfg['TOPO']['z0'])
    flatx = float(cfg['TOPO']['flatx'])
    flaty = float(cfg['TOPO']['flaty'])
    phi = float(cfg['TOPO']['phi'])
    dx = float(cfg['TOPO']['dx'])

    # initialise pyramid corners and center point
    points = np.asarray([[-1., -1., 0], [-1., 1., 0], [1., 1., 0], [1., -1, 0.], [0., 0., 1.]])
    dxPoints = abs(points[4,0] - points[1,0])

    # compute elevation of the apex point for given angle of pyramid facets
    zAlpha = dxPoints * np.tan(np.deg2rad(meanAlpha))
    points[4,2] = zAlpha
    dcoors = points * z0 /zAlpha

    # if desired rotate pyramid
    if cfg['TOPO'].getboolean('flagRot') == True:
        dcoorsRot = np.zeros((len(dcoors),3))
        for m in range(len(dcoorsRot)):
            dcoorsRot[m,0] = np.cos(np.deg2rad(phi)) * dcoors[m,0] - np.sin(np.deg2rad(phi)) * dcoors[m,1]
            dcoorsRot[m,1] = np.sin(np.deg2rad(phi)) * dcoors[m,0] + np.cos(np.deg2rad(phi)) * dcoors[m,1]
            dcoorsRot[m,2] = dcoors[m,2]
        dcoorsFin = dcoorsRot
    else:
        dcoorsFin = dcoors

    # split into horizontal and vertical coordinate points
    xyPoints = np.zeros((len(points),2))
    xyPoints[:,0] = dcoorsFin[:,0]
    xyPoints[:,1] = dcoorsFin[:,1]
    zPoints = dcoorsFin[:,2]

    # make meshgrid for final DEM
    xv = np.arange(-flatx+np.amin(dcoorsFin[:,0]), np.amax(dcoorsFin[:,0])+flatx, dx)
    yv = np.arange(-flaty+np.amin(dcoorsFin[:,1]), np.amax(dcoorsFin[:,1])+flaty, dx)
    nRows = len(yv)
    nCols = len(xv)
    x, y = np.meshgrid(xv, yv)

    # interpolate appex point information to meshgrid
    z = griddata(xyPoints, zPoints, (x, y), method='linear')
    zNan = np.isnan(z)
    z[zNan] = 0.0

    dX =  np.amax(dcoorsFin[:,0])+flatx - (-flatx+np.amin(dcoorsFin[:,0]))
    dY =  np.amax(dcoorsFin[:,1])+flaty - (-flaty+np.amin(dcoorsFin[:,1]))
    log.info('domain extent pyramid- inx: %f, in y: %f' % (dX, dY))

    return x, y, z


def writeDEM(cfg, z, outDir):
    """ Write topography information to file """
    nameExt = cfg['TOPO']['demType']
    nRows = z.shape[0]
    nCols = z.shape[1]

    # Read lower left center coordinates, cellsize and noDATA value
    xllcenter = float(cfg['DEMDATA']['xl'])
    yllcenter = float(cfg['DEMDATA']['yl'])
    cellsize = float(cfg['TOPO']['dx'])
    noDATA = float(cfg['DEMDATA']['nodata_value'])
    demName = cfg['DEMDATA']['demName']

    # Save elevation data to .asc file and add header lines
    z_mat = np.matrix(z)
    demFile = outDir / ('%s_%s_Topo.asc' % (demName, nameExt))
    with open(demFile, 'w') as f:
        f.write('nCols  %d\n' % (nCols))
        f.write('nRows  %d\n' % (nRows))
        f.write('xllcenter  %.02f\n' % (xllcenter))
        f.write('yllcenter %.02f\n' % (yllcenter))
        f.write('cellsize  %.02f\n' % (cellsize))
        f.write('nodata_value %.02f\n' % (noDATA))
        for line in z_mat:
            np.savetxt(f, line, fmt='%f')

    # Log info here
    log.info('DEM written to: %s/%s_%s_Topo.asc' % (outDir, demName, nameExt))


def generateTopo(cfg, avalancheDir):
    """ Compute coordinates of desired topography with given inputs """

    # Which DEM type
    demType = cfg['TOPO']['demType']

    log.info('DEM type is set to: %s' % demType)

    # Set Output directory
    outDir = pathlib.Path(avalancheDir, 'Inputs')
    if outDir.is_dir():
        log.info('The new DEM is saved to %s' % (outDir))
    else:
        log.error('Required folder structure: NameOfAvalanche/Inputs missing! \
                    Run runInitializeProject first!')

    # Call topography type
    if demType == 'FP':
        [x, y, z] = flatplane(cfg)

    elif demType == 'IP':
        [x, y, z] = inclinedplane(cfg)

    elif demType == 'PF':
        [x, y, z] = parabola(cfg)

    elif demType == 'HS':
        [x, y, z] = hockey(cfg)

    elif demType == 'BL':
        [x, y, z] = bowl(cfg)

    elif demType == 'HX':
        [x, y, z] = helix(cfg)

    elif demType == 'PY':
        [x, y, z] = pyramid(cfg)


    # If a drop shall be introduced
    if cfg['TOPO'].getboolean('drop'):
        z = addDrop(cfg, x, y, z)

    # Write DEM to file
    writeDEM(cfg, z, outDir)

    return(z, demType, outDir)
