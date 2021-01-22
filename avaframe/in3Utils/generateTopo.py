"""
  Create generic/idealised topographies
"""


# load modules
import logging
import numpy as np
from scipy.stats import norm
from scipy.interpolate import griddata
import os

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

    cFf = float(cfg['CHANNELS']['c_ff'])
    cRadius = float(cfg['CHANNELS']['c_radius'])

    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # Set surface elevation from slope and max. elevation
    zv = z0 - np.tan(np.radians(meanAlpha)) * x

    # If a channel shall be introduced
    if cfg['TOPO'].getboolean('channel'):
        # Compute cumulative distribution function and set horizontal extent of channel
        c_0 = norm.cdf(xv, 0, cFf)
        c_extent = cRadius
        yv = np.reshape(yv, (nRows , 1))

        # if location within horizontal extent of channel,
        # make half sphere shaped channel with radius given by channel horizontal extent
        mask = np.zeros(np.shape(yv))
        mask[np.where(abs(yv) < c_extent)] = 1
        if cfg['TOPO'].getboolean('topoconst'):
            zv = zv - c_extent*c_0*np.sqrt(np.abs(1. - (np.square(yv) / (c_extent**2))))*mask
        else:
            zv = zv + c_extent*c_0*(1. - np.sqrt(np.abs(1. - (np.square(yv) / (c_extent**2)))))*mask
        if not cfg['TOPO'].getboolean('topoconst'):
            mask = np.zeros(np.shape(yv))
            mask[np.where(abs(yv) >= c_extent)] = 1
            zv = zv + c_extent*c_0*mask

    # Log info here
    log.info('Inclined plane coordinates computed')

    return x, y, zv


def hockeysmooth(cfg):
    """
        Compute coordinates of an inclined plane with a flat foreland  defined by
        total fall height z0, angle to flat foreland (meanAlpha) and a radius (rCirc) to
        smooth the transition from inclined plane to flat foreland
    """

    # input parameters
    rCirc = float(cfg['TOPO']['rCirc'])
    meanAlpha = float(cfg['TOPO']['meanAlpha'])
    z0 = float(cfg['TOPO']['z0'])

    c_ff = float(cfg['CHANNELS']['c_ff'])
    c_radius = float(cfg['CHANNELS']['c_radius'])
    c_init = float(cfg['CHANNELS']['c_init'])
    c_mustart = float(cfg['CHANNELS']['c_mustart'])
    c_muendFP = float(cfg['CHANNELS']['c_muendFP'])

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
    x_circ = x1 + xc
    # for plotting
    d1 = np.tan(np.radians(beta)) * x1

    # Set surface elevation
    zv = np.zeros((nRows, nCols))
    mask = np.zeros(np.shape(x))
    mask[np.where(x < (x1 - yc))] = 1
    zv = zv + (z0 - np.tan(np.radians(meanAlpha)) * x)*mask

    mask = np.zeros(np.shape(x))
    mask[np.where(((x1 - yc) <= x) & (x <= (x1 + xc)))] = 1
    # rCirc + np.sqrt(rCirc**2 - (xv[m] - x_circ)**2)
    zv = zv + (rCirc - np.sqrt(np.abs(rCirc**2 - (x_circ - x)**2)))*mask

    # If a channel shall be introduced
    if cfg['TOPO'].getboolean('channel'):
        # Compute cumulative distribution function - c_1 for upper part (start)
        # of channel and c_2 for lower part (end) of channel
        c_1 = norm.cdf(xv, c_mustart * (x1), c_ff)
        c_2 = 1. - norm.cdf(xv, c_muendFP * (x1), c_ff)

        # combine both into one function separated at the the middle of
        #  the channel longprofile location
        mask = np.zeros(np.shape(xv))
        mask[np.where(xv < (x1 * (0.5 * (c_mustart + c_muendFP))))] = 1
        c_0 = c_1 * mask

        mask = np.zeros(np.shape(xv))
        mask[np.where(xv >= (x1 * (0.5 * (c_mustart + c_muendFP))))] = 1
        c_0 = c_0 + c_2 * mask

        # Is the channel of constant width or narrowing
        if cfg['TOPO'].getboolean('narrowing'):
            c_extent = c_init * (1 - c_0[:]) + (c_0[:] * c_radius)
        else:
            c_extent = np.zeros(np.shape(xv)) + c_radius

        # Set surface elevation
        mask = np.zeros(np.shape(y))
        mask[np.where(abs(y) < c_extent)] = 1
        # Add surface elevation modification introduced by channel
        if cfg['TOPO'].getboolean('topoconst'):
            zv = zv - c_extent*c_0*np.sqrt(np.abs(1. - (np.square(y) / (c_extent**2))))*mask
        else:
            zv = zv + c_extent*c_0*(1. - np.sqrt(np.abs(1. - (np.square(y) / (c_extent**2)))))*mask
        if not cfg['TOPO'].getboolean('topoconst'):
            # outside of the channel, add layer of channel depth
            mask = np.zeros(np.shape(y))
            mask[np.where(abs(y) >= c_extent)] = 1
            zv = zv + c_extent*c_0*mask

    # Log info here
    log.info('Hockeystick smooth coordinates computed')

    return x, y, zv


def hockey(cfg):
    """
        Compute coordinates of a parabolically-shaped slope with a flat foreland
        defined by total fall height C, angle (meanAlpha) or distance (fLen) to flat foreland
    """

    C = float(cfg['TOPO']['C'])
    c_ff = float(cfg['CHANNELS']['c_ff'])
    c_radius = float(cfg['CHANNELS']['c_radius'])
    c_init = float(cfg['CHANNELS']['c_init'])
    c_mustart = float(cfg['CHANNELS']['c_mustart'])
    c_muend = float(cfg['CHANNELS']['c_muend'])

    # Get grid definitons
    dx, xEnd, yEnd = getGridDefs(cfg)

    # Compute coordinate grid
    xv, yv, zv, x, y, nRows, nCols = computeCoordGrid(dx, xEnd, yEnd)

    # Get parabola Parameters
    [A, B, fLen] = getParabolaParams(cfg)

    # Set surface elevation
    zv = np.ones((nRows, nCols))
    zv = zv * ((-B**2) / (4. * A) + C)
    mask = np.zeros(np.shape(xv))
    mask[np.where(xv < fLen)] = 1

    zv = zv + (A * xv**2 + B * xv + C)*mask

    # If a channel shall be introduced
    if cfg['TOPO'].getboolean('channel'):
        c_1 = norm.cdf(xv, c_mustart * fLen, c_ff)
        c_2 = 1. - norm.cdf(xv, c_muend * fLen, c_ff)

        # combine both into one function separated at the the middle of
        #  the channel longprofile location
        mask = np.zeros(np.shape(xv))
        mask[np.where(xv < (fLen * (0.5 * (c_mustart + c_muend))))] = 1
        c_0 = c_1 * mask

        mask = np.zeros(np.shape(xv))
        mask[np.where(xv >= (fLen * (0.5 * (c_mustart + c_muend))))] = 1
        c_0 = c_0 + c_2 * mask

        # Is the channel of constant width or narrowing
        if cfg['TOPO'].getboolean('narrowing'):
            c_extent = c_init * (1 - c_0[:]) + (c_0[:] * c_radius)
        else:
            c_extent = np.zeros(nCols) + c_radius


        # Add surface elevation modification introduced by channel
        mask = np.zeros(np.shape(y))
        mask[np.where(abs(y) < c_extent)] = 1
        if cfg['TOPO'].getboolean('topoconst'):
            zv = zv - c_extent*c_0*np.sqrt(np.abs(1. - (np.square(y) / (c_extent**2))))*mask
        else:
            zv = zv + c_extent*c_0*(1. - np.sqrt(np.abs(1. - (np.square(y) / (c_extent**2)))))*mask
        if not cfg['TOPO'].getboolean('topoconst'):
            # outside of the channel, add layer of channel depth
            mask = np.zeros(np.shape(y))
            mask[np.where(abs(y) >= c_extent)] = 1
            zv = zv + c_extent*c_0*mask

    # Log info here
    log.info('Hockeystick coordinates computed')

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
    c_ff = float(cfg['CHANNELS']['c_ff'])
    c_radius = float(cfg['CHANNELS']['c_radius'])
    c_init = float(cfg['CHANNELS']['c_init'])
    c_mustart = float(cfg['CHANNELS']['c_mustart'])
    c_muend = float(cfg['CHANNELS']['c_muend'])

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
        c_0 = np.zeros(np.shape(x))
        mask = np.zeros(np.shape(x))
        mask[np.where((theta * rHelix) < (0.5 * (c_mustart + c_muend) * fLen))] = 1
        c_0 = c_0 + norm.cdf(theta * rHelix, c_mustart * fLen, c_ff)*mask
        mask = np.zeros(np.shape(x))
        mask[np.where((theta * rHelix) >= (0.5 * (c_mustart + c_muend) * fLen))] = 1
        c_0 = c_0 + (1. - norm.cdf(theta * rHelix, c_muend * fLen, c_ff))*mask
        # c_0 = np.ones(np.shape(zv))
        # If channel of constant width or becoming narrower in the middle
        if cfg['TOPO'].getboolean('narrowing'):
            c_extent = c_init * (1. - c_0) + c_0 * c_radius
        else:
            c_extent = c_radius

        if not cfg['TOPO'].getboolean('topoconst'):
            zv = zv + c_0 * c_extent

        # Inner and outer boundaries of the channel
        bound_in = rHelix - c_extent
        bound_ext = rHelix + c_extent

        # Set channel
        mask = np.zeros(np.shape(x))
        mask[np.where((radius >= rHelix) & (radius < bound_ext))] = 1
        radius1 = radius - rHelix
        zv = zv - c_extent*c_0*np.sqrt(np.abs(1. - (np.square(radius1) / np.square(c_extent))))*mask

        mask = np.zeros(np.shape(x))
        mask[np.where((radius < rHelix) & (radius > bound_in))] = 1
        radius2 = rHelix - radius
        zv = zv - c_extent*c_0*np.sqrt(np.abs(1. - (np.square(radius1) / np.square(c_extent))))*mask


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
    nameExt = cfg['TOPO']['DEM_type']
    nRows = z.shape[0]
    nCols = z.shape[1]

    # Read lower left corner coordinates, cellsize and noDATA value
    xllcenter = float(cfg['DEMDATA']['xl'])
    yllcenter = float(cfg['DEMDATA']['yl'])
    cellsize = float(cfg['TOPO']['dx'])
    noDATA = float(cfg['DEMDATA']['nodata_value'])
    dem_name = cfg['DEMDATA']['dem_name']

    # Save elevation data to .asc file and add header lines
    z_mat = np.matrix(z)
    with open(os.path.join(outDir, '%s_%s_Topo.asc' % (dem_name, nameExt)), 'w') as f:
        f.write('nCols  %d\n' % (nCols))
        f.write('nRows  %d\n' % (nRows))
        f.write('xllcenter  %.02f\n' % (xllcenter))
        f.write('yllcenter %.02f\n' % (yllcenter))
        f.write('cellsize  %d\n' % (cellsize))
        f.write('nodata_value %.02f\n' % (noDATA))
        for line in z_mat:
            np.savetxt(f, line, fmt='%f')

    # Log info here
    log.info('DEM written to: %s/%s_%s_Topo.asc' % (outDir, dem_name, nameExt))


def generateTopo(cfg, avalancheDir):
    """ Compute coordinates of desired topography with given inputs """

    # Which DEM type
    demType = cfg['TOPO']['DEM_type']

    log.info('DEM type is set to: %s' % demType)

    # Set Output directory
    outDir = os.path.join(avalancheDir, 'Inputs')
    if os.path.isdir(outDir):
        log.info('The new DEM is saved to %s' % (outDir))
    else:
        log.error('Required folder structure: NameOfAvalanche/Inputs missing! \
                    Run runInitializeProject first!')

    # Call topography type
    if demType == 'FP':
        [x, y, z] = flatplane(cfg)

    elif demType == 'IP':
        [x, y, z] = inclinedplane(cfg)

    elif demType == 'HS':
        [x, y, z] = hockey(cfg)

    elif demType == 'HS2':
        [x, y, z] = hockeysmooth(cfg)

    elif demType == 'BL':
        [x, y, z] = bowl(cfg)

    elif demType == 'HX':
        [x, y, z] = helix(cfg)

    elif demType == 'PY':
        [x, y, z] = pyramid(cfg)

    # Write DEM to file
    writeDEM(cfg, z, outDir)

    return(z, demType, outDir)
