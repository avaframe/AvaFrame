"""

  Create generic/idealised topographies - based on run_dem_creator.m JT Fischer (2011)

  This file is part of Avaframe.

"""


# load modules
import numpy as np
from scipy.stats import norm
import os
import configparser
import logging

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getParabolaParams(cfg):
    """ Compute parameters for parabola """

    # input parameters
    C = float(cfg['TOPO']['C'])
    fLens = float(cfg['TOPO']['f_lens'])
    meanAlpha = float(cfg['TOPO']['mean_alpha'])

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
    xEnd = float(cfg['TOPO']['x_end']) + dx
    yEnd = float(cfg['TOPO']['y_end']) + dx
    nRows = int(yEnd / dx)                    # number of rows
    nCols = int(xEnd / dx)                    # number of columns

    return dx, xEnd, yEnd, nRows, nCols


def flatplane(cfgTopo, ncols, nrows):
    """ Compute coordinates of flat plane topography """

    # input parameters
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx
    z_elev = float(cfgTopo['z_elev'])

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5 * y_end, 0.5 * y_end, dx)
    x, y = np.meshgrid(xv, yv)
    # Set elevation of surface
    zv = np.zeros((nrows, ncols)) + z_elev

    # Name extension for this type of topography
    name_ext = 'FP'

    # Log info here
    log.info('Flatplane coordinates computed')

    return x, y, zv, name_ext


def inclinedplane(cfgTopo, ncols, nrows, cfgChannel):
    """ Compute coordinates of inclined plane with given slope (mean_alpha)"""

    # input parameters
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx

    z0 = float(cfgTopo['z0'])
    mean_alpha = float(cfgTopo['mean_alpha'])
    c_ff = float(cfgChannel['c_ff'])
    c_radius = float(cfgChannel['c_radius'])

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5 * y_end, 0.5 * y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))

    # Set surface elevation from slope and max. elevation
    zv = z0 - np.tan(np.radians(mean_alpha)) * x

    # If a channel shall be introduced
    if cfgTopo.getboolean('channel'):
        # Compute cumulative distribution function and set horizontal extent of channel
        c_0 = norm.cdf(xv, 0, c_ff)
        c_extent = np.zeros(ncols) + c_radius

        # Introduce channel by cutting channel out as half sphere shaped feature
        # if flags['topoconst'] == 0 - add layer of channel depth everywhere and
        # cut out half-sphere shaped feature
        for m in range(ncols):
            for k in range(nrows):
                # if location within horizontal extent of channel,
                # make half sphere shaped channel with radius given by channel horizontal extent
                if abs(yv[k]) < c_extent[m]:
                    if cfgTopo.getboolean('topoconst'):
                        zv[k, m] = zv[k, m] - c_extent[m] * c_0[m] * \
                            np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2)))
                    else:
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                            (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))
                else:
                    if cfgTopo.getboolean('topoconst'):
                        # outside of the channel no modifcation
                        zv[k, m] = zv[k, m]
                    else:
                        # outside of the channel, add layer of channel depth
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m]

    # Name extension for this type of topography
    name_ext = 'IP'

    # Log info here
    log.info('Inclined plane coordinates computed')

    return x, y, zv, name_ext


def hockeysmooth(cfgTopo, ncols, nrows, cfgChannel):
    """
        Compute coordinates of an inclined plane with a flat foreland  defined by
        total fall height z0, angle to flat foreland (mean_alpha) and a radius (r_circ) to
        smooth the transition from inclined plane to flat foreland
    """

    # input parameters
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx
    r_circ = float(cfgTopo['r_circ'])
    mean_alpha = float(cfgTopo['mean_alpha'])
    z0 = float(cfgTopo['z0'])
    c_ff = float(cfgChannel['c_ff'])
    c_radius = float(cfgChannel['c_radius'])
    c_init = float(cfgChannel['c_init'])
    c_mustart = float(cfgChannel['c_mustart'])
    c_muendFP = float(cfgChannel['c_muendFP'])

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5 * y_end, 0.5 * y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))
    zv0 = np.zeros((nrows, ncols))
    # for plots
    cv = np.zeros((nrows, ncols))
    lv = np.zeros(len(xv))

    # Compute distance to flat foreland for given mean_alpha
    x1 = z0 / np.tan(np.radians(mean_alpha))
    if x1 >= x_end * 0.9:
        log.warning('Your domain (x_end) is to small or the slope angle (mean_alpha) to'
                        'shallow to produce a signifcant (>10 percent of domain, in your case:'
                        ' %.2f m) flat foreland!' % (0.1 * (x_end - dx)))

    # Compute circle parameters for smoothing the transition
    beta = (0.5 * (180. - (mean_alpha)))
    xc = r_circ / np.tan(np.radians(beta))
    yc = xc * np.cos(np.radians(mean_alpha))
    x_circ = x1 + xc
    # for plotting
    d1 = np.tan(np.radians(beta)) * x1

    # Set surface elevation
    for m in range(len(xv)):
        lv[m] = (r_circ / xc) * xv[m] - d1
        if xv[m] < x1 - yc:
            zv[:, m] = z0 - np.tan(np.radians(mean_alpha)) * xv[m]
            zv0[:, m] = z0 - np.tan(np.radians(mean_alpha)) * xv[m]
            lv[m] = np.nan
            cv[:, m] = np.nan
        elif x1 - yc <= xv[m] <= x1 + xc:
            r_circ + np.sqrt(r_circ**2 - (xv[m] - x_circ)**2)
            zv[:, m] = r_circ - np.sqrt(r_circ**2 - (x_circ - xv[m])**2)
            zv0[:, m] = r_circ - np.sqrt(r_circ**2 - (x_circ - xv[m])**2)
            cv[:, m] = r_circ - np.sqrt(r_circ**2 - (x_circ - xv[m])**2)
        else:
            zv[:, m] = 0.0
            zv0[:, m] = 0.0
            cv[:, m] = np.nan

    # If a channel shall be introduced
    if cfgTopo.getboolean('channel'):
        # Compute cumulative distribution function - c_1 for upper part (start)
        # of channel and c_2 for lower part (end) of channel
        c_1 = norm.cdf(xv, c_mustart * (x1), c_ff)
        c_2 = 1. - norm.cdf(xv, c_muendFP * (x1), c_ff)
        c_0 = np.zeros(ncols)

        # combine both into one function separated at the the middle of
        #  the channel longprofile location
        for l in range(ncols):
            if xv[l] < (x1 * (0.5 * (c_mustart + c_muendFP))):
                c_0[l] = c_1[l]
            else:
                c_0[l] = c_2[l]

        # Is the channel of constant width or narrowing
        if cfgTopo.getboolean('narrowing'):
            c_extent = c_init * (1 - c_0[:]) + (c_0[:] * c_radius)
        else:
            c_extent = np.zeros(ncols) + c_radius

        # Set surface elevation
        for m in range(ncols):
            for k in range(nrows):
                # Add surface elevation modification introduced by channel
                if cfgTopo.getboolean('channel'):
                    if abs(yv[k]) < c_extent[m]:
                        if cfgTopo.getboolean('topoconst'):
                            zv[k, m] = zv[k, m] - c_extent[m] * c_0[m] * \
                                np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2)))
                        else:
                            zv[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                                (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))
                    else:
                        if cfgTopo.getboolean('topoconst'):
                            # outside of the channel no modifcation
                            zv[k, m] = zv[k, m]
                        else:
                            # outside of the channel, add layer of channel depth
                            zv[k, m] = zv[k, m] + c_extent[m] * c_0[m]

    # Name extension for this type of topography
    name_ext = 'HS2'

    # Log info here
    log.info('Hockeystick smooth coordinates computed')

    return x, y, zv, name_ext


def hockey(cfgTopo, f_len, A, B, ncols, nrows, cfgChannel):
    """
        Compute coordinates of a parabolically-shaped slope with a flat foreland
        defined by total fall height C, angle (mean_alpha) or distance (f_len) to flat foreland
    """

    # input parameters
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx
    C = float(cfgTopo['C'])
    c_ff = float(cfgChannel['c_ff'])
    c_radius = float(cfgChannel['c_radius'])
    c_init = float(cfgChannel['c_init'])
    c_mustart = float(cfgChannel['c_mustart'])
    c_muend = float(cfgChannel['c_muend'])

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5 * y_end, 0.5 * y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))
    zv0 = np.zeros((nrows, ncols))
    zv1 = np.zeros((nrows, ncols))

    # If a channel shall be introduced

    # Get parabola Parameters
    [A, B, f_len] = getParabolaParams(cfgTopo)

    if cfgTopo.getboolean('channel'):
        c_1 = norm.cdf(xv, c_mustart * f_len, c_ff)
        c_2 = 1. - norm.cdf(xv, c_muend * f_len, c_ff)
        c_0 = np.zeros(ncols)

        for l in range(ncols):
            if xv[l] < (f_len * (0.5 * (c_mustart + c_muend))):
                c_0[l] = c_1[l]
            else:
                c_0[l] = c_2[l]

        # Is the channel of constant width or narrowing
        if cfgTopo.getboolean('narrowing'):
            c_extent = c_init * (1 - c_0[:]) + (c_0[:] * c_radius)
        else:
            c_extent = np.zeros(ncols) + c_radius

    # Set surface elevation
    for m in range(ncols):
        for k in range(nrows):
            if xv[m] < f_len:
                zv[k, m] = A * xv[m]**2 + B * xv[m] + C
                zv0[k, m] = A * xv[m]**2 + B * xv[m] + C
                zv1[k, m] = A * xv[m]**2 + B * xv[m] + C
            else:
                zv[k, m] = (-B**2) / (4. * A) + C
                zv0[k, m] = (-B**2) / (4. * A) + C
                zv1[k, m] = (-B**2) / (4. * A) + C

            # Add surface elevation modification introduced by channel
            if cfgTopo.getboolean('channel'):
                if abs(yv[k]) < c_extent[m]:
                    if cfgTopo.getboolean('topoconst'):
                        zv[k, m] = zv[k, m] - c_extent[m] * c_0[m] * \
                            np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2)))
                        zv1[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                            (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))
                    else:
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                            (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))

                else:
                    if cfgTopo.getboolean('topoconst'):
                        # outside of the channel no modifcation
                        zv[k, m] = zv[k, m]
                        zv1[k, m] = zv[k, m] + c_extent[m] * c_0[m]
                    else:
                        # outside of the channel, add layer of channel depth
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m]

    # Name extension for this type of topography
    name_ext = 'HS'

    # Log info here
    log.info('Hockeystick coordinates computed')

    return x, y, zv, name_ext


def bowl(cfgTopo, ncols, nrows):
    """ Compute coordinates of sphere with given radius (r_bowl) """

    # input parameters
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx
    r_bowl = float(cfgTopo['r_bowl'])

    # Compute coordinate grid
    xv = np.arange(-0.5 * x_end, 0.5 * x_end, dx)
    yv = np.arange(-0.5 * y_end, 0.5 * y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))

    # Set surface elevation
    for m in range(ncols):
        for k in range(nrows):
            radius = np.sqrt(xv[m]**2 + yv[k]**2)
            if radius <= r_bowl:
                zv[k, m] = r_bowl - r_bowl * np.sqrt(1 - (radius / r_bowl)**2)
            else:
                zv[k, m] = r_bowl

    # Name extension for this type of topography
    name_ext = 'BL'

    # Log info here
    log.info('Bowl coordinates computed')

    return x, y, zv, name_ext


def helix(cfgTopo, ncols, nrows, f_len, A, B, cfgChannel):
    """ Compute coordinates of helix-shaped topography with given radius (r_helix) """

    # input parameters
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx
    r_helix = float(cfgTopo['r_helix'])
    C = float(cfgTopo['C'])
    c_ff = float(cfgChannel['c_ff'])
    c_radius = float(cfgChannel['c_radius'])
    c_init = float(cfgChannel['c_init'])
    c_mustart = float(cfgChannel['c_mustart'])
    c_muend = float(cfgChannel['c_muend'])

    # Compute coordinate grid
    xv = np.arange(-0.5 * x_end, 0.5 * x_end, dx)
    yv = np.arange(-y_end, 0, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))

    # Set surface elevation

    # Get parabola Parameters
    [A, B, f_len] = getParabolaParams(cfgTopo)

    for m in range(ncols):
        for k in range(nrows):
            radius = np.sqrt(xv[m]**2 + yv[k]**2)
            theta = np.arctan2(yv[k], xv[m]) + np.pi
            if (theta * r_helix) < f_len:
                zv[k, m] = A * (theta * r_helix)**2 + B * (theta * r_helix) + C
            else:
                zv[k, m] = (-B**2) / (4. * A) + C

            # If channel is introduced to topography
            if cfgTopo.getboolean('channel'):
                if (theta * r_helix) < (0.5 * (c_mustart + c_muend) * f_len):
                    c_0 = norm.cdf(theta * r_helix, c_mustart * f_len, c_ff)
                else:
                    c_0 = 1. - norm.cdf(theta * r_helix, c_muend * f_len, c_ff)

                # If channel of constant width or becoming narrower in the middle
                if cfgTopo.getboolean('narrowing'):
                    c_extent = c_init * (1. - c_0) + c_0 * c_radius
                else:
                    c_extent = c_radius

                # Inner and outer boundaries of the channel
                bound_in = r_helix - c_extent
                bound_ext = r_helix + c_extent

                # Set channel
                if (radius >= r_helix) and (radius < bound_ext):
                    radius = radius - r_helix
                    if cfgTopo.getboolean('topoconst'):
                        zv[k, m] = zv[k, m] - c_0 * c_extent * \
                            np.sqrt(1. - (radius**2 / c_extent**2))
                    else:
                        zv[k, m] = zv[k, m] + c_0 * c_extent * \
                            (1. - np.sqrt(1. - (radius**2 / c_extent**2)))

                elif (radius < r_helix) and (radius > bound_in):
                    radius = r_helix - radius
                    if cfgTopo.getboolean('topoconst'):
                        zv[k, m] = zv[k, m] - c_0 * c_extent * \
                            np.sqrt(1. - (radius**2 / c_extent**2))
                    else:
                        zv[k, m] = zv[k, m] + c_0 * c_extent * \
                            (1. - np.sqrt(1. - (radius**2 / c_extent**2)))
                else:
                    if cfgTopo.getboolean('topoconst'):
                        zv[k, m] = zv[k, m]
                    else:
                        zv[k, m] = zv[k, m] + c_0 * c_extent

    # Name extension for this type of topography
    name_ext = 'HX'

    # Log info here
    log.info('Helix coordinates computed')

    return x, y, zv, name_ext


def writeDEM(z, name_ext, ncols, nrows, dx, cfgDEM, outDir):
    """ Write topography information to file """

    # Read lower left corner coordinates, cellsize and noDATA value
    xllcorner = float(cfgDEM['xl'])
    yllcorner = float(cfgDEM['yl'])
    cellsize = dx
    noDATA = float(cfgDEM['nodata_value'])
    dem_name = cfgDEM['dem_name']

    # Save elevation data to .asc file and add header lines
    z_mat = np.matrix(z)
    with open(os.path.join(outDir, '%s_%s_Topo.asc' % (dem_name, name_ext)), 'w') as f:
        f.write('ncols  %d\n' % (ncols))
        f.write('nrows  %d\n' % (nrows))
        f.write('xllcorner  %.02f\n' % (xllcorner))
        f.write('yllcorner %.02f\n' % (yllcorner))
        f.write('cellsize  %d\n' % (cellsize))
        f.write('nodata_value %.02f\n' % (noDATA))
        for line in z_mat:
            np.savetxt(f, line, fmt='%f')

    # Log info here
    log.info('DEM written to: %s/%s_%s_Topo.asc' % (outDir, dem_name, name_ext))


def generateTopo(cfg, avalancheDir):
    """ Compute coordinates of desired topography with given inputs """

    # Which DEM type
    demType = cfg['TOPO']['demType']

    log.info('DEM type is set to: %s' % demType)

    # Set Output directory
    outDir = os.path.join(avalancheDir, 'Inputs')
    if os.path.isdir(outDir):
        log.info('The new DEM is saved to %s' % (outDir))
    else:
        log.error('Required folder structure: NameOfAvalanche/Inputs missing! \
                    Run runInitializeProject first!')

    dx, xEnd, yEnd, nRows, nCols = getGridDefs(cfg)

    # Call topography type
    if demType == 'FP':
        [x, y, z, nameExt] = flatplane(cfg)

    elif demType == 'IP':
        [x, y, z, nameExt] = inclinedplane(cfg)

    elif demType == 'HS':
        [x, y, z, nameExt] = hockey(cfg)

    elif demType == 'HS2':
        [x, y, z, nameExt] = hockeysmooth(cfg)

    elif demType == 'BL':
        [x, y, z, nameExt] = bowl(cfgTopo)

    elif demType == 'HX':
        [x, y, z, nameExt] = helix(cfgTopo)

    # Write DEM to file
    writeDEM(z, nameExt, nCols, nRows, float(cfg['TOPO']['dx']), cfg['DEM'], outDir)

    return(z, nameExt, outDir)
