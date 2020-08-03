"""

  Create generic/idealised topographies - based on run_dem_creator.m JT Fischer (2011)

  This file is part of Avaframe.

"""


# load modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.stats import norm
import seaborn
import os
import configparser
import logging
logging.basicConfig(level=logging.INFO)


def getParams(cT):
    """ Compute parameters for parabola """

    # input parameters
    C = float(cT['C'])
    f_lens = float(cT['f_lens'])
    mean_alpha = float(cT['mean_alpha'])

    # If mean slope is given or distance to the start of the flat plane
    if mean_alpha != 0:
        f_len = C / np.tan(np.radians(mean_alpha))
        logging.info('flen computed from mean alpha: %.2f meters' % f_len)
    else:
        f_len = f_lens
        logging.info('flen directly set to: %.2f meters' % f_len)
    A = C / (f_len**2)
    B = (-C * 2.) / f_len

    return A, B, f_len


def flatplane(cG, ncols, nrows, z_elev):
    """ Compute coordinates of flat plane topography """

    # input parameters
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5*y_end, 0.5*y_end, dx)
    x, y = np.meshgrid(xv, yv)
    # Set elevation of surface
    zv = np.zeros((nrows, ncols)) + z_elev

    # Name extension for this type of topography
    name_ext = 'FP'

    # Log info here
    logging.info('Flatplane coordinates computed')

    return x, y, zv, name_ext


def inclinedplane(cG, ncols, nrows, cT, cC, cF):
    """ Compute coordinates of inclined plane with given slope (mean_alpha)"""

    # input parameters
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx
    z0 = float(cT['z0'])
    mean_alpha = float(cT['mean_alpha'])
    c_ff = float(cC['c_ff'])
    c_radius = float(cC['c_radius'])

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5*y_end, 0.5*y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))
    zv0 = np.zeros((nrows, ncols))

    # Set surface elevation from slope and max. elevation
    zv = z0 - np.tan(np.radians(mean_alpha)) * x

    # If a channel shall be introduced
    if float(float(cF['channel'])) == 1:
        # Compute cumulative distribution function and set horizontal extent of channel
        c_0 = norm.cdf(xv, 0, c_ff)
        c_extent = np.zeros(ncols) + c_radius

        # Introduce channel by cutting channel out as half sphere shaped feature
        # if flags['topoconst'] == 0 - add layer of channel depth everywhere and cut out half-sphere shaped feature
        for m in range(ncols):
            for k in range(nrows):
                # if location within horizontal extent of channel, make half sphere shaped channel with radius given by channel horizontal extent
                if abs(yv[k]) < c_extent[m]:
                    if float(cF['topoconst']) == 1:
                        zv[k, m] = zv[k, m] - c_extent[m] * c_0[m] * \
                            np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2)))
                    elif float(cF['topoconst']) == 0:
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                            (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))
                else:
                    if float(cF['topoconst']) == 1:
                        # outside of the channel no modifcation
                        zv[k, m] = zv[k, m]
                    elif float(cF['topoconst']) == 0:
                        # outside of the channel, add layer of channel depth
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m]

    # Name extension for this type of topography
    name_ext = 'IP'

    # Log info here
    logging.info('Inclined plane coordinates computed')

    return x, y, zv, name_ext


def hockeysmooth(cG, ncols, nrows, cT, cC, cF):
    """
        Compute coordinates of an inclined plane with a flat foreland
        defined by total fall height z0, angle to flat foreland (mean_alpha) and a radius (r_circ) to
        smooth the transition from inclined plane to flat foreland
    """

    # input parameters
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx
    r_circ = float(cT['r_circ'])
    mean_alpha = float(cT['mean_alpha'])
    z0 = float(cT['z0'])
    c_ff = float(cC['c_ff'])
    c_radius = float(cC['c_radius'])
    c_init = float(cC['c_init'])
    c_mustart = float(cC['c_mustart'])
    c_muendFP = float(cC['c_muendFP'])

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5*y_end, 0.5*y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))
    zv0 = np.zeros((nrows, ncols))
    # for plots
    cv = np.zeros((nrows, ncols))
    lv = np.zeros(len(xv))

    # Compute distance to flat foreland for given mean_alpha
    x1 = z0 / np.tan(np.radians(mean_alpha))
    if x1 >= x_end * 0.9:
        logging.warning('Your domain (x_end) is to small or the slope angle (mean_alpha) to shallow to produce a signifcant (>10 percent of domain, in your case: %.2f m) flat foreland!' % (0.1*(x_end-dx)))

    # Compute circle parameters for smoothing the transition
    beta = (0.5 * (180. - (mean_alpha)))
    xc = r_circ / np.tan(np.radians(beta))
    yc = xc * np.cos(np.radians(mean_alpha))
    zc = r_circ
    x_circ = x1 + xc
    # for plotting
    d1 = np.tan(np.radians(beta)) * x1

    # Set surface elevation
    for m in range(len(xv)):
        lv[m] = (r_circ / xc) * xv[m] - d1
        if xv[m] < x1-yc:
            zv[:, m] = z0 - np.tan(np.radians(mean_alpha)) * xv[m]
            zv0[:, m] = z0 - np.tan(np.radians(mean_alpha)) * xv[m]
            lv[m] = np.nan
            cv[:, m] = np.nan
        elif x1-yc <= xv[m] <= x1+xc:
            r_circ + np.sqrt(r_circ**2 - (xv[m] - x_circ)**2)
            zv[:, m] = r_circ - np.sqrt(r_circ**2 - (x_circ - xv[m])**2)
            zv0[:, m] = r_circ - np.sqrt(r_circ**2 - (x_circ - xv[m])**2)
            cv[:, m] = r_circ - np.sqrt(r_circ**2 - (x_circ - xv[m])**2)
        else:
            zv[:, m] = 0.0
            zv0[:, m] = 0.0
            cv[:, m] = np.nan

    # # Plot the smoothing transition
    # plt.plot(xv, zv[1,:])
    # plt.plot(xv, cv[1,:], 'r-')
    # plt.plot(xv, lv, 'g-')
    # plt.plot(x_circ, zc, 'g*')
    # plt.title('Long profile of hockey stick smooth topo')
    # plt.xlabel('Along valley distance [m]')
    # plt.ylabel('Surface elevation [m]')
    # plt.axis('equal')
    # plt.show()

    # If a channel shall be introduced
    if float(cF['channel']) == 1:
        # Compute cumulative distribution function - c_1 for upper part (start)
        # of channel and c_2 for lower part (end) of channel
        c_1 = norm.cdf(xv, c_mustart * (x1), c_ff)
        c_2 = 1. - norm.cdf(xv, c_muendFP * (x1), c_ff)
        c_0 = np.zeros(ncols)

        # combine both into one function separated at the the middle of the channel longprofile location
        for l in range(ncols):
            if xv[l] < (x1 * (0.5 * (c_mustart + c_muendFP))):
                c_0[l] = c_1[l]
            else:
                c_0[l] = c_2[l]

        # Is the channel of constant width or narrowing
        if float(cF['narrowing']) == 1:
            c_extent = c_init * (1 - c_0[:]) + (c_0[:] * c_radius)
        else:
            c_extent = np.zeros(ncols) + c_radius

        # Set surface elevation
        for m in range(ncols):
            for k in range(nrows):
                # Add surface elevation modification introduced by channel
                if float(cF['channel']) == 1:
                    if abs(yv[k]) < c_extent[m]:
                        if float(cF['topoconst']) == 1:
                            zv[k, m] = zv[k, m] - c_extent[m] * c_0[m] * \
                                np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2)))
                        elif float(cF['topoconst']) == 0:
                            zv[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                                (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))
                    else:
                        if float(cF['topoconst']) == 1:
                            # outside of the channel no modifcation
                            zv[k, m] = zv[k, m]
                        elif float(cF['topoconst']) == 0:
                            # outside of the channel, add layer of channel depth
                            zv[k, m] = zv[k, m] + c_extent[m] * c_0[m]

    # Name extension for this type of topography
    name_ext = 'HS2'

    # Log info here
    logging.info('Hockeystick smooth coordinates computed')

    return x, y, zv, name_ext


def hockey(cG, f_len, A, B, ncols, nrows, cT, cC, cF):
    """
        Compute coordinates of a parabolically-shaped slope with a flat foreland
        defined by total fall height C, angle (mean_alpha) or distance (f_len) to flat foreland
    """

    # input parameters
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx
    mean_alpha = float(cT['mean_alpha'])
    C = float(cT['C'])
    c_ff = float(cC['c_ff'])
    c_radius = float(cC['c_radius'])
    c_init = float(cC['c_init'])
    c_mustart = float(cC['c_mustart'])
    c_muend = float(cC['c_muend'])

    # Compute coordinate grid
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5*y_end, 0.5*y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))
    zv0 = np.zeros((nrows, ncols))
    zv1 = np.zeros((nrows, ncols))

    # If a channel shall be introduced
    if float(cF['channel']) == 1:
        c_1 = norm.cdf(xv, c_mustart * f_len, c_ff)
        c_2 = 1. - norm.cdf(xv, c_muend * f_len, c_ff)
        c_0 = np.zeros(ncols)

        for l in range(ncols):
            if xv[l] < (f_len * (0.5 * (c_mustart + c_muend))):
                c_0[l] = c_1[l]
            else:
                c_0[l] = c_2[l]

        # Is the channel of constant width or narrowing
        if float(cF['narrowing']) == 1:
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
                zv[k, m] = (-B**2) / (4.*A) + C
                zv0[k, m] = (-B**2) / (4.*A) + C
                zv1[k, m] = (-B**2) / (4.*A) + C

            # Add surface elevation modification introduced by channel
            if float(cF['channel']) == 1:
                if abs(yv[k]) < c_extent[m]:
                    if float(cF['topoconst']) == 1:
                        zv[k, m] = zv[k, m] - c_extent[m] * c_0[m] * \
                            np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2)))
                        zv1[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                            (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))
                    elif float(cF['topoconst']) == 0:
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m] * \
                            (1. - np.sqrt(1. - (yv[k]**2 / (c_extent[m]**2))))

                else:
                    if float(cF['topoconst']) == 1:
                        # outside of the channel no modifcation
                        zv[k, m] = zv[k, m]
                        zv1[k, m] = zv[k, m] + c_extent[m] * c_0[m]
                    elif float(cF['topoconst']) == 0:
                        # outside of the channel, add layer of channel depth
                        zv[k, m] = zv[k, m] + c_extent[m] * c_0[m]

    # Name extension for this type of topography
    name_ext = 'HS'

    # Log info here
    logging.info('Hockeystick coordinates computed')

    return x, y, zv, name_ext


def bowl(cG, ncols, nrows, r_bowl):
    """ Compute coordinates of sphere with given radius (r_bowl) """

    # input parameters
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx

    # Compute coordinate grid
    xv = np.arange(-0.5*x_end, 0.5*x_end, dx)
    yv = np.arange(-0.5*y_end, 0.5*y_end, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))

    # Set surface elevation
    for m in range(ncols):
        for k in range(nrows):
            radius = np.sqrt(xv[m]**2 + yv[k]**2)
            if radius <= r_bowl:
                zv[k, m] = r_bowl - r_bowl * np.sqrt(1 - (radius/r_bowl)**2)
            else:
                zv[k, m] = r_bowl

    # Name extension for this type of topography
    name_ext = 'BL'

    # Log info here
    logging.info('Bowl coordinates computed')

    return x, y, zv, name_ext


def helix(cG, ncols, nrows, cT, f_len, A, B, cC, cF):
    """ Compute coordinates of helix-shaped topography with given radius (r_helix) """

    # input parameters
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx
    r_helix = float(cT['r_helix'])
    C = float(cT['C'])
    c_ff = float(cC['c_ff'])
    c_radius = float(cC['c_radius'])
    c_init = float(cC['c_init'])
    c_mustart = float(cC['c_mustart'])
    c_muend = float(cC['c_muend'])

    # Compute coordinate grid
    xv = np.arange(-0.5*x_end, 0.5*x_end, dx)
    yv = np.arange(-y_end, 0, dx)
    x, y = np.meshgrid(xv, yv)
    zv = np.zeros((nrows, ncols))

    # Set surface elevation
    for m in range(ncols):
        for k in range(nrows):
            radius = np.sqrt(xv[m]**2 + yv[k]**2)
            theta = np.arctan2(yv[k], xv[m]) + np.pi
            if (theta*r_helix) < f_len:
                zv[k, m] = A * (theta*r_helix)**2 + B * (theta*r_helix) + C
            else:
                zv[k, m] = (-B**2) / (4.*A) + C

            # If channel is introduced to topography
            if float(cF['channel']) == 1:
                if (theta * r_helix) < (0.5 * (c_mustart + c_muend) * f_len):
                    c_0 = norm.cdf(theta*r_helix, c_mustart*f_len, c_ff)
                else:
                    c_0 = 1. - norm.cdf(theta*r_helix, c_muend*f_len, c_ff)

                # If channel of constant width or becoming narrower in the middle
                if float(cF['narrowing']) == 1:
                    c_extent = c_init * (1. - c_0) + c_0 * c_radius
                else:
                    c_extent = c_radius

                # Inner and outer boundaries of the channel
                bound_in = r_helix - c_extent
                bound_ext = r_helix + c_extent

                # Set channel
                if (radius >= r_helix) and (radius < bound_ext):
                    radius = radius - r_helix
                    if float(cF['topoconst']) == 1:
                        zv[k, m] = zv[k, m] - c_0 * c_extent * \
                            np.sqrt(1. - (radius**2 / c_extent**2))
                    elif float(cF['topoconst']) == 0:
                        zv[k, m] = zv[k, m] + c_0 * c_extent * \
                            (1. - np.sqrt(1. - (radius**2 / c_extent**2)))

                elif (radius < r_helix) and (radius > bound_in):
                    radius = r_helix - radius
                    if float(cF['topoconst']) == 1:
                        zv[k, m] = zv[k, m] - c_0 * c_extent * \
                            np.sqrt(1. - (radius**2 / c_extent**2))
                    elif float(cF['topoconst']) == 0:
                        zv[k, m] = zv[k, m] + c_0 * c_extent * \
                            (1. - np.sqrt(1. - (radius**2 / c_extent**2)))
                else:
                    if float(cF['topoconst']) == 1:
                        zv[k, m] = zv[k, m]
                    elif float(cF['topoconst']) == 0:
                        zv[k, m] = zv[k, m] + c_0 * c_extent

    # Name extension for this type of topography
    name_ext = 'HX'

    # Log info here
    logging.info('Helix coordinates computed')

    return x, y, zv, name_ext


def plotDEM(x, y, z, name_ext, cD, cG):
    """ Plot DEM with given information on the origin of the DEM """

    # input parameters
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx
    xl = float(cD['xl'])
    yl = float(cD['yl'])

    # Set coordinate grid with given origin
    xp = np.arange(xl, xl + x_end, dx)
    yp = np.arange(yl, yl + y_end, dx)
    X, Y = np.meshgrid(xp, yp)

    topoNames = {'IP': 'inclined Plane', 'FP': 'flat plane', 'HS': 'Hockeystick',
                 'HS2': 'Hockeystick smoothed', 'BL': 'bowl', 'HX': 'Helix'}

    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z, cmap=plt.cm.viridis,
                    linewidth=0, antialiased=False)

    ax.set_title('Generated DEM: %s' % (topoNames[name_ext]))
    ax.set_xlabel('along valley distance [m]')
    ax.set_ylabel('across valley distance [m]')
    ax.set_zlabel('surface elevation [m]')

    plt.show()


def writeDEM(z, name_ext, ncols, nrows, dx, cD):
    """ Write topography information to file """

    # Read lower left corner coordinates, cellsize and noDATA value
    xllcorner = float(cD['xl'])
    yllcorner = float(cD['yl'])
    cellsize = dx
    noDATA = float(cD['nodata_value'])
    dem_name = cD['dem_name']
    dpath = cD['path']

    # Save elevation data to .asc file and add header lines
    z_mat = np.matrix(z)
    with open(os.path.join(dpath, '%s_%s_Topo.asc' % (dem_name, name_ext)), 'w') as f:
        f.write('ncols  %d\n' % (ncols))
        f.write('nrows  %d\n' % (nrows))
        f.write('xllcorner  %.02f\n' % (xllcorner))
        f.write('yllcorner %.02f\n' % (yllcorner))
        f.write('cellsize  %d\n' % (cellsize))
        f.write('nodata_value %.02f\n' % (noDATA))
        for line in z_mat:
            np.savetxt(f, line, fmt='%f')

    # Log info here
    logging.info('DEM written to: %s_%s_Topo.asc' % (dem_name, name_ext))


def generateTopo():
    """ Compute coordinates of desired topography with given inputs """

    # Load all input Parameters
    cfg = configparser.ConfigParser()
    if os.path.isfile('avaframe/in3Utils/local_generateTopoCfg.ini') == True:
        cfg.read('avaframe/in3Utils/local_generateTopoCfg.ini')
    else:
        cfg.read('avaframe/in3Utils/generateTopoCfg.ini')
    cG = cfg['GENERAL']
    cT = cfg['TOPO']
    cF = cfg['FLAGS']
    cC = cfg['CHANNEL']
    cD = cfg['DEMDATA']

    # Which DEM type
    DEM_type = cT['DEM_type']

    logging.info('DEM type is set to: %s' % DEM_type)
    # determine number of rows and columns to define domain
    dx = float(cG['dx'])
    x_end = float(cG['x_end']) + dx
    y_end = float(cG['y_end']) + dx
    nrows = int(y_end / dx)                    # number of rows
    ncols = int(x_end / dx)                    # number of columns

    # Get parabola Parameters
    [A, B, f_len] = getParams(cT)

    # Call topography type
    if DEM_type == 'FP':
        [x, y, z, name_ext] = flatplane(cG, ncols, nrows, float(cT['z_elev']))

    elif DEM_type == 'IP':
        [x, y, z, name_ext] = inclinedplane(cG, ncols, nrows, cT, cC, cF)

    elif DEM_type == 'HS':
        [x, y, z, name_ext] = hockey(cG, f_len, A, B, ncols, nrows, cT, cC, cF)

    elif DEM_type == 'HS2':
        [x, y, z, name_ext] = hockeysmooth(cG, ncols, nrows, cT, cC, cF)

    elif DEM_type == 'BL':
        [x, y, z, name_ext] = bowl(cG, ncols, nrows, float(cT['r_bowl']))

    elif DEM_type == 'HX':
        [x, y, z, name_ext] = helix(cG, ncols, nrows, cT, f_len, A, B, cC, cF)

    # Plot new topogrpahy
    plotDEM(x, y, z, name_ext, cD, cG)

    # Write DEM to file
    writeDEM(z, name_ext, ncols, nrows, float(cG['dx']), cD)
