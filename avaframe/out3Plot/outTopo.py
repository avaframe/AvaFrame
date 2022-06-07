"""
    Simple plotting for DEMs
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np
import logging
from matplotlib.colors import LightSource
from matplotlib import cm

from scipy.interpolate import griddata
# local imports
from avaframe.in1Data import getInput
from avaframe.in3Utils import geoTrans

# create local logger
log = logging.getLogger(__name__)


def _generateDEMPlot(X, Y, z, title):
    """Generates 3d DEM plot, use this to style the plot"""

    plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')

    ls = LightSource(270, 45)
    # To use a custom hillshading mode, override the built-in shading and pass
    # in the rgb colors of the shaded surface calculated from "shade".
    rgb = ls.shade(z, cmap=cm.viridis, vert_exag=0.1, blend_mode='soft')
    surf = ax.plot_surface(X, Y, z, rstride=1, cstride=1, facecolors=rgb,
                           linewidth=0, antialiased=False, shade=False)

    # These are other options to plot in 3d in case another look is needed
    # ax.contour(X, Y, z, 20, linewidth=3, colors="g", linestyles="solid")
    # ax.plot_wireframe(X, Y, z, rstride=100, cstride=100,lw=1)
    # ax.plot_surface(X, Y, z, cmap=plt.cm.viridis,
    #                 linewidth=0, antialiased=False)

    ax.set_title('DEM: %s' % title)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('elevation (m)')

    return ax


def plotDEM3D(cfg, showPlot = False):
    """Plots the DEM from the avalancheDir in cfg alongside it

    Parameters
    ----------
    cfg : configparser object
      the main configuration
    showPlot : boolean
      If true shows the matplotlib plot

    """

    # get avalanche dir
    avalancheDir = cfg['MAIN']['avalancheDir']
    avaName = os.path.basename(avalancheDir)
    log.info('Plotting DEM in : %s', avalancheDir)

    # read DEM
    dem = getInput.readDEM(avalancheDir)

    # get DEM Path
    demPath = getInput.getDEMPath(avalancheDir)

    header = dem['header']
    xl = header['xllcenter']
    yl = header['yllcenter']
    ncols = header['ncols']
    nrows = header['nrows']
    dx = header['cellsize']
    z = dem['rasterData']

    # this line is needed for plot_surface to be able to handle the nans
    z[np.isnan(z)] = np.nanmin(z)

    # Set coordinate grid with given origin
    X, Y = geoTrans.makeCoordinateGrid(xl, yl, dx, ncols, nrows)

    ax = _generateDEMPlot(X, Y, z, avaName)

    # Save figure to file
    outName = os.path.splitext(demPath)[0] + '_plot.png'
    log.info('Saving plot to: %s', outName)
    plt.savefig(outName)

    # If flag is set, plot figure
    if showPlot:
        plt.show()

    plt.close()


def plotGeneratedDEM(z, nameExt, cfg, outDir, cfgMain):
    """ Plot DEM with given information on the origin of the DEM """

    cfgTopo = cfg['TOPO']
    cfgDEM = cfg['DEMDATA']

    # input parameters
    dx = float(cfgTopo['dx'])
    xl = float(cfgDEM['xl'])
    yl = float(cfgDEM['yl'])
    demName = cfgDEM['demName']

    # Set coordinate grid with given origin
    nrows, ncols = z.shape
    X, Y = geoTrans.makeCoordinateGrid(xl, yl, dx, ncols, nrows)

    topoNames = {'IP': 'inclined Plane', 'FP': 'flat plane', 'PF': 'parabola flat', 'TPF': 'triple parabola flat',
                 'HS': 'Hockeystick smoothed', 'BL': 'bowl', 'HX': 'Helix', 'PY': 'Pyramid'}

    ax = _generateDEMPlot(X, Y, z, topoNames[nameExt])

    # Save figure to file
    outName = os.path.join(outDir, '%s_%s_plot' % (demName, nameExt))

    log.info('Saving plot to: %s', outName)

    plt.savefig(outName)

    # If flag is set, plot figure
    if cfgMain['FLAGS'].getboolean('showPlot'):
        plt.show()

    plt.close()


def plotReleasePoints(xv, yv, xyPoints, demType):

    plt.figure()
    plt.plot(xv, np.zeros(len(xv))+yv[0], 'k-')
    plt.plot(xv, np.zeros(len(xv))+yv[-1], 'k-')
    plt.plot(np.zeros(len(yv))+xv[0], yv, 'k-')
    plt.plot(np.zeros(len(yv))+xv[-1], yv, 'k-')
    plt.plot(xyPoints[:, 0], xyPoints[:, 1], 'r*')
    plt.title('Domain and release area of %s - projected' % demType)
    plt.xlabel('along valley [m]')
    plt.ylabel('across valley [m]')

    plt.show()
