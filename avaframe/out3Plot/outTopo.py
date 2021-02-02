"""
    Simple plotting for DEMs
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np
import logging

# local imports
from avaframe.in1Data import getInput

# create local logger
log = logging.getLogger(__name__)

def _generateDEMPlot(X, Y, z,):
    """Generates 3d DEM plot, use this to style the plot"""
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z, cmap=plt.cm.viridis,
                    linewidth=0, antialiased=False)
    return ax

def plotDEM(cfg):
    """Plots the DEM from the avalancheDir in cfg alongside it

    Parameters
    ----------
    cfg : configparser object
      the main configuration

    Returns
    -------
    plt : matplotlib object
      plot object in case displaying is needed
    """

    # get avalanche dir
    avalancheDir = cfg['MAIN']['avalancheDir']
    log.info('Plotting DEM in : %s', avalancheDir)

    # read DEM
    dem = getInput.readDEM(avalancheDir)

    header = dem['header']
    xl = header.xllcenter
    yl = header.yllcenter
    dx = header.cellsize
    z = dem['rasterData']

    # Set coordinate grid with given origin
    X,Y = _setCoordinateGrid(xl,yl,dx,z)

def plotGeneratedDEM(z, nameExt, cfg, outDir):
    """ Plot DEM with given information on the origin of the DEM """

    cfgTopo = cfg['TOPO']
    cfgDEM = cfg['DEMDATA']

    # input parameters
    dx = float(cfgTopo['dx'])
    xl = float(cfgDEM['xl'])
    yl = float(cfgDEM['yl'])
    demName = cfgDEM['demName']

    # Set coordinate grid with given origin
    X,Y = _setCoordinateGrid(xl,yl,dx,z)

    topoNames = {'IP': 'inclined Plane', 'FP': 'flat plane', 'HS': 'Hockeystick',
                 'HS2': 'Hockeystick smoothed', 'BL': 'bowl', 'HX': 'Helix', 'PY': 'Pyramid'}

    ax = _generateDEMPlot(X, Y, z)

    ax.set_title('Generated DEM: %s' % (topoNames[nameExt]))
    ax.set_xlabel('along valley distance [m]')
    ax.set_ylabel('across valley distance [m]')
    ax.set_zlabel('surface elevation [m]')

    # Save figure to file
    outName = os.path.join(outDir, '%s_%s_plot' % (demName, nameExt))

    log.info('Saving plot to: %s', outName)

    plt.savefig(outName)

    # If flag is set, plot figure
    if cfgDEM.getboolean('showplot'):
        plt.show()

    plt.close('all')


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


def _setCoordinateGrid(xl,yl,dx,z):
    """get a Coordinate Grid for plotting"""

    xEnd = z.shape[1] * dx
    yEnd = z.shape[0] * dx

    xp = np.arange(xl, xl + xEnd, dx)
    yp = np.arange(yl, yl + yEnd, dx)
    X, Y = np.meshgrid(xp, yp)
    return(X,Y)
