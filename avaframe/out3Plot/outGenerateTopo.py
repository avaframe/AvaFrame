"""
    Simple plotting for idealised/generic DEM results

    This file is part of Avaframe.
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np
import logging

# create local logger
log = logging.getLogger(__name__)

def plotGeneratedDEM(z, name_ext, cfg, outDir):
    """ Plot DEM with given information on the origin of the DEM """

    cfgTopo = cfg['TOPO']
    cfgDEM = cfg['DEMDATA']

    # input parameters
    dx = float(cfgTopo['dx'])
    xl = float(cfgDEM['xl'])
    yl = float(cfgDEM['yl'])
    dem_name = cfgDEM['dem_name']
    xEnd = z.shape[1] * dx
    yEnd = z.shape[0] * dx

    # Set coordinate grid with given origin
    xp = np.arange(xl, xl + xEnd, dx)
    yp = np.arange(yl, yl + yEnd, dx)
    X, Y = np.meshgrid(xp, yp)

    topoNames = {'IP': 'inclined Plane', 'FP': 'flat plane', 'HS': 'Hockeystick',
                 'HS2': 'Hockeystick smoothed', 'BL': 'bowl', 'HX': 'Helix', 'PY': 'Pyramid'}

    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z, cmap=plt.cm.viridis,
                    linewidth=0, antialiased=False)

    ax.set_title('Generated DEM: %s' % (topoNames[name_ext]))
    ax.set_xlabel('along valley distance [m]')
    ax.set_ylabel('across valley distance [m]')
    ax.set_zlabel('surface elevation [m]')

    # Save figure to file
    outName = os.path.join(outDir, '%s_%s_plot' % (dem_name, name_ext))

    log.info('Saving plot to: %s', outName)

    plt.savefig(outName)

    # If flag is set, plot figure
    if cfgDEM.getboolean('showplot'):
        plt.show()

    plt.close('all')


def plotReleasePoints(xv, yv, xyPoints, DEM_type):

    plt.figure()
    plt.plot(xv, np.zeros(len(xv))+yv[0], 'k-')
    plt.plot(xv, np.zeros(len(xv))+yv[-1], 'k-')
    plt.plot(np.zeros(len(yv))+xv[0], yv, 'k-')
    plt.plot(np.zeros(len(yv))+xv[-1], yv, 'k-')
    plt.plot(xyPoints[:, 0], xyPoints[:, 1], 'r*')
    plt.title('Domain and release area of %s - projected' % DEM_type)
    plt.xlabel('along valley [m]')
    plt.ylabel('across valley [m]')

    plt.show()
