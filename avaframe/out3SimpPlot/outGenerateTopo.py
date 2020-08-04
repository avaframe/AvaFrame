"""
    Simple plotting for idealised/generic DEM results

    This file is part of Avaframe.
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np
import configparser


def plotDEM(z, name_ext):
    """ Plot DEM with given information on the origin of the DEM """

    # Load all input Parameters
    cfg = configparser.ConfigParser()
    if os.path.isfile('avaframe/in3Utils/local_generateTopoCfg.ini'):
        cfg.read('avaframe/in3Utils/local_generateTopoCfg.ini')
    else:
        cfg.read('avaframe/in3Utils/generateTopoCfg.ini')
    cfgTopo = cfg['TOPO']
    cfgDEM = cfg['DEMDATA']

    # input parameters
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx
    xl = float(cfgDEM['xl'])
    yl = float(cfgDEM['yl'])
    dem_name = cfgDEM['dem_name']
    dpath = cfgDEM['path']

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

    # Save figure to file
    plt.savefig(os.path.join(dpath, '%s_%s_plot.png' % (dem_name, name_ext)))

    # If flag is set, plot figure
    if int(cfgDEM['showplot']) == 1:
        plt.show()
