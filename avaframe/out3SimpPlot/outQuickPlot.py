"""

This is a simple function for a quick plot of datasets and option for comparison
between two datasets of identical shape. Also plots cross and longprofiles.

set desired values in

This file is part of Avaframe.

"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from avaframe.in3Utils import fileHandlerUtils as fU
import numpy as np
import os
import logging
import shutil
import glob

# Local imports
import avaframe.in3Utils.ascUtils as IOf
from avaframe.out3SimpPlot.plotSettings import *

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def quickPlot(avaDir, suffix, cfg, simName):
    """ Plot two raster datasets of identical dimension:

        Inputs:

        avaDir          avalanche directory
        suffix          result parameter abbreviation (e.g. 'ppr')
        com1DFAOutput   folder where results to be plotted are located
        simName         entres or null for simulation type
        cfgR            configuration for plots

        Outputs:

        figure 1: plot raster data for dataset1, dataset2 and their difference
        figure 2: plot cross and longprofiles for both datasets (ny_loc and nx_loc define location of profiles)
        -plots are saved to Outputs/out3SimpPlot
    """

    # Create required directories
    workDir = os.path.join(avaDir, 'Work', 'out3SimpPlot')
    fU.makeADir(workDir)
    outDir = os.path.join(avaDir, 'Outputs', 'out3SimpPlot')
    fU.makeADir(outDir)

    # Setup input from com1DFA
    fU.getDFAData(avaDir, workDir, suffix)

    # Get data from reference run
    fU.getRefData(avaDir, workDir, suffix)

    # prepare data
    data = fU.makeSimDict(workDir)
    cellSize = data['cellSize'][0]

    # Count the number of release areas
    relAreas = set(data['releaseArea'])

    for rel in relAreas:
        # get list of indices of files that are of correct simulation type and result paramete
        indSuffix = []
        for m in range(len(data['files'])):
            if data['resType'][m] == suffix and data['simType'][m] == simName and data['releaseArea'][m] == rel:
                indSuffix.append(m)

        # Load data
        data1 = np.loadtxt(data['files'][indSuffix[0]], skiprows=6)
        data2 = np.loadtxt(data['files'][indSuffix[1]], skiprows=6)
        ny = data1.shape[0]
        nx = data1.shape[1]
        Ly = ny*cellSize
        Lx = nx*cellSize
        log.info('dataset1: %s' % data['files'][indSuffix[0]])
        log.info('dataset2: %s' % data['files'][indSuffix[1]])

        # Location of Profiles
        ny_loc = int(nx * 0.5)
        nx_loc = int(ny * 0.5)

        # Plot data
        # Figure 1 shows the result parameter data
        fig = plt.figure(figsize=(figW*3, figH), dpi=figReso)
        ax1 = fig.add_subplot(131)
        cmap = cmapGB
        im1 = plt.imshow(data1, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
        fig.colorbar(im1, ax=ax1)
        ax1.set_aspect('auto')
        ax1.set_title('%s' % data['names'][indSuffix[0]])
        ax1.set_xlabel('$x\;[m]$')
        ax1.set_ylabel('$y\;[m]$')
        ax2 = fig.add_subplot(132)
        im2 = plt.imshow(data2, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
        fig.colorbar(im2, ax=ax2)
        ax2.set_aspect('auto')
        ax2.set_xlabel('$x\;[m]$')
        ax2.set_title('%s' % data['names'][indSuffix[1]])
        ax3 = fig.add_subplot(133)
        cmap = cmapdiv
        im3 = plt.imshow(data1-data2, cmap=cmap,
                         extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
        fig.colorbar(im3, ax=ax3)
        ax3.set_aspect('auto')
        ax3.set_xlabel('$x\;[m]$')
        ax3.set_title('Difference ref-sim')
        fig.tight_layout()
        fig.savefig(os.path.join(outDir, 'refDfa_%s_%s.png' % (rel, suffix)))

        # Fgiure 2 cross and lonprofile
        fig, ax = plt.subplots(ncols=2, figsize=(figW*2, figH), dpi=figReso)
        ax[0].plot(data1[:, ny_loc], 'k', linewidth=lw, label='Reference')
        ax[0].plot(data2[:, ny_loc], 'b--', label='Simulation')
        ax[0].set_xlabel('Location across track [nrows]')
        ax[0].set_ylabel('Result parameter %s' % suffix, fontsize=fs)
        ax[0].set_title('Cross profile at y =  %d' % ny_loc)
        ax[1].plot(data1[nx_loc, :], 'k', linewidth=lw, label='Reference')
        ax[1].plot(data2[nx_loc, :], 'b--', label='Simulation')
        ax[1].set_xlabel('Location along track [ncols]')
        ax[1].set_ylabel('Result parameter %s' % suffix, fontsize=fs)
        ax[1].set_title('Long profile at x =  %d' % nx_loc)

        ax[0].legend()
        ax[1].legend()
        fig.tight_layout()
        fig.savefig(os.path.join(outDir, 'refDfaProfiles_%s_%s' % (rel, suffix)))

        log.info('Figures saved to: %s' % outDir)

        if cfg['FLAGS'].getboolean('showPlot'):
            plt.show()
