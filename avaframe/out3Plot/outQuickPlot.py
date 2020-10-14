"""

This function creates a quick plot of datasets and option for comparison
between two datasets of identical shape. Also plots cross and longprofiles.

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
from avaframe.out3Plot.plotUtils import *

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def generatePlot(dataDict, avaName, outDir, cfg, plotList):
    """ Create plots of two ascii datasets that shall be compared  """

    # Extract info for plotting
    data1 = dataDict['data1']
    data2 = dataDict['data2']
    name1 = dataDict['name1']
    name2 = dataDict['name2']
    cellSize = dataDict['cellSize']
    if dataDict['compareType'] == 'compToRef':
        simName = dataDict['simName'] + '_' + dataDict['suffix']
        unit = dataDict['unit']
    else:
        simName = 'compare'
        unit = ''

    # Set dimensions of plots
    ny = data1.shape[0]
    nx = data1.shape[1]
    Ly = ny*cellSize
    Lx = nx*cellSize

    # Location of Profiles
    ny_loc = int(nx * 0.5)
    nx_loc = int(ny * 0.5)

    # Difference between datasets
    dataDiff = data1 - data2
    diffMax = np.amax(dataDiff)
    diffMin = np.amin(dataDiff)
    diffMean = np.mean(dataDiff)

    # Location of box
    nybox = int(nx * 0.2)
    nxbox = int(ny * 0.2)

    # Plot data
    # Figure 1 shows the result parameter data
    fig = plt.figure(figsize=(figW*3, figH))
    suptitle = fig.suptitle(avaName, fontsize=14, color='0.5')
    ax1 = fig.add_subplot(131)
    cmap = cmapGB
    im1 = plt.imshow(data1, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    fig.colorbar(im1, ax=ax1)
    ax1.set_aspect('auto')
    title = str('%s' % name1)
    ax1.set_title(title)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')

    ax2 = fig.add_subplot(132)
    im2 = plt.imshow(data2, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    fig.colorbar(im2, ax=ax2)
    ax2.set_aspect('auto')
    ax2.set_xlabel('x [m]')
    title = str('%s' % name2)
    ax2.set_title(title)

    ax3 = fig.add_subplot(133)
    cmap = cmapdiv
    im3 = plt.imshow(dataDiff, cmap=cmap,
                     extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    fig.colorbar(im3, ax=ax3)
    ax3.text(nybox, nxbox, 'Mean: %.2f %s\n Max: %.2f %s\n Min: %.2f %s' %
            (diffMean, unit, diffMax, unit, diffMin, unit),
             bbox=dict(boxstyle="square", ec='white', fc='white'),
             horizontalalignment='left', verticalalignment='bottom')
    ax3.set_aspect('auto')
    ax3.set_xlabel('x [m]')
    ax3.set_title('Difference ref-sim')
    fig.savefig(os.path.join(outDir, 'Diff_%s_%s.%s' % (avaName, simName, outputFormat)))

    # Fgiure 2 cross and lonprofile
    fig, ax = plt.subplots(ncols=2, figsize=(figW*2, figH))
    suptitle = fig.suptitle(avaName, fontsize=14, color='0.5')
    ax[0].plot(data1[:, ny_loc], 'k', label='Reference')
    ax[0].plot(data2[:, ny_loc], 'b--', label='Simulation')
    ax[0].set_xlabel('Location across track [nrows]')
    ax[0].set_ylabel('Result parameter')
    ax[0].set_title('Cross profile at y =  %d' % ny_loc)
    ax[1].plot(data1[nx_loc, :], 'k', label='Reference')
    ax[1].plot(data2[nx_loc, :], 'b--', label='Simulation')
    ax[1].set_xlabel('Location along track [ncols]')
    ax[1].set_ylabel('Result parameter')
    ax[1].set_title('Long profile at x =  %d' % nx_loc)
    ax[0].legend()
    ax[1].legend()
    fig.savefig(os.path.join(outDir, 'Profiles_%s_%s.%s' % (avaName, simName, outputFormat)))

    log.info('Figures saved to: %s' % outDir)

    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()

    plotList.append(os.path.join(outDir, 'Diff_%s_%s.%s' % (avaName, simName, outputFormat)))

    return plotList


def quickPlot(avaDir, suffix, val, parameter, cfg, cfgPlot):
    """ Plot simulation result and compare to reference solution (two raster datasets of identical dimension):

        Inputs:

        avaDir          avalanche directory
        suffix          result parameter abbreviation (e.g. 'ppr')
        com1DFAOutput   folder where results to be plotted are located
        simName         entres or null for simulation type
        cfgPlot         configuration settings for plots
        cfg             global configuration settings

        Outputs:

        figure 1: plot raster data for dataset1, dataset2 and their difference
        figure 2: plot cross and longprofiles for both datasets (ny_loc and nx_loc define location of profiles)
        -plots are saved to Outputs/out3Plot
    """

    # Create required directories
    workDir = os.path.join(avaDir, 'Work', 'out3Plot')
    fU.makeADir(workDir)
    outDir = os.path.join(avaDir, 'Outputs', 'out3Plot')
    fU.makeADir(outDir)

    # Initialise plotList
    plotList = []

    # Setup input from com1DFA
    fU.getDFAData(avaDir, workDir, suffix)

    # Get data from reference run
    fU.getRefData(avaDir, workDir, suffix)

    # prepare data
    if parameter == 'Mu' or parameter == 'RelTh':
        data = fU.makeSimDict(workDir, parameter, avaDir)
    else:
        data = fU.makeSimDict(workDir, '', avaDir)

    cellSize = data['cellSize'][0]
    unit = cfgPlot['PLOT']['unit%s' % suffix]

    # Count the number of release areas
    relAreas = set(data['releaseArea'])

    for rel in relAreas:
        # get list of indices of files that are of correct simulation type and result paramete
        indSuffix = [-9999, -9999]
        for m in range(len(data['files'])):
            if data['resType'][m] == suffix and data['releaseArea'][m] == rel and data[parameter][m] == val:
                if data['modelType'][m] == 'dfa':
                    indSuffix[0] = m
                elif data['modelType'][m] == cfgPlot['PLOT']['refModel']:
                    indSuffix[1] = m

        # Load data
        data1 = np.loadtxt(data['files'][indSuffix[0]], skiprows=6)
        data2 = np.loadtxt(data['files'][indSuffix[1]], skiprows=6)
        log.info('dataset1: %s' % data['files'][indSuffix[0]])
        log.info('dataset2: %s' % data['files'][indSuffix[1]])

        # Get name of Avalanche
        avaName = data['avaName'][indSuffix[0]]

        # Create dataDict to be passed to generatePlot
        dataDict = {'data1': data1, ' data2': data2, 'name1': data['names'][indSuffix[0]],
                    'name2': data['names'][indSuffix[1]], 'compareType': 'compToRef',
                    'simName': data['simName'][indSuffix[0]], 'suffix': suffix, 'cellSize': cellSize, 'unit': unit}

        # Create Plots
        plotList = generatePlot(dataDict, avaName, outDir, cfg, plotList)

        return plotList


def quickPlotSimple(avaDir, inputDir, cfg):
    """ Plot two raster datasets of identical dimension:

        Inputs:

        avaDir          avalanche directory
        inputDir        where input data can be found

        Outputs:

        figure 1: plot raster data for dataset1, dataset2 and their difference
        figure 2: plot cross and longprofiles for both datasets (ny_loc and nx_loc define location of profiles)
        -plots are saved to Outputs/out3SimpPlot
    """

    outDir = os.path.join(avaDir, 'Outputs', 'out3Plot')
    fU.makeADir(outDir)

    # Get name of Avalanche
    avaName = os.path.basename(avaDir)

    # Load input datasets from input directory
    datafiles = glob.glob(inputDir+os.sep + '*.asc')
    datafiles.extend(glob.glob(inputDir+os.sep + '*.txt'))

    name1 = os.path.basename(datafiles[0])
    name2 = os.path.basename(datafiles[1])
    log.info('input dataset #1 is %s' % name1)
    log.info('input dataset #2 is %s' % name2)

    # Load data
    data1 = np.loadtxt(datafiles[0], skiprows=6)
    data2 = np.loadtxt(datafiles[1], skiprows=6)
    header = IOf.readASCheader(datafiles[0])
    cellSize = header.cellsize

    # Create dataDict to be passed to generatePlot
    dataDict = {'data1': data1, 'data2': data2, 'name1': name1,
                'name2': name2, 'compareType': '', 'cellSize': cellSize}

    # Create Plots
    generatePlot(dataDict, avaName, outDir, cfg)
