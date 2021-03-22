"""

This function creates a quick plot of datasets and option for comparison
between two datasets of identical shape. Also plots cross and longprofiles.

This file is part of Avaframe.

"""

import matplotlib.pyplot as plt
from avaframe.in3Utils import fileHandlerUtils as fU
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import os
import logging
import glob

# Local imports
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.statsPlots as statsPlots
import avaframe.out3Plot.makePalette as makePalette
import avaframe.out3Plot.plotUtils as pU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def generatePlot(dataDict, avaName, outDir, cfg, plotDict):
    """ Create plots of two ascii datasets that shall be compared

        Parameters
        ----------
        dataDict : dict
            dictionary with info on both datasets to be plotted
        avaName : str
            name of avalanche
        outDir : str
            path to dictionary where plots shall be saved to
        cfg : dict
            main configuration settings
        plotDict : dict
            dictionary with information about plots, for example release area

        Returns
        -------
        plotDict : dict
            updated plot dictionary with info about e.g. min, mean, max of difference between datasets
    """

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
    ny = data2.shape[0]
    nx = data2.shape[1]
    Ly = ny*cellSize
    Lx = nx*cellSize

    # Location of Profiles
    ny_loc = int(nx * 0.5)
    nx_loc = int(ny * 0.5)

    # Difference between datasets
    dataDiff = np.nan_to_num(data1 - data2)
    diffMax = np.amax(dataDiff)
    diffMin = np.amin(dataDiff)
    diffMean = np.mean(dataDiff)

    # Location of box
    nybox = int(nx * 0.2)
    nxbox = int(ny * 0.2)

    # Plot data
    # Figure 1 shows the result parameter data
    fig = plt.figure(figsize=(pU.figW*3, pU.figH))
    suptitle = fig.suptitle(avaName, fontsize=14, color='0.5')
    ax1 = fig.add_subplot(131)
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        pU.cmapPres, np.amin(data1), np.amax(data1), continuous=pU.contCmap)

    im1 = plt.imshow(data1, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    pU.addColorBar(im1, ax1, ticks, unit)

    ax1.set_aspect('auto')
    title = str('%s' % name1)
    ax1.set_title(title)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')

    ax2 = fig.add_subplot(132)
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        pU.cmapPres, np.amin(data2), np.amax(data2), continuous=pU.contCmap)

    im2 = plt.imshow(data2, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    pU.addColorBar(im2, ax2, ticks, unit)

    ax2.set_aspect('auto')
    ax2.set_xlabel('x [m]')
    title = str('%s' % name2)
    ax2.set_title(title)

    ax3 = fig.add_subplot(133)
    cmap = pU.cmapdiv
    elev_max = np.max(np.abs(dataDiff))
    im3 = plt.imshow(dataDiff, cmap=cmap, clim=(-elev_max, elev_max),
                     extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    fig.colorbar(im3, ax=ax3)
    ax3.text(nybox, nxbox, 'Mean: %.2f %s\n Max: %.2f %s\n Min: %.2f %s' %
             (diffMean, unit, diffMax, unit, diffMin, unit),
             bbox=dict(boxstyle="square", ec='white', fc='white'),
             horizontalalignment='left', verticalalignment='bottom')
    ax3.set_aspect('auto')
    ax3.set_xlabel('x [m]')
    ax3.set_title('Difference ref-sim')

    # for difference histogramm - remove dataDiff == 0 values from array
    dataDiffPlot = dataDiff[dataDiff != 0]
    axin2 = ax3.inset_axes([0.75, 0.1, 0.25, 0.25])
    axin2.patch.set_alpha(0.0)
    axin2.hist(dataDiffPlot, bins=30)
    axin2.get_yaxis().set_ticks([])

    fig.savefig(os.path.join(outDir, 'Diff_%s_%s.%s' % (avaName, simName, pU.outputFormat)))

    # Fgiure 2 cross and lonprofile
    fig, ax = plt.subplots(ncols=2, figsize=(pU.figW*2, pU.figH))
    suptitle = fig.suptitle(avaName, fontsize=14, color='0.5')
    ax[0].plot(data1[:, ny_loc], 'k', label='Reference')
    ax[0].plot(data2[:, ny_loc], 'b--', label='Simulation')
    ax[0].set_xlabel('Location across track [nrows]')
    ax[0].set_ylabel('Result parameter')
    ax[0].set_title('Cross profile at x =  %d' % ny_loc)
    ax[1].plot(data1[nx_loc, :], 'k', label='Reference')
    ax[1].plot(data2[nx_loc, :], 'b--', label='Simulation')
    ax[1].set_xlabel('Location along track [ncols]')
    ax[1].set_ylabel('Result parameter')
    ax[1].set_title('Long profile at y =  %d' % nx_loc)
    ax[0].legend()
    ax[1].legend()
    fig.savefig(os.path.join(outDir, 'Profiles_%s_%s.%s' % (avaName, simName, pU.outputFormat)))

    log.info('Figures saved to: %s' % outDir)

    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()

    plotDict['plots'].append(os.path.join(outDir, 'Diff_%s_%s.%s' % (avaName, simName, pU.outputFormat)))
    plotDict['difference'].append(diffMax)
    plotDict['difference'].append(diffMean)
    plotDict['difference'].append(diffMin)
    plotDict['stats'].append(np.amax(data2))
    plotDict['stats'].append(np.amin(data2))

    plt.close('all')

    return plotDict


def quickPlot(avaDir, testDir, suffix, val, parameter, cfg, cfgPlot, rel='', simType='null', comModule='com1DFA'):
    """ Plot simulation result and compare to reference solution (two raster datasets of identical dimension) and save to
        Outputs/out3Plot within avalanche directoy

        figure 1: plot raster data for dataset1, dataset2 and their difference
        figure 2: plot cross and longprofiles for both datasets (ny_loc and nx_loc define location of profiles)
        -plots are saved to Outputs/out3Plot

        Parameters
        ----------
        avaDir : str
            path to avalanche directory
        suffix : str
            result parameter abbreviation (e.g. 'ppr')
        val : str
            value of parameter
        parameter : str
            parameter that is used to filter simulation results within folder, for example, symType, parameter variation, etc.
        cfg : dict
            global configuration settings
        cfgPlot : dict
            configuration settings for plots, required for flag if plots shall be shown or only saved
        rel : str
            optional - name of release area scenarios
        simType : str
            optional - simulation type null or entres

        Returns
        -------
        plotList : list
            list of plot dictionaries (path to plots, min, mean and max difference
            between plotted datasets, max and mean value of reference dataset )

    """

    # Create required directories
    workDir = os.path.join(avaDir, 'Work', 'out3Plot')
    fU.makeADir(workDir)
    outDir = os.path.join(avaDir, 'Outputs', 'out3Plot')
    fU.makeADir(outDir)

    # Initialise plotDictList
    plotList = []

    # Setup input from com1DFA
    fU.getDFAData(avaDir, workDir, suffix, comModule=comModule)

    # Get data from reference run
    fU.getRefData(testDir, workDir, suffix)

    # prepare data
    if parameter == 'Mu' or parameter == 'RelTh':
        data = fU.makeSimDict(workDir, parameter, avaDir)
    else:
        data = fU.makeSimDict(workDir, '', avaDir)

    cellSize = data['cellSize'][0]
    unit = pU.cfgPlotUtils['unit%s' % suffix]

    # check if release Area and simType area provided
    if rel != '':
        relAreas = [rel]
    else:
        # Count the number of release areas
        relAreas = set(data['releaseArea'])
    if parameter == 'simType':
        simType = val

    for rel in relAreas:

        # Initialise plotList
        plotDict = {'relArea' : rel, 'plots': [], 'difference': [], 'stats': []}

        # get list of indices of files that are of correct simulation type and result paramete
        indSuffix = [-9999, -9999]
        for m in range(len(data['files'])):
            if data['resType'][m] == suffix and data['releaseArea'][m] == rel and data[parameter][m] == val and data['simType'][m] == simType:
                if data['modelType'][m] == 'dfa':
                    indSuffix[0] = m
                elif data['modelType'][m] == cfgPlot['PLOT']['refModel']:
                    indSuffix[1] = m

        # Load data
        raster = IOf.readRaster(data['files'][indSuffix[0]])
        rasterRef = IOf.readRaster(data['files'][indSuffix[1]])
        data1, data2 = geoTrans.resizeData(raster, rasterRef)
        log.debug('dataset1: %s' % data['files'][indSuffix[0]])
        log.debug('dataset2: %s' % data['files'][indSuffix[1]])

        # Get name of Avalanche
        avaName = data['avaName'][indSuffix[0]]

        # Create dataDict to be passed to generatePlot
        dataDict = {'data1': data1, 'data2': data2, 'name1': data['names'][indSuffix[0]],
                    'name2': data['names'][indSuffix[1]], 'compareType': 'compToRef',
                    'simName': data['simName'][indSuffix[0]], 'suffix': suffix, 'cellSize': cellSize, 'unit': unit}

        # Create Plots
        plotDictNew = generatePlot(dataDict, avaName, outDir, cfg, plotDict)
        plotList.append(plotDictNew)

    return plotList


def quickPlotSimple(avaDir, inputDir, cfg):
    """ Plot two raster datasets of identical dimension and difference between two datasets

        figure 1: plot raster data for dataset1, dataset2 and their difference
        figure 2: plot cross and longprofiles for both datasets (ny_loc and nx_loc define location of profiles)
        -plots are saved to Outputs/out3Plot

        Parameters
        ----------
        avaDir : str
            path to avalanche directory
        inputDir : str
            path to directory of input data (only 2 raster files allowed)

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
    raster = IOf.readRaster(datafiles[0])
    rasterRef = IOf.readRaster(datafiles[1])
    data1, data2 = geoTrans.resizeData(raster, rasterRef)
    header = IOf.readASCheader(datafiles[0])
    cellSize = header.cellsize

    # Create dataDict to be passed to generatePlot
    dataDict = {'data1': data1, 'data2': data2, 'name1': name1,
                'name2': name2, 'compareType': '', 'cellSize': cellSize}

    # Initialise plotList
    plotDict = {'plots': [], 'difference': [], 'stats': []}

    # Create Plots
    plotList = generatePlot(dataDict, avaName, outDir, cfg, plotDict)


def quickPlotOne(inputDir, datafile, cfg, locVal, axis, resType=''):
    """ Plots one raster dataset and a cross profile

        figure 1: plot raster data for dataset and profile
        -plot is saved to Outputs/out3Plot

        Parameters
        ----------
        inputDir : str
            path to directory of input data (takes first dataset)
        datafile : str
            path to data file
        cfg : dict
            configuration including flags for plotting
        locVal : float
            location of cross profile
        resType : str
            result parameter type e.g. 'pfd' - optional

    """

    outDir = os.path.join(inputDir, 'out3Plot')
    fU.makeADir(outDir)

    name1 = os.path.basename(datafile)
    log.info('input dataset #1 is %s' % name1)

    # Load data
    data1 = np.loadtxt(datafile, skiprows=6)
    header = IOf.readASCheader(datafile)
    cellSize = header.cellsize

    # Create dataDict to be passed to generatePlot
    dataDict = {'data1': data1, 'name1': name1,
                'cellSize': cellSize}

    # Initialise plotList
    plotDict = {'plots': [], 'location': locVal, 'resType': resType, 'axis': axis}

    # Create Plots
    plotList = generateOnePlot(dataDict, outDir, cfg, plotDict)



def generateOnePlot(dataDict, outDir, cfg, plotDict):
    """ Create plots of ascii dataset

        Parameters
        ----------
        dataDict : dict
            dictionary with info on both datasets to be plotted
        avaName : str
            name of avalanche
        outDir : str
            path to dictionary where plots shall be saved to
        cfg : dict
            main configuration settings
        plotDict : dict
            dictionary with information about plots, for example release area

        Returns
        -------
        plotDict : dict
            updated plot dictionary with info about e.g. min, mean, max of difference between datasets
    """

    # Extract info for plotting
    data1 = dataDict['data1']
    name1 = dataDict['name1']
    cellSize = dataDict['cellSize']
    simName = 'Analyse'
    if plotDict['resType'] != '':
        unit = pU.cfgPlotUtils['unit%s' % plotDict['resType']]
        nameRes = pU.cfgPlotUtils['name%s' % plotDict['resType']]
    else:
        unit = ''
        nameRes = 'Result parameter'

    # Set dimensions of plots
    ny = data1.shape[0]
    nx = data1.shape[1]
    Ly = ny*cellSize
    Lx = nx*cellSize
    axis = plotDict['axis']

    # Location of Profiles
    location = float(plotDict['location'])
    if axis == 'x':
        nx_loc = int(location / cellSize)
    elif axis == 'y':
        ny_loc = int(location / cellSize)
    else:
        log.error('Not an axis, please provide axis of profile')

    # Plot data
    # Figure 1 shows the result parameter data
    fig = plt.figure(figsize=(pU.figW*2, pU.figH))
    suptitle = fig.suptitle(name1, fontsize=14, color='0.5')
    ax1 = fig.add_subplot(121)
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        pU.cmapPres, np.amin(data1), np.amax(data1), continuous=pU.contCmap)

    im1 = plt.imshow(data1, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    pU.addColorBar(im1, ax1, ticks, unit)

    ax1.set_aspect('auto')
    title = str('%s' % name1)
    ax1.set_title(title)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')

    ax3 = fig.add_subplot(122)
    if axis == 'x':
        ax3.plot(data1[:, nx_loc], 'k', label='Reference')
    else:
        ax3.plot(data1[ny_loc, :], 'k', label='Reference')

    ax3.set_xlabel('Location across track [nrows]')
    ax3.set_ylabel('%s [%s]' % (nameRes, unit))
    if axis == 'x':
        ax3.set_title('Profile at x ~ %d [%s] (%d)' % (location, unit, nx_loc))
    else:
        ax3.set_title('Profile at y ~ %d [%s] (%d)' % (location, unit, ny_loc))

    fig.savefig(os.path.join(outDir, 'Profiles_%s.%s' % (name1, pU.outputFormat)))


    log.info('Figures saved to: %s' % outDir)

    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()

    plotDict['plots'].append(os.path.join(outDir, 'Profiles_%s.%s' % (name1, pU.outputFormat)))

    plt.close('all')

    return plotDict
