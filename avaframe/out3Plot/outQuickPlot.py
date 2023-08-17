"""
Functions to plot 2D avalanche simulation results as well as comparison plots
between two datasets of identical shape.

"""

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pathlib
import logging
from cmcrameri import cm as cmapCrameri
import matplotlib as mpl

# Local imports
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.out3Plot import statsPlots as sPlot
import avaframe.out3Plot.plotUtils as pU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def generatePlot(dataDict, avaName, outDir, cfg, plotDict, crossProfile=True):
    """ Create comparison plots of two ascii datasets

        This function creates two plots, one plot with four panels with, first dataset, second dataset,
        the absolute difference of the two datasets and the absolute difference capped to a smaller range
        of differences (ppr: +- 100kPa, pft: +-1m, pfv:+- 10ms-1).
        The difference plots also include an insert showing the histogram and the cumulative density function
        of the differences. The second plot shows a cross- and along profile cut of the two datasets.
        The folder and simulation name of the datasets has to be passed to the function.

        Parameters
        ----------
        dataDict : dict
            dictionary with info on both datasets to be plotted:
                name1: str
                    string with the name of the first data set
                name2: str
                    string with the name of the second data set
                data1: 2D numpy arrays
                    raster of the first data sets
                data2: 2D numpy arrays
                    raster of the first second sets
                cellSize: float
                    cell size
                suffix: str
                    optional information about the data type compared ('ppr', 'pft', 'pfv', 'P', 'FV', 'FT', 'Vx'...)
        avaName : str
            name of avalanche
        outDir : pathlib path
            path to dictionary where plots shall be saved to
        cfg : configParser
            cfg['FLAGS'].getboolean('showPlot')
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
    if 'suffix' in dataDict:
        simName = dataDict['simName'] + '_' + dataDict['suffix']
        unit = pU.cfgPlotUtils['unit' + dataDict['suffix']]
        cmapType = pU.colorMaps[dataDict['suffix']]
    else:
        simName = 'compare'
        unit = ''
        cmapType = pU.cmapNN

    # Set dimensions of plots
    ny = data2.shape[0]
    nx = data2.shape[1]
    Ly = ny*cellSize
    Lx = nx*cellSize

    # Location of Profiles
    ny_loc = int(nx * 0.5)
    nx_loc = int(ny * 0.5)

    # Difference between datasets
    dataDiff = data1 - data2
    dataDiff = np.where((data1 == 0) & (data2 == 0), np.nan, dataDiff)
    dataExtend = np.where((data1 == 0) | (data2 == 0), 0, 1)
    diffMax = np.nanmax(dataDiff)
    diffMin = np.nanmin(dataDiff)
    diffMean = np.nanmean(dataDiff)

    minVal = min(np.nanmin(data1), np.nanmin(data2))
    maxVal = max(np.nanmax(data1), np.nanmax(data2))

    # constrain data to where there is data
    rowsMinPlot, rowsMaxPlot, colsMinPlot, colsMaxPlot = pU.constrainPlotsToData(dataExtend, cellSize, extentOption=True)

    # Location of box
    nybox = 0.05
    nxbox = 0.05
    # Plot data
    # Figure 1 shows the result parameter data
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))
    suptitle = fig.suptitle(avaName, fontsize=14, color='0.5')
    ax1 = fig.add_subplot(221)
    cmap, _, ticks, norm = pU.makeColorMap(cmapType, minVal, maxVal, continuous=pU.contCmap)

    cmap.set_bad('w')
    data1P = ma.masked_where(data1 == 0.0, data1)
    im1 = plt.imshow(data1P, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower',
                     aspect=nx/ny, norm=norm)
    ax1.set_xlim([colsMinPlot, colsMaxPlot])
    ax1.set_ylim([rowsMinPlot, rowsMaxPlot])
    pU.addColorBar(im1, ax1, ticks, unit)

    ax1.set_aspect('auto')
    title = str('%s - simulation' % name1)
    ax1.set_title(title)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')

    ax2 = fig.add_subplot(222)
    cmap, _, ticks, norm = pU.makeColorMap(cmapType, minVal, maxVal, continuous=pU.contCmap)

    cmap.set_bad('w')
    data2P = ma.masked_where(data2 == 0.0, data2)
    im2 = plt.imshow(data2P, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower',
                     aspect=nx/ny, norm=norm)
    pU.addColorBar(im2, ax2, ticks, unit)
    ax2.set_xlim([colsMinPlot, colsMaxPlot])
    ax2.set_ylim([rowsMinPlot, rowsMaxPlot])
    ax2.set_aspect('auto')
    ax2.set_xlabel('x [m]')
    title = str('%s - reference' % name2)
    ax2.set_title(title)

    ax3 = fig.add_subplot(223)
    cmap = pU.cmapdiv
    elevMax = np.nanmax(np.abs(dataDiff))
    im3 = plt.imshow(dataDiff, cmap=cmap, clim=(-elevMax, elevMax),
                     extent=[0, Lx, 0, Ly],
                     origin='lower', aspect=nx/ny)
    ax3.set_xlim([colsMinPlot, colsMaxPlot])
    ax3.set_ylim([rowsMinPlot, rowsMaxPlot])
    pU.addColorBar(im3, ax3, None, unit)
    ax3.text(nybox, nxbox, 'Mean: %.2e %s\n Max: %.2e %s\n Min: %.2e %s' %
             (diffMean, unit, diffMax, unit, diffMin, unit),
             horizontalalignment='left', verticalalignment='bottom', transform=ax3.transAxes)
    ax3.set_aspect('auto')
    ax3.set_xlabel('x [m]')
    ax3.set_title('Difference simulation-reference')

    # for difference histogramm - remove dataDiff == 0 values from array
    dataDiffPlot = dataDiff[np.isnan(dataDiff) == False]

    ax4 = fig.add_subplot(224)
    cmap = pU.cmapdiv

    if 'suffix' in dataDict:
        elevMax = pU.cfgPlotUtils.getfloat('elevMax' + dataDict['suffix'])
        ax4.set_title('Difference capped at max difference in %s: +-%.2g %s' % (dataDict['suffix'], elevMax, unit))

    else:
        cutVal = 0.5
        elevMax = cutVal * elevMax
        ax4.set_title('Difference capped at %.1f times max difference: +-%.2f' % (cutVal, elevMax))

    # for difference histogramm - remove dataDiff == 0 values from array
    dataDiffZoom = np.where((dataDiffPlot < -elevMax) | (dataDiffPlot > elevMax), np.nan, dataDiffPlot)
    diffMaxZoom = np.nanmax(dataDiffZoom)
    diffMinZoom = np.nanmin(dataDiffZoom)
    diffMeanZoom = np.nanmean(dataDiffZoom)

    im4 = plt.imshow(dataDiff, cmap=cmap, clim=(-elevMax, elevMax),
                     extent=[0, Lx, 0, Ly],
                     origin='lower', aspect=nx/ny)
    ax4.set_xlim([colsMinPlot, colsMaxPlot])
    ax4.set_ylim([rowsMinPlot, rowsMaxPlot])
    pU.addColorBar(im4, ax4, None, unit, extend='both')
    ax4.set_aspect('auto')
    ax4.set_xlabel('x [m]')

    # if difference is zero dont insert CDF plots
    indDiff = dataDiffPlot > 0
    if indDiff.any():
        axin3 = ax3.inset_axes([0.6, 0.1, 0.4, 0.25])
        axin3.patch.set_alpha(0.0)

        axin4 = ax4.inset_axes([0.6, 0.1, 0.4, 0.25])
        axin4.patch.set_alpha(0.0)

        centiles = sPlot.plotHistCDFDiff(dataDiffPlot, axin4, axin3)
        ax4.text(nybox, nxbox, '95%% centile: %.2e %s\n 99%% centile: %.2e %s' %
                 (centiles[0], unit, centiles[1], unit),
                 horizontalalignment='left', verticalalignment='bottom', transform=ax4.transAxes)

    saveNameDiff = outDir / ('Diff_%s_%s.%s' % (avaName, simName, pU.outputFormat))

    # save and or show figure
    plotPath = pU.saveAndOrPlot({'pathResult': outDir}, saveNameDiff.stem, fig)

    if crossProfile:
        # Fgiure 2 cross and lonprofile
        fig1, ax = plt.subplots(ncols=2, nrows=2, figsize=(pU.figW*2, pU.figH*2))
        suptitle = fig1.suptitle(avaName, fontsize=14, color='0.5')
        ax[0,0].plot(data1[:, ny_loc], 'k', label='Reference')
        ax[0,0].plot(data2[:, ny_loc], 'b--', label='Simulation')
        ax[0,0].set_xlabel('Location across track [nrows]')
        ax[0,0].set_ylabel('Result parameter')
        ax[0,0].set_title('Cross profile at x =  %d' % ny_loc)
        ax[0,1].plot(data1[nx_loc, :], 'k', label='Reference')
        ax[0,1].plot(data2[nx_loc, :], 'b--', label='Simulation')
        ax[0,1].set_xlabel('Location along track [ncols]')
        ax[0,1].set_ylabel('Result parameter')
        ax[0,1].set_title('Long profile at y =  %d' % nx_loc)
        ax[0,0].legend()
        ax[0,1].legend()
        ax[1,0].imshow(data1P, cmap=cmap, extent=[0, nx, 0, ny], origin='lower', aspect=nx/ny)
        ax[1,1].imshow(data2P, cmap=cmap, extent=[0, nx, 0, ny], origin='lower', aspect=nx/ny)
        ax[1,1].axhline(ny - nx_loc, label='profile y=%d' % nx_loc)
        ax[1,0].axvline(ny_loc, label='profile x=%d' % ny_loc)
        ax[1,0].legend()
        ax[1,1].legend()
        saveNameProfile = outDir / ('Profiles_%s_%s.%s' % (avaName, simName, pU.outputFormat))
        # save and or show figure
        plotPath = pU.saveAndOrPlot({'pathResult': outDir}, saveNameProfile.stem, fig1)

    log.info('Figures saved to: %s' % outDir)

    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()

    plotDict['plots'].append(saveNameDiff)
    plotDict['difference'].append(diffMax)
    plotDict['difference'].append(diffMean)
    plotDict['difference'].append(diffMin)
    # stats is the max and min value of the reference
    plotDict['stats'].append(np.amax(data2))
    plotDict['stats'].append(np.amin(data2))
    if 'differenceZoom' in plotDict:
        plotDict['differenceZoom'].append(diffMaxZoom)
        plotDict['differenceZoom'].append(diffMeanZoom)
        plotDict['differenceZoom'].append(diffMinZoom)

    return plotDict


def quickPlotBench(avaDir, simNameRef, simNameComp, refDir, compDir, cfg, suffix):
    """ Plot simulation result and compare to reference solution
        (two raster datasets of identical dimension) and save to
        Outputs/out3Plot within avalanche directoy.

        figure 1, plot raster data for dataset1, dataset2 and their difference, their difference limited to
        specified range, including a histogram and the cumulative density function of the differences

        figure 2, plot cross and longprofiles for both datasets (ny_loc and nx_loc define location of profiles)

        plots are saved to Outputs/out3Plot

        Parameters
        ----------
        avaDir : str or pathlib path
            path to avalanche directory
        simNameRef: str
            name of reference simulation
        simNameComp: str
            name of comparison simulation
        refDir: str or pathlib path
            path to reference file
        compDir: str or pathlib path
            path to comparison file
        cfg : dict
            global configuration settings
        suffix: str
            result type

        Returns
        -------
        plotList : dict
            plot dictionaries (path to plots, min, mean and max difference
            between plotted datasets, max and mean value of reference dataset )

    """

    # Create required directories
    avaDir = fU.checkPathlib(avaDir)
    outDir = avaDir / 'Outputs' / 'out3Plot'
    fU.makeADir(outDir)

    # Initialise plotDictList
    plotList = []

    # Initialise plotList
    plotDict = {'plots': [], 'difference': [], 'stats': []}

    refDir = fU.checkPathlib(refDir)
    compDir = fU.checkPathlib(compDir)
    simRefFile = refDir / (simNameRef + '_' + suffix + '.asc')
    simCompFile = compDir / (simNameComp + '_' + suffix + '.asc')

    if not simRefFile.is_file() or not simCompFile.is_file():
        log.error('File for result type: %s not found' % suffix)

    # Load data
    raster = IOf.readRaster(simCompFile, noDataToNan=True)
    rasterRef = IOf.readRaster(simRefFile, noDataToNan=True)
    dataComp, dataRef = geoTrans.resizeData(raster, rasterRef)
    log.debug('dataset1: %s' % simCompFile)
    log.debug('dataset2: %s' % simRefFile)

    cellSize = rasterRef['header']['cellsize']
    unit = pU.cfgPlotUtils['unit%s' % suffix]

    # Get name of Avalanche
    avaName = avaDir.stem
    # Create dataDict to be passed to generatePlot
    dataDict = {'data1': dataComp, 'data2': dataRef, 'name1': simNameComp + '_' + suffix,
                'name2': simNameRef + '_' + suffix, 'compareType': 'compToRef',
                'simName': simNameComp, 'suffix': suffix, 'cellSize': cellSize, 'unit': unit}
    # Create Plots
    plotDictNew = generatePlot(dataDict, avaName, outDir, cfg, plotDict)

    return plotDictNew


def quickPlotSimple(avaDir, inputDir, cfg):
    """ Plot two raster datasets of identical dimension and difference between two datasets

        figure 1: plot raster data for dataset1, dataset2 and their difference
        figure 2: plot cross and longprofiles for both datasets (ny_loc and nx_loc define location of profiles)
        -plots are saved to Outputs/out3Plot

        Be aware: files are being sorted after getting them from the directory!
        (Important for the differences)

        Parameters
        ----------
        avaDir : str or pathlib path
            path to avalanche directory
        inputDir : str or pathlib path
            path to directory of input data (only 2 raster files allowed)
        cfg: configParser object
            global configuration settings

    """

    avaDir = fU.checkPathlib(avaDir)
    outDir = avaDir / 'Outputs' / 'out3Plot'
    fU.makeADir(outDir)

    # Get name of Avalanche
    avaName = avaDir.stem

    # Load input datasets from input directory
    inputDir = fU.checkPathlib(inputDir)
    datafiles = sorted(list(inputDir.glob('*.asc')))
    datafiles.extend(list(inputDir.glob('*.txt')))

    name1 = datafiles[0].name
    name2 = datafiles[1].name
    log.info('input dataset #1 is %s' % name1)
    log.info('input dataset #2 is %s' % name2)

    # Load data
    raster = IOf.readRaster(datafiles[0], noDataToNan=True)
    rasterRef = IOf.readRaster(datafiles[1], noDataToNan=True)
    data1, data2 = geoTrans.resizeData(raster, rasterRef)
    header = IOf.readASCheader(datafiles[0])
    cellSize = header['cellsize']

    # Create dataDict to be passed to generatePlot
    dataDict = {'data1': data1, 'data2': data2, 'name1': name1,
                'name2': name2, 'compareType': '', 'cellSize': cellSize}

    # Initialise plotList
    plotDict = {'plots': [], 'difference': [], 'stats': []}

    # Create Plots
    plotDictNew = generatePlot(dataDict, avaName, outDir, cfg, plotDict)

    return plotDictNew


def quickPlotOne(avaDir, datafile, cfg, locVal, axis, resType=''):
    """ Plots one raster dataset and a cross profile

        figure 1: plot raster data for dataset and profile
        -plot is saved to Outputs/out3Plot

        Parameters
        ----------
        avaDir : str or pathlib Path
            path to avalanche directory
        datafile : str or pathlib path
            path to data file
        cfg : dict
            configuration including flags for plotting
        locVal : float
            location of cross profile
        resType : str
            result parameter type e.g. 'pft' - optional

    """

    avaDir = fU.checkPathlib(avaDir)
    outDir = avaDir / 'Outputs' / 'out3Plot'
    fU.makeADir(outDir)

    datafile = fU.checkPathlib(datafile)
    name1 = datafile.stem
    log.info('input dataset #1 is %s' % name1)

    # Load data
    raster = IOf.readRaster(datafile, noDataToNan=True)
    data1 = raster['rasterData']
    header = IOf.readASCheader(datafile)
    cellSize = header['cellsize']

    # Create dataDict to be passed to generatePlot
    dataDict = {'data1': data1, 'name1': name1,
                'cellSize': cellSize}

    # Initialise plotList
    plotDict = {'plots': [], 'location': locVal, 'resType': resType, 'axis': axis}

    # Create Plots
    plotList = generateOnePlot(dataDict, outDir, cfg, plotDict)

    return plotDict


def generateOnePlot(dataDict, outDir, cfg, plotDict):
    """ Create plot of ascii dataset

        Parameters
        ----------
        dataDict : dict
            dictionary with info of the dataset to be plotted
        outDir : pathlib path
            path to dictionary where plots shall be saved to
        cfg : dict
            main configuration settings
        plotDict : dict
            dictionary with information about plots, for example release area...

        Returns
        -------
        plotDict : dict
            updated plot dictionary with path to plot
    """

    # Extract info for plotting
    data1 = dataDict['data1']
    name1 = dataDict['name1']
    cellSize = dataDict['cellSize']
    simName = 'Analyse'
    if plotDict['resType'] != '':
        unit = pU.cfgPlotUtils['unit%s' % plotDict['resType']]
        nameRes = pU.cfgPlotUtils['name%s' % plotDict['resType']]
        cmapType = pU.colorMaps[plotDict['resType']]
    else:
        unit = ''
        nameRes = 'Result parameter'
        cmapType = pU.cmapNN

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

    cmap, _, ticks, norm = pU.makeColorMap(cmapType, np.nanmin(data1), np.nanmax(data1), continuous=pU.contCmap)
    cmap.set_bad('w')
    data1P = ma.masked_where(data1 == 0.0, data1)
    im1 = plt.imshow(data1P, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny, norm=norm)
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

    saveNameProfile = outDir / ('Profiles_%s.%s' % (name1, pU.outputFormat))
    # save and or show figure
    plotPath = pU.saveAndOrPlot({'pathResult': outDir}, saveNameProfile.stem, fig)

    log.info('Figures saved to: %s' % outDir)
    plotDict['plots'].append(saveNameProfile)

    return plotDict


def plotContours(contourDict, resType, thresholdValue, pathDict, addLegend=True):
    """ plot contour lines of all transformed fields

        Parameters
        -----------
        contourDict: dict
            dictionary with contourline coordinates
        resType: str
            result type
        thresholdValue: float
            value for contour level
        pathDict: dict
            dictionary with info on project name, ...
        addLegend: bool
            if True add legend to plot
    """

    unit = pU.cfgPlotUtils['unit' + resType]
    fig = plt.figure(figsize=(pU.figW*2, pU.figH))

    # show flow path
    ax1 = fig.add_subplot(121)
    ax1.set_title('%s %s %s contour lines' % (resType, thresholdValue, unit))
    ax1.set_ylabel('x [m]')
    ax1.set_xlabel('y [m]')
    pU.putAvaNameOnPlot(ax1, pathDict['avaDir'])

    # create color map
    nSims = len(contourDict)
    norm = mpl.colors.Normalize(vmin=0, vmax=nSims)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=cmapCrameri.hawaii)
    cmap.set_array([])

    # loop over all sims
    for ind, simName in enumerate(contourDict):
        keyLabel = list(contourDict[simName])[0]
        for key in contourDict[simName]:
            if key == keyLabel and 'line' in key:
                ax1.plot(contourDict[simName][key]['x'], contourDict[simName][key]['y'],
                    c=cmap.to_rgba(ind), label=simName)
            elif 'line' in key:
                ax1.plot(contourDict[simName][key]['x'], contourDict[simName][key]['y'],
                    c=cmap.to_rgba(ind))

    if addLegend:
        ax1.legend()

    # save and or plot
    outFileName = pathDict['plotScenario'] + '_ContourLines'
    pU.saveAndOrPlot(pathDict, outFileName, fig)
