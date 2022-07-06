"""
Functions to plot 2D simulation results: plot of all peak files at once

"""

import os
import logging
import numpy as np
from matplotlib import pyplot as plt
import pathlib

import avaframe.out3Plot.plotUtils as pU
import avaframe.in1Data.getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)


def plotAllPeakFields(avaDir, cfgFLAGS, modName, demData=''):
    """ Plot all peak fields and return dictionary with paths to plots
        with DEM in background

        Parameters
        ----------
        avaDir : str
            path to avalanche directoy
        cfgFLAGS : str
            general configuration, required to define if plots saved to reports directoy
        modName : str
            name of module that has been used to produce data to be plotted
        demData: dictionary
            optional - if not the dem in the avaDir/Inputs folder has been used but a different one

        Returns
        -------
        plotDict : dict
            dictionary with info on plots, like path to plot
        """

    # Load all infos on simulations
    avaDir = pathlib.Path(avaDir)
    inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    if demData == '':
        demFile = gI.getDEMPath(avaDir)
        demData = IOf.readRaster(demFile, noDataToNan=True)
        demDataField = demData['rasterData']
    else:
        # check if noDataValue is found and if replace with nans for plotting
        demDataField = np.where(demData['rasterData'] == demData['header']['noDataValue'], np.nan, demData['rasterData'])
    demField = demDataField

    # Output directory
    if cfgFLAGS.getboolean('ReportDir'):
        outDir = avaDir / 'Outputs' / modName / 'reports'
        fU.makeADir(outDir)
    else:
        outDir = avaDir / 'Outputs' / 'out1Peak'
        fU.makeADir(outDir)

    # Initialise plot dictionary with simulation names
    plotDict = {}
    for sName in peakFilesDF['simName']:
        plotDict[sName] = {}

    # Loop through peakFiles and generate plot
    for m in range(len(peakFilesDF['names'])):

        # Load names and paths of peakFiles
        name = peakFilesDF['names'][m]
        fileName = peakFilesDF['files'][m]
        avaName = peakFilesDF['avaName'][m]
        resType = peakFilesDF['resType'][m]
        log.debug('now plot %s:' % (fileName))

        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))

        # only produce a plot if it does not already exists
        # make sure to remove the outpu folder if you want to regenerate the plot
        # this enables to append simulations to an already existing output without regenerating all plits
        if not plotName.is_file():
            # Load data
            raster = IOf.readRaster(fileName, noDataToNan=True)
            data = raster['rasterData']

            # constrain data to where there is data
            cellSize = peakFilesDF['cellSize'][m]
            rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(data, cellSize)
            dataConstrained = data[rowsMin:rowsMax+1, colsMin:colsMax+1]
            demConstrained = demField[rowsMin:rowsMax+1, colsMin:colsMax+1]

            data = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)
            unit = pU.cfgPlotUtils['unit%s' % resType]

            # Set extent of peak file
            ny = data.shape[0]
            nx = data.shape[1]
            Ly = ny*cellSize
            Lx = nx*cellSize

            # Figure  shows the result parameter data
            fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
            # choose colormap
            cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.amin(data), np.amax(data), continuous=pU.contCmap)
            cmap.set_bad(alpha=0)
            # uncomment this to set the under value for discrete cmap transparent
            # cmap.set_under(alpha=0)
            rowsMinPlot = rowsMin*cellSize
            rowsMaxPlot = (rowsMax+1)*cellSize
            colsMinPlot = colsMin*cellSize
            colsMaxPlot = (colsMax+1)*cellSize
            extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]

            # add DEM hillshade with contour lines
            ls, CS = pU.addHillShadeContours(ax, demConstrained, cellSize, extent)

            # add peak field data
            im1 = ax.imshow(data, cmap=cmap, norm=norm, extent=extent, origin='lower', aspect='equal', zorder=4)
            pU.addColorBar(im1, ax, ticks, unit)

            # add title, labels and ava Info
            title = str('%s' % name)
            ax.set_title(title)
            ax.set_xlabel('x [m]')
            ax.set_ylabel('y [m]')
            pU.putAvaNameOnPlot(ax, avaDir)

            if cfgFLAGS.getboolean('showPlot'):
                plt.show()
            fig.savefig(plotName)
            plotPath = pathlib.Path.cwd() / plotName
            plotDict[peakFilesDF['simName'][m]].update({peakFilesDF['resType'][m]: plotPath})
            plt.close(fig)

    return plotDict


def plotAllFields(avaDir, inputDir, outDir, unit='', constrainData=True):
    """ Plot all fields within given directory and save to outDir

        Parameters
        ----------
        avaDir : str
            path to avalanche directoy
        inputDir : str
            path to input directoy
        outDir : str
            path to directoy where plots shall be saved to
        unit: str
            unit of result type

        """

    # Load all infos on simulations
    inputDir = pathlib.Path(inputDir)
    outDir = pathlib.Path(outDir)
    if outDir.is_dir() is False:
        # create out dir if not already existing
        outDir.mkdir()
    peakFiles = list(inputDir.glob('*.asc'))

    # Loop through peakFiles and generate plot
    for filename in peakFiles:

        # Load data
        raster = IOf.readRaster(filename)
        data = raster['rasterData']
        data = np.ma.masked_where(data == 0.0, data)
        name = filename.stem

        # get header info for file writing
        header = raster['header']
        cellSize = header['cellsize']

        # Set extent of peak file
        ny = data.shape[0]
        nx = data.shape[1]
        Ly = ny*cellSize
        Lx = nx*cellSize

        # Figure  shows the result parameter data
        fig = plt.figure(figsize=(pU.figW, pU.figH))
        fig, ax = plt.subplots()
        # choose colormap
        cmap, _, ticks, norm  = pU.makeColorMap(pU.cmapNN, np.amin(data), np.amax(data), continuous=pU.contCmap)
        cmap.set_bad('w')

        if constrainData:
            rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(data, cellSize)
            dataConstrained = data[rowsMin:rowsMax+1, colsMin:colsMax+1]
            data = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)
            rowsMinPlot = rowsMin*cellSize
            rowsMaxPlot = (rowsMax+1)*cellSize
            colsMinPlot = colsMin*cellSize
            colsMaxPlot = (colsMax+1)*cellSize
            im1 = ax.imshow(data, cmap=cmap, norm=norm, extent=[colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot], origin='lower', aspect=nx/ny)
        else:
            im1 = ax.imshow(data, cmap=cmap, norm=norm, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)

        pU.addColorBar(im1, ax, ticks, unit)

        title = str('%s' % name)
        ax.set_title(title)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')

        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))

        pU.putAvaNameOnPlot(ax, avaDir)

        fig.savefig(plotName)
        plt.close(fig)
