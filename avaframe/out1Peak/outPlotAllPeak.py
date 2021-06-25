"""
This is a simple function for a quick plot of all peak files at once

This file is part of Avaframe.
"""

import os
import logging
import numpy as np
from matplotlib import pyplot as plt
import numpy.ma as ma
import glob

import avaframe.out3Plot.plotUtils as pU
import avaframe.in1Data.getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.makePalette as makePalette

# create local logger
log = logging.getLogger(__name__)


def plotAllPeakFields(avaDir, cfg, cfgFLAGS, modName):
    """ Plot all peak fields and return dictionary with paths to plots

        Parameters
        ----------
        avaDir : str
            path to avalanche directoy
        cfg : dict
            configuration used to perform simulations
        cfgFLAGS : str
            general configuration, required to define if plots saved to reports directoy
        modName : str
            name of module that has been used to produce data to be plotted

        Returns
        -------
        plotDict : dict
            dictionary with info on plots, like path to plot
        """

    # Load all infos on simulations
    inputDir = os.path.join(avaDir, 'Outputs', modName, 'peakFiles')
    peakFiles, _ = fU.makeSimDict(inputDir, '', avaDir)

    demFile = gI.getDEMPath(avaDir)
    demData = IOf.readRaster(demFile)
    demField = demData['rasterData']

    # Output directory
    if cfgFLAGS.getboolean('ReportDir'):
        outDir = os.path.join(avaDir, 'Outputs', modName, 'reports')
        fU.makeADir(outDir)
    else:
        outDir = os.path.join(avaDir, 'Outputs', 'out1Peak')
        fU.makeADir(outDir)

    # Initialise plot dictionary with simulation names
    plotDict = {}
    for sName in peakFiles['simName']:
        plotDict[sName] = {}

    # Loop through peakFiles and generate plot
    for m in range(len(peakFiles['names'])):

        # Load names and paths of peakFiles
        name = peakFiles['names'][m]
        fileName = peakFiles['files'][m]
        avaName = peakFiles['avaName'][m]
        resType = peakFiles['resType'][m]
        log.debug('now plot %s:' % (fileName))

        # Load data
        raster = IOf.readRaster(fileName)
        data = raster['rasterData']

        # constrain data to where there is data
        cellSize = peakFiles['cellSize'][m]
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
        fig = plt.figure(figsize=(pU.figW, pU.figH))
        fig, ax = plt.subplots()
        # choose colormap
        cmap, _, _, norm, ticks = makePalette.makeColorMap(
            pU.colorMaps[resType], np.amin(data), np.amax(data), continuous=pU.contCmap)
        cmap.set_bad(alpha=0)
        rowsMinPlot = rowsMin*cellSize
        rowsMaxPlot = (rowsMax+1)*cellSize
        colsMinPlot = colsMin*cellSize
        colsMaxPlot = (colsMax+1)*cellSize
        im0 = ax.imshow(demConstrained, cmap='Greys', extent=[colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot], origin='lower', aspect='equal')
        im1 = ax.imshow(data, cmap=cmap, norm=norm, extent=[colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot], origin='lower', aspect='equal')
        pU.addColorBar(im1, ax, ticks, unit)

        title = str('%s' % name)
        ax.set_title(title)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')

        plotName = os.path.join(outDir, '%s.%s' % (name, pU.outputFormat))

        pU.putAvaNameOnPlot(ax, avaDir)

        fig.savefig(plotName)
        if cfgFLAGS.getboolean('showPlot'):
            plt.show()
        plotPath = os.path.join(os.getcwd(), plotName)
        plotDict[peakFiles['simName'][m]].update({peakFiles['resType'][m]: plotPath})
        plt.close('all')

    return plotDict


def plotAllFields(avaDir, inputDir, outDir, unit=''):
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
    peakFiles = glob.glob(inputDir+os.sep + '*.asc')

    # create out dir if not allready existing
    fU.makeADir(outDir)

    # Loop through peakFiles and generate plot
    for filename in peakFiles:

        # Load data
        raster = IOf.readRaster(filename)
        data = raster['rasterData']
        data = np.ma.masked_where(data == 0.0, data)
        name = os.path.splitext(os.path.basename(filename))[0]

        # get header info for file writing
        header = raster['header']
        cellSize = header.cellsize

        # Set extent of peak file
        ny = data.shape[0]
        nx = data.shape[1]
        Ly = ny*cellSize
        Lx = nx*cellSize

        # Figure  shows the result parameter data
        fig = plt.figure(figsize=(pU.figW, pU.figH))
        fig, ax = plt.subplots()
        # choose colormap
        cmap, _, _, norm, ticks = makePalette.makeColorMap(
            pU.cmapPres, np.amin(data), np.amax(data), continuous=pU.contCmap)
        cmap.set_bad('w')
        im1 = ax.imshow(data, cmap=cmap, norm=norm, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
        pU.addColorBar(im1, ax, ticks, unit)

        title = str('%s' % name)
        ax.set_title(title)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')

        plotName = os.path.join(outDir, '%s.%s' % (name, pU.outputFormat))

        pU.putAvaNameOnPlot(ax, avaDir)

        fig.savefig(plotName)
        plt.close('all')
