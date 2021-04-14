"""
This is a simple function for a quick plot of all peak files at once

This file is part of Avaframe.
"""

import os
import logging
import numpy as np
from matplotlib import pyplot as plt
import glob

import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.makePalette as makePalette

# create local logger
log = logging.getLogger(__name__)


def plotAllPeakFields(avaDir, cfg, cfgFLAGS, modName='com1DFA'):
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
            name of module that has been used to produce data to be plotted - default com1DFA

        Returns
        -------
        plotDict : dict
            dictionary with info on plots, like path to plot
        """

    # Load all infos on simulations
    inputDir = os.path.join(avaDir, 'Outputs', modName, 'peakFiles')
    peakFiles = fU.makeSimDict(inputDir, '', avaDir)

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
        log.debug('now plot %s:' % (fileName))

        # Load data
        raster = IOf.readRaster(fileName)
        data = raster['rasterData']
        data = np.ma.masked_where(data == 0.0, data)
        unit = pU.cfgPlotUtils['unit%s' % peakFiles['resType'][m]]

        # Set extent of peak file
        cellSize = peakFiles['cellSize'][m]
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
        im1 = ax.imshow(data, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
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


def plotAllFields(avaDir, inputDir, outDir, cfg):
    """ Plot all fields within given directory and save to outDir

        Parameters
        ----------
        avaDir : str
            path to avalanche directoy
        inputDir : str
            path to input directoy
        outDir : str
            path to directoy where plots shall be saved to
        cfg : dict
            configuration settings

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
        unit = pU.cfgPlotUtils['unit%s' % cfg['GENERAL']['peakVar']]

        # Figure  shows the result parameter data
        fig = plt.figure(figsize=(pU.figW, pU.figH))
        fig, ax = plt.subplots()
        # choose colormap
        cmap, _, _, norm, ticks = makePalette.makeColorMap(
            pU.cmapPres, np.amin(data), np.amax(data), continuous=pU.contCmap)
        cmap.set_bad('w')
        im1 = ax.imshow(data, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
        pU.addColorBar(im1, ax, ticks, unit)

        title = str('%s' % name)
        ax.set_title(title)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')

        plotName = os.path.join(outDir, '%s.%s' % (name, pU.outputFormat))

        pU.putAvaNameOnPlot(ax, avaDir)

        fig.savefig(plotName)
        plt.close('all')
