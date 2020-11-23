"""

This is a simple function for a quick plot of all peak files at once

This file is part of Avaframe.

"""

import os
import numpy as np
from matplotlib import pyplot as plt
import glob

from avaframe.out3Plot.plotUtils import *
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.makePalette as makePalette


def plotAllPeakFields(avaDir, cfg, cfgFLAGS, modName='com1DFA'):
    """ Plot all peak fields and return dictionary with paths to plots """

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
        data = np.loadtxt(fileName, skiprows=6)
        unit = cfg['REPORT']['unit%s' % peakFiles['resType'][m]]

        # Set extent of peak file
        cellSize = peakFiles['cellSize'][m]
        ny = data.shape[0]
        nx = data.shape[1]
        Ly = ny*cellSize
        Lx = nx*cellSize

        # Figure  shows the result parameter data
        fig = plt.figure(figsize=(figW, figH))
        fig, ax = plt.subplots()
        # choose colormap
        cmap, _, _, norm, ticks = makePalette.makeColorMap(
            cmapPres, np.amin(data), np.amax(data), continuous=contCmap)

        im1 = ax.imshow(data, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
        addColorBar(im1, ax, ticks, unit)

        title = str('%s' % name)
        ax.set_title(title)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')

        plotName = os.path.join(outDir, '%s.%s' % (name, outputFormat))

        plotUtils.putAvaNameOnPlot(ax, avaDir)

        fig.savefig(plotName)
        plotPath = os.path.join(os.getcwd(), plotName)
        plotDict[peakFiles['simName'][m]].update({peakFiles['resType'][m]: plotPath})

    return plotDict


def plotAllFields(avaDir, inputDir, outDir, cfg=''):
    """ Plot all fields within given directory and save to outDir"""

    # Load all infos on simulations
    peakFiles = glob.glob(inputDir+os.sep + '*.asc')

    # create out dir if not allready existing
    fU.makeADir(outDir)

    # Loop through peakFiles and generate plot
    for m in range(len(peakFiles)):

        # Load data
        data = np.loadtxt(peakFiles[m], skiprows=6)
        name = os.path.splitext(os.path.basename(peakFiles[m]))[0]

        # get header info for file writing
        header = IOf.readASCheader(peakFiles[m])
        cellSize = header.cellsize

        # Set extent of peak file
        ny = data.shape[0]
        nx = data.shape[1]
        Ly = ny*cellSize
        Lx = nx*cellSize
        unit = cfg['unit']

        # Figure  shows the result parameter data
        fig = plt.figure(figsize=(figW, figH))
        fig, ax = plt.subplots()
        # choose colormap
        cmap, _, _, norm, ticks = makePalette.makeColorMap(
            cmapPres, np.amin(data), np.amax(data), continuous=contCmap)

        im1 = ax.imshow(data, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
        addColorBar(im1, ax, ticks, unit)

        title = str('%s' % name)
        ax.set_title(title)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')

        plotName = os.path.join(outDir, '%s.%s' % (name, outputFormat))

        plotUtils.putAvaNameOnPlot(ax, avaDir)

        fig.savefig(plotName)
