"""

This is a simple function for a quick plot of peak files to be used in reports

This file is part of Avaframe.

"""

import os
import numpy as np
import copy
from matplotlib import pyplot as plt
import seaborn as sns
import avaframe.in3Utils.ascUtils as IOf
from avaframe.out3Plot.plotUtils import *
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.out3Plot.makePalette as makePalette

def plotPeakField(avaDir, simName, resType, unit):
    """ Plot peak field """

    # input dir where peak files are located
    inputDir = os.path.join(avaDir, 'Outputs', 'com1DFA', 'peakFiles')
    name = simName + '_' + resType
    fullName = name + '.asc'
    fileName = os.path.join(inputDir, fullName)
    log.info('now plot %s :' % fileName)

    # Load data
    data = np.loadtxt(fileName, skiprows=6)

    # Output directory
    outDir = os.path.join(avaDir, 'Work', 'log2Report')
    fU.makeADir(outDir)

    # Set extent of peak file
    header = IOf.readASCheader(fileName)
    cellSize = header.cellsize
    ny = data.shape[0]
    nx = data.shape[1]
    Ly = ny*cellSize
    Lx = nx*cellSize

    # Figure 1 shows the result parameter data
    fig = plt.figure(figsize=(figW, figH))
    fig, ax1 = plt.subplots()
    # choose colormap
    cmapDict = {'ppr' : cmapPres, 'pfd' : cmapDepth, 'pv' : cmapSpeed}
    maskedArray = np.ma.masked_where(data == 0, data)
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        cmapPres, 0.0, np.nanmax(maskedArray), continuous=contCmap)

    im1 = ax1.imshow(data, cmap=cmap, extent=[0, Lx, 0, Ly], origin='lower', aspect=nx/ny)
    addColorBar(im1, ax1, ticks, unit)
    title = str('%s' % fullName)
    ax1.set_title(title)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    plotName = os.path.join(outDir, '%s.%s' % (name, outputFormat))
    fig.savefig(plotName)

    return plotName
