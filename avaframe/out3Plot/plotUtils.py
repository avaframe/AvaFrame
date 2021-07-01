"""
    Plot settings for output figures

    This file is part of Avaframe.
"""

import os
import warnings
import seaborn as sns
import copy
import numpy as np
import matplotlib
import datetime
import pathlib
from matplotlib.image import NonUniformImage
from matplotlib import pyplot as plt
import logging
from cmcrameri import cm as cmapCameri

from avaframe.in3Utils import cfgUtils
from avaframe.out3Plot import makePalette
from avaframe.out3Plot import plotUtils


# create local logger
log = logging.getLogger(__name__)

# Load all input Parameters from config file
# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(plotUtils)
cfgPlotUtils = cfg['UNITS']
cfg = cfg['MAIN']

# define seaborn style and color maps
sns.set(font_scale=1)
sns.set_style("ticks", {'axes.linewidth': 1, 'axes.edgecolor': 'black',
                        'font.family': [cfg['fontFamily']]})


# define figure dimentions
figW = float(cfg['figW'])
figH = float(cfg['figH'])
# define lines and marker properties
markers = cfg['markerStyle']
ms = float(cfg['markerSize'])
matplotlib.rcParams['lines.linewidth'] = float(cfg['lineWidth'])
matplotlib.rcParams['lines.markersize'] = ms
# font size
fs = float(cfg['fontSize'])
matplotlib.rcParams['figure.titlesize'] = cfg['titleSize']
matplotlib.rcParams['figure.dpi'] = float(cfg['figResolution'])
matplotlib.rcParams['figure.autolayout'] = True

matplotlib.rcParams['axes.labelsize'] = cfg['labelSize']
matplotlib.rcParams['axes.linewidth'] = 1.0
matplotlib.rcParams['axes.edgecolor'] = 'lightgrey'
matplotlib.rcParams['axes.labelcolor'] = 'grey'
matplotlib.rcParams['xtick.color'] = 'grey'
matplotlib.rcParams['xtick.major.width'] = 1.0
matplotlib.rcParams['xtick.major.size'] = 3
matplotlib.rcParams['xtick.labelsize'] = cfg['tickLabelSize']
matplotlib.rcParams['ytick.color'] = 'grey'
matplotlib.rcParams['ytick.major.width'] = 1.0
matplotlib.rcParams['ytick.major.size'] = 3
matplotlib.rcParams['ytick.labelsize'] = cfg['tickLabelSize']


# set output extension {png, ps, pdf, svg}
outputFormat = cfg['savefigFormat']
matplotlib.rcParams['savefig.format'] = outputFormat
matplotlib.rcParams['savefig.bbox'] = 'tight'

matplotlib.rcParams['legend.edgecolor'] = 'None'
matplotlib.rcParams['text.usetex'] = cfg.getboolean('usetex')

matplotlib.rcParams['grid.color'] = 'whitesmoke'
matplotlib.rcParams['grid.linestyle'] = ':'
matplotlib.rcParams['grid.linewidth'] = 0.3


############################
# Color maps
############################
# hell white/green to dark blue
cmapGB = copy.copy(sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True))
cmapGB.set_bad(color='k')

cmapReds = copy.copy(matplotlib.cm.Reds)
cmapReds.set_bad(color='k')

cmapBlues = copy.copy(matplotlib.cm.Blues)
cmapBlues.set_bad(color='k')

cmapGreys = copy.copy(matplotlib.cm.get_cmap("Greys"))

cmapPlasma = copy.copy(matplotlib.cm.plasma)
cmapPlasma.set_bad(color='k')

cmapViridis = copy.copy(matplotlib.cm.viridis)
cmapViridis.set_bad(color='k')

# divergent color map
cmapdiv = copy.copy(matplotlib.cm.RdBu_r)  # sns.color_palette("RdBu_r")


# custom colomaps
# cmap based on avaframe logo colors
colorAvaframe = ['#0EF8EA', '#12E4E6', '#28D0DF', '#3CBCD5', '#4AA8C9',
                 '#5595B9', '#5C82A8', '#5F6F95', '#5E5E81', '#5A4D6C', '#523E58', '#483045']
cmapAvaframe = makePalette.get_continuous_cmap(colorAvaframe)
cmapAvaframe.set_bad(color='k')
# add a continuous version
cmapAvaframeCont = makePalette.get_continuous_cmap(colorAvaframe, continuous=True)


# for the choice of the colormaps, check https://www.fabiocrameri.ch/colourmaps/
# and http://hclwizard.org:3000/hclwizard/
# multi sequential colormap for pressure
levP = [1.0, 10.0, 25.0, 50.0, 1000.0]
ticksP = [0., 1.0, 10.0, 25.0, 50.0]
# Hawaii color map
colorsP = ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"]
cmapP = cmapCameri.hawaii.reversed()

# multi sequential colormap for flow depth
levD = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 50.0]
ticksD = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
# Lajolla color map
colorsD = ["#FCFFC9", "#EBCE7B", "#DE9529", "#BE5A32", "#7F2B3F", "#1D0B14"]
cmapD = cmapCameri.lajolla

# multi sequential colormap for speed
levS = [1, 5, 10, 15, 20, 25, 30, 100]
ticksS = [1, 5, 10, 15, 20, 25, 30]
# Batflow color map
colorsS = ['#FFCEF4', '#FFA7A8', '#C19A1B', '#578B21', '#007054', '#004960', '#201158']
cmapS = cmapCameri.batlow.reversed()

# colormap used if no resType provided
cmapNNcmap = cmapCameri.imola.reversed()

# colormap for probabilities
levProb = [0.0, 0.25, 0.5, 0.75, 0.9999999999, 10]
ticksProb = [0, 0.25, 0.50, 0.75, 1.]
# lapaz color map
colorsProb = ['#FEF1F1', '#B2AB96', '#5B8BA3', '#2D5393', '#1A0C64']
cmapProbmap = cmapCameri.lapaz.reversed()

###############################################
# Set colormaps to use
###############################################
# ploting with a descrete (contCmap = continuousCmap = False) or continuous colormap (contCmap = True)?
# if continuous, only the cmap argument in the cmapDictionnary maters
# replace it with the wanted colormap
contCmap = cfg.getboolean('contCmap')
# for pressure
cmapPres = {}
cmapPres['cmap'] = cmapP
cmapPres['colors'] = colorsP
cmapPres['lev'] = levP
cmapPres['ticks'] = ticksP


cmapDepth = {}
cmapDepth['cmap'] = cmapD
cmapDepth['colors'] = colorsD
cmapDepth['lev'] = levD
cmapDepth['ticks'] = ticksD

cmapSpeed = {}
cmapSpeed['cmap'] = cmapS
cmapSpeed['colors'] = colorsS
cmapSpeed['lev'] = levS
cmapSpeed['ticks'] = ticksS


cmapProb = {}
cmapProb['cmap'] = cmapProbmap
cmapProb['colors'] = colorsProb
cmapProb['lev'] = levProb
cmapProb['ticks'] = ticksProb

colorMaps = {'ppr' : cmapPres, 'pfv' : cmapSpeed, 'pfd' : cmapDepth, 'PR' : cmapPres,
             'FV' : cmapSpeed, 'FD' : cmapDepth, 'prob': cmapProb}

cmapNN = {}
cmapNN['cmap'] = cmapNNcmap

cmapDEM = cmapGreys

cmapDEM2 = {}
cmapDEM2['cmap'] = cmapGreys
cmapDEM2['colors'] = colorsS
cmapDEM2['lev'] = None
cmapDEM2['ticks'] = None

cmapAimec = cmapAvaframe

cmapVar = {}
cmapVar['cmap'] = cmapAvaframeCont
cmapVar['colors'] = colorAvaframe
cmapVar['lev'] = None
cmapVar['ticks'] = None

cmapGreys1 = {}
cmapGreys1['cmap'] = cmapGreys
###################################
# shortcut plot functions
###################################
def NonUnifIm(ax, x, y, z, xlab, ylab, **kwargs):
    im = NonUniformImage(ax, **kwargs)
    im.set_data(x, y, z)
    # im.set_clim(vmin=vmin, vmax=vmax)
    ref = ax.images.append(im)
    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min(), y.max()])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ref, im


def saveAndOrPlot(cfgPath, cfgFlags, outFileName, fig):
    """
    Receive a plot handle and config and check whether to save and or plot
    closes it afterwards
    """

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()

    if cfgFlags.getboolean('savePlot'):
        outname = os.path.join(cfgPath['pathResult'], 'pics', outFileName)
        if not os.path.exists(os.path.dirname(outname)):
            os.makedirs(os.path.dirname(outname))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.savefig(outname)

    plt.close(fig)

    return


def constrainPlotsToData(inputData, cellSize, extentOption=False):
    """ constrain inut raster dataset to where there is data plus buffer zone

        Parameters
        -----------
        inputData: numpy array
            raster data to be plotted
        cellSize: float
            cellsize of raster data
        extentOption: bool
            if True rows and columns limits converted to acutal extent in meters

        Returns
        --------
        rowsMin, rowsMax: int or float
            indices where there is data in x direction (or min/max extent in meters)
        colsMin, colsMax: int or float
            indices where there is data in y direction (or min/max extent in meters)
        dataConstrained: numpy array
            constrained array where there is data
        """

    ind = np.where(inputData > 0)
    if len(ind[0]) > 0:
        plotBuffer = int(cfg.getfloat('plotBuffer') / cellSize)
        rowsMin = max(np.amin(ind[0])-plotBuffer, 0)
        rowsMax = min(np.amax(ind[0])+plotBuffer, inputData.shape[0])
        colsMin = max(np.amin(ind[1])-plotBuffer, 0)
        colsMax = min(np.amax(ind[1])+plotBuffer, inputData.shape[1])
    else:
        rowsMin = 0
        rowsMax = inputData.shape[0]
        colsMin = 0
        colsMax = inputData.shape[1]

    if extentOption:
        dataConstrained = inputData[rowsMin:rowsMax+1, colsMin:colsMax+1]
        rowsMinPlot = rowsMin*cellSize
        rowsMaxPlot = (rowsMax+1)*cellSize
        colsMinPlot = colsMin*cellSize
        colsMaxPlot = (colsMax+1)*cellSize

        return rowsMinPlot, rowsMaxPlot, colsMinPlot, colsMaxPlot, dataConstrained
    else:
        return rowsMin, rowsMax, colsMin, colsMax


def addColorBar(im, ax2, ticks, myUnit, title='', extend='neither', pad=0.05):
    '''
    Adds, styles and labels a colorbar to the given image and axes
    '''
    cbar = ax2.figure.colorbar(im, ax=ax2, ticks=ticks, extend=extend, pad=pad)
    cbar.outline.set_visible(False)
    cbar.ax.set_title('[' + myUnit + ']')
    if title != '':
        cbar.set_label(title)


def putAvaNameOnPlot(ax, avaDir):
    '''
    Puts the date and avalanche name (or a list of ava names) in the lower left corner of the given
    matplotlib axes
    '''

    # if avaDir is just a single avaDir or a list of avaDirs
    if isinstance(avaDir, str):
        avaName = pathlib.PurePath(avaDir).name
        infoText = datetime.datetime.now().strftime("%d.%m.%y") + \
               '; ' + str(avaName)
    else:
        infoText = datetime.datetime.now().strftime("%d.%m.%y")
        for ava in avaDir:
            avaName = pathlib.PurePath(ava).name
            infoText = infoText + ';' + str(avaName)

    plt.text(0, 0, infoText, fontsize=8, verticalalignment='bottom',
             horizontalalignment='left', transform=ax.transAxes,
             color='0.6')
