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
from matplotlib import colors as mplCol
import logging
from cmcrameri import cm as cmapCameri

from avaframe.in3Utils import cfgUtils
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


# define settings for colormaps creation
discreteLevels = cfg.getint('discreteLevels')
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
cmapdiv = cmapCameri.broc

# custom colomaps
# cmap based on avaframe logo colors
colorAvaframe = ['#0EF8EA', '#12E4E6', '#28D0DF', '#3CBCD5', '#4AA8C9',
                 '#5595B9', '#5C82A8', '#5F6F95', '#5E5E81', '#5A4D6C', '#523E58', '#483045']
cmapAvaframe = mplCol.ListedColormap(colorAvaframe)
cmapAvaframe.set_bad(color='k')
# add a continuous version
cmapAvaframeCont = mplCol.LinearSegmentedColormap.from_list('cmapAvaframeCont', colorAvaframe, N=256)


# for the choice of the colormaps, check https://www.fabiocrameri.ch/colourmaps/
# and http://hclwizard.org:3000/hclwizard/
# multi sequential colormap for pressure
levP = [1.0, 10.0, 25.0, 50.0]
# Hawaii color map
colorsP = ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"]
cmapP = cmapCameri.hawaii.reversed()

# multi sequential colormap for flow depth
levD = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
# Lajolla color map
colorsD = ["#FCFFC9", "#EBCE7B", "#DE9529", "#BE5A32", "#7F2B3F", "#1D0B14"]
cmapD = cmapCameri.lajolla

# multi sequential colormap for speed
levS = [1, 5, 10, 15, 20, 25, 30]
# Batflow color map
colorsS = ['#FFCEF4', '#FFA7A8', '#C19A1B', '#578B21', '#007054', '#004960', '#201158']
cmapS = cmapCameri.batlow.reversed()

# colormap used if no resType provided
cmapNN = cmapCameri.imola.reversed()

# colormap for probabilities
levProb = [0, 0.25, 0.50, 0.75, 1.]
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
cmapPres['levels'] = levP


cmapDepth = {}
cmapDepth['cmap'] = cmapD
cmapDepth['colors'] = colorsD
cmapDepth['levels'] = levD

cmapSpeed = {}
cmapSpeed['cmap'] = cmapS
cmapSpeed['colors'] = colorsS
cmapSpeed['levels'] = levS


cmapProb = {}
cmapProb['cmap'] = cmapProbmap
cmapProb['colors'] = colorsProb
cmapProb['levels'] = levProb

colorMaps = {'ppr' : cmapPres, 'pfv' : cmapSpeed, 'pfd' : cmapDepth, 'PR' : cmapPres,
             'FV' : cmapSpeed, 'FD' : cmapDepth, 'prob': cmapProb}

cmapDEM = cmapGreys

cmapDEM2 = {}
cmapDEM2['cmap'] = cmapGreys
cmapDEM2['colors'] = colorsS
cmapDEM2['levels'] = None

cmapAimec = cmapAvaframe

cmapVar = {}
# cmapVar['cmap'] = cmapAvaframeCont
cmapVar['colors'] = colorAvaframe
cmapVar['levels'] = None

cmapGreys1 = {}
cmapGreys1['cmap'] = cmapGreys


def makeColorMap(colormapDict, levMin, levMax, continuous=False):
    """Get colormap for plot

    get the colormap, norm, levels... for ploting results


    1) colormapDict = dict (details see Parameters) - depending on the continuous flag and what is
    given in the colormapDict dictionary, the cmap is created and the levelsNew and
    colorsNew are created.
    2) colormapDict = matplotlib cmap - the continuous flag is ignored and a continuous cmap is
    returned with levelsNew set to None.


    Parameters
    ----------
        colormapDict: dict or cmap
            if colormap= dict :
                cmap: matplotlib colormap
                    continuous colormap
                colors: list
                    list of hex or rgb or dec colors
                levels: list
                    levels corresponding to the colors (same number of levels as colors)
        levMin : float
            minimum value of the colormap
        levMax: float
            maximum value of the colormap
        continuous: boolean
            True for continuous cmaps

    Returns
    -------
        cmap: matplotlib colormap
            new color map
        colorsNew: list
            new list of hex or rgb or dec colors
        levelsNew: list
            new levels corresponding to the colors (one more level than colors)
        norm: matplotlib norm
            norm associated to the levels and the colors (includes the vmin and vmax for
            continuous colormaps)
    """

    if type(colormapDict) is matplotlib.colors.LinearSegmentedColormap:
        cmap = colormapDict
        colorsNew = None
        norm = mplCol.Normalize(vmin=levMin, vmax=levMax, clip=False)
        levelsNew = None
    elif continuous:
        # make a continuous color map
        # check if cmap is provided
        if 'cmap' in colormapDict.keys():
            cmap = colormapDict['cmap']
            colorsNew = None
        # check if list of colors is provided
        elif 'colors' in colormapDict.keys():
            colorsNew = colormapDict['colors']
            cmap = mplCol.LinearSegmentedColormap.from_list('myCmap', colorsNew, N=256)
        # Houston ve have a problem
        else:
            message = 'You need a `colors` list or a `cmap` to be able to create the colormap'
            log.error(message)
            raise FileNotFoundError(message)

        norm = mplCol.Normalize(vmin=levMin, vmax=levMax, clip=False)
        levelsNew = None
    else:
        # make a discrete color map
        if 'levels' in colormapDict:
            levels = colormapDict['levels']
        else:
            if 'colors' in colormapDict:
                levels = list(np.linspace(levMin, levMax, len(colormapDict['colors'])+1))
                log.warning('No `levels` list is provided to generate a discrete colormap, \
                            creating %d levels ranging from %.2f to %.2f' %
                            (len(colormapDict['colors']), levMin, levMax))
            else:
                levels = list(np.linspace(levMin, levMax, discreteLevels))
                log.warning('No `levels` list is provided to generate a discrete colormap, \
                            creating %d levels ranging from %.2f to %.2f' %
                            (discreteLevels, levMin, levMax))
        try:
            indEnd = np.where(np.asarray(levels) >= levMax)[0][0]
        except IndexError:
            indEnd = len(levels)
        levelsNew = levels[:indEnd]
        levelsNew.append(levMax)

        # check if list of colors is provided
        if 'colors' in colormapDict.keys():
            if indEnd > len(colormapDict['colors']):
                message = 'Number of levels is not allowed to exceed number of colors'
                log.error(message)
                raise AssertionError(message)
            colors = colormapDict['colors']
            colorsNew = colors[:indEnd]
            colorsNew.append(colors[indEnd-1])
            cmap = mplCol.ListedColormap(colorsNew)
            cmap.set_bad(color='w')
            cmap.set_under(color='w')
        # check if a cmap is provided
        elif 'cmap' in colormapDict.keys():
            cmap = colormapDict['cmap']
            colorsNew = cmap(np.asarray(levelsNew[:-1])/levelsNew[-1])
        # Houston ve have a problem
        else:
            message = 'A `colors` list or a `cmap` is required to create the colormap'
            log.error(message)
            raise FileNotFoundError(message)

        if len(levelsNew) == 1:
            levelsNew = [0] + levelsNew
        norm = mplCol.BoundaryNorm(levelsNew, cmap.N)

    return cmap, colorsNew, levelsNew, norm


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


def saveAndOrPlot(pathDict, cfgFlags, outFileName, fig):
    """
    Receive a plot handle and config and check whether to save and or plot
    closes it afterwards
    """

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()

    if cfgFlags.getboolean('savePlot'):
        outname = os.path.join(pathDict['pathResult'], 'pics', outFileName)
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
    cbar = ax2.figure.colorbar(im, ax=ax2, ticks=ticks, extend=extend, pad=pad, shrink=0.9)
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
    if isinstance(avaDir, str) or isinstance(avaDir, pathlib.Path):
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

    return infoText
