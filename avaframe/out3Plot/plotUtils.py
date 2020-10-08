"""
    Plot settings for output figures

    This file is part of Avaframe.
"""

import os
import seaborn as sns
import copy
import matplotlib
import datetime
import pathlib
from matplotlib.image import NonUniformImage
from matplotlib import pyplot as plt
import cmocean
import logging

from avaframe.in3Utils import cfgUtils
from avaframe.out3Plot import makePalette
from avaframe.out3Plot import plotUtils


# create local logger
log = logging.getLogger(__name__)

# Load all input Parameters from config file
# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(plotUtils)
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

cmapIce = copy.copy(cmocean.cm.ice)
cmapIce.set_bad(color='k')

cmapDense = copy.copy(cmocean.cm.dense)
cmapDense.set_bad(color='k')

# divergent color map
cmapdiv = copy.copy(matplotlib.cm.RdBu_r)  # sns.color_palette("RdBu_r")


# custom colomaps
# cmap based on avaframe logo colors
colorAvaframe = ['#0EF8EA', '#12E4E6', '#28D0DF', '#3CBCD5', '#4AA8C9',
                 '#5595B9', '#5C82A8', '#5F6F95', '#5E5E81', '#5A4D6C', '#523E58', '#483045']
cmapAvaframe = makePalette.get_continuous_cmap(colorAvaframe)
cmapAvaframe.set_bad(color='k')


# multi sequential colormap for pressure
levP = [0., 1.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0, 500.0, 1000.0]
ticksP = [0., 1.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0, 500.0]
colorsP = ['#F0FF94', '#DFF244', '#E4C82F', '#D77411', '#C5491E',
           '#C51F2E', '#A30A54', '#232D5F', '#102D5B']
cmapP = makePalette.get_continuous_cmap(colorsP, continuous=True)
# multi sequential colormap for flow depth
levD = [0., 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 50.0]
ticksD = [0., 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
colorsD = ['#E1F88F', '#D7C410', '#D58B15', '#C51F2E', '#A30A54', '#232D5F', '#232D5F']
cmapD = makePalette.get_continuous_cmap(colorsD, continuous=True)

# multi sequential colormap for speed
levS = [0., 1, 5, 10, 15, 20, 25, 30, 35, 100]
ticksS = [0., 1, 5, 10, 15, 20, 25, 30, 35]
colorsS = ['#F0FF94', '#DFF244', '#E4C82F', '#D77411', '#C5491E', '#BC3334', '#A70753', '#5B2967', '#102D5B']
# '#851D62' or '#A70753' and '#5B2967'
cmapS = makePalette.get_continuous_cmap(colorsS, continuous=True)

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


cmapDEM = cmapGreys

cmapAimec = cmapAvaframe


###################################
# shortcut plot functions
###################################
def NonUnifIm(ax, x, y, z, xlab, ylab, **kwargs):
    im = NonUniformImage(ax, **kwargs)
    im.set_data(x, y, z)
    ref = ax.images.append(im)
    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min(), y.max()])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ref, im


def saveAndOrPlot(cfgPath, cfgFlags, outFileName, fig):
    """
    Receive a plot handle and config and check whether to save and or plot
    """

    if cfgFlags.getboolean('savePlot'):
        outname = os.path.join(cfgPath['pathResult'], 'pics', outFileName)
        if not os.path.exists(os.path.dirname(outname)):
            os.makedirs(os.path.dirname(outname))
        fig.savefig(outname)

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()

    plt.close(fig)

    return


def addColorBar(im, ax2, ticks, myUnit):
    '''
    Adds, styles and labels a colorbar to the given image and axes
    '''
    cbar = ax2.figure.colorbar(im, ax=ax2, ticks=ticks)
    cbar.outline.set_visible(False)
    cbar.ax.set_title(myUnit)


def putAvaNameOnPlot(ax, avaDir):
    '''
    Puts the date and avalanche name in the lower left corner of the given
    matplotlib axes
    '''
    avaName = pathlib.PurePath(avaDir).name

    infoText = datetime.datetime.now().strftime("%d.%m.%y") + \
               '; ' + str(avaName)
    plt.text(0, 0, infoText, fontsize=8, verticalalignment='bottom',
             horizontalalignment='left', transform=ax.transAxes,
             color='0.6')
