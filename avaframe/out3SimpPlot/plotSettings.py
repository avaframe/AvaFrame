"""
    Plot settings for output figures

    This file is part of Avaframe.
"""

import seaborn as sns
from matplotlib import cm
import copy
import matplotlib
import cmocean
import copy

from avaframe.out3SimpPlot.makePalette import *


# define seaborn style and color maps
sns.set(font_scale=1)
sns.set_style("ticks", {'axes.linewidth': 1, 'axes.edgecolor': 'black',  'font.family': ['sans-serif']})
# print(sns.axes_style())


# define figure dimentions
figW = 6
figH = 6
# define lines and marker properties
lw = 1.5
ms = 5
markers = 'o'
matplotlib.rcParams['lines.linewidth'] = lw
matplotlib.rcParams['lines.markersize'] = ms
# font size
fs = 12
matplotlib.rcParams['figure.titlesize'] = 'xx-large'
matplotlib.rcParams['axes.labelsize'] = 'x-large'
# set output extension {png, ps, pdf, svg}
matplotlib.rcParams["savefig.format"] = 'png'
# define figure resolution (dpi)
matplotlib.rcParams['figure.dpi'] = 150

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['figure.autolayout'] = True


############################
###### Color maps ##########
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
cmapAvaframe = get_continuous_cmap(colorAvaframe)
cmapAvaframe.set_bad(color='k')


# multi sequential colormap for pressure
levP = [0., 1.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0, 500.0, 1000.0]
ticksP = [0., 1.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0, 500.0]
# colorsP = ['#E1F88F', '#D7C410', '#D58B15', '#CE5F1E', '#C51F2E', '#A30A54', '#232D5F']
colorsP = ['#F0FF94', '#DFF244', '#E4C82F', '#D77411', '#C5491E', '#BC3334', '#A70753', '#5B2967', '#102D5B']
#C84618 or #CE5F1E
cmapP  = get_continuous_cmap(colorsP, continuous=True)

# multi sequential colormap for flow depth
levD = [0., 0.5, 1.0, 2.0, 3.0, 5.0, 50.0]
ticksD = [0., 0.5, 1.0, 2.0, 3.0, 5.0]
colorsD  = ['#E1F88F', '#D7C410', '#D58B15', '#C51F2E', '#A30A54', '#232D5F']
cmapD  = get_continuous_cmap(colorsD, continuous=True)

# multi sequential colormap for speed
levS = [0., 1, 5, 10, 15, 20, 25, 30, 35, 100]
ticksS = [0., 1, 5, 10, 15, 20, 25, 30, 35]
colorsS  = ['#F0FF94', '#DFF244', '#E4C82F', '#D77411', '#C5491E', '#BC3334', '#A70753', '#5B2967', '#102D5B']
# '#851D62' or '#A70753' and '#5B2967'
cmapS  = get_continuous_cmap(colorsS, continuous=True)

###############################################
############ Set colormaps to use #############
###############################################
# ploting with a descrete or continuous colormap?
# if continuous, only the cmap argument in the cmapDictionnary maters
# replace it with the wanted colormap
contCmap = True
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
