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
sns.set_style("ticks", {'axes.linewidth': 1, 'axes.edgecolor':'black',  'font.family': ['serif']})
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
cmapdiv = copy.copy(matplotlib.cm.RdBu_r) #sns.color_palette("RdBu_r")


# custom colomaps
# cmap based on avaframe logo colors
colorAvaframe = ['#0EF8EA', '#12E4E6', '#28D0DF', '#3CBCD5', '#4AA8C9', '#5595B9', '#5C82A8', '#5F6F95', '#5E5E81', '#5A4D6C', '#523E58', '#483045']
cmapAvaframe = get_continuous_cmap(colorAvaframe)
cmapAvaframe.set_bad(color='k')


# multi sequential colormap for pressure
levP  = [0., 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
        5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 25.0, 30.0, 35.0,
        40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0]
ticksP=[0, 1, 3, 5, 10, 20, 40, 60, 100, 150, 200]
threshold = [1, 3, 5, 10]
h = [140, 180, 250, 300, 350]
colorsP, cmapP, normP = createColorMap(lev, threshold, h, c=[10, 80], l=[10, 80], power=[1, 1], test=True)


###############################################
############ Set colormaps to use #############
###############################################
# for pressure
cmapPres =  cmapP
colorsPres = colorsP
levPres = levP
ticksPres = ticksP


cmapDepth = cmapBlues
cmapSpeed = cmapReds
cmapDEM = cmapGreys
cmapAimec = cmapAvaframe
