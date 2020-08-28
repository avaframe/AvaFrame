"""
    Plot settings for output figures

    This file is part of Avaframe.
"""

import seaborn as sns
from matplotlib import cm
import matplotlib
import cmocean

# define figure dimentions
figW = 6
figH = 6
# define figure resolution
figReso = 150
# define lines and marker properties
lw = 2
ms = 10
markers = ['+', 'o', 'x', '*', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '.',
           '^', 'v', '>', '<', 'p', 'h', '.']
# font size
fs = 12

# define seaborn style and color maps
sns.set(font_scale = 1)
sns.set_style("dark", {'axes.edgecolor': '.8', 'xtick.bottom': True, 'ytick.left': True})
# sns.set_style("dark", {'xtick.bottom': True,
#  'ytick.left': True})

# hell white/green to dark blue
cmapGB = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)
cmapGB.set_bad(color='k')
# hell pink to dark purple
cmapPP = sns.cubehelix_palette(8, as_cmap=True)
cmapPP.set_bad(color='k')

cmapReds = matplotlib.cm.Reds
cmapReds.set_bad(color='k')
cmapBlues = matplotlib.cm.Blues
cmapBlues.set_bad(color='k')

cmapjet = matplotlib.cm.jet
cmapjet.set_bad(color='k')
cmapPlasma = matplotlib.cm.plasma
cmapPlasma.set_bad(color='k')
cmapViridis = matplotlib.cm.viridis
cmapViridis.set_bad(color='k')

cmapIce = cmocean.cm.ice
cmapIce.set_bad(color='k')

cmapDense = cmocean.cm.dense
cmapDense.set_bad(color='k')

# divergent color map
cmapdiv = matplotlib.cm.RdBu_r #sns.color_palette("RdBu_r")


cmapPres = cmapViridis
cmapDepth = cmapBlues
cmapSpeed = cmapReds
