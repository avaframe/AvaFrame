import seaborn as sns
from matplotlib import cm
import copy
import matplotlib

# define figure dimentions
figW = 12
figH = 8
# define figure resolution
figReso = 150
# define lines and marker properties
lw = 2
ms = 10
markers = ['+', 'o', 'x', '*', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '.',
           '^', 'v', '>', '<', 'p', 'h', '.']
# font size
fs = 10

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

cmapjet = copy.copy(matplotlib.cm.jet)
cmapjet.set_bad(color='k')
cmapPlasma = copy.copy(matplotlib.cm.plasma)
cmapPlasma.set_bad(color='k')
cmapViridis = copy.copy(matplotlib.cm.viridis)
cmapPlasma.set_bad(color='k')

# divergent color map
cmapdiv = sns.color_palette("RdBu_r")


cmap1 = cmapViridis
cmap2 = cmapViridis
