import logging
import numpy as np
from scipy import sparse
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.ascUtils as IOf
from avaframe.out3Plot.plotUtils import *

debugPlot = True


def getRelease(dem, releaseLine):
    header = dem['header']
    # adim and center dem and release line
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    cellsize = header.cellsize
    xCoord = (releaseLine['x'] - xllc) / cellsize
    yCoord = (releaseLine['y'] - yllc) / cellsize

    relRaster = geoTrans.poly2maskSimple(xCoord, yCoord, ncols, nrows)

    if debugPlot:
        fig, ax = plt.subplots(figsize=(figW, figH))
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        im = plt.imshow(relRaster, cmap, origin='lower')
        ax.plot(xCoord, yCoord, 'k')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        plt.show()

    return relRaster


def initializeRelease(relRaster):
    indX, indY = np.nonzero(relRaster)
    for indx, indy in zip(indX, indY):
        
