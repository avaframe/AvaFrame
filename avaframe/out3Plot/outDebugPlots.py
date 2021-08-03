import numpy as np
import copy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


import avaframe.out3Plot.plotUtils as pU


def plotPartIni(particles, dem):
    header = dem['header']
    x = np.arange(header.ncols) * header.cellsize
    y = np.arange(header.nrows) * header.cellsize
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    ref0, im = pU.NonUnifIm(ax, x, y, dem['areaRaster'], 'x [m]', 'y [m]',
                            extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, norm=None)

    ax.plot(particles['x'], particles['y'], 'or', linestyle='None')
    pU.addColorBar(im, ax, None, 'm²')
    plt.show()


def plotAreaDebug(dem, avapath, Raster):
    ncols = dem['header'].ncols
    nrows = dem['header'].nrows
    cellsize = dem['header'].cellsize
    x = np.arange(ncols) * cellsize
    y = np.arange(nrows) * cellsize
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Release area')
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    ref0, im = pU.NonUnifIm(ax, x, y, Raster, 'x [m]', 'y [m]',
                            extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, norm=None)
    ax.plot(avapath['x'] * cellsize, avapath['y'] * cellsize, 'r', label='release polyline')
    plt.legend()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax)
    plt.show()


def plotRemovePart(xCoord0, yCoord0, header, X, Y, Mask, mask):
    x = np.arange(header.ncols) * header.cellsize
    y = np.arange(header.nrows) * header.cellsize
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Release area')
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    ref0, im = pU.NonUnifIm(ax, x, y, Mask, 'x [m]', 'y [m]',
                            extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, norm=None)
    ax.plot(xCoord0 * header.cellsize, yCoord0 * header.cellsize, 'r', label='release polyline')
    ax.plot(X[mask] * header.cellsize, Y[mask] * header.cellsize, '.b')
    plt.legend()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax)
    plt.show()


def plotPartAfterRemove(points, xCoord0, yCoord0, mask):
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Release area')
    ax.plot(xCoord0, yCoord0, 'r', label='release polyline')
    ax.plot(points['x'], points['y'], '.b')
    ax.plot(points['x'][mask], points['y'][mask], '.g')
    plt.legend()
    plt.show()
