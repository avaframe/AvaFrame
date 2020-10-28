import logging
import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.ascUtils as IOf
from avaframe.out3Plot.plotUtils import *
from avaframe.DFAkernel.setParam import *

debugPlot = False


def getRelease(demHeader, releaseLine):
    # adim and center dem and release line
    ncols = demHeader.ncols
    nrows = demHeader.nrows
    xllc = demHeader.xllcenter
    yllc = demHeader.yllcenter
    csz = demHeader.cellsize
    xCoord = (releaseLine['x'] - xllc) / csz
    yCoord = (releaseLine['y'] - yllc) / csz
    # get the raster corresponding to the release polygon
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


def initializeRelease(relRaster, dem):
    header = dem['header']
    csz = header.cellsize
    S = csz * csz
    # initialize
    partPerCell = np.zeros(np.shape(relRaster), dtype=np.int64)
    Npart = 0
    Xpart = np.empty(0)
    Ypart = np.empty(0)
    Mpart = np.empty(0)
    Dpart = np.empty(0)
    InCell = np.empty((0, 2), int)
    # find all non empty cells
    indY, indX = np.nonzero(relRaster)
    # loop on non empty cells
    for indx, indy in zip(indX, indY):
        # number of particles for this cell
        depth = relRaster[indy][indx]
        V = S * depth
        mass = V * rho
        nPart = np.ceil(mass / massPerPart).astype('int')
        Npart += nPart
        partPerCell[indy][indx] = nPart
        mPart = mass / nPart
        xpart = csz * (np.random.rand(nPart) - 0.5 + indx)
        ypart = csz * (np.random.rand(nPart) - 0.5 + indy)
        Xpart = np.append(Xpart, xpart)
        Ypart = np.append(Ypart, ypart)
        Mpart = np.append(Mpart, mPart * np.ones(nPart))
        Dpart = np.append(Dpart, depth * np.ones(nPart))
        InCell = np.append(InCell, np.tile(np.array([indx, indy]), (nPart, 1)), axis=0)

    # create dictionnary to store particles properties
    particles = {}
    particles['Npart'] = Npart
    particles['mTot'] = np.sum(Mpart)
    particles['x'] = Xpart
    particles['y'] = Ypart
    # adding z component
    particles, _, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')
    particles['m'] = Mpart
    particles['d'] = Dpart
    particles['InCell'] = InCell
    particles['ux'] = np.zeros(np.shape(Xpart))
    particles['uy'] = np.zeros(np.shape(Xpart))
    particles['uz'] = np.zeros(np.shape(Xpart))

    if debugPlot:
        ncols = header.ncols
        nrows = header.nrows
        x = np.arange(ncols) * csz
        y = np.arange(nrows) * csz
        fig, ax = plt.subplots(figsize=(figW, figH))
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        ref0, im = NonUnifIm(ax, x, y, relRaster, 'x [m]', 'y [m]',
                             extent=[x.min(), x.max(), y.min(), y.max()],
                             cmap=cmap, norm=None)
        ax.plot(Xpart, Ypart, 'or', linestyle='None')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        plt.show()

    return particles


def computeBodyForce():
    1


def getNormalVect(z, csz):
    n, m = np.shape(z)
    # first and last row, first and last column are inacurate
    # normal calculation with 4 triangles
    Nx = np.ones((n, m))
    # (Zl - Zr) * csz
    Nx[1:n-1, 1:m-1] = (z[1:n-1, 0:m-2] - z[1:n-1, 2:m]) / csz
    Ny = np.ones((n, m))
    # (Zd - Zu) * csz
    Ny[1:n-1, 1:m-1] = (z[0:n-2, 1:m-1] - z[2:n, 1:m-1]) / csz
    Nz = 2 * np.ones((n, m))

    # # normal calculation with 6 triangles
    # Nx = np.ones((n, m))
    # # (2*(Zl - Zr) - Zur + Zdl + Zu - Zd) * csz
    # Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m])
    #                     - z[2:n, 2:m] + z[0:n-2, 0:m-2]
    #                     + z[2:n, 1:m-1] - z[0:n-2, 1:m-1]) / csz
    # Ny = np.ones((n, m))
    # # (2*(Zd - Zu) + Zur + Zdl - Zu - Zl) * csz
    # Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1])
    #                     + z[2:n, 2:m] + z[0:n-2, 0:m-2]
    #                     - z[2:n, 1:m-1] - z[1:n-1, 0:m-2]) / csz
    # Nz = 6 * np.ones((n, m))

    Nx, Ny, Nz = normalize(Nx, Ny, Nz)

    return Nx, Ny, Nz


def normalize(x, y, z):
    norm = np.sqrt(x*x + y*y + z*z)
    xn = x / norm
    yn = y / norm
    zn = z / norm
    return xn, yn, zn
