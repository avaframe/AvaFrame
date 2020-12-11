# import pyximport
# pyximport.install()

import logging
import time
import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as geoTrans
from avaframe.out3Plot.plotUtils import *
import avaframe.com1DFAPy.timeDiscretizations as tD
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.SPHfunctions as SPH
import avaframe.com1DFAPy.frictionLaws as fricLaws

from SPHfunctionsCython import *
# import avaframe.in2Trans.shpConversion as shpConv
# import avaframe.in2Trans.ascUtils as IOf
# from avaframe.DFAkernel.setParam import *

# create local logger
log = logging.getLogger(__name__)

slopeAnglex = 40*math.pi/180
slopeAngley = 20*math.pi/180


def createDEM(Lx, Ly, csz):
    # define grid
    ncols = np.int(Lx/csz + 1)
    nrows = np.int(Ly/csz + 1)
    header = IOf.cASCheader()
    header.ncols = ncols
    header.nrows = nrows
    header.xllcenter = 0
    header.yllcenter = 0
    header.cellsize = csz
    dem = {}
    dem['header'] = header
    X = np.linspace(0, Lx, ncols)
    Y = np.linspace(0, Ly, nrows)
    XX, YY = np.meshgrid(X, Y)
    ZZ, _, _, _ = Sfunction(XX, YY, Lx, Ly)
    dem['rasterData'] = ZZ
    return dem


def createParticles(dem, lam, rho):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    Lx = ncols*csz-1
    Ly = nrows*csz-1
    dx = 1/lam
    dy = 1/lam
    # define particles
    nx = np.int(Lx * lam)
    ny = np.int(Ly * lam)
    Npart = nx*ny
    x = np.linspace(2*csz, Lx-2*csz, nx+1, endpoint=False)[1:]
    y = np.linspace(2*csz, Ly-2*csz, ny+1, endpoint=False)[1:]
    Xpart, Ypart = np.meshgrid(x, y)
    Xpart = Xpart.flatten()
    Ypart = Ypart.flatten()
    # Xpart = Xpart + (np.random.rand(Npart) - 0.5) * dx * coef
    # Ypart = Ypart + (np.random.rand(Npart) - 0.5) * dy * coef
    # adding z component
    Zpart, sx, sy, area = Sfunction(Xpart, Ypart, Lx, Ly)
    Hpart, _, _, _ = Hfunction(Xpart, Ypart, Zpart)
    Mpart = Hpart * dx * dy * area * rho
    # create dictionnary to store particles properties
    particles = {}
    particles['Npart'] = Npart
    particles['mTot'] = np.sum(Mpart)
    particles['x'] = Xpart
    particles['y'] = Ypart
    particles['z'] = Zpart
    particles['s'] = np.zeros(np.shape(Ypart))
    particles['ux'] = np.zeros(np.shape(Ypart))
    particles['uy'] = np.zeros(np.shape(Ypart))
    particles['uz'] = np.zeros(np.shape(Ypart))

    particles['m'] = Mpart
    particles['h'] = Hpart
    particles['kineticEne'] = 0
    particles['potentialEne'] = 0

    PP = np.zeros((nrows, ncols))
    fields = {}
    fields['pv'] = PP
    fields['ppr'] = PP
    fields['pfd'] = PP
    fields['V'] = PP
    fields['P'] = PP
    fields['FD'] = PP

    particles = SPH.getNeighboursVect(particles, dem)

    return particles, fields


def Hfunction(x, y, z):
    # h = x*x/1000 + 1
    GHx = 2*x*y/5000
    GHy = x*x/5000
    h = x*x*y/5000 + 1
    # GHx = np.ones(np.shape(x))/10
    # h = np.ones(np.shape(x))
    # GHx = np.zeros(np.shape(x))
    # h = y/50 + 1
    # GHy = np.ones(np.shape(x))/50
    # GHy = np.zeros(np.shape(x))
    GHz = np.zeros(np.shape(x))
    return h, GHx, GHy, GHz


# Choose the surface shape you want to use
def Sfunction(x, y, Lx, Ly):
    # horizontal plane
    # Z = np.ones(np.shape(x))
    # area = 1

    # plane
    Z = 0*x*np.tan(slopeAnglex) + 0*y*np.tan(slopeAngley)
    sx = 0*np.tan(slopeAnglex)*np.ones(np.shape(x))
    sy = 0*np.tan(slopeAngley)*np.ones(np.shape(x))
    area = np.sqrt(1 + 0*(np.tan(slopeAnglex))*(np.tan(slopeAnglex)) + 0*(np.tan(slopeAngley))*(np.tan(slopeAngley)))*np.ones(np.shape(x))
    # quadratic surface
    # Z = 1/(Lx*Lx) * x*x + 2/(Ly*Ly) * y*y
    # area = np.sqrt(1 + (2/(Lx*Lx)*x)*(2/(Lx*Lx)*x) + (2*2/(Ly*Ly)*y)*(2*2/(Ly*Ly)*y))
    # sx = 2/(Lx*Lx)*x
    # sy = 2*2/(Ly*Ly)*y
    # other function
    # Z = Lx*np.exp(1/(Lx*Lx) * x*x) + 1/(Ly) * y*y
    # sx = (Lx*2/(Lx*Lx)*x)*np.exp(1/(Lx*Lx) * x*x)
    # sy = (2/(Ly)*y)
    # area = np.sqrt(1 + (Lx*2/(Lx*Lx)*x)*np.exp(1/(Lx*Lx) * x*x)*(Lx*2/(Lx*Lx)*x)*np.exp(1/(Lx*Lx) * x*x) + (2/(Ly)*y)*(2/(Ly)*y))
    return Z, sx, sy, area
