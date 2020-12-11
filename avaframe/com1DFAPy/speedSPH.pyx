import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.ascUtils as IOf
from avaframe.out3Plot.plotUtils import *
import avaframe.com1DFAPy.timeDiscretizations as tD
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.com1DFAPy.SPHfunctions as SPH
from avaframe.in3Utils import cfgUtils
cfg = cfgUtils.getModuleConfig(com1DFA)['GENERAL']
cfgFull = cfgUtils.getModuleConfig(com1DFA)

########################################################################
# CHOOSE YOUR SETUP
##########################################################################
# Choose the snow depth you want to use (h function)

def coputeGrad(Npart, particles, Nx, Ny, Nz, NX, NY):
    Npart = particles['Npart']
    # initialize
    GHX = np.zeros(Npart)
    GHY = np.zeros(Npart)
    GHZ = np.zeros(Npart)
    # loop on particles
    # TcpuSPH = 0
    # Tcpuadd = 0
    for j in range(Npart):
        mass = particles['m'][j]
        # adding lateral force (SPH component)
        # startTime = time.time()

        x = particles['x'][j]
        y = particles['y'][j]
        nx, ny, nz = DFAtls.getNormal(x, y, Nx, Ny, Nz, csz)
        gradhX, gradhY,  gradhZ, _ = SPH.calcGradHSPHVect(particles, j, NX, NY, csz, nx, ny, nz)
        # tcpuSPH = time.time() - startTime
        # TcpuSPH = TcpuSPH + tcpuSPH
        # startTime = time.time()
        GHX[j] = GHX[j] - gradhX / rho
        GHY[j] = GHY[j] - gradhY / rho
        GHZ[j] = GHZ[j] - gradhZ / rho
        # tcpuadd = time.time() - startTime
        # Tcpuadd = Tcpuadd + tcpuadd

    return GHX, GHY, GHZ


def coputeGradcython(int Npart, double[:] m, double[:] X, double[:] Y, double[:, :] Nx, double[:, :] Ny, double[:, :] Nz, particles, int NX, int NY):
    cdef int N = Npart
    cdef double[:] GHX = np.zeros(N, dtype=np.float)
    cdef double[:] GHY = np.zeros(N, dtype=np.float)
    cdef double[:] GHZ = np.zeros(N, dtype=np.float)
    cdef double[:] mass = m
    cdef double[:] x = X
    cdef double[:] y = Y
    cdef double gradhX
    cdef double gradhY
    cdef double gradhZ
    cdef int j
    cdef double xx
    cdef double yy
    # loop on particles
    # TcpuSPH = 0
    # Tcpuadd = 0
    for j in range(N):
      xx = x[j]
      yy = y[j]

      nx, ny, nz = DFAtls.getNormal(xx, yy, Nx, Ny, Nz, csz)
      gradhX, gradhY,  gradhZ, _ = SPH.calcGradHSPHVect(particles, j, NX, NY, csz, nx, ny, nz)
      # tcpuSPH = time.time() - startTime
      # TcpuSPH = TcpuSPH + tcpuSPH
      # startTime = time.time()
      GHX[j] = GHX[j] - gradhX / rho
      GHY[j] = GHY[j] - gradhY / rho
      GHZ[j] = GHZ[j] - gradhZ / rho
      # tcpuadd = time.time() - startTime
      # Tcpuadd = Tcpuadd + tcpuadd
    return GHX, GHY, GHZ



cdef double rho = 200
cdef double csz = 5
