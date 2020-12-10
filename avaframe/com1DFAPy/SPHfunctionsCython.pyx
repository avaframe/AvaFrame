"""
    function related to SPH calculations in com1DFA
"""

# Load modules
import copy
import logging
import math
import numpy as np
from libc.math cimport sqrt

# Local imports
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.SPHfunctions as SPH


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)

# how to compute SPH gradient:
# 1) like in SAMOS AT just project on the plane
# 2) project on the plane, compute the gradient in the local coord sys related
# to the local plane and flow direction (orthogonal coord sys)
# enables to choose earth pressure coefficients
# 3) project on the plane, compute the gradient in the local coord sys related
# to the local plane (tau1, tau2, n) non orthogonal coord sys
cdef int SPHoption = 2

cdef double rho = 200
cdef double csz = 5

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

def coputeGradcython(int Npart, double[:] m, double[:] X, double[:] Y, double[:] Z, double[:, :] Nx, double[:, :] Ny, double[:, :] Nz, particles, int NX, int NY):
    cdef int N = Npart
    cdef double[:] GHX = np.zeros(N, dtype=np.float)
    cdef double[:] GHY = np.zeros(N, dtype=np.float)
    cdef double[:] GHZ = np.zeros(N, dtype=np.float)
    cdef double[:] mass = m
    cdef double[:] x = X
    cdef double[:] y = Y
    cdef double[:] z = Z
    cdef double[:] indPartInCell = particles['indPartInCell']
    cdef double[:] partInCell = particles['partInCell']
    cdef double[:] indx = particles['InCell'][:, 0]
    cdef double[:] indy = particles['InCell'][:, 1]
    cdef double gradhX, gradhY, gradhZ
    cdef int j
    cdef double xx, yy, zz
    # loop on particles
    # TcpuSPH = 0
    # Tcpuadd = 0
    for j in range(N):
      xx = x[j]
      yy = y[j]
      # nx, ny, nz = DFAtls.getNormal(xx, yy, Nx, Ny, Nz, csz)
      gradhX, gradhY,  gradhZ, _ = calcGradHSPH(x, y, z, m, indPartInCell, partInCell, indx[j], indy[j], j, NX, NY, csz)
      # tcpuSPH = time.time() - startTime
      # TcpuSPH = TcpuSPH + tcpuSPH
      # startTime = time.time()
      GHX[j] = GHX[j] - gradhX / rho
      GHY[j] = GHY[j] - gradhY / rho
      GHZ[j] = GHZ[j] - gradhZ / rho
      # tcpuadd = time.time() - startTime
      # Tcpuadd = Tcpuadd + tcpuadd
    return GHX, GHY, GHZ


def calcGradHSPH(double[:] X, double[:] Y, double[:] Z, double[:] M, double[:] indPartInCell, double[:] partInCell, int indx, int indy, int j, int ncols, int nrows, double csz):
    """ Compute gradient of Flow Depth using SPH (for loop implementation)

    Parameters
    ----------
    particles : dict
    j: int
        index of particle under consideration
    ncols: int
        number of columns of the DEM
    nrows: int
        number of rows of the DEM
    csz  : float
        cellsize of the DEM

    Returns
    -------
    gradhX: float
        x coordinate of the gradient of the flow depth at particle j location
    gradhY: float
        x coordinate of the gradient of the flow depth at particle j location
    gradhZ: float
        x coordinate of the gradient of the flow depth at particle j location
    L: 1D numpy array
        index of particles within the kernel function radius
    """
    # SPH kernel
    # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
    cdef double rKernel = csz
    cdef double facKernel = 10.0 / (3.1415 * rKernel*rKernel*rKernel*rKernel*rKernel)
    cdef double dfacKernel = -3.0 * facKernel
    cdef double x = X[j]
    cdef double y = Y[j]
    cdef double z = Z[j]
    cdef double dx, dy, dz, r, hr, dwdr, massl
    # With loop
    cdef double gradhX = 0
    cdef double gradhY = 0
    cdef double gradhZ = 0
    # startTime = time.time()
    L = np.empty((0), dtype=int)
    # check if we are on the bottom ot top row!!!
    cdef int lInd = -1
    cdef int rInd = 2
    cdef int ic, n, p, l, imax, imin, iPstart, iPend
    if indy == 0:
        lInd = 0
    if indy == nrows - 1:
        rInd = 1
    for n in range(lInd, rInd):
        ic = (indx - 1) + ncols * (indy + n)
        # make sure not to take particles from the other edge
        imax = max(ic, ncols * (indy + n))
        imin = min(ic+3, ncols * (indy + n + 1))
        iPstart = int(indPartInCell[imax])
        iPend = int(indPartInCell[imin])
        # loop on all particles in neighbour boxes
        for p in range(iPstart, iPend):
            # index of particle in neighbour box
            l = int(partInCell[p])
            if j != l:
                L = np.append(L, l)
                dx = X[l] - x
                dy = Y[l] - y
                dz = Z[l] - z
                r = norm(dx, dy, dz)
                if r < 0.001 * rKernel:
                    # impose a minimum distance between particles
                    r = 0.001 * rKernel
                if r < rKernel:
                    hr = rKernel - r
                    dwdr = dfacKernel * hr * hr
                    massl = M[l]
                    gradhX = gradhX + massl * dwdr * dx / r
                    gradhY = gradhY + massl * dwdr * dy / r
                    gradhZ = gradhZ + massl * dwdr * dz / r
    # print((time.time() - startTime)/1)

    return gradhX, gradhY,  gradhZ, L

def norm(double x, double y, double z):
  """ Compute the Euclidean norm of the vector (x, y, z).

  (x, y, z) can be numpy arrays.

  Parameters
  ----------
      x: numpy array
          x component of the vector
      y: numpy array
          y component of the vector
      z: numpy array
          z component of the vector

  Returns
  -------
      norme: numpy array
          norm of the vector
  """
  cdef double norme
  norme = sqrt(x*x + y*y + z*z)
  return norme
