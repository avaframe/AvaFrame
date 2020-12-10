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

def coputeGradcython(int Npart, double[:] m, double[:] X, double[:] Y, double[:] Z, double[:, :] Nx, double[:, :] Ny, double[:, :] Nz, long[:] indPartInCell, long[:] partInCell, long[:] indX, long[:] indY, int ncols, int nrows):
    cdef int N = X.shape[0]
    cdef double rKernel = csz
    cdef double facKernel = 10.0 / (3.1415 * rKernel*rKernel*rKernel*rKernel*rKernel)
    cdef double dfacKernel = -3.0 * facKernel
    cdef double[:] GHX = np.zeros(N)
    cdef double[:] GHY = np.zeros(N)
    cdef double[:] GHZ = np.zeros(N)
    cdef double[:] mass = m
    cdef double gradhX, gradhY, gradhZ, uMag
    cdef int j
    cdef double xx, yy, zz, ux, uy, uz, uxOrtho, uyOrtho, uzOrtho
    cdef double dx, dy, dz, r, hr, dwdr, massl
    cdef int lInd = -1
    cdef int rInd = 2
    cdef long indx
    cdef long indy
    cdef int ic, n, p, l, imax, imin, iPstart, iPend
    # With loop
    # loop on particles
    # TcpuSPH = 0
    # Tcpuadd = 0
    for j in range(N):
      xx = X[j]
      yy = Y[j]
      zz = Z[j]
      gradhX = 0
      gradhY = 0
      gradhZ = 0
      indx = indX[j]
      indy = indY[j]
      ux = UX[j]
      uy = UY[j]
      uz = UZ[j]
      uMag = norm(ux, uy, uz)
      if uMag < 0.1:
          # ax = 1
          # ay = 0
          ux = 1
          uy = 0
          uz = -(1*nx + 0*ny) / nz
          # uz = -(ax*nx + ay*ny) / nz
          ux, uy, uz = normalize(ux, uy, uz)
      else:
          # TODO check if direction is non zero, if it is define another u1 direction
          ux, uy, uz = normalize(ux, uy, uz)

      uxOrtho, uyOrtho, uzOrtho = croosProd(nx, ny, nz, ux, uy, uz)
      uxOrtho, uyOrtho, uzOrtho = normalize(uxOrtho, uyOrtho, uzOrtho)

      v1 = np.array([[ux, uy, uz]])
      v2 = np.array([[uxOrtho, uyOrtho, uzOrtho]])
      v3 = np.array([[nx, ny, nz]])
      # build the transformation matrix from (e1, e2, e3) to (ex, ey, ez)
      M = np.concatenate((v1.T, v2.T), axis=1)
      M = np.concatenate((M, v3.T), axis=1)
      # compute the transformation matrix from (ex, ey, ez) to (e1, e2, e3)
      MM1 = M.T  # because M is orthogonal, it inverse is its transpose !!! np.linalg.inv(M).T

      # now take into accout the fact that we are on the surface so the r3 or x3
      # component is not independent from the 2 other ones!!
      ux = ux - nx*uz/nz
      uy = uy - nx*uzOrtho/nz
      uxOrtho = uxOrtho - ny*uz/nz
      uyOrtho = uyOrtho - ny*uzOrtho/nz
      # print(M)

      # buil the matrix that transforms the gradient in (r1, r2, r3) to the one in (x1, x2, x3)
      MMGrad = MM1.T

      # SPH kernel
      # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
      # startTime = time.time()
      # L = np.empty((0), dtype=int)
      # check if we are on the bottom ot top row!!!
      if indy == 0:
          lInd = 0
      if indy == nrows - 1:
          rInd = 1
      for n in range(lInd, rInd):
          ic = (indx - 1) + ncols * (indy + n)
          # make sure not to take particles from the other edge
          imax = max(ic, ncols * (indy + n))
          imin = min(ic+3, ncols * (indy + n + 1))
          iPstart = indPartInCell[imax]
          iPend = indPartInCell[imin]
          # loop on all particles in neighbour boxes
          for p in range(iPstart, iPend):
              # index of particle in neighbour box
              l = int(partInCell[p])
              if j != l:
                  # L = np.append(L, l)
                  dx = X[l] - xx
                  dy = Y[l] - yy
                  dz = Z[l] - zz
                  # get coordinates in local coord system
                  r1 = scalProd(dx, dy, dz, ux, uy, uz)
                  r2 = scalProd(dx, dy, dz, uxOrtho, uyOrtho, uzOrtho)
                  # impse r3=0 even if the particle is not exactly on the tengent plane
                  # get norm of r = xj - xl
                  r = norm(r1, r2, 0)
                  if r < 0.001 * rKernel:
                      # impose a minimum distance between particles
                      r1 = 0.0001 * rKernel * r1
                      r2 = 0.0001 * rKernel * r2
                      r = 0.0001 * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      dwdr = dfacKernel * hr * hr
                      massl = mass[l]
                      mdwdrr = massl * dwdr / r
                      gradhX = gradhX + massl * dwdr * dx / r
                      gradhY = gradhY + massl * dwdr * dy / r
                      gradhZ = gradhZ + massl * dwdr * dz / r
                      K1, K2 = 1, 1
                      G1 = mdwdrr * K1*r1
                      G2 = mdwdrr * K2*r2
                      G3 = 0

                      g1 = nx/(nz)
                      g2 = ny/(nz)

                      GX1 = MMGrad[0, 0]*G1 + MMGrad[0, 1]*G2
                      GY1 = MMGrad[1, 0]*G1 + MMGrad[1, 1]*G2
                      GZ1 = (- g1*GX1 - g2*GY1)

                      GX = np.sum(GX1)
                      GY = np.sum(GY1)
                      GZ = np.sum(GZ1)
      # tcpuSPH = time.time() - startTime
      # TcpuSPH = TcpuSPH + tcpuSPH
      # startTime = time.time()
      GHX[j] = GHX[j] - gradhX / rho
      GHY[j] = GHY[j] - gradhY / rho
      GHZ[j] = GHZ[j] - gradhZ / rho
      # tcpuadd = time.time() - startTime
      # Tcpuadd = Tcpuadd + tcpuadd
    return GHX, GHY, GHZ
#
#
# def calcGradHSPH(double[:] X, double[:] Y, double[:] Z, double[:] M, long[:] indPartInCell, long[:] partInCell, int indx, int indy, int j, int ncols, int nrows, double csz):
#     """ Compute gradient of Flow Depth using SPH (for loop implementation)
#
#     Parameters
#     ----------
#     particles : dict
#     j: int
#         index of particle under consideration
#     ncols: int
#         number of columns of the DEM
#     nrows: int
#         number of rows of the DEM
#     csz  : float
#         cellsize of the DEM
#
#     Returns
#     -------
#     gradhX: float
#         x coordinate of the gradient of the flow depth at particle j location
#     gradhY: float
#         x coordinate of the gradient of the flow depth at particle j location
#     gradhZ: float
#         x coordinate of the gradient of the flow depth at particle j location
#     L: 1D numpy array
#         index of particles within the kernel function radius
#     """
#     # SPH kernel
#     # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
#     cdef double rKernel = csz
#     cdef double facKernel = 10.0 / (3.1415 * rKernel*rKernel*rKernel*rKernel*rKernel)
#     cdef double dfacKernel = -3.0 * facKernel
#     cdef double x = X[j]
#     cdef double y = Y[j]
#     cdef double z = Z[j]
#     cdef double dx, dy, dz, r, hr, dwdr, massl
#     # With loop
#     cdef double gradhX = 0
#     cdef double gradhY = 0
#     cdef double gradhZ = 0
#     # startTime = time.time()
#     L = np.empty((0), dtype=int)
#     # check if we are on the bottom ot top row!!!
#     cdef int lInd = -1
#     cdef int rInd = 2
#     cdef int ic, n, p, l, imax, imin, iPstart, iPend
#     if indy == 0:
#         lInd = 0
#     if indy == nrows - 1:
#         rInd = 1
#     for n in range(lInd, rInd):
#         ic = (indx - 1) + ncols * (indy + n)
#         # make sure not to take particles from the other edge
#         imax = max(ic, ncols * (indy + n))
#         imin = min(ic+3, ncols * (indy + n + 1))
#         iPstart = int(indPartInCell[imax])
#         iPend = int(indPartInCell[imin])
#         # loop on all particles in neighbour boxes
#         for p in range(iPstart, iPend):
#             # index of particle in neighbour box
#             l = int(partInCell[p])
#             if j != l:
#                 L = np.append(L, l)
#                 dx = X[l] - x
#                 dy = Y[l] - y
#                 dz = Z[l] - z
#                 r = norm(dx, dy, dz)
#                 if r < 0.001 * rKernel:
#                     # impose a minimum distance between particles
#                     r = 0.001 * rKernel
#                 if r < rKernel:
#                     hr = rKernel - r
#                     dwdr = dfacKernel * hr * hr
#                     massl = M[l]
#                     gradhX = gradhX + massl * dwdr * dx / r
#                     gradhY = gradhY + massl * dwdr * dy / r
#                     gradhZ = gradhZ + massl * dwdr * dz / r
#     # print((time.time() - startTime)/1)
#
#     return gradhX, gradhY,  gradhZ, L

cdef double norm(double x, double y, double z):
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
  return sqrt(x*x + y*y + z*z)

  cdef normalize(double x, double y, double z):
      """ Normalize vector (x, y, z) for the Euclidean norm.

      (x, y, z) can be np arrays.

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
          xn: numpy array
              x component of the normalized vector
          yn: numpy array
              y component of the normalized vector
          zn: numpy array
              z component of the normalized vector
      """
      # TODO : avoid error message when input vector is zero and make sure
      # to return zero
      double norme = norm(x, y, z)
      double xn = x / norme
      # xn = np.where(np.isnan(xn), 0, xn)
      double yn = y / norme
      # yn = np.where(np.isnan(yn), 0, yn)
      double zn = z / norme
      # zn = np.where(np.isnan(zn), 0, zn)

      return xn, yn, zn


  cdef croosProd(double ux, double uy, double uz, double vx, double vy, double vz):
      """ Compute cross product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
      """
      double wx = uy*vz - uz*vy
      double wy = uz*vx - ux*vz
      double wz = ux*vy - uy*vx

      return wx, wy, wz


  cdef double scalProd(double ux, double uy, double uz, double vx, double vy, double vz):
      """ Compute scalar product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
      """
      return ux*vx + uy*vy + uz*vz
