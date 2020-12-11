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
cdef double gravAcc = 9.81

def computeGrad(Npart, particles, Nx, Ny, Nz, NX, NY):
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

def computeGradcython(particles, header, double[:] Nx, double[:] Ny, double[:] Nz, long[:] indX, long[:] indY):
    cdef double[:] mass = particles['m']
    cdef double[:] X = particles['x']
    cdef double[:] Y = particles['y']
    cdef double[:] Z = particles['z']
    cdef double[:] UX = particles['ux']
    cdef double[:] UY = particles['uy']
    cdef double[:] UZ = particles['uz']
    cdef long[:] indPartInCell = particles['indPartInCell']
    cdef long[:] partInCell = particles['partInCell']
    cdef int N = X.shape[0]
    cdef int nrows = header.nrows
    cdef int ncols = header.ncols
    cdef double rKernel = csz
    cdef double facKernel = 10.0 / (3.1415 * rKernel*rKernel*rKernel*rKernel*rKernel)
    cdef double dfacKernel = -3.0 * facKernel
    cdef double[:] GHX = np.zeros(N)
    cdef double[:] GHY = np.zeros(N)
    cdef double[:] GHZ = np.zeros(N)
    cdef double K1 = 1
    cdef double K2 = 1
    cdef double gradhX, gradhY, gradhZ, uMag, g1, g2, nx, ny, nz, G1, G2, G3
    cdef int j
    cdef double xx, yy, zz, ux, uy, uz, vx, vy, wx, wy, uxOrtho, uyOrtho, uzOrtho
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
      nx = Nx[j]
      ny = Ny[j]
      nz = Nz[j]
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

      # now take into accout the fact that we are on the surface so the r3 or x3
      # component is not independent from the 2 other ones!!
      vx = ux - nx*uz/nz
      vy = uy - ny*uz/nz
      wx = uxOrtho - nx*uzOrtho/nz
      wy = uyOrtho - ny*uzOrtho/nz

      # SPH kernel
      # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
      # startTime = time.time()
      # L = np.empty((0), dtype=int)
      # check if we are on the bottom ot top row!!!
      lInd = -1
      rInd = 2
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
              l = partInCell[p]
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
                  # r = norm(dx, dy, dz)
                  if r < 0.0001 * rKernel:
                      # impose a minimum distance between particles
                      r1 = 0.0001 * rKernel * r1
                      r2 = 0.0001 * rKernel * r2
                      r = 0.0001 * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      dwdr = dfacKernel * hr * hr
                      massl = mass[l]
                      mdwdrr = massl * dwdr / r
                      G1 = mdwdrr * K1*r1
                      G2 = mdwdrr * K2*r2
                      G3 = 0

                      g1 = nx/(nz)
                      g2 = ny/(nz)

                      gradhX = gradhX + vx*G1 + wx*G2
                      gradhY = gradhY + vy*G1 + wy*G2
                      gradhZ = gradhZ + (- g1*(vx*G1 + wx*G2) - g2*(vy*G1 + wy*G2))

                      # gradhX = gradhX + mdwdrr * dx
                      # gradhY = gradhY + mdwdrr * dy
                      # gradhZ = gradhZ + mdwdrr * dz
      # tcpuSPH = time.time() - startTime
      # TcpuSPH = TcpuSPH + tcpuSPH
      # startTime = time.time()
      # GHX[j] = GHX[j] - gradhX / rho
      # GHY[j] = GHY[j] - gradhY / rho
      # GHZ[j] = GHZ[j] - gradhZ / rho

      GHX[j] = GHX[j] + gradhX / rho* mass[j] * gravAcc
      GHY[j] = GHY[j] + gradhY / rho* mass[j] * gravAcc
      GHZ[j] = GHZ[j] + gradhZ / rho* mass[j] * gravAcc
      # tcpuadd = time.time() - startTime
      # Tcpuadd = Tcpuadd + tcpuadd
    return GHX, GHY, GHZ


def computeFDcython(particles, header, double[:] Nx, double[:] Ny, double[:] Nz, long[:] indX, long[:] indY):
  cdef double[:] mass = particles['m']
  cdef double[:] X = particles['x']
  cdef double[:] Y = particles['y']
  cdef double[:] Z = particles['z']
  cdef double[:] UX = particles['ux']
  cdef double[:] UY = particles['uy']
  cdef double[:] UZ = particles['uz']
  cdef long[:] indPartInCell = particles['indPartInCell']
  cdef long[:] partInCell = particles['partInCell']
  cdef int N = X.shape[0]
  cdef int nrows = header.nrows
  cdef int ncols = header.ncols
  cdef double rKernel = csz
  cdef double facKernel = 10.0 / (3.1415 * rKernel*rKernel*rKernel*rKernel*rKernel)
  cdef double[:] H = np.zeros(N)
  cdef double dn, h, nx, ny, nz
  cdef int j
  cdef double xx, yy, zz
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
    h = 0
    indx = indX[j]
    indy = indY[j]
    nx = Nx[j]
    ny = Ny[j]
    nz = Nz[j]

    # SPH kernel
    # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
    # startTime = time.time()
    # L = np.empty((0), dtype=int)
    # check if we are on the bottom ot top row!!!
    lInd = -1
    rInd = 2
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
            l = partInCell[p]
            if j != l:
                # L = np.append(L, l)
                dx = X[l] - xx
                dy = Y[l] - yy
                dz = Z[l] - zz
                # get coordinates in local coord system
                dn = nx*dx + ny*dy + nz*dz
                dx = dx - dn*nx
                dy = dy - dn*ny
                dz = dz - dn*nz
                # impse r3=0 even if the particle is not exactly on the tengent plane
                # get norm of r = xj - xl
                r = norm(dx, dy, dz)
                # r = norm(dx, dy, dz)
                if r < 0.0001 * rKernel:
                    # impose a minimum distance between particles
                    dx = 0.0001 * rKernel * dx
                    dy = 0.0001 * rKernel * dy
                    dz = 0.0001 * rKernel * dz
                    r = 0.0001 * rKernel
                if r < rKernel:
                    hr = rKernel - r
                    w = facKernel * hr * hr * hr
                    massl = mass[l]
                    dh = massl * w
                    h = h + dh

                    # gradhX = gradhX + mdwdrr * dx
                    # gradhY = gradhY + mdwdrr * dy
                    # gradhZ = gradhZ + mdwdrr * dz
    # tcpuSPH = time.time() - startTime
    # TcpuSPH = TcpuSPH + tcpuSPH
    # startTime = time.time()
    H[j] = H[j] + h / rho
    # tcpuadd = time.time() - startTime
    # Tcpuadd = Tcpuadd + tcpuadd
  return H


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

cpdef normalize(double x, double y, double z):
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
  cdef double norme
  norme = norm(x, y, z)
  cdef double xn = x / norme
  # xn = np.where(np.isnan(xn), 0, xn)
  cdef double yn = y / norme
  # yn = np.where(np.isnan(yn), 0, yn)
  cdef double zn = z / norme
  # zn = np.where(np.isnan(zn), 0, zn)
  return xn, yn, zn


cpdef croosProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute cross product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  cdef double wx = uy*vz - uz*vy
  cdef double wy = uz*vx - ux*vz
  cdef double wz = ux*vy - uy*vx

  return wx, wy, wz


cdef double scalProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute scalar product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  return ux*vx + uy*vy + uz*vz
