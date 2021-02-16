"""
    function related to SPH calculations in com1DFA

    to build: go to repository containing this file and run:
    python setup.py build_ext --inplace
"""

# Load modules
import copy
import logging
import math
import numpy as np
cimport numpy as np
from libc cimport math as math
# from libc.math cimport log as ln
cimport cython

# Local imports
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.in3Utils.geoTrans as geoTrans


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
# cdef int SPHOption = 1


ctypedef double dtypef_t
ctypedef long dtypel_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def pointsToRasterC(x, y, z, Z0, csz=1, xllc=0, yllc=0):
    """ Interpolate from unstructured points to grid

    Interpolate unstructured points on a structured grid using bilinear interpolation
    The (x, y) points have to be on the extend of the DEM!!

    Parameters
      ----------
      x: 1D numpy array
          x coordinate of the points
      y: 1D numpy array
          y coordinate of the points
      z: 1D numpy array
          quantity to interpolate associated to (x, y) points
      Z0: 2D numpy array
          initial ratser for interpolated result
      csz : float
          cellsize
      xllc : float
          x coord of the lower left center
      yllc : float
          y coord of the lower left center

      Returns
      -------

      Z1: 2D numpy array
          interpolated results
    """
    n, m = np.shape(Z0)
    cdef int nrow = int(n)
    cdef int ncol = int(m)
    cdef int Lx0, Ly0, Lx1, Ly1
    cdef double Lx, Ly, xx, yy, zz
    cdef double xllc0 = xllc
    cdef double yllc0 = yllc
    cdef double csz0 = csz
    cdef double[:] XX = x
    cdef double[:] YY = y
    cdef double[:] ZZ = z
    cdef double[:] Z = Z0.flatten()
    cdef int Npart = len(x)
    cdef int j, ic

    for j in range(Npart):
      xx = XX[j]
      yy = YY[j]
      zz = ZZ[j]
      # find coordinates in normalized ref (origin (0,0) and cellsize 1)
      Lx = (xx - xllc0) / csz0
      Ly = (yy - yllc0) / csz0

      # find coordinates of the 4 nearest cornes on the raster
      Lx0 = <int>Lx
      Ly0 = <int>Ly
      Lx1 = Lx0 + 1
      Ly1 = Ly0 + 1
      # prepare for bilinear interpolation
      dx = Lx - Lx0
      dy = Ly - Ly0

      # add the component of the points value to the 4 neighbour grid points
      # start with the lower left
      f11 = zz*(1-dx)*(1-dy)
      ic = Lx0 + ncol * Ly0
      Z[ic] = Z[ic] + f11
      # lower right
      f21 = zz*dx*(1-dy)
      ic = Lx1 + ncol * Ly0
      Z[ic] = Z[ic] + f21
      # uper left
      f12 = zz*(1-dx)*dy
      ic = Lx0 + ncol * Ly1
      Z[ic] = Z[ic] + f12
      # and uper right
      f22 = zz*dx*dy
      ic = Lx1 + ncol * Ly1
      Z[ic] = Z[ic] + f22

    Z1 = np.reshape(np.asarray(Z), (np.shape(Z0)))

    return Z1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def computeForceC(cfg, particles, fields, dem, Ment, Cres, dT):
  """ compute forces acting on the particles (without the SPH component)

  Cython implementation implementation

  Parameters
  ----------
  cfg: configparser
      configuration for DFA simulation
  particles : dict
      particles dictionary at t
  dem : dict
      dictionary with dem information
  Ment : 2D numpy array
      entrained mass raster
  Cres : 2D numpy array
      resistance raster
  dT : float
      time step

  Returns
  -------
  force : dict
      force dictionary
  """
  cdef double Rs0 = cfg.getfloat('Rs0')
  cdef double kappa = cfg.getfloat('kappa')
  cdef double B = cfg.getfloat('B')
  cdef double R = cfg.getfloat('R')
  cdef double rho = cfg.getfloat('rho')
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef int frictType = cfg.getint('frictType')
  cdef int interpOption = cfg.getint('interpOption')
  cdef double subgridMixingFactor = cfg.getfloat('subgridMixingFactor')
  cdef double dt = dT
  cdef double mu = cfg.getfloat('mu')
  cdef int Npart = particles['Npart']
  cdef double csz = dem['header'].cellsize
  cdef double[:, :] Nx = dem['Nx']
  cdef double[:, :] Ny = dem['Ny']
  cdef double[:, :] Nz = dem['Nz']

  cdef double[:] Fnormal = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceX = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceY = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceZ = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceFrict = np.zeros(Npart, dtype=np.float64)
  # cdef double[:] forceFrictY = np.zeros(Npart, dtype=np.float64)
  # cdef double[:] forceFrictZ = np.zeros(Npart, dtype=np.float64)
  cdef double[:] dM = np.zeros(Npart, dtype=np.float64)

  cdef double[:] mass = particles['m']
  cdef double[:] H = particles['h']
  cdef double[:] X = particles['x']
  cdef double[:] Y = particles['y']
  cdef double[:] UX = particles['ux']
  cdef double[:] UY = particles['uy']
  cdef double[:] UZ = particles['uz']
  cdef double[:, :] VX = fields['Vx']
  cdef double[:, :] VY = fields['Vy']
  cdef double[:, :] VZ = fields['Vz']
  cdef long[:] IndCellX = particles['indX']
  cdef long[:] IndCellY = particles['indY']
  cdef long indCellX, indCellY
  cdef double A, uMag, m, h, x, y, z, ux, uy, uz, nx, ny, nz
  cdef double vMeanx, vMeany, vMeanz, vMeanNorm, dvX, dvY, dvZ
  cdef double uxDir, uyDir, uzDir, nxEnd, nyEnd, nzEnd, nxAvg, nyAvg, nzAvg
  cdef double gravAccNorm, accNormCurv, effAccNorm, gravAccTangX, gravAccTangY, gravAccTangZ, forceBotTang, sigmaB, tau
  cdef int j
  force = {}
  # loop on particles
  for j in range(Npart):
      m = mass[j]
      x = X[j]
      y = Y[j]
      h = H[j]
      ux = UX[j]
      uy = UY[j]
      uz = UZ[j]
      indCellX = IndCellX[j]
      indCellY = IndCellY[j]
      # deduce area
      A = m / (h * rho)
      # get normal at the particle location
      nx, ny, nz = getVector(x, y, Nx, Ny, Nz, csz, interpOption)
      nx, ny, nz = normalize(nx, ny, nz)
      # get velocity magnitude and direction
      uMag = norm(ux, uy, uz)
      if uMag>0:
        uxDir, uyDir, uzDir = normalize(ux, uy, uz)
      # else:
      #   ux = 1
      #   uy = 0
      #   uz = -(1*nx + 0*ny) / nz
      #   uxDir, uyDir, uzDir = normalize(ux, uy, uz)

      # add artificial viscosity
      vMeanx, vMeany, vMeanz = getVector(x, y, VX, VY, VZ, csz, interpOption)
      # normal component of the velocity
      vMeanNorm = scalProd(vMeanx, vMeany, vMeanz, nx, ny, nz)
      vMeanx = vMeanx - vMeanNorm * nx
      vMeany = vMeany - vMeanNorm * ny
      vMeanz = vMeanz - vMeanNorm * nz
      dvX = vMeanx - ux
      dvY = vMeany - uy
      dvZ = vMeanz - uz
      dvMag = norm(dvX, dvY, dvZ)
      Alat = 2.0 * math.sqrt((m * h) / rho)
      fDrag = (subgridMixingFactor * 0.5 * rho * dvMag * Alat * dt) / m

      # get normal at the particle estimated end location
      xEnd = x + dt * ux
      yEnd = y + dt * uy
      nxEnd, nyEnd, nzEnd = getVector(xEnd, yEnd, Nx, Ny, Nz, csz, interpOption)
      nxEnd, nyEnd, nzEnd = normalize(nxEnd, nyEnd, nzEnd)
      # get average of those normals
      nxAvg = nx + nxEnd
      nyAvg = ny + nyEnd
      nzAvg = nz + nzEnd
      nxAvg, nyAvg, nzAvg = normalize(nxAvg, nyAvg, nzAvg)

      # acceleration due to curvature
      accNormCurv = (ux*(nxEnd-nx) + uy*(nyEnd-ny) + uz*(nzEnd-nz)) / dt
      # normal component of the acceleration of gravity
      gravAccNorm = - gravAcc * nzAvg
      effAccNorm = gravAccNorm + accNormCurv
      if(effAccNorm < 0.0):
          Fnormal[j] = m * effAccNorm

      # body forces (tangential component of acceleration of gravity)
      gravAccTangX = - gravAccNorm * nxAvg
      gravAccTangY = - gravAccNorm * nyAvg
      gravAccTangZ = -gravAcc - gravAccNorm * nzAvg
      # adding gravity force contribution
      forceX[j] = forceX[j] + gravAccTangX * m
      forceY[j] = forceY[j] + gravAccTangY * m
      forceZ[j] = forceZ[j] + gravAccTangZ * m

      # Calculating bottom shear and normal stress
      if(effAccNorm > 0.0):
          # if fluid detatched
          # log.info('fluid detatched for particle %s', j)
          tau = 0.0
      else:
          # bottom normal stress sigmaB
          sigmaB = - effAccNorm * rho * h
          if frictType == 1:
            # SamosAT friction type (bottom shear stress)
            tau = SamosATfric(rho, Rs0, mu, kappa, B, R, uMag, sigmaB, h)
          elif frictType == 2:
            # coulomb friction type (bottom shear stress)
            tau = mu * sigmaB

      # adding bottom shear resistance contribution
      forceBotTang = - A * tau
      forceFrict[j] = forceFrict[j] - forceBotTang
      # forceFrictX[j] = forceFrictX[j] + forceBotTang * uxDir
      # forceFrictY[j] = forceFrictY[j] + forceBotTang * uyDir
      # forceFrictZ[j] = forceFrictZ[j] + forceBotTang * uzDir

      # update velocity with artificial viscosity - implicit method
      ux = ux + fDrag * vMeanx
      uy = uy + fDrag * vMeany
      uz = uz + fDrag * vMeanz
      ux = ux / (1.0 + fDrag)
      uy = uy / (1.0 + fDrag)
      uz = uz / (1.0 + fDrag)
      UX[j] = ux
      UY[j] = uy
      UZ[j] = uz

  # save results
  force['dM'] = np.asarray(dM)
  force['forceX'] = np.asarray(forceX)
  force['forceY'] = np.asarray(forceY)
  force['forceZ'] = np.asarray(forceZ)

  force['forceFrict'] = np.asarray(forceFrict)
  particles['ux'] = np.asarray(UX)
  particles['uy'] = np.asarray(UY)
  particles['uz'] = np.asarray(UZ)

  # force['forceFrictY'] = np.asarray(forceFrictY)
  # force['forceFrictZ'] = np.asarray(forceFrictZ)
  return particles, force


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def updatePositionC(cfg, particles, dem, force):
  """ update particle position using euler forward scheme

  Cython implementation

  Parameters
  ----------
  cfg: configparser
      configuration for DFA simulation
  particles : dict
      particles dictionary at t
  dem : dict
      dictionary with dem information
  force : dict
      force dictionary
  Returns
  -------
  particles : dict
      particles dictionary at t + dt
  """
  DT = cfg.getfloat('dt')
  cdef double dt = DT
  cdef double stopCrit = cfg.getfloat('stopCrit')
  log.debug('dt used now is %f' % DT)
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef double rho = cfg.getfloat('rho')
  cdef int interpOption = cfg.getint('interpOption')
  cdef double csz = dem['header'].cellsize
  cdef double mu = cfg.getfloat('mu')
  cdef int Npart = particles['Npart']
  cdef double[:, :] Nx = dem['Nx']
  cdef double[:, :] Ny = dem['Ny']
  cdef double[:, :] Nz = dem['Nz']
  cdef double[:, :] Z = dem['rasterData']
  # initializeinter = np.zeros(N, dtype=np.float64)
  cdef double[:] MNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] XNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] YNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] ZNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] SNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] UXNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] UYNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] UZNew = np.zeros(Npart, dtype=np.float64)

  cdef double[:] mass = particles['m']
  cdef double[:] H = particles['h']
  cdef double[:] S = particles['s']
  cdef double[:] X = particles['x']
  cdef double[:] Y = particles['y']
  cdef double[:] UX = particles['ux']
  cdef double[:] UY = particles['uy']
  cdef double[:] UZ = particles['uz']
  cdef double[:] forceX = force['forceX']
  cdef double[:] forceY = force['forceY']
  cdef double[:] forceZ = force['forceZ']
  cdef double[:] forceFrict = force['forceFrict']
  cdef double[:] forceSPHX = force['forceSPHX']
  cdef double[:] forceSPHY = force['forceSPHY']
  cdef double[:] forceSPHZ = force['forceSPHZ']
  cdef double[:] dM = force['dM']
  cdef double TotkinEne = particles['kineticEne']
  cdef double TotpotEne = particles['potentialEne']
  cdef double peakKinEne = particles['peakKinEne']
  cdef double TotkinEneNew = 0
  cdef double TotpotEneNew = 0
  cdef double m, h, x, y, z, s, ux, uy, uz, nx, ny, nz, dtStop
  cdef double xDir, yDir, zDir, ForceDriveX, ForceDriveY, ForceDriveZ, zeroCrossing
  cdef double mNew, xNew, yNew, zNew, uxNew, uyNew, uzNew, sNew, uN, uMag, uMagNew
  cdef int j
  # loop on particles
  for j in range(Npart):
    m = mass[j]
    x = X[j]
    y = Y[j]
    h = H[j]
    ux = UX[j]
    uy = UY[j]
    uz = UZ[j]
    s = S[j]

    # Force magnitude (without friction)
    ForceDriveX = forceX[j] + forceSPHX[j]
    ForceDriveY = forceY[j] + forceSPHY[j]
    ForceDriveZ = forceZ[j] + forceSPHZ[j]

    # velocity magnitude
    uMag = norm(ux, uy, uz)

    # procede to time integration
    # operator splitting
    # estimate new velocity due to driving force
    uxNew = ux + ForceDriveX * dt / m
    uyNew = uy + ForceDriveY * dt / m
    uzNew = uz + ForceDriveZ * dt / m
    uMagNew = norm(uxNew, uyNew, uzNew)
    # will friction force stop the particle
    if uMagNew<dt*forceFrict[j]/m:
      # stop the particle
      uxNew = 0
      uyNew = 0
      uzNew = 0
      # particle stops after
      if uMag<=0:
        dtStop = 0
      else:
        dtStop = m * uMagNew / (dt * forceFrict[j])
    else:
      # add friction force in the opposite direction of the motion
      xDir, yDir, zDir = normalize(uxNew, uyNew, uzNew)
      uxNew = uxNew - xDir * forceFrict[j] * dt / m
      uyNew = uyNew - yDir * forceFrict[j] * dt / m
      uzNew = uzNew - zDir * forceFrict[j] * dt / m
      dtStop = dt

    # update mass
    mNew = m + dM[j]
    # update position
    xNew = x + dtStop * 0.5 * (ux + uxNew)
    yNew = y + dtStop * 0.5 * (uy + uyNew)
    sNew = s + math.sqrt((xNew-x)*(xNew-x) + (yNew-y)*(yNew-y))
    # make sure particle is on the mesh (recompute the z component)
    zNew = getScalar(xNew, yNew, Z, csz, interpOption)
    nx, ny, nz = getVector(xNew, yNew, Nx, Ny, Nz, csz, interpOption)
    nx, ny, nz = normalize(nx, ny, nz)
    # velocity magnitude
    uMag = norm(uxNew, uyNew, uzNew)
    # normal component of the velocity
    # uN = scalProd(uxNew, uyNew, uzNew, nx, ny, nz)
    uN = uxNew*nx + uyNew*ny + uzNew*nz
    # remove normal component of the velocity
    uxNew = uxNew - uN * nx
    uyNew = uyNew - uN * ny
    uzNew = uzNew - uN * nz
    TotkinEneNew = TotkinEneNew + 0.5 * m * uMag * uMag
    TotpotEneNew = TotpotEneNew + mNew * gravAcc * zNew
    XNew[j] = xNew
    YNew[j] = yNew
    ZNew[j] = zNew
    UXNew[j] = uxNew
    UYNew[j] = uyNew
    UZNew[j] = uzNew
    SNew[j] = sNew
    MNew[j] = mNew
  particles['ux'] = np.asarray(UXNew)
  particles['uy'] = np.asarray(UYNew)
  particles['uz'] = np.asarray(UZNew)
  particles['s'] = np.asarray(SNew)
  particles['m'] = np.asarray(MNew)
  particles['mTot'] = np.sum(particles['m'])
  particles['x'] = np.asarray(XNew)
  particles['y'] = np.asarray(YNew)
  particles['z'] = np.asarray(ZNew)
  particles['kineticEne'] = TotkinEneNew
  particles['potentialEne'] = TotpotEneNew
  if peakKinEne < TotkinEneNew:
    particles['peakKinEne'] = TotkinEneNew
  if TotkinEneNew <= stopCrit*peakKinEne:
    particles['iterate'] = False

  # make sure particle is on the mesh (recompute the z component)
  # particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')
  #################################################################
  # this is dangerous!!!!!!!!!!!!!!
  ###############################################################
  # remove particles that are not located on the mesh any more
  particles = com1DFA.removeOutPart(cfg, particles, dem)
  return particles

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def updateFieldsC(cfg, particles, dem, fields):
  """ update fields and particles fow depth

  Cython implementation

 Parameters
 ----------
 cfg: configparser
     configuration for DFA simulation
 particles : dict
     particles dictionary
 dem : dict
     dictionary with dem information
 fields : dict
     fields dictionary
 Returns
 -------

 particles : dict
     particles dictionary
 fields : dict
     fields dictionary
 """
  cdef double rho = cfg.getfloat('rho')
  cdef int interpOption = cfg.getint('interpOption')
  header = dem['header']
  CSZ = dem['header'].cellsize
  cdef double[:, :]A = dem['Area']
  ncols = header.ncols
  nrows = header.nrows
  cdef double[:, :] MassBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] VBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] PBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] FDBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] MomBilinearX = np.zeros((nrows, ncols))
  cdef double[:, :] MomBilinearY = np.zeros((nrows, ncols))
  cdef double[:, :] MomBilinearZ = np.zeros((nrows, ncols))
  cdef double[:, :] VXBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] VYBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] VZBilinear = np.zeros((nrows, ncols))

  cdef double[:] mass = particles['m']
  cdef double[:] X = particles['x']
  cdef double[:] Y = particles['y']
  cdef double[:] UX = particles['ux']
  cdef double[:] UY = particles['uy']
  cdef double[:] UZ = particles['uz']
  cdef double[:, :] PFV = fields['pfv']
  cdef double[:, :] PP = fields['ppr']
  cdef double[:, :] PFD = fields['pfd']
  cdef long[:] IndX = particles['indX']
  cdef long[:] IndY = particles['indY']
  cdef int nrow = int(nrows)
  cdef int ncol = int(ncols)
  cdef int Lx0, Ly0, Lx1, Ly1
  cdef double f00, f01, f10, f11
  cdef double m, h, x, y, z, s, ux, uy, uz, nx, ny, nz, hbb, hLim, aPart
  cdef double xllc = 0
  cdef double yllc = 0
  cdef double csz = CSZ
  cdef int Npart = np.size(particles['x'])
  cdef double[:] hBB = np.zeros((Npart))
  cdef int j, i
  cdef long indx, indy

  for j in range(Npart):
    x = X[j]
    y = Y[j]
    ux = UX[j]
    uy = UY[j]
    uz = UZ[j]
    m = mass[j]
    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    # find coordinates of the 4 nearest cornes on the raster
    # prepare for bilinear interpolation
    Lx0, Lx1, Ly0, Ly1, f00, f10, f01, f11 = getWeights(x, y, csz, interpOption)

    # add the component of the points value to the 4 neighbour grid points
    # start with the lower left
    MassBilinear[Ly0, Lx0] = MassBilinear[Ly0, Lx0] + m * f00
    FDBilinear[Ly0, Lx0] = FDBilinear[Ly0, Lx0] + m / (A[Ly0, Lx0] * rho) * f00
    MomBilinearX[Ly0, Lx0] = MomBilinearX[Ly0, Lx0] + m * ux * f00
    MomBilinearY[Ly0, Lx0] = MomBilinearY[Ly0, Lx0] + m * uy * f00
    MomBilinearZ[Ly0, Lx0] = MomBilinearZ[Ly0, Lx0] + m * uz * f00
    # lower right
    MassBilinear[Ly0, Lx1] = MassBilinear[Ly0, Lx1] + m * f10
    FDBilinear[Ly0, Lx1] = FDBilinear[Ly0, Lx1] + m / (A[Ly0, Lx1] * rho) * f10
    MomBilinearX[Ly0, Lx1] = MomBilinearX[Ly0, Lx1] + m * ux * f10
    MomBilinearY[Ly0, Lx1] = MomBilinearY[Ly0, Lx1] + m * uy * f10
    MomBilinearZ[Ly0, Lx1] = MomBilinearZ[Ly0, Lx1] + m * uz * f10
    # uper left
    MassBilinear[Ly1, Lx0] = MassBilinear[Ly1, Lx0] + m * f01
    FDBilinear[Ly1, Lx0] = FDBilinear[Ly1, Lx0] + m / (A[Ly1, Lx0] * rho) * f01
    MomBilinearX[Ly1, Lx0] = MomBilinearX[Ly1, Lx0] + m * ux * f01
    MomBilinearY[Ly1, Lx0] = MomBilinearY[Ly1, Lx0] + m * uy * f01
    MomBilinearZ[Ly1, Lx0] = MomBilinearZ[Ly1, Lx0] + m * uz * f01
    # and uper right
    MassBilinear[Ly1, Lx1] = MassBilinear[Ly1, Lx1] + m * f11
    FDBilinear[Ly1, Lx1] = FDBilinear[Ly1, Lx1] + m / (A[Ly1, Lx1] * rho) * f11
    MomBilinearX[Ly1, Lx1] = MomBilinearX[Ly1, Lx1] + m * ux * f11
    MomBilinearY[Ly1, Lx1] = MomBilinearY[Ly1, Lx1] + m * uy * f11
    MomBilinearZ[Ly1, Lx1] = MomBilinearZ[Ly1, Lx1] + m * uz * f11

  for i in range(ncol):
    for j in range(nrow):
      if MassBilinear[j, i] > 0:
        VXBilinear[j, i] = MomBilinearX[j, i]/MassBilinear[j, i]
        VYBilinear[j, i] = MomBilinearY[j, i]/MassBilinear[j, i]
        VZBilinear[j, i] = MomBilinearZ[j, i]/MassBilinear[j, i]
        VBilinear[j, i] = norm(VXBilinear[j, i], VYBilinear[j, i], VZBilinear[j, i])
        PBilinear[j, i] = VBilinear[j, i] * VBilinear[j, i] * rho
      if VBilinear[j, i] > PFV[j, i]:
        PFV[j, i] = VBilinear[j, i]
      if PBilinear[j, i] > PP[j, i]:
        PP[j, i] = PBilinear[j, i]
      if FDBilinear[j, i] > PFD[j, i]:
        PFD[j, i] = FDBilinear[j, i]


  fields['FV'] = np.asarray(VBilinear)
  fields['Vx'] = np.asarray(VXBilinear)
  fields['Vy'] = np.asarray(VYBilinear)
  fields['Vz'] = np.asarray(VZBilinear)
  fields['P'] = np.asarray(PBilinear)
  fields['FD'] = np.asarray(FDBilinear)
  fields['pfv'] = np.asarray(PFV)
  fields['ppr'] = np.asarray(PP)
  fields['pfd'] = np.asarray(PFD)


  for j in range(Npart):
    x = X[j]
    y = Y[j]
    hbb = getScalar(x, y, FDBilinear, csz, interpOption)
    indx = IndX[j]
    indy = IndY[j]
    # aPart = A[indy, indx]
    # hLim = mass[j]/(rho*aPart)
    # if hbb< hLim:
    #   hbb = hLim
    hBB[j] = hbb

  particles['hBilinearBilinear'] = np.asarray(hBB)
  particles['h'] = np.asarray(hBB)

  # remove particles that have a too small height
  # particles = removeSmallPart(hmin, particles, dem)

  return particles, fields

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def getNeighboursC(particles, dem):
    """ Locate particles in cell for SPH computation (for loop implementation)

    Ĺocate each particle in a grid cell and build the indPartInCell and
    partInCell arrays. See issue #200 and documentation for details

    Parameters
    ----------
    particles : dict
    dem : dict

    Returns
    -------
    particles : dict
      updated particles dictionary with indPartInCell and partInCell arrays
    """
    # get grid information
    header = dem['header']
    cdef int ncols = header.ncols
    cdef int nrows = header.nrows
    cdef float csz = header.cellsize
    # get particle location
    cdef int Npart = particles['Npart']
    cdef int j
    cdef double[:] x = particles['x']
    cdef double[:] y = particles['y']

    # initialize outputs
    cdef int Ncells = ncols*nrows
    cdef long[:] indPartInCell = np.zeros(Ncells + 1).astype('int')
    cdef long[:] indPartInCell2 = np.zeros(Ncells + 1).astype('int')
    cdef long[:] partInCell = np.zeros(Npart).astype('int')
    cdef long[:] indX = np.zeros(Npart).astype('int')
    cdef long[:] indY = np.zeros(Npart).astype('int')
    cdef long[:] InCell = np.zeros(Npart).astype('int')
    # Count number of particles in each cell
    cdef long indx, indy, ic
    for j in range(Npart):
      indx = int((x[j] + csz/2) / csz)
      indy = int((y[j] + csz/2) / csz)
      # get index of cell containing the particle
      ic = indx + ncols * indy
      indPartInCell[ic+1] = indPartInCell[ic+1] + 1
    for j in range(Ncells):
      indPartInCell[j+1] = indPartInCell[j] + indPartInCell[j+1]
      indPartInCell2[j+1] = indPartInCell[j+1]

    # make the list of which particles are in which cell
    for j in range(Npart):
        indx = int((x[j] + csz/2) / csz)
        indy = int((y[j] + csz/2) / csz)
        ic = int(indx + ncols * indy)
        partInCell[indPartInCell2[ic]] = j
        indPartInCell2[ic] = indPartInCell2[ic] + 1
        indX[j] = indx
        indY[j] = indy
        InCell[j] = ic

    particles['indX'] = np.asarray(indX)
    particles['indY'] = np.asarray(indY)
    particles['InCell'] = np.asarray(InCell)
    particles['indPartInCell'] = np.asarray(indPartInCell)
    particles['partInCell'] = np.asarray(partInCell)

    return particles


def computeForceSPHC(cfg, particles, force, dem, SPHOption=2, gradient=0):
  """ Prepare data for C computation of lateral forces (SPH component)
  acting on the particles (SPH component)

  Parameters
  ----------
  cfg: configparser
      configuration for DFA simulation
  particles : dict
      particles dictionary at t
  force : dict
      force dictionary
  dem : dict
      dictionary with dem information
  Returns
  -------
  particles : dict
      particles dictionary at t
  force : dict
      force dictionary
  """
  # Load required parameters
  Npart = particles['Npart']
  header = dem['header']
  nrows = dem['header'].nrows
  ncols = dem['header'].ncols
  csz = dem['header'].cellsize
  Nx = dem['Nx']
  Ny = dem['Ny']
  Nz = dem['Nz']

  indX = particles['indX'].astype('int')
  indY = particles['indY'].astype('int')
  forceSPHX, forceSPHY, forceSPHZ = computeGradC(cfg, particles, header, Nx, Ny, Nz, indX, indY, SPHOption, gradient)
  forceSPHX = np.asarray(forceSPHX)
  forceSPHY = np.asarray(forceSPHY)
  forceSPHZ = np.asarray(forceSPHZ)
  # log.info(('cpu time SPH = %s s' % (TcpuSPH / Npart)))
  # log.info(('cpu time SPH add = %s s' % (Tcpuadd / Npart)))

  force['forceSPHX'] = forceSPHX
  force['forceSPHY'] = forceSPHY
  force['forceSPHZ'] = forceSPHZ
  # particles['GHX'] = GHX
  # particles['GHY'] = GHY
  # particles['GHZ'] = GHZ

  return particles, force

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def computeGradC(cfg, particles, header, double[:, :] Nx, double[:, :] Ny,
                 double[:, :] Nz, long[:] indX, long[:] indY, SPHOption, gradient):
  """ compute lateral forces acting on the particles (SPH component)

  Cython implementation

  Parameters
  ----------
  cfg: configparser
      configuration for DFA simulation
  particles : dict
      particles dictionary at t
  header : dict
      header dictionary
  Nx : 2D numpy array
      x component of the normal vector of the DEM
  Ny : 2D numpy array
      y component of the normal vector of the DEM
  Nz : 2D numpy array
      z component of the normal vector of the DEM
  indX : 1D numpy array
      column index of the location of the particles
  indY : 1D numpy array
      row index of the location of the particles
  gradient : int
    Return the gradient (if 1) or the force associated (if 0, default)
  Returns
  -------
  GHX : 1D numpy array
      x component of the lateral force
  GHY : 1D numpy array
      y component of the lateral force
  GHZ : 1D numpy array
      z component of the lateral force
  """
  cdef double rho = cfg.getfloat('rho')
  cdef double minRKern = cfg.getfloat('minRKern')
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef int interpOption = cfg.getint('interpOption')
  cdef double gravAcc3
  cdef double csz = header.cellsize
  cdef double[:] mass = particles['m']
  cdef double[:] h = particles['h']
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
  # SPH kernel
  # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
  cdef double rKernel = csz
  cdef double facKernel = 10.0 / (3.1415 * rKernel * rKernel * rKernel * rKernel * rKernel)
  cdef double dfacKernel = - 3.0 * facKernel
  cdef double[:] GHX = np.zeros(N, dtype=np.float64)
  cdef double[:] GHY = np.zeros(N, dtype=np.float64)
  cdef double[:] GHZ = np.zeros(N, dtype=np.float64)
  cdef double K1 = 1
  cdef double K2 = 1
  cdef double gradhX, gradhY, gradhZ, uMag, nx, ny, nz, G1, G2, mdwdrr
  cdef double g1, g2, g11, g12, g22, g33
  cdef double m11, m12, m22, GG1, GG2
  cdef double xx, yy, zz, ux, uy, uz, vx, vy, wx, wy, uxOrtho, uyOrtho, uzOrtho
  cdef double dx, dy, dz, dn, r, hr, dwdr
  cdef int lInd, rInd
  cdef long indx, indy
  cdef int j, ic, n, p, l, imax, imin, iPstart, iPend
  cdef int SPHoption = SPHOption
  cdef int grad = gradient
  # L = np.empty((0), dtype=int)
  # indL = np.zeros((N+1), dtype=int)
  # loop on particles
  for j in range(N):
    xx = X[j]
    yy = Y[j]
    zz = Z[j]
    gradhX = 0
    gradhY = 0
    gradhZ = 0
    G1 = 0
    G2 = 0
    m11 = 0
    m12 = 0
    m22 = 0
    indx = indX[j]
    indy = indY[j]
    ux = UX[j]
    uy = UY[j]
    uz = UZ[j]
    nx, ny, nz = getVector(xx, yy, Nx, Ny, Nz, csz, interpOption)
    nx, ny, nz = normalize(nx, ny, nz)
    gravAcc3 = scalProd(nx, ny, nz, 0, 0, gravAcc)
    uMag = norm(ux, uy, uz)
    if uMag < 0.1:
        ux = 1
        uy = 0
        uz = -(1*nx + 0*ny) / nz
        ux, uy, uz = normalize(ux, uy, uz)
        K1 = 1
        K2 = 1
    else:
        ux, uy, uz = normalize(ux, uy, uz)

    uxOrtho, uyOrtho, uzOrtho = croosProd(nx, ny, nz, ux, uy, uz)
    uxOrtho, uyOrtho, uzOrtho = normalize(uxOrtho, uyOrtho, uzOrtho)

    g1 = nx/(nz)
    g2 = ny/(nz)

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
                # indL[j+1] = indL[j+1] + 1
                # L = np.append(L, l)
                dx = X[l] - xx
                dy = Y[l] - yy
                dz = Z[l] - zz
                if SPHoption == 1:
                  dz = 0
                  # get norm of r = xj - xl
                  r = norm(dx, dy, dz)
                  if r < minRKern * rKernel:
                      # impose a minimum distance between particles
                      # dx = minRKern * rKernel * dx
                      # dy = minRKern * rKernel * dy
                      r = minRKern * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      dwdr = dfacKernel * hr * hr
                      mdwdrr = mass[l] * dwdr / r
                      gradhX = gradhX + mdwdrr*dx
                      gradhY = gradhY + mdwdrr*dy
                      gradhZ = gradhZ + mdwdrr*dz
                      gravAcc3 = gravAcc

                if SPHoption == 2:
                  # get coordinates in local coord system
                  r1 = scalProd(dx, dy, dz, ux, uy, uz)
                  r2 = scalProd(dx, dy, dz, uxOrtho, uyOrtho, uzOrtho)
                  # impse r3=0 even if the particle is not exactly on the tengent plane
                  # get norm of r = xj - xl
                  r = norm(r1, r2, 0)
                  if r < minRKern * rKernel:
                      # impose a minimum distance between particles
                      r1 = minRKern * rKernel * r1
                      r2 = minRKern * rKernel * r2
                      r = minRKern * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      dwdr = dfacKernel * hr * hr
                      mdwdrr = mass[l] * dwdr / r
                      G1 = mdwdrr * K1*r1
                      G2 = mdwdrr * K2*r2

                      gradhX = gradhX + ux*G1 + uxOrtho*G2
                      gradhY = gradhY + uy*G1 + uyOrtho*G2
                      gradhZ = gradhZ + (- g1*(ux*G1 + uxOrtho*G2) - g2*(uy*G1 + uyOrtho*G2))

                elif SPHoption == 3:
                  # constant exact gradient correction
                  # get coordinates in local coord system
                  r1 = scalProd(dx, dy, dz, ux, uy, uz)
                  r2 = scalProd(dx, dy, dz, uxOrtho, uyOrtho, uzOrtho)
                  # impse r3=0 even if the particle is not exactly on the tengent plane
                  # get norm of r = xj - xl
                  r = norm(r1, r2, 0)
                  if r < minRKern * rKernel:
                      # impose a minimum distance between particles
                      r1 = minRKern * rKernel * r1
                      r2 = minRKern * rKernel * r2
                      r = minRKern * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      dwdr = dfacKernel * hr * hr
                      mdwdrr = mass[l] * (1 - h[j]/h[l]) * dwdr / r
                      G1 = mdwdrr * K1*r1
                      G2 = mdwdrr * K2*r2

                      gradhX = gradhX + ux*G1 + uxOrtho*G2
                      gradhY = gradhY + uy*G1 + uyOrtho*G2
                      gradhZ = gradhZ + (- g1*(ux*G1 + uxOrtho*G2) - g2*(uy*G1 + uyOrtho*G2))

                if SPHoption == 4:
                  # linear exact gradient correction
                  # get coordinates in local coord system
                  r1 = scalProd(dx, dy, dz, ux, uy, uz)
                  r2 = scalProd(dx, dy, dz, uxOrtho, uyOrtho, uzOrtho)
                  # impse r3=0 even if the particle is not exactly on the tengent plane
                  # get norm of r = xj - xl
                  r = norm(r1, r2, 0)
                  if r < minRKern * rKernel:
                      # impose a minimum distance between particles
                      r1 = minRKern * rKernel * r1
                      r2 = minRKern * rKernel * r2
                      r = minRKern * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      dwdr = dfacKernel * hr * hr
                      mdwdrr = mass[l] * (1 - h[j]/h[l]) * dwdr / r
                      m11 = m11 + mass[l] / h[l] * dwdr / r * r1 * r1 / rho
                      m12 = m12 + mass[l] / h[l] * dwdr / r * r1 * r2 / rho
                      m22 = m22 + mass[l] / h[l] * dwdr / r * r2 * r2 / rho
                      G1 = G1 + mdwdrr * K1*r1
                      G2 = G2 + mdwdrr * K2*r2

                if SPHoption == 5:
                  # integral method
                  # get coordinates in local coord system
                  r1 = scalProd(dx, dy, dz, ux, uy, uz)
                  r2 = scalProd(dx, dy, dz, uxOrtho, uyOrtho, uzOrtho)
                  # impse r3=0 even if the particle is not exactly on the tengent plane
                  # get norm of r = xj - xl
                  r = norm(r1, r2, 0)
                  if r < minRKern * rKernel:
                      # impose a minimum distance between particles
                      r1 = minRKern * rKernel * r1
                      r2 = minRKern * rKernel * r2
                      r = minRKern * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      w = facKernel * hr * hr * hr
                      m11 = m11 + mass[l] / h[l] * w * r1 * r1 / rho
                      m12 = m12 + mass[l] / h[l] * w * r1 * r2 / rho
                      m22 = m22 + mass[l] / h[l] * w * r2 * r2 / rho
                      G1 = G1 + mass[l] * (1 - h[j]/h[l]) * w * K1*r1
                      G2 = G2 + mass[l] * (1 - h[j]/h[l]) * w * K2*r2

    if grad == 1:
      if SPHoption >= 4:
        # check that M is invertible
        if (m11*m22-m12*m12)>0:
          GG1 = 1/(m11*m22-m12*m12)*(m22*G1-m12*G2)
          GG2 = 1/(m11*m22-m12*m12)*(m11*G2-m12*G1)
        else:
          GG1 = G1
          GG2 = G2
        gradhX = ux*GG1 + uxOrtho*GG2
        gradhY = uy*GG1 + uyOrtho*GG2
        gradhZ = (- g1*(ux*GG1 + uxOrtho*GG2) - g2*(uy*GG1 + uyOrtho*GG2))
        GHX[j] = GHX[j] + gradhX / rho
        GHY[j] = GHY[j] + gradhY / rho
        GHZ[j] = GHZ[j] + gradhZ / rho
      else:
        GHX[j] = GHX[j] - gradhX / rho
        GHY[j] = GHY[j] - gradhY / rho
        GHZ[j] = GHZ[j] - gradhZ / rho
    else:
      if SPHoption >= 4:
          # check that M is invertible
        if (m11*m22-m12*m12)>0:
          GG1 = 1/(m11*m22-m12*m12)*(m22*G1-m12*G2)
          GG2 = 1/(m11*m22-m12*m12)*(m11*G2-m12*G1)
        else:
          GG1 = G1
          GG2 = G2
        gradhX = ux*GG1 + uxOrtho*GG2
        gradhY = uy*GG1 + uyOrtho*GG2
        gradhZ = (- g1*(ux*GG1 + uxOrtho*GG2) - g2*(uy*GG1 + uyOrtho*GG2))
        GHX[j] = GHX[j] - gradhX / rho * mass[j] * gravAcc3
        GHY[j] = GHY[j] - gradhY / rho * mass[j] * gravAcc3
        GHZ[j] = GHZ[j] - gradhZ / rho * mass[j] * gravAcc3
      else:
        GHX[j] = GHX[j] + gradhX / rho* mass[j] * gravAcc3
        GHY[j] = GHY[j] + gradhY / rho* mass[j] * gravAcc3
        GHZ[j] = GHZ[j] + gradhZ / rho* mass[j] * gravAcc3
  return GHX, GHY, GHZ# , L, indL


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def computeFDC(cfg, particles, header, double[:, :] Nx, double[:, :] Ny, double[:, :] Nz, long[:] indX, long[:] indY):
  """ compute flow depth at particle location (SPH flow depth)

  Cython implementation

  Parameters
  ----------
  cfg: configparser
      configuration for DFA simulation
  particles : dict
      particles dictionary at t
  header : dict
      header dictionary
  Nx : 2D numpy array
      x component of the normal vector of the DEM
  Ny : 2D numpy array
      y component of the normal vector of the DEM
  Nz : 2D numpy array
      z component of the normal vector of the DEM
  indX : 1D numpy array
      column index of the location of the particles
  indY : 1D numpy array
      row index of the location of the particles
  Returns
  -------
  H : 1D numpy array
      flow depth
  """
  cdef double rho = cfg.getfloat('rho')
  cdef double minRKern = cfg.getfloat('minRKern')
  cdef int interpOption = cfg.getint('interpOption')
  cdef double csz = header.cellsize
  cdef double[:] mass = particles['m']
  cdef double[:] X = particles['x']
  cdef double[:] Y = particles['y']
  cdef double[:] Z = particles['z']
  cdef double[:] H1 = particles['h']
  cdef double[:] UX = particles['ux']
  cdef double[:] UY = particles['uy']
  cdef double[:] UZ = particles['uz']
  # cdef double[:] gradx = particles['gradx']
  # cdef double[:] grady = particles['grady']
  # cdef double[:] gradz = particles['gradz']
  cdef long[:] indPartInCell = particles['indPartInCell']
  cdef long[:] partInCell = particles['partInCell']
  cdef int N = X.shape[0]
  cdef int nrows = header.nrows
  cdef int ncols = header.ncols
  # SPH kernel
  # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
  cdef double rKernel = csz
  cdef double facKernel = 10.0 / (3.1415 * rKernel*rKernel*rKernel*rKernel*rKernel)
  cdef double[:] H = np.zeros(N)
  cdef double[:] W = np.zeros(N)
  cdef double[:] C = np.zeros(N)
  cdef double dn, h, nx, ny, nz, gx, gy, gz
  cdef int j
  cdef double xx, yy, zz
  cdef double dx, dy, dz, r, hr, dwdr, massl, hl, ww
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
    ww = 0
    gx = 0
    gy = 0
    gz = 0
    indx = indX[j]
    indy = indY[j]
    nx, ny, nz = getVector(xx, yy, Nx, Ny, Nz, csz, interpOption)
    nx, ny, nz = normalize(nx, ny, nz)

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
                if r < minRKern * rKernel:
                    # impose a minimum distance between particles
                    dx = minRKern * rKernel * dx
                    dy = minRKern * rKernel * dy
                    dz = minRKern * rKernel * dz
                    r = minRKern * rKernel
                if r < rKernel:
                    hr = rKernel - r
                    w = facKernel * hr * hr * hr
                    massl = mass[l]
                    hl = H1[l]
                    dh = massl * w
                    # gx = gx + massl * w / hl * dx
                    # gy = gy + massl * w / hl * dy
                    # gz = gz + massl * w / hl * dz
                    ww = ww + massl * w / hl
                    h = h + dh

                    # gradhX = gradhX + mdwdrr * dx
                    # gradhY = gradhY + mdwdrr * dy
                    # gradhZ = gradhZ + mdwdrr * dz
    # tcpuSPH = time.time() - startTime
    # TcpuSPH = TcpuSPH + tcpuSPH
    # startTime = time.time()
    H[j] = H[j] + h / rho
    # C[j] = C[j] + (gx*gradx[l] + gy*grady[l] + gz*gradz[l]) / rho
    W[j] = W[j] + ww / rho
    # if H[j]>0:
    #   H[j] = H[j]/ww
    # tcpuadd = time.time() - startTime
    # Tcpuadd = Tcpuadd + tcpuadd
  return H, C, W


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
  return math.sqrt(x*x + y*y + z*z)

def normpy(x, y, z): # <-- small wrapper to expose norm() to Python
    return norm(x, y, z)

cdef double norm2(double x, double y, double z):
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
  return x*x + y*y + z*z

def norm2py(x, y, z): # <-- small wrapper to expose norm2() to Python
    return norm2(x, y, z)


@cython.cdivision(True)
cdef (double, double, double) normalize(double x, double y, double z):
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
  if norme>0:
    x = x / norme
    # xn = np.where(np.isnan(xn), 0, xn)
    y = y / norme
    # yn = np.where(np.isnan(yn), 0, yn)
    z = z / norme
    # zn = np.where(np.isnan(zn), 0, zn)
  return x, y, z


def normalizepy(x, y, z): # <-- small wrapper to expose normalize() to Python
    # cdef double xx = x
    # cdef double yy = y
    # cdef double zz = z
    return normalize(x, y, z)


@cython.cdivision(True)
cdef (double, double, double) croosProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute cross product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  cdef double wx = uy * vz - uz * vy
  cdef double wy = uz * vx - ux * vz
  cdef double wz = ux * vy - uy * vx
  return wx, wy, wz


def croosProdpy(x, y, z, u, v, w): # <-- small wrapper to expose croosProd() to Python
    return croosProd(x, y, z, u, v, w)


cdef double scalProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute scalar product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  return ux*vx + uy*vy + uz*vz

def scalProdpy(x, y, z, u, v, w): # <-- small wrapper to expose scalProd() to Python
    return scalProd(x, y, z, u, v, w)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef (int, int, int, int, double, double, double, double) getWeights(double x, double y, double csz, int interpOption):
  """ Prepare weight for interpolation from grid to single point location

  3 Options available : -0: nearest neighbour interpolation
                        -1: equal weights interpolation
                        -2: bilinear interpolation

  Parameters
  ----------
      x: float
          location in the x location of desiered interpolation
      y: float
          location in the y location of desiered interpolation
      csz: float
          cellsize of the grid
      interpOption: int
          -0: nearest neighbour interpolation
          -1: equal weights interpolation
          -2: bilinear interpolation

  Returns
  -------
      Lx0, Ly0, Lx1, Ly1: int
          colomn and row indices for interpolation
      f00, f10, f01, f11: float
          corresponding weights
  """
  cdef int Lx0, Ly0, Lx1, Ly1
  cdef double f00, f10, f01, f11
  cdef double Lx, Ly
  cdef double xllc = 0.
  cdef double yllc = 0.

  # find coordinates in normalized ref (origin (0,0) and cellsize 1)
  Lx = (x - xllc) / csz
  Ly = (y - yllc) / csz
  # find coordinates of the 4 nearest cornes on the raster
  Lx0 = <int>Lx
  Ly0 = <int>Ly
  Lx1 = Lx0 + 1
  Ly1 = Ly0 + 1
  # prepare for bilinear interpolation
  dx = Lx - Lx0
  dy = Ly - Ly0

  if interpOption == 0:
    dx = 1.*math.round(dx)
    dy = 1.*math.round(dy)
  elif interpOption == 1:
    dx = 1./2.
    dy = 1./2.
    
  # lower left
  f00 = (1-dx)*(1-dy)
  # lower right
  f10 = dx*(1-dy)
  # uper left
  f01 = (1-dx)*dy
  # and uper right
  f11 = dx*dy

  return Lx0, Lx1, Ly0, Ly1, f00, f10, f01, f11
  # return Lxy, w


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double getScalar(double x, double y, double[:, :] V, double csz, int interpOption):
  """ Interpolate vector field from grid to single point location

  Originaly created to get the normal vector at location (x,y) given the
  normal vector field on the grid. Grid has its origin in (0,0).
  Can be used to interpolate any vector field.
  Interpolation using a bilinear interpolation

  Parameters
  ----------
      x: float
          location in the x location of desiered interpolation
      y: float
          location in the y location of desiered interpolation
      V: 2D numpy array
          scalar field at the grid nodes
      csz: float
          cellsize of the grid
      interpOption: int
          -0: nearest neighbour interpolation
          -1: equal weights interpolation
          -2: bilinear interpolation

  Returns
  -------
      v: float
          interpolated scalar at position (x, y)
  """
  cdef int Lx0, Ly0, Lx1, Ly1
  cdef double f00, f01, f10, f11
  Lx0, Lx1, Ly0, Ly1, f00, f10, f01, f11 = getWeights(x, y, csz, interpOption)
  cdef double v = (V[Ly0, Lx0]*f00 +
                   V[Ly0, Lx1]*f10 +
                   V[Ly1, Lx0]*f01 +
                   V[Ly1, Lx1]*f11)


  return v


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef (double, double, double) getVector(double x, double y, double[:, :] Nx, double[:, :] Ny, double[:, :] Nz, double csz, int interpOption):
  """ Interpolate vector field from grid to single point location

  Originaly created to get the normal vector at location (x,y) given the
  normal vector field on the grid. Grid has its origin in (0,0).
  Can be used to interpolate any vector field.
  Interpolation using a bilinear interpolation

  Parameters
  ----------
      x: float
          location in the x location of desiered interpolation
      y: float
          location in the y location of desiered interpolation
      Nx: 2D numpy array
          x component of the vector field at the grid nodes
      Ny: 2D numpy array
          y component of the vector field at the grid nodes
      Nz: 2D numpy array
          z component of the vector field at the grid nodes
      csz: float
          cellsize of the grid
      interpOption: int
          -0: nearest neighbour interpolation
          -1: equal weights interpolation
          -2: bilinear interpolation

  Returns
  -------
      nx: float
          x component of the interpolated vector field at position (x, y)
      ny: float
          y component of the interpolated vector field at position (x, y)
      nz: float
          z component of the interpolated vector field at position (x, y)
  """
  cdef int Lx0, Ly0, Lx1, Ly1
  cdef double f00, f01, f10, f11
  Lx0, Lx1, Ly0, Ly1, f00, f10, f01, f11 = getWeights(x, y, csz, interpOption)
  cdef double nx = (Nx[Ly0, Lx0]*f00 +
                    Nx[Ly0, Lx1]*f10 +
                    Nx[Ly1, Lx0]*f01 +
                    Nx[Ly1, Lx1]*f11)
  cdef double ny = (Ny[Ly0, Lx0]*f00 +
                    Ny[Ly0, Lx1]*f10 +
                    Ny[Ly1, Lx0]*f01 +
                    Ny[Ly1, Lx1]*f11)
  cdef double nz = (Nz[Ly0, Lx0]*f00 +
                    Nz[Ly0, Lx1]*f10 +
                    Nz[Ly1, Lx0]*f01 +
                    Nz[Ly1, Lx1]*f11)
  return nx, ny, nz

@cython.cdivision(True)
cdef double SamosATfric(double rho, double Rs0, double mu, double kappa, double B, double R, double v, double p, double h):
  cdef double Rs = rho * v * v / (p + 0.001)
  cdef double div = h / R
  if div < 1.0:
    div = 1.0
  div = math.log(div) / kappa + B
  cdef double tau = p * mu * (1.0 + Rs0 / (Rs0 + Rs)) + rho * v * v / (div * div)
  return tau
