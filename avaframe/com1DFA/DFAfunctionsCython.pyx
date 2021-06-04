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
import avaframe.com1DFA.DFAtools as DFAtls
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
ctypedef int dtypel_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def pointsToRasterC(double[:] xArray, double[:] yArray, double[:] zArray, Z0, double csz=1, double xllc=0, double yllc=0):
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
    cdef double Lx, Ly, x, y, z
    cdef double[:] zRaster = Z0.flatten()
    cdef int Npart = len(xArray)
    cdef int j, ic

    for j in range(Npart):
      x = xArray[j]
      y = yArray[j]
      z = zArray[j]
      # find coordinates in normalized ref (origin (0,0) and cellsize 1)
      Lx = (x - xllc) / csz
      Ly = (y - yllc) / csz

      # find coordinates of the 4 nearest cornes on the raster
      Lx0 = <int>math.floor(Lx)
      Ly0 = <int>math.floor(Ly)
      Lx1 = Lx0 + 1
      Ly1 = Ly0 + 1
      # prepare for bilinear interpolation
      dx = Lx - Lx0
      dy = Ly - Ly0

      # add the component of the points value to the 4 neighbour grid points
      # start with the lower left
      f11 = z*(1-dx)*(1-dy)
      ic = Lx0 + ncol * Ly0
      zRaster[ic] = zRaster[ic] + f11
      # lower right
      f21 = z*dx*(1-dy)
      ic = Lx1 + ncol * Ly0
      zRaster[ic] = zRaster[ic] + f21
      # uper left
      f12 = z*(1-dx)*dy
      ic = Lx0 + ncol * Ly1
      zRaster[ic] = zRaster[ic] + f12
      # and uper right
      f22 = z*dx*dy
      ic = Lx1 + ncol * Ly1
      zRaster[ic] = zRaster[ic] + f22

    Z1 = np.reshape(np.asarray(zRaster), (np.shape(Z0)))

    return Z1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def computeForceC(cfg, particles, fields, dem, dT, int frictType):
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
  dT : float
      time step

  Returns
  -------
  force : dict
      force dictionary
  """
  # read input parameters
  cdef double Rs0 = cfg.getfloat('Rs0')
  cdef double kappa = cfg.getfloat('kappa')
  cdef double B = cfg.getfloat('B')
  cdef double R = cfg.getfloat('R')
  cdef double entEroEnergy = cfg.getfloat('entEroEnergy')
  cdef double entShearResistance = cfg.getfloat('entShearResistance')
  cdef double entDefResistance = cfg.getfloat('entDefResistance')
  cdef double rho = cfg.getfloat('rho')
  cdef double rhoEnt = cfg.getfloat('rhoEnt')
  cdef double hRes = cfg.getfloat('hRes')
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef double curvAcceleration = cfg.getfloat('curvAcceleration')
  cdef double velMagMin = cfg.getfloat('velMagMin')
  cdef int interpOption = cfg.getint('interpOption')
  cdef int explicitFriction = cfg.getint('explicitFriction')
  cdef double subgridMixingFactor = cfg.getfloat('subgridMixingFactor')
  cdef double dt = dT
  cdef double mu = cfg.getfloat('mu')
  cdef int Npart = particles['Npart']
  cdef double csz = dem['header'].cellsize
  cdef double[:, :] ZDEM = dem['rasterData']
  cdef double[:, :] nxArray = dem['Nx']
  cdef double[:, :] nyArray = dem['Ny']
  cdef double[:, :] nzArray = dem['Nz']
  cdef double[:, :] areaRatser = dem['areaRaster']
  # read particles and fields
  cdef double[:] mass = particles['m']
  cdef double[:] hArray = particles['h']
  cdef double[:] xArray = particles['x']
  cdef double[:] yArray = particles['y']
  cdef double[:] zArray = particles['z']
  cdef double[:] uxArray = particles['ux']
  cdef double[:] uyArray = particles['uy']
  cdef double[:] uzArray = particles['uz']
  cdef double[:, :] VX = fields['Vx']
  cdef double[:, :] VY = fields['Vy']
  cdef double[:, :] VZ = fields['Vz']
  cdef double[:, :] entrMassRaster = fields['entrMassRaster']
  cdef double[:, :] cResRaster = fields['cResRaster']
  cdef int[:] indXDEM = particles['indXDEM']
  cdef int[:] indYDEM = particles['indYDEM']
  # initialize outputs
  cdef double[:] Fnormal = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceX = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceY = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceZ = np.zeros(Npart, dtype=np.float64)
  cdef double[:] forceFrict = np.zeros(Npart, dtype=np.float64)
  cdef double[:] dM = np.zeros(Npart, dtype=np.float64)
  # declare intermediate step variables
  cdef int indCellX, indCellY
  cdef double areaPart, areaCell, araEntrPart, cResCell, cResPart, uMag, m, dm, h, entrMassCell, dEnergyEntr, dis
  cdef double vMeanx, vMeany, vMeanz, vMeanNorm, dvX, dvY, dvZ
  cdef double x, y, z, xEnd, yEnd, zEnd, ux, uy, uz, uxDir, uyDir, uzDir
  cdef double nx, ny, nz, nxEnd, nyEnd, nzEnd, nxAvg, nyAvg, nzAvg
  cdef double gravAccNorm, accNormCurv, effAccNorm, gravAccTangX, gravAccTangY, gravAccTangZ, forceBotTang, sigmaB, tau
  cdef int Lxy[4]
  cdef double w[4]
  cdef int j
  force = {}
  # loop on particles
  for j in range(Npart):
      m = mass[j]
      x = xArray[j]
      y = yArray[j]
      z = zArray[j]
      h = hArray[j]
      ux = uxArray[j]
      uy = uyArray[j]
      uz = uzArray[j]
      indCellX = indXDEM[j]
      indCellY = indYDEM[j]
      # deduce area
      areaPart = m / (h * rho)
      # get normal at the particle location
      Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, csz, interpOption)
      nx, ny, nz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
      nx, ny, nz = normalize(nx, ny, nz)

      # add artificial viscosity
      vMeanx, vMeany, vMeanz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], VX, VY, VZ)
      # compute normal component of the velocity
      vMeanNorm = scalProd(vMeanx, vMeany, vMeanz, nx, ny, nz)
      # remove normal component (make sure vMean is in the tangent plane)
      vMeanx = vMeanx - vMeanNorm * nx
      vMeany = vMeany - vMeanNorm * ny
      vMeanz = vMeanz - vMeanNorm * nz
      # compute particle to field velocity difference
      dvX = vMeanx - ux
      dvY = vMeany - uy
      dvZ = vMeanz - uz
      dvMag = norm(dvX, dvY, dvZ)
      Alat = 2.0 * math.sqrt((m * h) / rho)
      fDrag = (subgridMixingFactor * 0.5 * rho * dvMag * Alat * dt) / m

      # update velocity with artificial viscosity - implicit method
      ux = ux + fDrag * vMeanx
      uy = uy + fDrag * vMeany
      uz = uz + fDrag * vMeanz
      ux = ux / (1.0 + fDrag)
      uy = uy / (1.0 + fDrag)
      uz = uz / (1.0 + fDrag)

      # get normal at the particle estimated end location
      xEnd = x + dt * ux
      yEnd = y + dt * uy
      zEnd = z + dt * uz
      # this point is not necessarily on the surface, project it on the surface
      xEnd, yEnd, Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = normalProjectionIteratrive(xEnd, yEnd, zEnd, ZDEM, nxArray, nyArray, nzArray, csz, interpOption)
      # get the normal at this location
      nxEnd, nyEnd, nzEnd = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
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
      # turn off curvature with the  curvAcceleration coefficient
      effAccNorm = gravAccNorm + curvAcceleration * accNormCurv
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
      # get new velocity magnitude (leave 0 if uMag is 0)
      # this is important because uMag is first used to compute tau
      uMag = norm(ux, uy, uz)
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
          else:
            tau = 0.0

      # adding bottom shear resistance contribution
      # make sure uMag is not 0
      if uMag<velMagMin:
        uMag = velMagMin
      forceBotTang = - areaPart * tau
      if explicitFriction>0:
        # explicit formulation
        forceFrict[j] = forceFrict[j] - forceBotTang
      else:
        forceFrict[j] = forceFrict[j] - forceBotTang/uMag

      # compute entrained mass
      entrMassCell = entrMassRaster[indCellY, indCellX]
      dm, areaEntrPart = computeEntMassAndForce(dt, entrMassCell, areaPart, uMag, tau, entEroEnergy, rhoEnt)
      # update velocity
      ux = ux * m / (m + dm)
      uy = uy * m / (m + dm)
      uz = uz * m / (m + dm)
      # update mass
      m = m + dm
      mass[j] = m
      dM[j] = dm

      # speed loss due to energy loss due to entrained mass
      dEnergyEntr = areaEntrPart * entShearResistance + dm * entDefResistance
      dis = 1.0 - dEnergyEntr / (0.5 * m * (uMag*uMag + velMagMin))
      if dis < 0.0:
        dis = 0.0
      # update velocity
      ux = ux * dis
      uy = uy * dis
      uz = uz * dis

      # adding resistance force due to obstacles
      cResCell = cResRaster[indCellY][indCellX]
      cResPart = computeResForce(hRes, h, areaPart, rho, cResCell, uMag, explicitFriction)
      forceFrict[j] = forceFrict[j] - cResPart

      uxArray[j] = ux
      uyArray[j] = uy
      uzArray[j] = uz

  # save results
  force['dM'] = np.asarray(dM)
  force['forceX'] = np.asarray(forceX)
  force['forceY'] = np.asarray(forceY)
  force['forceZ'] = np.asarray(forceZ)
  force['forceFrict'] = np.asarray(forceFrict)
  particles['ux'] = np.asarray(uxArray)
  particles['uy'] = np.asarray(uyArray)
  particles['uz'] = np.asarray(uzArray)
  particles['m'] = np.asarray(mass)

  # update mass available for entrainement
  # TODO: this allows to entrain more mass then available...
  for j in range(Npart):
    indCellX = indXDEM[j]
    indCellY = indYDEM[j]
    entrMassCell = entrMassRaster[indCellY, indCellX]
    areaCell = areaRatser[indCellY, indCellX]
    dm = dM[j]
    # update surface entrainment mass available
    entrMassCell = entrMassCell - dm/areaCell
    if entrMassCell < 0:
      entrMassCell = 0
    entrMassRaster[indCellY, indCellX] = entrMassCell
  fields['entrMassRaster'] = np.asarray(entrMassRaster)

  return particles, force, fields


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef (double, double) computeEntMassAndForce(double dt, double entrMassCell, double areaPart, double uMag, double tau, double entEroEnergy, double rhoEnt):
  """ compute force component due to entrained mass

  Parameters
  ----------
  entrMassCell : float
      available mass for entrainement
  areaPart : float
      particle area
  uMag : float
      particle speed (velocity magnitude)
  tau : float
      bottom shear stress

  Returns
  -------
  dm : float
      entrained mass
  areaEntrPart : float
      Area for entrainement energy loss computation
  """
  cdef double width, ABotSwiped, areaEntrPart
  # compute entrained mass
  cdef double dm = 0
  if entrMassCell > 0:
      # either erosion or ploughing but not both

      if(entEroEnergy > 0):
          # erosion: erode according to shear and erosion energy
          dm = areaPart * tau * uMag * dt / entEroEnergy
          areaEntrPart = areaPart
      else:
          # ploughing in at avalanche front: erode full area weight
          # mass available in the cell [kg/m²]
          # width of the particle
          width = math.sqrt(areaPart)
          # bottom area covered by the particle during dt
          ABotSwiped = width * uMag * dt
          dm = entrMassCell * ABotSwiped
          areaEntrPart = entrMassCell / rhoEnt

  return dm, areaEntrPart


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double computeResForce(double hRes, double h, double areaPart, double rho, double cResCell, double uMag, int explicitFriction):
  """ compute force component due to resistance

  Parameters
  ----------
  hRes: float
      resistance height
  h : float
      particle flow depth
  areaPart : float
      particle area
  rho : float
      snow density
  cResCell : float
      resisance coefficient of cell
  uMag : float
      particle speed (velocity magnitude)
  explicitFriction: int
    if 1 add resistance with an explicit method. Implicit otherwise

  Returns
  -------
  cResPart : float
      resistance component for particle
  """
  cdef double hResEff = hRes
  cdef double cRecResPart
  if(h < hRes):
      hResEff = h
  if explicitFriction>0:
    # explicit formulation
    cRecResPart = - rho * areaPart * hResEff * cResCell * uMag * uMag
  else:
    cRecResPart = - rho * areaPart * hResEff * cResCell * uMag
  return cRecResPart


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def updatePositionC(cfg, particles, dem, force, DT):
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

  cdef double dt = DT
  cdef double stopCrit = cfg.getfloat('stopCrit')
  cdef double uFlowingThreshold = cfg.getfloat('uFlowingThreshold')
  log.debug('dt used now is %f' % DT)
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef double velMagMin = cfg.getfloat('velMagMin')
  cdef double rho = cfg.getfloat('rho')
  cdef int interpOption = cfg.getint('interpOption')
  cdef int explicitFriction = cfg.getint('explicitFriction')
  cdef double csz = dem['header'].cellsize
  cdef int Npart = particles['Npart']
  cdef double[:, :] nxArray = dem['Nx']
  cdef double[:, :] nyArray = dem['Ny']
  cdef double[:, :] nzArray = dem['Nz']
  cdef double[:, :] ZDEM = dem['rasterData']
  # initializeinter = np.zeros(N, dtype=np.float64)
  cdef double[:] mass = particles['m']
  cdef double[:] S = particles['s']
  cdef double[:] L = particles['l']
  cdef double[:] xArray = particles['x']
  cdef double[:] yArray = particles['y']
  cdef double[:] zArray = particles['z']
  cdef double[:] uxArray = particles['ux']
  cdef double[:] uyArray = particles['uy']
  cdef double[:] uzArray = particles['uz']
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
  cdef double m, h, x, y, z, s, l, ux, uy, uz, nx, ny, nz, dtStop
  cdef double ForceDriveX, ForceDriveY, ForceDriveZ, zeroCrossing
  cdef double mNew, xNew, yNew, zNew, uxNew, uyNew, uzNew, sNew, lNew, uN, uMag, uMagNew
  cdef double[:] mNewArray = np.zeros(Npart, dtype=np.float64)
  cdef double[:] xNewArray = np.zeros(Npart, dtype=np.float64)
  cdef double[:] yNewArray = np.zeros(Npart, dtype=np.float64)
  cdef double[:] zNewArray = np.zeros(Npart, dtype=np.float64)
  cdef double[:] sNewArray = np.zeros(Npart, dtype=np.float64)
  cdef double[:] uxArrayNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] uyArrayNew = np.zeros(Npart, dtype=np.float64)
  cdef double[:] uzArrayNew = np.zeros(Npart, dtype=np.float64)
  cdef int Lxy[4]
  cdef double w[4]
  cdef int LxyNew[4]
  cdef double wNew[4]
  cdef double massEntrained = 0, massFlowing = 0
  cdef int j
  # loop on particles
  for j in range(Npart):
    m = mass[j]
    x = xArray[j]
    y = yArray[j]
    z = zArray[j]
    ux = uxArray[j]
    uy = uyArray[j]
    uz = uzArray[j]
    s = S[j]
    l = L[j]

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

    # take friction force into account
    uxNew, uyNew, uzNew, dtStop = account4FrictionForce(uxNew, uyNew, uzNew, m, dt, forceFrict[j], uMag, explicitFriction)

    # update mass (already done un computeForceC)
    mNew = m
    massEntrained = massEntrained + dM[j]
    # update position
    xNew = x + dtStop * 0.5 * (ux + uxNew)
    yNew = y + dtStop * 0.5 * (uy + uyNew)
    zNew = z + dtStop * 0.5 * (uz + uzNew)
    # make sure particle is on the mesh (normal reprojection!!)
    xNew, yNew, LxyNew[0], LxyNew[1], LxyNew[2], LxyNew[3], wNew[0], wNew[1], wNew[2], wNew[3] = normalProjectionIteratrive(xNew, yNew, zNew, ZDEM, nxArray, nyArray, nzArray, csz, interpOption)
    zNew = getScalar(LxyNew[0], LxyNew[1], LxyNew[2], LxyNew[3], wNew[0], wNew[1], wNew[2], wNew[3], ZDEM)
    nxNew, nyNew, nzNew = getVector(LxyNew[0], LxyNew[1], LxyNew[2], LxyNew[3], wNew[0], wNew[1], wNew[2], wNew[3], nxArray, nyArray, nzArray)
    nxNew, nyNew, nzNew = normalize(nxNew, nyNew, nzNew)
    # compute the distance traveled by the particle
    lNew = l + norm((xNew-x), (yNew-y), (zNew-z))
    # compute average normal between the beginning and end of the time step
    Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, csz, interpOption)
    nx, ny, nz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
    nx, ny, nz = normalize(nx, ny, nz)
    nx, ny, nz = normalize(nx + nxNew, ny + nyNew, nz + nzNew)
    # compute the horizontal distance traveled by the particle (corrected with
    # the angle difference between the slope and the normal)
    sNew = s + nzNew*norm((xNew-x), (yNew-y), (zNew-z))
    # velocity magnitude
    uMag = norm(uxNew, uyNew, uzNew)
    # normal component of the velocity
    # uN = scalProd(uxNew, uyNew, uzNew, nx, ny, nz)
    uN = uxNew*nxNew + uyNew*nyNew + uzNew*nzNew
    # remove normal component of the velocity
    uxNew = uxNew - uN * nxNew
    uyNew = uyNew - uN * nyNew
    uzNew = uzNew - uN * nzNew

    # velocity magnitude new
    uMagNew = norm(uxNew, uyNew, uzNew)

    if uMag > 0.0:
      # ensure that velocitity magnitude stays the same also after reprojection onto terrain
      uxNew = uxNew * uMag / (uMagNew + velMagMin)
      uyNew = uyNew * uMag / (uMagNew + velMagMin)
      uzNew = uzNew * uMag / (uMagNew + velMagMin)

    # prepare for stopping criterion
    if uMag > uFlowingThreshold:
      # if velocity is bigger then threshold add to flowing mass
      massFlowing = massFlowing + mNew

    TotkinEneNew = TotkinEneNew + 0.5 * m * uMag * uMag
    TotpotEneNew = TotpotEneNew + mNew * gravAcc * zNew

    xNewArray[j] = xNew
    yNewArray[j] = yNew
    zNewArray[j] = zNew
    uxArrayNew[j] = uxNew
    uyArrayNew[j] = uyNew
    uzArrayNew[j] = uzNew
    sNewArray[j] = sNew
    mNewArray[j] = mNew
  particles['ux'] = np.asarray(uxArrayNew)
  particles['uy'] = np.asarray(uyArrayNew)
  particles['uz'] = np.asarray(uzArrayNew)
  particles['s'] = np.asarray(sNewArray)
  particles['m'] = np.asarray(mNewArray)
  particles['mTot'] = np.sum(particles['m'])
  particles['x'] = np.asarray(xNewArray)
  particles['y'] = np.asarray(yNewArray)
  particles['z'] = np.asarray(zNewArray)
  particles['massEntrained'] = np.asarray(massEntrained)
  log.debug('Entrained DFA mass: %s kg', np.asarray(massEntrained))
  particles['kineticEne'] = TotkinEneNew
  particles['potentialEne'] = TotpotEneNew
  if peakKinEne < TotkinEneNew:
    particles['peakKinEne'] = TotkinEneNew
  if particles['peakMassFlowing'] < massFlowing:
    particles['peakMassFlowing'] = massFlowing

  if cfg['stopCritType'] == 'massFlowing':
    value = massFlowing
    peakValue = particles['peakMassFlowing']
    stop = (value <= stopCrit*peakValue) and peakValue > 0
  elif cfg['stopCritType'] == 'kinEnergy':
    value = TotkinEneNew
    peakValue = particles['peakKinEne']
    stop = value <= stopCrit*peakValue
  else:
    message = 'stopCritType %s is not defined' % cfg['stopCritType']
    log.error(message)
    raise AssertionError(message)

  if stop:
    particles['iterate'] = False
    log.info('stopping because of %s stopCriterion.' % (cfg['stopCritType']))

  # make sure particle is on the mesh (recompute the z component)
  # particles, _ = geoTrans.projectOnRaster(dem, particles, interp='bilinear')
  #################################################################
  # this is dangerous!!!!!!!!!!!!!!
  ###############################################################
  # remove particles that are not located on the mesh any more
  particles = DFAtls.removeOutPart(cfg, particles, dem, dt)

  # split particles with too much mass
  particles = DFAtls.splitPart(cfg, particles, dem)
  return particles


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef (double, double, double, double) account4FrictionForce(double ux, double uy, double uz, double m, double dt, double forceFrict, double uMag, int explicitFriction):
  """ update velocity with friction force

  Parameters
  ----------
  ux: float
      x velocity
  uy: float
      y velocity
  uz: float
      z velocity
  m : float
      particle area
  dt : float
      snow density
  forceFrict : float
      resisance coefficient of cell
  uMag : float
      particle speed (velocity magnitude)
  explicitFriction: int
    if 1 add resistance with an explicit method. Implicit otherwise

  Returns
  -------
  cResPart : float
      resistance component for particle
  """
  cdef double xDir, yDir, zDir, uxNew, uyNew, uzNew, uMagNew, dtStop

  if explicitFriction>0:
    # explicite method
    uMagNew = norm(ux, uy, uz)
    # will friction force stop the particle
    if uMagNew<dt*forceFrict/m:
      # stop the particle
      uxNew = 0
      uyNew = 0
      uzNew = 0
      # particle stops after
      if uMag<=0:
        dtStop = 0
      else:
        dtStop = m * uMagNew / (dt * forceFrict)
    else:
      # add friction force in the opposite direction of the motion
      xDir, yDir, zDir = normalize(ux, uy, uz)
      uxNew = ux - xDir * forceFrict * dt / m
      uyNew = uy - yDir * forceFrict * dt / m
      uzNew = uz - zDir * forceFrict * dt / m
      dtStop = dt
  else:
    # implicite method
    uxNew = ux / (1.0 + dt * forceFrict / m)
    uyNew = uy / (1.0 + dt * forceFrict / m)
    uzNew = uz / (1.0 + dt * forceFrict / m)
    dtStop = dt

  return uxNew, uyNew, uzNew, dtStop


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
  cdef double[:, :] areaRaster = dem['areaRaster']
  ncols = header.ncols
  nrows = header.nrows
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
  cdef double[:] xArray = particles['x']
  cdef double[:] yArray = particles['y']
  cdef double[:] uxArray = particles['ux']
  cdef double[:] uyArray = particles['uy']
  cdef double[:] uzArray = particles['uz']
  cdef double[:, :] PFV = fields['pfv']
  cdef double[:, :] PP = fields['ppr']
  cdef double[:, :] PFD = fields['pfd']
  cdef int nrow = int(nrows)
  cdef int ncol = int(ncols)
  cdef int Lxy[4]
  cdef double w[4]
  cdef double m, h, x, y, z, s, ux, uy, uz, nx, ny, nz, hbb, hLim, areaPart
  cdef double xllc = 0
  cdef double yllc = 0
  cdef double csz = CSZ
  cdef int Npart = np.size(particles['x'])
  cdef double[:] hBB = np.zeros((Npart))
  cdef int j, i
  cdef int indx, indy

  for j in range(Npart):
    x = xArray[j]
    y = yArray[j]
    ux = uxArray[j]
    uy = uyArray[j]
    uz = uzArray[j]
    m = mass[j]
    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    # find coordinates of the 4 nearest cornes on the raster
    # prepare for bilinear interpolation
    Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, csz, interpOption)
    # add the component of the points value to the 4 neighbour grid points
    for i in range(4):
      indx = Lxy[i%2] # 0 1 0 1
      indy = Lxy[i/2+2] # 2 2 3 3
      FDBilinear[indy, indx] = FDBilinear[indy, indx] + m * w[i]
      MomBilinearX[indy, indx] = MomBilinearX[indy, indx] + m * ux * w[i]
      MomBilinearY[indy, indx] = MomBilinearY[indy, indx] + m * uy * w[i]
      MomBilinearZ[indy, indx] = MomBilinearZ[indy, indx] + m * uz * w[i]

  for i in range(ncol):
    for j in range(nrow):
      m = FDBilinear[j, i]
      if m > 0:
        # TODO: here we devide by the area of the vertex, would it not make
        # more sense to devide by the area of the cell in the previous loop?
        FDBilinear[j, i] = m / (areaRaster[j, i] * rho)
        VXBilinear[j, i] = MomBilinearX[j, i]/m
        VYBilinear[j, i] = MomBilinearY[j, i]/m
        VZBilinear[j, i] = MomBilinearZ[j, i]/m
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
    x = xArray[j]
    y = yArray[j]
    Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, csz, interpOption)
    hbb = getScalar(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], FDBilinear)
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
    """ Locate particles on DEM and neighbour search grid

    Ĺocate each particle in a grid cell (both DEM and neighbour search grid) and build the
    indPartInCell and partInCell arrays for SPH neighbourSearch and
    InCell, IndX and IndY arrays location on the DEM grid.
    See issue #200 and documentation for details

    Parameters
    ----------
    particles : dict
    dem : dict
      dem dict with neighbour search grid header (information about neighbour search grid)

    Returns
    -------
    particles : dict
      updated particles dictionary with
      neighbours info indPartInCell and partInCell arrays for neighbour search and
      with InCell, IndX and IndY arrays location on the DEM grid
    """
    # get DEM grid information
    header = dem['header']
    cdef int nColsDEM = header.ncols
    cdef int nRowsDEM = header.nrows
    cdef float cszDEM = header.cellsize
    # get neighbour search grid information
    headerNeighbourGrid = dem['headerNeighbourGrid']
    cdef int nColsNeighbourGrid = headerNeighbourGrid.ncols
    cdef int nRowsNeighbourGrid = headerNeighbourGrid.nrows
    cdef float cszNeighbourGrid = headerNeighbourGrid.cellsize
    # get particle location
    cdef int Npart = particles['Npart']
    cdef int j
    cdef double[:] xArray = particles['x']
    cdef double[:] yArray = particles['y']

    # initialize outputs
    cdef int nCellsNeighbourGrid = (nColsNeighbourGrid-1)*(nRowsNeighbourGrid-1)
    cdef int[:] indPartInCell = np.zeros(nCellsNeighbourGrid + 1).astype('intc')
    cdef int[:] indPartInCell2 = np.zeros(nCellsNeighbourGrid + 1).astype('intc')
    cdef int[:] partInCell = np.zeros(Npart).astype('intc')
    cdef int[:] indXDEM = np.zeros(Npart).astype('intc')
    cdef int[:] indYDEM = np.zeros(Npart).astype('intc')
    cdef int[:] inCellDEM = np.zeros(Npart).astype('intc')
    # Count number of particles in each SPH grid cell
    cdef int indx, indy, ic
    for j in range(Npart):
      indx = int(xArray[j] / cszNeighbourGrid)
      indy = int(yArray[j] / cszNeighbourGrid)
      # get index of cell containing the particle
      ic = indx + (nColsNeighbourGrid-1) * indy
      indPartInCell[ic+1] = indPartInCell[ic+1] + 1
    for j in range(nCellsNeighbourGrid):
      indPartInCell[j+1] = indPartInCell[j] + indPartInCell[j+1]
      indPartInCell2[j+1] = indPartInCell[j+1]

    # make the list of which particles are in which cell
    for j in range(Npart):
        indx = int(xArray[j] / cszNeighbourGrid)
        indy = int(yArray[j] / cszNeighbourGrid)
        ic = indx + (nColsNeighbourGrid-1) * indy
        partInCell[indPartInCell2[ic+1]-1] = j
        indPartInCell2[ic+1] = indPartInCell2[ic+1] - 1
        indXDEM[j] = int(xArray[j] / cszDEM)
        indYDEM[j] = int(yArray[j] / cszDEM)
        # get index of cell containing the particle
        inCellDEM[j] = indXDEM[j] + (nColsDEM-1) * indYDEM[j]

    particles['inCellDEM'] = np.asarray(inCellDEM)
    particles['indXDEM'] = np.asarray(indXDEM)
    particles['indYDEM'] = np.asarray(indYDEM)
    particles['indPartInCell'] = np.asarray(indPartInCell)
    particles['partInCell'] = np.asarray(partInCell)

    return particles


def computeForceSPHC(cfg, particles, force, dem, gradient=0):
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
  # get SPH grid information
  headerNeighbourGrid = dem['headerNeighbourGrid']
  # get DEM header and normal field
  header = dem['header']
  cszNormal = header.cellsize
  nxArray = dem['Nx']
  nyArray = dem['Ny']
  nzArray = dem['Nz']

  forceSPHX, forceSPHY, forceSPHZ = computeGradC(cfg, particles, headerNeighbourGrid, cszNormal, nxArray, nyArray, nzArray, gradient)
  forceSPHX = np.asarray(forceSPHX)
  forceSPHY = np.asarray(forceSPHY)
  forceSPHZ = np.asarray(forceSPHZ)

  force['forceSPHX'] = forceSPHX
  force['forceSPHY'] = forceSPHY
  force['forceSPHZ'] = forceSPHZ

  return particles, force

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def computeGradC(cfg, particles, headerNeighbourGrid, double cszNormal, double[:, :] nxArray, double[:, :] nyArray,
                 double[:, :] nzArray, gradient):
  """ compute lateral forces acting on the particles (SPH component)

  Cython implementation

  Parameters
  ----------
  cfg: configparser
      configuration for DFA simulation
  particles : dict
      particles dictionary at t
  headerNeighbourGrid : dict
      neighbour search header dictionary (information about SPH grid)
  cszNormal : double
      cell size of normal vector Grid
  nxArray : 2D numpy array
      x component of the normal vector of the DEM
  nyArray : 2D numpy array
      y component of the normal vector of the DEM
  nzArray : 2D numpy array
      z component of the normal vector of the DEM
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
  # get all inputs
  # configuration parameters
  cdef double rho = cfg.getfloat('rho')
  cdef double minRKern = cfg.getfloat('minRKern')
  cdef double velMagMin = cfg.getfloat('velMagMin')
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef int interpOption = cfg.getint('interpOption')
  cdef int SPHoption = cfg.getint('sphOption')
  # neighbour search grid information and neighbour information
  cdef double cszNeighbourGrid = headerNeighbourGrid.cellsize
  cdef int nRowsNeighbourGrid = headerNeighbourGrid.nrows
  cdef int nColsNeighbourGrid = headerNeighbourGrid.ncols
  cdef int[:] indPartInCell = particles['indPartInCell']
  cdef int[:] partInCell = particles['partInCell']
  # SPH kernel
  # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
  cdef double rKernel = cszNeighbourGrid
  cdef double facKernel = 10.0 / (math.pi * rKernel * rKernel * rKernel * rKernel * rKernel)
  cdef double dfacKernel = - 3.0 * facKernel
  # particle information
  cdef double[:] mass = particles['m']
  cdef double[:] h = particles['h']
  cdef double[:] xArray = particles['x']
  cdef double[:] yArray = particles['y']
  cdef double[:] zArray = particles['z']
  cdef double[:] uxArray = particles['ux']
  cdef double[:] uyArray = particles['uy']
  cdef double[:] uzArray = particles['uz']
  cdef int N = xArray.shape[0]

  # initialize variables and outputs
  cdef double[:] GHX = np.zeros(N, dtype=np.float64)
  cdef double[:] GHY = np.zeros(N, dtype=np.float64)
  cdef double[:] GHZ = np.zeros(N, dtype=np.float64)
  cdef double K1 = 1
  cdef double K2 = 1
  cdef double gravAcc3
  cdef double gradhX, gradhY, gradhZ, uMag, nx, ny, nz, G1, G2, mdwdrr
  cdef double g1, g2, g11, g12, g22, g33
  cdef double m11, m12, m22, GG1, GG2
  cdef double x, y, z, ux, uy, uz, vx, vy, wx, wy, uxOrtho, uyOrtho, uzOrtho
  cdef double dx, dy, dz, dn, r, hr, dwdr, wKernel
  cdef int Lxy[4]
  cdef double w[4]
  cdef int lInd, rInd
  cdef int indx, indy
  cdef int j, ic, n, p, l, imax, imin, iPstart, iPend
  cdef int grad = gradient

  # loop on particles
  for j in range(N):
    gradhX = 0
    gradhY = 0
    gradhZ = 0
    G1 = 0
    G2 = 0
    m11 = 0
    m12 = 0
    m22 = 0
    x = xArray[j]
    y = yArray[j]
    z = zArray[j]
    ux = uxArray[j]
    uy = uyArray[j]
    uz = uzArray[j]
    # locate particle in SPH grid
    indx = int(x / cszNeighbourGrid)
    indy = int(y / cszNeighbourGrid)

    if SPHoption > 1:
      # get normal vector
      Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, cszNormal, interpOption)
      nx, ny, nz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
      nx, ny, nz = normalize(nx, ny, nz)
      # projection of gravity on normal vector
      gravAcc3 = scalProd(nx, ny, nz, 0, 0, gravAcc)
      uMag = norm(ux, uy, uz)
      if uMag < velMagMin:
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

    # check if we are on the bottom ot top row!!!
    lInd = -1
    rInd = 2
    if indy == 0:
        lInd = 0
    if indy == nRowsNeighbourGrid - 1:
        rInd = 1
    for n in range(lInd, rInd):
        ic = (indx - 1) + (nColsNeighbourGrid-1) * (indy + n)
        # make sure not to take particles from the other edge
        imax = max(ic, (nColsNeighbourGrid-1) * (indy + n))
        imin = min(ic+3, (nColsNeighbourGrid-1) * (indy + n + 1))
        iPstart = indPartInCell[imax]
        iPend = indPartInCell[imin]
        # loop on all particles in neighbour boxes
        for p in range(iPstart, iPend):
            # index of particle in neighbour box
            l = partInCell[p]
            if j != l:
                dx = xArray[l] - x
                dy = yArray[l] - y
                dz = zArray[l] - z
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
                      wKernel = facKernel * hr * hr * hr
                      m11 = m11 + mass[l] / h[l] * wKernel * r1 * r1 / rho
                      m12 = m12 + mass[l] / h[l] * wKernel * r1 * r2 / rho
                      m22 = m22 + mass[l] / h[l] * wKernel * r2 * r2 / rho
                      G1 = G1 + mass[l] * (1 - h[j]/h[l]) * wKernel * K1*r1
                      G2 = G2 + mass[l] * (1 - h[j]/h[l]) * wKernel * K2*r2

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

  return GHX, GHY, GHZ


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def computeFDC(cfg, particles, header, double[:, :] nxArray, double[:, :] nyArray, double[:, :] nzArray, int[:] indX, int[:] indY):
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
  cdef double[:] xArray = particles['x']
  cdef double[:] yArray = particles['y']
  cdef double[:] zArray = particles['z']
  cdef double[:] H1 = particles['h']
  cdef double[:] uxArray = particles['ux']
  cdef double[:] uyArray = particles['uy']
  cdef double[:] uzArray = particles['uz']
  # cdef double[:] gradx = particles['gradx']
  # cdef double[:] grady = particles['grady']
  # cdef double[:] gradz = particles['gradz']
  cdef int[:] indPartInCell = particles['indPartInCell']
  cdef int[:] partInCell = particles['partInCell']
  cdef int N = xArray.shape[0]
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
  cdef double x, y, z
  cdef double dx, dy, dz, r, hr, dwdr, massl, hl, ww
  cdef int Lxy[4]
  cdef double w[4]
  cdef int lInd = -1
  cdef int rInd = 2
  cdef int indx
  cdef int indy
  cdef int ic, n, p, l, imax, imin, iPstart, iPend
  # With loop
  # loop on particles
  # TcpuSPH = 0
  # Tcpuadd = 0
  for j in range(N):
    x = xArray[j]
    y = yArray[j]
    z = zArray[j]
    h = 0
    ww = 0
    gx = 0
    gy = 0
    gz = 0
    indx = indX[j]
    indy = indY[j]
    Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, csz, interpOption)
    nx, ny, nz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
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
                dx = xArray[l] - x
                dy = yArray[l] - y
                dz = zArray[l] - z
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
                    wKernel = facKernel * hr * hr * hr
                    massl = mass[l]
                    hl = H1[l]
                    dh = massl * wKernel
                    # gx = gx + massl * w / hl * dx
                    # gy = gy + massl * w / hl * dy
                    # gz = gz + massl * w / hl * dz
                    ww = ww + massl * wKernel / hl
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
    return np.asarray(norm(x, y, z))

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
    return np.asarray(norm2(x, y, z))


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
  x, y, z = normalize(x, y, z)
  return np.asarray(x), np.asarray(y), np.asarray(z)


@cython.cdivision(True)
cdef (double, double, double) croosProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute cross product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  cdef double wx = uy * vz - uz * vy
  cdef double wy = uz * vx - ux * vz
  cdef double wz = ux * vy - uy * vx
  return wx, wy, wz


def crossProdpy(x, y, z, u, v, w): # <-- small wrapper to expose croosProd() to Python
  wx, wy, wz = croosProd(x, y, z, u, v, w)
  return np.asarray(wx), np.asarray(wy), np.asarray(wz)


cdef double scalProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute scalar product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  return ux*vx + uy*vy + uz*vz

def scalProdpy(x, y, z, u, v, w): # <-- small wrapper to expose scalProd() to Python
  s = scalProd(x, y, z, u, v, w)
  return np.asarray(s)


@cython.cdivision(True)
cdef (int, int, int, int, double, double, double, double) getWeights(double x, double y, double csz, int interpOption):#, int *Lxy, double *w):
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
  cdef int Lxy[4]
  cdef double w[4]
  cdef double Lx, Ly
  cdef double xllc = 0.
  cdef double yllc = 0.

  # find coordinates in normalized ref (origin (0,0) and cellsize 1)
  Lx = (x - xllc) / csz
  Ly = (y - yllc) / csz
  # find coordinates of the 4 nearest cornes on the raster
  Lxy[0] = <int>math.floor(Lx)
  Lxy[2] = <int>math.floor(Ly)
  Lxy[1] = Lxy[0] + 1
  Lxy[3] = Lxy[2] + 1
  # prepare for bilinear interpolation
  dx = Lx - Lxy[0]
  dy = Ly - Lxy[2]

  if interpOption == 0:
    dx = 1.*math.round(dx)
    dy = 1.*math.round(dy)
  elif interpOption == 1:
    dx = 1./2.
    dy = 1./2.

  # lower left
  w[0] = (1-dx)*(1-dy)
  # lower right
  w[1] = dx*(1-dy)
  # uper left
  w[2] = (1-dx)*dy
  # and uper right
  w[3] = dx*dy

  return Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3]


def getWeightspy(x, y, csz, interpOption): # <-- small wrapper to expose getWeightspy() to Python
  cdef int Lxy[4]
  cdef double w[4]
  Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, csz, interpOption)
  return np.asarray(Lxy[0]), np.asarray(Lxy[1]), np.asarray(Lxy[2]), np.asarray(Lxy[3]), np.asarray(w[0]), np.asarray(w[1]), np.asarray(w[2]), np.asarray(w[3])


@cython.wraparound(False)
@cython.cdivision(True)
cdef (double, double, int, int, int, int, double, double, double, double) normalProjectionIteratrive(
  double xOld, double yOld, double zOld, double[:,:] ZDEM, double[:,:] nxArray, double[:,:] nyArray, double[:,:] nzArray, double csz, int interpOption):
  """ Find the orthogonal projection of a point on a mesh

  Iterative method to find the projection of a point on a surface defined by its mesh
  Parameters
  ----------
      xOld: float
          x coordinate of the point to project
      yOld: float
          y coordinate of the point to project
      zOld: float
          z coordinate of the point to project
      ZDEM: 2D array
          z component of the DEM field at the grid nodes
      Nx: 2D array
          x component of the normal vector field at the grid nodes
      Ny: 2D array
          y component of the normal vector field at the grid nodes
      Nz: 2D array
          z component of the normal vector field at the grid nodes
      csz: float
          cellsize of the grid
      interpOption: int
          -0: nearest neighbour interpolation
          -1: equal weights interpolation
          -2: bilinear interpolation
    Returns
    -------
    xNew: float
        x coordinate of the projected point
    yNew: float
        y coordinate of the projected point
    Lxy: int[4]
        colomn and row indices for interpolation
    w: float[4]
        corresponding weights
    """
  cdef double zTemp, zn, xNew, yNew, zNew
  cdef double nx, ny, nz
  cdef int reprojectionIterations = 5
  cdef int Lxi, Lyi, Lxj, Lyj
  cdef int Lxy[4]
  cdef double w[4]
  xNew = xOld
  yNew = yOld
  zNew = zOld
  Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, csz, interpOption)
  # get the cell location of the point
  Lxi = Lxy[0]
  Lyi = Lxy[2]
  # iterate
  while reprojectionIterations > 0:
    reprojectionIterations = reprojectionIterations - 1
    # vertical projection of the point
    zTemp = getScalar(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], ZDEM)
    # normal vector at this vertical projection location
    # if we take the normal at the new particle position)
    # nx, ny, nz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], Nx, Ny, Nz)
    # What Peter does: get the normal of the cell center (but still not exactly giving the same results)
    nx, ny, nz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], 0.25, 0.25, 0.25, 0.25, nxArray, nyArray, nzArray)
    nx, ny, nz = normalize(nx, ny, nz)
    zn = (xNew-xNew) * nx + (yNew-yNew) * ny + (zTemp-zNew) * nz
    # correct position with the normal part
    xNew = xNew + zn * nx
    yNew = yNew + zn * ny
    zNew = zNew + zn * nz
    # What I think is better for normal projection
    # # Normal component of the vector between the initial and projected point
    # zn = (xNew-xOld) * nx + (yNew-yOld) * ny + (zTemp-zNew) * nz
    # # correct position with the normal part
    # xNew = xOld + zn * nx
    # yNew = yOld + zn * ny
    # zNew = zOld + zn * nz
    # get cell
    Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, csz, interpOption)
    Lxj = Lxy[0]
    Lyj = Lxy[2]
    # are we in the same cell?
    if Lxi==Lxj and Lyi==Lyj:
      return xNew, yNew, Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3]
    Lxj = Lxi
    Lyj = Lyi
  return xNew, yNew, Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3]



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def projOnRaster(double[:] xArray, double[:] yArray, double[:, :] vArray, double csz, int interpOption):
  """ Interpolate vector field from grid to points
  """
  cdef int N = xArray.shape[0]
  cdef double x, y
  cdef int Lxy[4]
  cdef double w[4]
  cdef int j
  cdef double[:] v = np.zeros(N)
  for j in range(N):
    x = xArray[j]
    y = yArray[j]
    Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3] = getWeights(x, y, csz, interpOption)
    v[j] = getScalar(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], vArray)

  return np.asarray(v)


@cython.boundscheck(False)
cdef double getScalar(int Lxy0, int Lxy1, int Lxy2, int Lxy3, double w0, double w1, double w2, double w3, double[:, :] V):
  """ Interpolate vector field from grid to single point location

  Originaly created to get the normal vector at location (x,y) given the
  normal vector field on the grid. Grid has its origin in (0,0).
  Can be used to interpolate any vector field.
  Interpolation using a bilinear interpolation

  Parameters
  ----------
      Lxy: int[4]
          colomn and row indices for interpolation
      w: float[4]
          corresponding weights
      V: 2D numpy array
          scalar field at the grid nodes

  Returns
  -------
      v: float
          interpolated scalar at position (x, y)
  """
  cdef double v = (V[Lxy2, Lxy0]*w0 +
                   V[Lxy2, Lxy1]*w1 +
                   V[Lxy3, Lxy0]*w2 +
                   V[Lxy3, Lxy1]*w3)


  return v


@cython.boundscheck(False)
cdef (double, double, double) getVector(
  int Lxy0, int Lxy1, int Lxy2, int Lxy3, double w0, double w1, double w2, double w3,
  double[:, :] Nx, double[:, :] Ny, double[:, :] Nz):
  """ Interpolate vector field from grid to single point location

  Originaly created to get the normal vector at location (x,y) given the
  normal vector field on the grid. Grid has its origin in (0,0).
  Can be used to interpolate any vector field.
  Interpolation using a bilinear interpolation

  Parameters
  ----------
      Lxy: int[4]
          colomn and row indices for interpolation
      w: float[4]
          corresponding weights
          location in the y location of desiered interpolation
      Nx: 2D numpy array
          x component of the vector field at the grid nodes
      Ny: 2D numpy array
          y component of the vector field at the grid nodes
      Nz: 2D numpy array
          z component of the vector field at the grid nodes

  Returns
  -------
      nx: float
          x component of the interpolated vector field at position (x, y)
      ny: float
          y component of the interpolated vector field at position (x, y)
      nz: float
          z component of the interpolated vector field at position (x, y)
  """
  cdef double nx = (Nx[Lxy2, Lxy0]*w0 +
                   Nx[Lxy2, Lxy1]*w1 +
                   Nx[Lxy3, Lxy0]*w2 +
                   Nx[Lxy3, Lxy1]*w3)
  cdef double ny = (Ny[Lxy2, Lxy0]*w0 +
                   Ny[Lxy2, Lxy1]*w1 +
                   Ny[Lxy3, Lxy0]*w2 +
                   Ny[Lxy3, Lxy1]*w3)
  cdef double nz = (Nz[Lxy2, Lxy0]*w0 +
                   Nz[Lxy2, Lxy1]*w1 +
                   Nz[Lxy3, Lxy0]*w2 +
                   Nz[Lxy3, Lxy1]*w3)
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
