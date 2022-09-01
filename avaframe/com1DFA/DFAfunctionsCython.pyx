#!python
# cython: boundscheck=False, wraparound=False, cdivision=True
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

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.com1DFA.particleTools as particleTools
import avaframe.in3Utils.geoTrans as geoTrans


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


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
    cdef int nPart = len(xArray)
    cdef int k, ic

    for k in range(nPart):
      x = xArray[k]
      y = yArray[k]
      z = zArray[k]
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


def computeForceC(cfg, particles, fields, dem, int frictType):
  """ compute forces acting on the particles (without the SPH component)

  Cython implementation implementation

  Parameters
  ----------
  cfg: configparser
      configuration for DFA simulation
  particles : dict
      particles dictionary at t
  fields: dict
    fields dictionary
  dem : dict
      dictionary with dem information
  frictType: int
    identifier for friction law to be used

  Returns
  -------
  force : dict
      force dictionary
  """
  # read input parameters
  cdef double tau0 = cfg.getfloat('tau0')
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
  cdef double xsi = cfg.getfloat('xsi')
  cdef double curvAccInFriction = cfg.getfloat('curvAccInFriction')
  cdef double curvAccInTangent = cfg.getfloat('curvAccInTangent')
  cdef int curvAccInGradient = cfg.getint('curvAccInGradient')
  cdef double velMagMin = cfg.getfloat('velMagMin')
  cdef int interpOption = cfg.getint('interpOption')
  cdef int explicitFriction = cfg.getint('explicitFriction')
  cdef int reprojMethod = cfg.getint('reprojMethodForce')
  cdef int reprojectionIterations = cfg.getint('reprojectionIterations')
  cdef double thresholdProjection = cfg.getfloat('thresholdProjection')
  cdef double subgridMixingFactor = cfg.getfloat('subgridMixingFactor')
  cdef int viscOption = cfg.getint('viscOption')
  cdef double dt = particles['dt']
  cdef double mu = cfg.getfloat('mu')
  cdef int nPart = particles['nPart']
  cdef double csz = dem['header']['cellsize']
  cdef int nrows = dem['header']['nrows']
  cdef int ncols = dem['header']['ncols']
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
  cdef double[:] forceX = np.zeros(nPart, dtype=np.float64)
  cdef double[:] forceY = np.zeros(nPart, dtype=np.float64)
  cdef double[:] forceZ = np.zeros(nPart, dtype=np.float64)
  cdef double[:] forceFrict = np.zeros(nPart, dtype=np.float64)
  cdef double[:] dM = np.zeros(nPart, dtype=np.float64)
  cdef double[:] gEff = np.zeros(nPart, dtype=np.float64)
  cdef double[:] curvAcc = np.zeros(nPart, dtype=np.float64)
  # declare intermediate step variables
  cdef int indCellX, indCellY
  cdef double areaPart, areaCell, araEntrPart, cResCell, cResPart, uMag, m, dm, h, entrMassCell, dEnergyEntr, dis
  cdef double x, y, z, xEnd, yEnd, zEnd, ux, uy, uz, uxDir, uyDir, uzDir
  cdef double nx, ny, nz, nxEnd, nyEnd, nzEnd, nxAvg, nyAvg, nzAvg
  cdef double gravAccNorm, accNormCurv, effAccNorm, gravAccTangX, gravAccTangY, gravAccTangZ, forceBotTang, sigmaB, tau
  # variables for interpolation
  cdef int Lx0, Ly0, LxEnd0, LyEnd0, iCell, iCellEnd
  cdef double w[4]
  cdef double wEnd[4]
  cdef int k
  force = {}
  # loop on particles
  for k in range(nPart):
      m = mass[k]
      x = xArray[k]
      y = yArray[k]
      z = zArray[k]
      h = hArray[k]
      ux = uxArray[k]
      uy = uyArray[k]
      uz = uzArray[k]
      indCellX = indXDEM[k]
      indCellY = indYDEM[k]
      # deduce area
      areaPart = m / (h * rho)

      # get cell and weights
      Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = getCellAndWeights(x, y, ncols, nrows, csz, interpOption)

      # get normal at the particle location
      nx, ny, nz = getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
      nx, ny, nz = normalize(nx, ny, nz)

      if viscOption == 1:
        # add artificial viscosity
        ux, uy, uz = addArtificialViscosity(m, h, dt, rho, ux, uy, uz, subgridMixingFactor, Lx0, Ly0,
                                            w[0], w[1], w[2], w[3], VX, VY, VZ, nx, ny, nz)

      # get normal at the particle estimated end location
      xEnd = x + dt * ux
      yEnd = y + dt * uy
      zEnd = z + dt * uz
      # this point is not necessarily on the surface, project it on the surface
      if reprojMethod == 0:
        # Project vertically on the dem
        iCellEnd = getCells(xEnd, yEnd, ncols, nrows, csz)
        if iCellEnd >= 0:
          LxEnd0, LyEnd0, wEnd[0], wEnd[1], wEnd[2], wEnd[3] = getWeights(xEnd, yEnd, iCellEnd, csz, ncols, interpOption)
      elif reprojMethod == 1:
        # project trying to keep the travelled distance constant
        xEnd, yEnd, zEnd, iCellEnd, LxEnd0, LyEnd0, wEnd[0], wEnd[1], wEnd[2], wEnd[3] = distConservProjectionIteratrive(
          x, y, z, ZDEM, nxArray, nyArray, nzArray, xEnd, yEnd, zEnd, csz, ncols, nrows, interpOption,
          reprojectionIterations, thresholdProjection)
      elif reprojMethod == 2:
        # project using samos method
        xEnd, yEnd, iCellEnd, LxEnd0, LyEnd0, wEnd[0], wEnd[1], wEnd[2], wEnd[3] = samosProjectionIteratrive(
          xEnd, yEnd, zEnd, ZDEM, nxArray, nyArray, nzArray, csz, ncols, nrows, interpOption, reprojectionIterations)

      if iCellEnd < 0:
        # if not on the DEM take x, y as end point
        LxEnd0 = Lx0
        LyEnd0 = Ly0
        wEnd = w

      # get the normal at this location
      nxEnd, nyEnd, nzEnd = getVector(LxEnd0, LyEnd0, wEnd[0], wEnd[1], wEnd[2], wEnd[3], nxArray, nyArray, nzArray)
      nxEnd, nyEnd, nzEnd = normalize(nxEnd, nyEnd, nzEnd)
      # get average of those normals
      nxAvg = nx + nxEnd
      nyAvg = ny + nyEnd
      nzAvg = nz + nzEnd
      nxAvg, nyAvg, nzAvg = normalize(nxAvg, nyAvg, nzAvg)

      # acceleration due to curvature
      accNormCurv = (ux*(nxEnd-nx) + uy*(nyEnd-ny) + uz*(nzEnd-nz)) / dt
      # normal component of the acceleration of gravity
      gravAccNorm = - gravAcc * nzAvg  # use nzAvg because we want the average gNorm on the time step
      # turn off curvature with the  curvAccInFriction coefficient
      effAccNorm = gravAccNorm + curvAccInFriction * accNormCurv
      # save the normal component of the gravity (plus the curvature term if desiered)
      # this is needed to compute the pressure gradient
      # save the norm of the gEff, because we know in which direction it goes (minus the normal vector)
      if curvAccInGradient == 1:
        # use gravity plus curvature (if the effAccNorm > 0, material is detatched, then gEff = 0)
        if(effAccNorm <= 0.0):
          gEff[k] = -effAccNorm
        else:
          gEff[k] = 0
      else:
        # only use gravity
        gEff[k] = -gravAccNorm

      # body forces (tangential component of acceleration of gravity with the term due to curvature)
      gravAccTangX = - (gravAccNorm + curvAccInTangent * accNormCurv) * nxAvg
      gravAccTangY = - (gravAccNorm + curvAccInTangent * accNormCurv) * nyAvg
      gravAccTangZ = - gravAcc - (gravAccNorm + curvAccInTangent * accNormCurv) * nzAvg
      # adding gravity force contribution
      forceX[k] = forceX[k] + gravAccTangX * m
      forceY[k] = forceY[k] + gravAccTangY * m
      forceZ[k] = forceZ[k] + gravAccTangZ * m

      # Calculating bottom shear and normal stress
      # get new velocity magnitude (leave 0 if uMag is 0)
      # this is important because uMag is first used to compute tau
      uMag = norm(ux, uy, uz)
      if(effAccNorm > 0.0):
          # if fluid detatched
          # log.info('fluid detatched for particle %s' % j)
          tau = 0.0
      else:
          # bottom normal stress sigmaB
          sigmaB = - effAccNorm * rho * h
          if frictType == 1:
            # SamosAT friction type (bottom shear stress)
            tau = SamosATfric(rho, tau0, Rs0, mu, kappa, B, R, uMag, sigmaB, h)
          elif frictType == 2:
            # coulomb friction type (bottom shear stress)
            tau = mu * sigmaB
          elif frictType == 3:
            # voellmy friction type
            tau = mu * sigmaB + rho * uMag * uMag * gravAcc / xsi
          else:
            tau = 0.0

      # adding bottom shear resistance contribution
      # norm of the friction force >=0 (if 0 -> detatchment)
      forceBotTang = - areaPart * tau
      if explicitFriction == 1:
        # explicit formulation
        forceFrict[k] = forceFrict[k] - forceBotTang
      elif explicitFriction == 0:
        # make sure uMag is not 0
        if uMag<velMagMin:
          uMag = velMagMin
        forceFrict[k] = forceFrict[k] - forceBotTang/uMag

      # compute entrained mass
      entrMassCell = entrMassRaster[indCellY, indCellX]
      dm, areaEntrPart = computeEntMassAndForce(dt, entrMassCell, areaPart, uMag, tau, entEroEnergy, rhoEnt)
      # update velocity
      ux = ux * m / (m + dm)
      uy = uy * m / (m + dm)
      uz = uz * m / (m + dm)
      # update mass
      m = m + dm
      mass[k] = m
      dM[k] = dm

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
      forceFrict[k] = forceFrict[k] - cResPart

      uxArray[k] = ux
      uyArray[k] = uy
      uzArray[k] = uz

  # save results
  force['dM'] = np.asarray(dM)
  force['forceX'] = np.asarray(forceX)
  force['forceY'] = np.asarray(forceY)
  force['forceZ'] = np.asarray(forceZ)
  force['forceFrict'] = np.asarray(forceFrict)
  particles['gEff'] = np.asarray(gEff)
  particles['curvAcc'] = np.asarray(curvAcc)
  particles['ux'] = np.asarray(uxArray)
  particles['uy'] = np.asarray(uyArray)
  particles['uz'] = np.asarray(uzArray)
  particles['m'] = np.asarray(mass)

  # update mass available for entrainement
  # TODO: this allows to entrain more mass then available...
  for k in range(nPart):
    indCellX = indXDEM[k]
    indCellY = indYDEM[k]
    entrMassCell = entrMassRaster[indCellY, indCellX]
    areaCell = areaRatser[indCellY, indCellX]
    dm = dM[k]
    # update surface entrainment mass available
    entrMassCell = entrMassCell - dm/areaCell
    if entrMassCell < 0:
      entrMassCell = 0
    entrMassRaster[indCellY, indCellX] = entrMassCell
  fields['entrMassRaster'] = np.asarray(entrMassRaster)

  return particles, force, fields


cpdef (double, double) computeEntMassAndForce(double dt, double entrMassCell,
                                              double areaPart, double uMag,
                                              double tau, double entEroEnergy,
                                              double rhoEnt):
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


cpdef double computeResForce(double hRes, double h, double areaPart, double rho,
                             double cResCell, double uMag, int explicitFriction):
  """ compute force component due to resistance

  Parameters
  ----------
  hRes: float
      resistance height
  h : float
      particle flow thickness
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
  if explicitFriction == 1:
    # explicit formulation
    cRecResPart = - rho * areaPart * hResEff * cResCell * uMag * uMag
  elif explicitFriction == 0:
    cRecResPart = - rho * areaPart * hResEff * cResCell * uMag
  return cRecResPart


cdef (double, double, double) addArtificialViscosity(double m, double h, double dt, double rho,
                                                     double ux, double uy, double uz, double subgridMixingFactor,
                                                     int Lx0, int Ly0, double w0, double w1, double w2, double w3,
                                                     double[:, :] VX, double[:, :] VY, double[:, :] VZ,
                                                     double nx, double ny, double nz):
  """ add artificial viscosity

  Add the artificial viscosity in an implicit way and this before adding the other forces.

  Parameters
  ----------
  m : float
      mass of the particle
  h : float
      flow thickness of the particle
  dt : float
      time step
  rho : float
      density
  subgridMixingFactor : float
      viscosity coefficient
  Lx0: int
      colomn of the nearest lower left cell
  Ly0: int
      row of the nearest lower left cell
  w: float[4]
      corresponding weights
      location in the y location of desiered interpolation
  VX: 2D numpy array
      x component of the velocity vector field at the grid nodes
  VY: 2D numpy array
      y component of the velocity vector field at the grid nodes
  VZ: 2D numpy array
      z component of the velocity vector field at the grid nodes
  nx: float
      x component of the interpolated vector field at particle position
  ny: float
      y component of the interpolated vector field at particle position
  nz: float
      z component of the interpolated vector field at particle position

  Returns
  -------
  ux: float
    x component of the velocity uptated with the viscous force
  uy: float
    y component of the uptated with the viscous force
  uz: float
    z component of the uptated with the viscous force
  """
  cdef double vMeanx, vMeany, vMeanz, vMeanNorm, dvX, dvY, dvZ
  vMeanx, vMeany, vMeanz = getVector(Lx0, Ly0, w0, w1, w2, w3, VX, VY, VZ)
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
  return ux, uy, uz


def updatePositionC(cfg, particles, dem, force, int typeStop=0):
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
  typeStop: int
    0 if standard stopping criterion, if 1 stopping criterion based on SPHforce - used for iniStep
  Returns
  -------
  particles : dict
      particles dictionary at t + dt
  """

  # read input parameters
  cdef double stopCrit = cfg.getfloat('stopCrit')
  cdef double stopCritIni = cfg.getfloat('stopCritIni')
  cdef double stopCritIniSmall = cfg.getfloat('stopCritIniSmall')
  cdef double uFlowingThreshold = cfg.getfloat('uFlowingThreshold')
  cdef double dt = particles['dt']
  log.debug('dt used now is %f' % dt)
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef double velMagMin = cfg.getfloat('velMagMin')
  cdef double rho = cfg.getfloat('rho')
  cdef int interpOption = cfg.getint('interpOption')
  cdef int explicitFriction = cfg.getint('explicitFriction')
  cdef int reprojMethod = cfg.getint('reprojMethodPosition')
  cdef int reprojectionIterations = cfg.getint('reprojectionIterations')
  cdef double thresholdProjection = cfg.getfloat('thresholdProjection')
  cdef double centeredPosition = cfg.getfloat('centeredPosition')
  cdef double csz = dem['header']['cellsize']
  cdef int nrows = dem['header']['nrows']
  cdef int ncols = dem['header']['ncols']
  cdef int nPart = particles['nPart']
  cdef double[:, :] nxArray = dem['Nx']
  cdef double[:, :] nyArray = dem['Ny']
  cdef double[:, :] nzArray = dem['Nz']
  cdef double[:, :] ZDEM = dem['rasterData']
  cdef np.ndarray[np.uint8_t, ndim = 1, cast=True] outOfDEM
  outOfDEM = np.array(dem['outOfDEM'], dtype=bool)
  # read particles and fields
  cdef double[:] mass = particles['m']
  cdef double[:] idFixed = particles['idFixed']
  cdef double[:] sArray = particles['s']
  cdef double[:] sCorArray = particles['sCor']
  cdef double[:] lArray = particles['l']
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
  cdef double peakForceSPH = particles['peakForceSPH']
  cdef double totKForceSPH = particles['forceSPHIni']
  # initialize outputs
  cdef double TotkinEneNew = 0
  cdef double TotpotEneNew = 0
  cdef double totForceSPHNew = 0
  cdef double[:] mNewArray = np.zeros(nPart, dtype=np.float64)
  cdef double[:] xNewArray = np.zeros(nPart, dtype=np.float64)
  cdef double[:] yNewArray = np.zeros(nPart, dtype=np.float64)
  cdef double[:] zNewArray = np.zeros(nPart, dtype=np.float64)
  cdef double[:] sNewArray = np.zeros(nPart, dtype=np.float64)
  cdef double[:] sCorNewArray = np.zeros(nPart, dtype=np.float64)
  cdef double[:] lNewArray = np.zeros(nPart, dtype=np.float64)
  cdef double[:] uxArrayNew = np.zeros(nPart, dtype=np.float64)
  cdef double[:] uyArrayNew = np.zeros(nPart, dtype=np.float64)
  cdef double[:] uzArrayNew = np.zeros(nPart, dtype=np.float64)
  cdef int[:] keepParticle = np.ones(nPart, dtype=np.int32)
  # declare intermediate step variables
  cdef double m, h, x, y, z, sCor, s, l, ux, uy, uz, nx, ny, nz, dtStop, idfixed
  cdef double mNew, xNew, yNew, zNew, uxNew, uyNew, uzNew
  cdef double sCorNew, sNew, lNew, ds, dl, uN, uMag, uMagNew
  cdef double ForceDriveX, ForceDriveY, ForceDriveZ
  cdef double massEntrained = 0, massFlowing = 0
  cdef int k
  cdef int nRemove = 0
  # variables for interpolation
  cdef int Lx0, Ly0, LxNew0, LyNew0, iCell, iCellNew
  cdef double w[4]
  cdef double wNew[4]
  # loop on particles
  for k in range(nPart):
    m = mass[k]
    x = xArray[k]
    y = yArray[k]
    z = zArray[k]
    ux = uxArray[k]
    uy = uyArray[k]
    uz = uzArray[k]
    s = sArray[k]
    sCor = sCorArray[k]
    l = lArray[k]
    idfixed = idFixed[k]

    # Force magnitude (without friction)
    ForceDriveX = forceX[k] + forceSPHX[k]
    ForceDriveY = forceY[k] + forceSPHY[k]
    ForceDriveZ = forceZ[k] + forceSPHZ[k]

    # velocity magnitude
    uMag = norm(ux, uy, uz)

    # procede to time integration
    # operator splitting
    # estimate new velocity due to driving force
    uxNew = ux + ForceDriveX * dt / m
    uyNew = uy + ForceDriveY * dt / m
    uzNew = uz + ForceDriveZ * dt / m

    # take friction force into account
    if typeStop != 1:
      uxNew, uyNew, uzNew, dtStop = account4FrictionForce(uxNew, uyNew, uzNew, m, dt, forceFrict[k], uMag, explicitFriction)
    else:
      dtStop = dt

    # update mass (already done un computeForceC)
    mNew = m
    massEntrained = massEntrained + dM[k]
    # update position
    if centeredPosition:
      xNew = x + dtStop * 0.5 * (ux + uxNew)
      yNew = y + dtStop * 0.5 * (uy + uyNew)
      zNew = z + dtStop * 0.5 * (uz + uzNew)
    else:
      xNew = x + dtStop * uxNew
      yNew = y + dtStop * uyNew
      zNew = z + dtStop * uzNew
    # make sure particle is on the mesh (normal reprojection!!)

    if reprojMethod == 0:
      # Project vertically on the dem
      iCellNew = getCells(xNew, yNew, ncols, nrows, csz)
      if iCellNew >= 0:
        Lx0, Ly0, wNew[0], wNew[1], wNew[2], wNew[3] = getWeights(xNew, yNew, iCellNew, csz, ncols, interpOption)
        zNew = getScalar(Lx0, Ly0, wNew[0], wNew[1], wNew[2], wNew[3], ZDEM)

    elif reprojMethod == 1:
      # project trying to keep the travelled distance constant
      xNew, yNew, zNew, iCellNew, LxNew0, LyNew0, wNew[0], wNew[1], wNew[2], wNew[3] = distConservProjectionIteratrive(
        x, y, z, ZDEM, nxArray, nyArray, nzArray, xNew, yNew, zNew, csz, ncols, nrows, interpOption,
        reprojectionIterations, thresholdProjection)
    elif reprojMethod == 2:
      # project using samos method
      xNew, yNew, iCellNew, LxNew0, LyNew0, wNew[0], wNew[1], wNew[2], wNew[3] = samosProjectionIteratrive(
        xNew, yNew, zNew, ZDEM, nxArray, nyArray, nzArray, csz, ncols, nrows, interpOption, reprojectionIterations)
      zNew = getScalar(LxNew0, LyNew0, wNew[0], wNew[1], wNew[2], wNew[3], ZDEM)

    if iCellNew < 0:
      # if the particle is not on the DEM, memorize it and remove it at the next update
      keepParticle[k] = 0
      LxNew0 = 0
      LyNew0 = 0
      wNew = [0, 0, 0, 0]
      nRemove = nRemove + 1
    elif outOfDEM[iCellNew]:
      # if the particle is on the DEM but in a noData area,
      # memorize it and remove it at the next update
      keepParticle[k] = 0
      nRemove = nRemove + 1

    # get cell and weights at old position
    Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = getCellAndWeights(x, y, ncols, nrows, csz, interpOption)
    # get normal at the old particle location
    nx, ny, nz = getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
    # get normal at the new particle location
    nxNew, nyNew, nzNew = getVector(LxNew0, LyNew0, wNew[0], wNew[1], wNew[2], wNew[3], nxArray, nyArray, nzArray)
    # get average normal between old and new position
    nx, ny, nz = normalize(nx+nxNew, ny+nyNew, nz+nzNew)
    # normalize new normal
    nxNew, nyNew, nzNew = normalize(nxNew, nyNew, nzNew)
    # compute the distance traveled by the particle
    dl = norm((xNew-x), (yNew-y), (zNew-z))
    lNew = l + dl
    # compute the horizontal distance traveled by the particle
    ds = norm((xNew-x), (yNew-y), 0)
    sNew = s + ds
    # compute the horizontal distance traveled by the particle (corrected with
    # the angle difference between the slope and the normal)
    sCorNew = sCor + nz*dl
    # velocity magnitude
    uMag = norm(uxNew, uyNew, uzNew)
    # normal component of the velocity
    # uN = scalProd(uxNew, uyNew, uzNew, nxNew, nyNew, nzNew)
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
    totForceSPHNew = totForceSPHNew + mNew * norm(forceSPHX[k], forceSPHY[k], forceSPHZ[k])

    if idfixed == 1:
      # idfixed = 1 if particles belong to 'fixed' boundary particles - so zero velocity and fixed position
      xNewArray[k] = x
      yNewArray[k] = y
      zNewArray[k] = z
      uxArrayNew[k] = ux
      uyArrayNew[k] = uy
      uzArrayNew[k] = uz
      sNewArray[k] = s
      sCorNewArray[k] = s
      mNewArray[k] = m
    else:
      # idfixed = 0 particles belong to the actual releae area
      xNewArray[k] = xNew
      yNewArray[k] = yNew
      zNewArray[k] = zNew
      uxArrayNew[k] = uxNew
      uyArrayNew[k] = uyNew
      uzArrayNew[k] = uzNew
      sNewArray[k] = sNew
      sCorNewArray[k] = sCorNew
      mNewArray[k] = mNew

  particles['ux'] = np.asarray(uxArrayNew)
  particles['uy'] = np.asarray(uyArrayNew)
  particles['uz'] = np.asarray(uzArrayNew)
  particles['l'] = np.asarray(lNewArray)
  particles['s'] = np.asarray(sNewArray)
  particles['sCor'] = np.asarray(sCorNewArray)
  particles['m'] = np.asarray(mNewArray)
  particles['mTot'] = np.sum(particles['m'])
  particles['x'] = np.asarray(xNewArray)
  particles['y'] = np.asarray(yNewArray)
  particles['z'] = np.asarray(zNewArray)
  particles['massEntrained'] = massEntrained
  log.debug('Entrained DFA mass: %s kg', np.asarray(massEntrained))
  particles['kineticEne'] = TotkinEneNew
  particles['potentialEne'] = TotpotEneNew

  if typeStop == 1:
    # typeStop = 1 for initialisation step where particles are redistributed to reduce SPH force
    # here stop criterion based on SPHForce
    value = totForceSPHNew
    peakValue = particles['peakForceSPH']
    oldValue = particles['forceSPHIni']
    particles['forceSPHIni'] = totForceSPHNew
    if oldValue == 0.0:
      stop = False
    elif value < 1.:
      stop = oldValue/value < stopCritIniSmall
    else:
      stop = value <= stopCritIni*peakValue
      log.debug('SPHFORCE value %f and stop value %f' % (totForceSPHNew, stopCritIni*peakValue))
    if peakForceSPH < totForceSPHNew:
      particles['peakForceSPH'] = totForceSPHNew
  else:
    # avalanche computation stop criterion based on kinetic energy or massFlowing
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
    if typeStop == 1:
      log.debug('stopping initial particle distribution')
    else:
      log.info('stopping because of %s stopCriterion.' % (cfg['stopCritType']))

  # remove particles that are not located on the mesh any more
  if nRemove > 0:
    mask = np.array(np.asarray(keepParticle), dtype=bool)
    particles = particleTools.removePart(particles, mask, nRemove, 'because they exited the domain')

  return particles


cpdef (double, double, double, double) account4FrictionForce(double ux, double uy,
                                                            double uz, double m,
                                                            double dt, double forceFrict,
                                                            double uMag, int explicitFriction):
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
      time step
  forceFrict : float
      friction force (actually a force for the explicit method, force/vel for implicit)
  uMag : float
      particle speed (velocity magnitude)
  explicitFriction: int
    if 1 add resistance with an explicit method. Implicit otherwise

  Returns
  -------
  uxNew: float
      x velocity
  uyNew: float
      y velocity
  uzNew: float
      z velocity
  dtStop : float
      time step (smaller then dt if the particle stops during this time step)
  """
  cdef double xDir, yDir, zDir, uxNew, uyNew, uzNew, uMagNew, dtStop

  if explicitFriction == 1:
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
        dtStop = m * uMagNew / (forceFrict)
    else:
      # add friction force in the opposite direction of the motion
      xDir, yDir, zDir = normalize(ux, uy, uz)
      uxNew = ux - xDir * forceFrict * dt / m
      uyNew = uy - yDir * forceFrict * dt / m
      uzNew = uz - zDir * forceFrict * dt / m
      dtStop = dt
  elif explicitFriction == 0:
    # implicite method
    uxNew = ux / (1.0 + dt * forceFrict / m)
    uyNew = uy / (1.0 + dt * forceFrict / m)
    uzNew = uz / (1.0 + dt * forceFrict / m)
    dtStop = dt

  return uxNew, uyNew, uzNew, dtStop


def updateFieldsC(cfg, particles, dem, fields):
  """ update fields and particles flow thickness

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
  # read input parameters
  cdef double rho = cfg.getfloat('rho')
  cdef int interpOption = cfg.getint('interpOption')
  header = dem['header']
  cdef int nrows = header['nrows']
  cdef int ncols = header['ncols']
  cdef double xllc = 0
  cdef double yllc = 0
  cdef double csz = header['cellsize']
  cdef int nPart = np.size(particles['x'])
  cdef double[:, :] areaRaster = dem['areaRaster']
  # read particles and fields
  cdef double[:] mass = particles['m']
  cdef double[:] xArray = particles['x']
  cdef double[:] yArray = particles['y']
  cdef double[:] uxArray = particles['ux']
  cdef double[:] uyArray = particles['uy']
  cdef double[:] uzArray = particles['uz']
  cdef double[:] travelAngleArray = particles['travelAngle']
  cdef bint computeTA = fields['computeTA']
  cdef bint computeKE = fields['computeKE']
  cdef bint computeP = fields['computeP']
  cdef double[:, :] PFV = fields['pfv']
  cdef double[:, :] PP = fields['ppr']
  cdef double[:, :] PFT = fields['pft']
  cdef double[:, :] PTA = fields['pta']
  cdef double[:, :] PKE = fields['pke']
  # initialize outputs
  cdef double[:, :] MassBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] VBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] PBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] FTBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] MomBilinearX = np.zeros((nrows, ncols))
  cdef double[:, :] MomBilinearY = np.zeros((nrows, ncols))
  cdef double[:, :] MomBilinearZ = np.zeros((nrows, ncols))
  cdef double[:, :] VXBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] VYBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] VZBilinear = np.zeros((nrows, ncols))
  cdef double[:, :] kineticEnergy = np.zeros((nrows, ncols))
  cdef double[:, :] travelAngleField = np.zeros((nrows, ncols))
  # declare intermediate step variables
  cdef double[:] hBB = np.zeros((nPart))
  cdef double m, h, x, y, z, s, ux, uy, uz, nx, ny, nz, hbb, hLim, areaPart, travelAngle
  cdef int k, i
  cdef int indx, indy
  cdef int ind1[4]
  cdef int ind2[4]
  ind1[:] = [0, 1, 0, 1]
  ind2[:] = [0, 0, 1, 1]
  # variables for interpolation
  cdef int Lx0, Ly0, iCell
  cdef double w[4]
  cdef double mwi

  for k in range(nPart):
    x = xArray[k]
    y = yArray[k]
    ux = uxArray[k]
    uy = uyArray[k]
    uz = uzArray[k]
    m = mass[k]
    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    # find coordinates of the 4 nearest cornes on the raster
    # prepare for bilinear interpolation
    Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = getCellAndWeights(x, y, ncols, nrows, csz, interpOption)
    # for the travel angle we simply do a nearest interpolation
    indx = <int>math.round(x / csz)
    indy = <int>math.round(y / csz)
    if computeTA:
      travelAngle = travelAngleArray[k]
      travelAngleField[indy, indx] = max(travelAngleField[indy, indx], travelAngle)
    # add the component of the points value to the 4 neighbour grid points
    # TODO : check if giving the arrays [0 1 0 1].. is faster
    for i in range(4):
      indx = Lx0 + ind1[i]
      indy = Ly0 + ind2[i]
      mwi = m * w[i]
      MassBilinear[indy, indx] = MassBilinear[indy, indx] + mwi
      MomBilinearX[indy, indx] = MomBilinearX[indy, indx] + mwi * ux
      MomBilinearY[indy, indx] = MomBilinearY[indy, indx] + mwi * uy
      MomBilinearZ[indy, indx] = MomBilinearZ[indy, indx] + mwi * uz

  for i in range(ncols):
    for j in range(nrows):
      m = MassBilinear[j, i]
      if m > 0:
        # TODO: here we devide by the area of the vertex, would it not make
        # more sense to devide by the area of the cell in the previous loop?
        FTBilinear[j, i] = m / (areaRaster[j, i] * rho)
        VXBilinear[j, i] = MomBilinearX[j, i]/m
        VYBilinear[j, i] = MomBilinearY[j, i]/m
        VZBilinear[j, i] = MomBilinearZ[j, i]/m
        VBilinear[j, i] = norm(VXBilinear[j, i], VYBilinear[j, i], VZBilinear[j, i])
        if VBilinear[j, i] > PFV[j, i]:
          PFV[j, i] = VBilinear[j, i]
        if FTBilinear[j, i] > PFT[j, i]:
          PFT[j, i] = FTBilinear[j, i]
        if computeP:
          PBilinear[j, i] = computePressure(VBilinear[j, i], rho)
          if PBilinear[j, i] > PP[j, i]:
            PP[j, i] = PBilinear[j, i]
        if computeTA:
          if travelAngleField[j, i] > PTA[j, i]:
            PTA[j, i] = travelAngleField[j, i]
        if computeKE:
          # in J/cell (this is not normalized yet and depends on the cell size used for the computation)
          kineticEnergy[j, i] = 0.5*m*VBilinear[j, i]*VBilinear[j, i]
          if kineticEnergy[j, i] > PKE[j, i]:
            PKE[j, i] = kineticEnergy[j, i]

  fields['FM'] = np.asarray(MassBilinear)
  fields['FV'] = np.asarray(VBilinear)
  fields['Vx'] = np.asarray(VXBilinear)
  fields['Vy'] = np.asarray(VYBilinear)
  fields['Vz'] = np.asarray(VZBilinear)
  fields['FT'] = np.asarray(FTBilinear)
  fields['pfv'] = np.asarray(PFV)
  fields['pft'] = np.asarray(PFT)
  if computeP:
    fields['ppr'] = np.asarray(PP)
    fields['P'] = np.asarray(PBilinear)
  if computeTA:
    fields['TA'] = np.asarray(travelAngleField)
    fields['pta'] = np.asarray(PTA)
  if computeKE:
    fields['pke'] = np.asarray(PKE)


  for k in range(nPart):
    x = xArray[k]
    y = yArray[k]
    Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = getCellAndWeights(x, y, ncols, nrows, csz, interpOption)
    hbb = getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], FTBilinear)
    hBB[k] = hbb

  particles['h'] = np.asarray(hBB)

  return particles, fields


cpdef double computePressure(double v, double rho):
  """Compute pressure using the p = rho*v² equation

  Parameters
  ----------
  v : float
      velocity magnitude
  rho : float
      snow density

  Returns
  -------
  p : float
      pressure computed using p = rho*v²
  """
  cdef double p
  p = rho * v * v
  return p


def computeTravelAngleC(particles, zPartArray0):
  """Compute the travel angle associated to the particles

  Parameters
  ----------
  particles : dict
      particles dictionary at t
  zPartArray0 : dict
      z coordinate of particles at t=0s

  Returns
  -------
  particles : dict
      particles dictionary updated with the travel angle
  """
  cdef int[:] parentIDArray = particles['parentID'].astype('intc')
  cdef int nPart = particles['nPart']
  cdef double[:] zArray = particles['z']
  cdef double[:] sArray = particles['s']
  cdef double[:] z0Array = zPartArray0
  cdef double[:] gammaArray = np.zeros(nPart)
  cdef int parentID, j
  cdef double tanGamma, gamma, s, z, z0
  # get particle location
  # first compute travel angle for each particle
  for k in range(nPart):
    # get parent Id in order to  get the first z position
    parentID = parentIDArray[k]
    # get z0
    z0 = z0Array[parentID]
    # compute tan of the travel angle
    z = zArray[k]
    s = sArray[k]
    if s > 0:
      tanGamma = (z0 - z) / s
    else:
      tanGamma = 0
    # get the travel angle
    gamma = math.atan(tanGamma) * 180 / math.pi
    gammaArray[k] = gamma
  particles['travelAngle'] = np.asarray(gammaArray)
  return particles


def getNeighborsC(particles, dem):
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
    cdef int nColsDEM = header['ncols']
    cdef int nRowsDEM = header['nrows']
    cdef float cszDEM = header['cellsize']
    # get neighbour search grid information
    headerNeighbourGrid = dem['headerNeighbourGrid']
    cdef int nColsNeighbourGrid = headerNeighbourGrid['ncols']
    cdef int nRowsNeighbourGrid = headerNeighbourGrid['nrows']
    cdef float cszNeighbourGrid = headerNeighbourGrid['cellsize']
    # get particle location
    cdef int nPart = particles['nPart']
    cdef int k
    cdef double[:] xArray = particles['x']
    cdef double[:] yArray = particles['y']

    # initialize outputs
    cdef int nCellsNeighbourGrid = nColsNeighbourGrid*nRowsNeighbourGrid
    cdef int[:] indPartInCell = np.zeros(nCellsNeighbourGrid + 1).astype('intc')
    cdef int[:] indPartInCell2 = np.zeros(nCellsNeighbourGrid + 1).astype('intc')
    cdef int[:] partInCell = np.zeros(nPart).astype('intc')
    cdef int[:] indXDEM = np.zeros(nPart).astype('intc')
    cdef int[:] indYDEM = np.zeros(nPart).astype('intc')
    cdef int[:] inCellDEM = np.zeros(nPart).astype('intc')
    # Count number of particles in each SPH grid cell
    cdef int indx, indy, ic
    for k in range(nPart):
      indx = <int>math.round(xArray[k] / cszNeighbourGrid)
      indy = <int>math.round(yArray[k] / cszNeighbourGrid)
      # get index of cell containing the particle
      ic = indx + nColsNeighbourGrid * indy
      indPartInCell[ic+1] = indPartInCell[ic+1] + 1
    for ic in range(nCellsNeighbourGrid):
      indPartInCell[ic+1] = indPartInCell[ic] + indPartInCell[ic+1]
      indPartInCell2[ic+1] = indPartInCell[ic+1]

    # make the list of which particles are in which cell
    for k in range(nPart):
        indx = <int>math.round(xArray[k] / cszNeighbourGrid)
        indy = <int>math.round(yArray[k] / cszNeighbourGrid)
        ic = indx + nColsNeighbourGrid * indy
        partInCell[indPartInCell2[ic+1]-1] = k
        indPartInCell2[ic+1] = indPartInCell2[ic+1] - 1
        indXDEM[k] = <int>math.round(xArray[k] / cszDEM)
        indYDEM[k] = <int>math.round(yArray[k] / cszDEM)
        # get index of cell containing the particle
        inCellDEM[k] = indXDEM[k] + nColsDEM * indYDEM[k]

    particles['inCellDEM'] = np.asarray(inCellDEM)
    particles['indXDEM'] = np.asarray(indXDEM)
    particles['indYDEM'] = np.asarray(indYDEM)
    particles['indPartInCell'] = np.asarray(indPartInCell)
    particles['partInCell'] = np.asarray(partInCell)

    return particles


def computeForceSPHC(cfg, particles, force, dem, int sphOption, gradient=0):
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
  sphOption: int
      which sphOption should be use
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
  headerNormalGrid = dem['header']
  nxArray = dem['Nx']
  nyArray = dem['Ny']
  nzArray = dem['Nz']

  forceSPHX, forceSPHY, forceSPHZ = computeGradC(cfg, particles, headerNeighbourGrid, headerNormalGrid, nxArray, nyArray, nzArray, gradient, sphOption)
  forceSPHX = np.asarray(forceSPHX)
  forceSPHY = np.asarray(forceSPHY)
  forceSPHZ = np.asarray(forceSPHZ)

  force['forceSPHX'] = forceSPHX
  force['forceSPHY'] = forceSPHY
  force['forceSPHZ'] = forceSPHZ

  return particles, force


def computeGradC(cfg, particles, headerNeighbourGrid, headerNormalGrid, double[:, :] nxArray, double[:, :] nyArray,
                 double[:, :] nzArray, gradient, int SPHoption):
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
  headerNormalGrid : double
      normal grid header dictionary (information about the DEM grid)
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
  cdef int viscOption = cfg.getint('viscOption')

  # grid normal raster information
  cdef double cszNormal = headerNormalGrid['cellsize']
  cdef int nRowsNormal = headerNormalGrid['nrows']
  cdef int nColsNormal = headerNormalGrid['ncols']
  # neighbour search grid information and neighbour information
  cdef double cszNeighbourGrid = headerNeighbourGrid['cellsize']
  cdef int nRowsNeighbourGrid = headerNeighbourGrid['nrows']
  cdef int nColsNeighbourGrid = headerNeighbourGrid['ncols']
  cdef int[:] indPartInCell = particles['indPartInCell']
  cdef int[:] partInCell = particles['partInCell']
  # SPH kernel
  # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
  cdef double rKernel = cszNeighbourGrid
  cdef double facKernel = 10.0 / (math.pi * rKernel * rKernel * rKernel * rKernel * rKernel)
  cdef double dfacKernel = - 3.0 * facKernel
  # particle information
  cdef double[:] gEff = particles['gEff']
  cdef double[:] mass = particles['m']
  cdef double[:] hArray = particles['h']
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
  cdef double x, y, z, ux, uy, uz, vx, vy, wx, wy, uxOrtho, uyOrtho, uzOrtho
  cdef double dx, dy, dz, dux, duy, duz, dn, r, hr, dwdr, wKernel
  cdef int Lx0, Ly0, iCell
  cdef double w[4]
  cdef int lInd, rInd
  cdef int indx, indy
  cdef int k, ic, n, p, l, imax, imin, iPstart, iPend
  cdef int grad = gradient
  # artificial viscosity parameters
  cdef double dwdrr, area
  cdef double pikl, flux
  cdef double mk, ml, hk, hl, ck, cl, lambdakl

  # loop on particles
  for k in range(N):
    gradhX = 0
    gradhY = 0
    gradhZ = 0
    pikl = 0
    G1 = 0
    G2 = 0
    m11 = 0
    m12 = 0
    m22 = 0
    x = xArray[k]
    y = yArray[k]
    z = zArray[k]
    ux = uxArray[k]
    uy = uyArray[k]
    uz = uzArray[k]
    hk = hArray[k]
    mk = mass[k]

    # locate particle in SPH grid
    indx = <int>math.round(x / cszNeighbourGrid)
    indy = <int>math.round(y / cszNeighbourGrid)

    if SPHoption > 1:
        # get normal vector
        Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = getCellAndWeights(x, y, nColsNormal, nRowsNormal, cszNormal, interpOption)
        nx, ny, nz = getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
        nx, ny, nz = normalize(nx, ny, nz)
        # projection of gravity on normal vector or use the effective gravity (gravity + curvature acceleration part)
        # This was done I computeForce and is passed over here
        gravAcc3 = gEff[k]
        if SPHoption > 2:
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

                uxOrtho, uyOrtho, uzOrtho = crossProd(nx, ny, nz, ux, uy, uz)
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
        ic = (indx - 1) + nColsNeighbourGrid * (indy + n)
        # make sure not to take particles from the other edge
        imax = max(ic, nColsNeighbourGrid * (indy + n))
        imin = min(ic+3, nColsNeighbourGrid * (indy + n + 1))
        iPstart = indPartInCell[imax]
        iPend = indPartInCell[imin]
        # loop on all particles in neighbour boxes
        for p in range(iPstart, iPend):
            # index of particle in neighbour box
            l = partInCell[p]
            if k != l:
                dx = xArray[l] - x
                dy = yArray[l] - y
                dz = zArray[l] - z

#----------------------------SPHOPTION = 1--------------------------------------
                # SamosAT style (no reprojecion on the surface, dz = 0 and gz is used)
                if SPHoption == 1:
                    dz = 0
                    # get norm of r = xj - xl
                    r = norm(dx, dy, dz)
                    if r < minRKern * rKernel:
                        # impose a minimum distance between particles
                        r = minRKern * rKernel
                    if r < rKernel:
                        hr = rKernel - r
                        dwdr = dfacKernel * hr * hr
                        mdwdrr = mass[l] * dwdr / r
                        gradhX = gradhX + mdwdrr*dx
                        gradhY = gradhY + mdwdrr*dy
                        gradhZ = gradhZ + mdwdrr*dz
                        gravAcc3 = gravAcc

#----------------------------SPHOPTION = 2--------------------------------------
                # Compute the gradient in the cartesian coord system with reprojecion on the surface
                # and dz != 0 and g3 are used
                # It is here also possible to add the artificial viscosity comming from the Ata theory
                elif SPHoption == 2:
                    # like option 1 but with dz!=0
                    # get norm of r = xk - xl
                    r = norm(dx, dy, dz)
                    if r < minRKern * rKernel:
                        # impose a minimum distance between particles
                        dx = minRKern * rKernel * dx
                        dy = minRKern * rKernel * dy
                        dz = minRKern * rKernel * dz
                        r = minRKern * rKernel
                    if r < rKernel:
                        hl = hArray[l]
                        hr = rKernel - r
                        dwdrr = dfacKernel * hr * hr / r
                        ml = mass[l]
                        area = ml/(rho*hl)
                        flux = gravAcc3*hl
                        if viscOption == 2:
                            # ATA artificial viscosity
                            dux = uxArray[l] - ux
                            duy = uyArray[l] - uy
                            duz = uzArray[l] - uz
                            ck = math.sqrt(gravAcc3*hk)
                            cl = math.sqrt(gravAcc3*hl)
                            lamdbakl = (ck+cl)/2
                            pikl = - lamdbakl * scalProd(dux, duy, duz, dx, dy, dz) / r
                        # SPH gradient computation - standard SPH formulation
                        gradhX = gradhX + (flux + pikl)*dwdrr*dx*area
                        gradhY = gradhY + (flux + pikl)*dwdrr*dy*area
                        gradhZ = gradhZ + (flux + pikl)*dwdrr*dz*area


#----------------------------SPHOPTION = 3--------------------------------------
                # Compute the gradient in the local coord system (will allow us to add the earth pressure coef)
                # and this time reprojecion on the surface, dz != 0 and g3 are used
                if SPHoption == 3:
                  # get coordinates in local coord system
                  r1 = scalProd(dx, dy, dz, ux, uy, uz)
                  r2 = scalProd(dx, dy, dz, uxOrtho, uyOrtho, uzOrtho)
                  # impse r3=0 even if the particle is not exactly on the tengent plane
                  # get norm of r = xk - xl
                  r = norm(r1, r2, 0)
                  if r < minRKern * rKernel:
                      # impose a minimum distance between particles
                      r1 = minRKern * rKernel * r1
                      r2 = minRKern * rKernel * r2
                      r = minRKern * rKernel
                  if r < rKernel:
                      hr = rKernel - r
                      dwdr = dfacKernel * hr * hr
                      ml = mass[l]
                      mdwdrr = ml * dwdr / r
                      G1 = mdwdrr * K1*r1
                      G2 = mdwdrr * K2*r2

                      gradhX = gradhX + ux*G1 + uxOrtho*G2
                      gradhY = gradhY + uy*G1 + uyOrtho*G2
                      gradhZ = gradhZ + (- g1*(ux*G1 + uxOrtho*G2) - g2*(uy*G1 + uyOrtho*G2))


    if grad == 1:
      if SPHoption == 2:
        GHX[k] = GHX[k] + gradhX / gravAcc3
        GHY[k] = GHY[k] + gradhY / gravAcc3
        GHZ[k] = GHZ[k] + gradhZ / gravAcc3
      else:
        GHX[k] = GHX[k] - gradhX / rho
        GHY[k] = GHY[k] - gradhY / rho
        GHZ[k] = GHZ[k] - gradhZ / rho
    elif grad == 0:
      if SPHoption == 2:
        # grad here is g*grad(h)
        GHX[k] = GHX[k] + gradhX * mk
        GHY[k] = GHY[k] + gradhY * mk
        GHZ[k] = GHZ[k] + gradhZ * mk
      else :
        GHX[k] = GHX[k] + gradhX*gravAcc3 / rho* mk
        GHY[k] = GHY[k] + gradhY*gravAcc3 / rho* mk
        GHZ[k] = GHZ[k] + gradhZ*gravAcc3 / rho* mk


  return GHX, GHY, GHZ


cpdef double norm(double x, double y, double z):
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

cpdef double norm2(double x, double y, double z):
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


cpdef (double, double, double) normalize(double x, double y, double z):
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


cpdef (double, double, double) crossProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute cross product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  cdef double wx = uy * vz - uz * vy
  cdef double wy = uz * vx - ux * vz
  cdef double wz = ux * vy - uy * vx
  return wx, wy, wz


cpdef double scalProd(double ux, double uy, double uz, double vx, double vy, double vz):
  """ Compute scalar product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
  """
  return ux*vx + uy*vy + uz*vz


cpdef (int) getCells(double x, double y, int ncols, int nrows, double csz):
  """ Locate point on grid.

  Parameters
  ----------
      x: float
          location in the x location of desiered interpolation
      y: float
          location in the y location of desiered interpolation
      ncols: int
          number of columns
      nrows: int
          number of rows
      csz: float
          cellsize of the grid

  Returns
  -------
      iCell: int
          index of the nearest lower left cell
  """
  cdef double Lx, Ly
  cdef int Lx0, Ly0
  cdef double xllc = 0.
  cdef double yllc = 0.
  cdef int iCell

  # find coordinates in normalized ref (origin (0,0) and cellsize 1)
  Lx = (x - xllc) / csz
  Ly = (y - yllc) / csz
  # find coordinates of the lower left cell
  Lx0 = <int>math.floor(Lx)
  Ly0 = <int>math.floor(Ly)
  iCell = Ly0*ncols + Lx0
  if (Lx0<0) | (Ly0<0) | (Lx0+1>=ncols) | (Ly0+1>=nrows):
    # check whether we are in the domain or not
    return -1

  return iCell


cpdef (int, int, double, double, double, double) getWeights(double x, double y, int iCell, double csz, int ncols, int interpOption):
  """ Get weight for interpolation from grid to single point location

  3 Options available : -0: nearest neighbour interpolation
                        -1: equal weights interpolation
                        -2: bilinear interpolation

  Parameters
  ----------
    x: float
        location in the x location of desiered interpolation
    y: float
        location in the y location of desiered interpolation
    iCell: int
        index of the nearest lower left cell
    csz: float
        cellsize of the grid
    ncols: int
      number of columns
    interpOption: int
        -0: nearest neighbour interpolation
        -1: equal weights interpolation
        -2: bilinear interpolation

  Returns
  -------
      Lx0, Ly0: ints
          lower left cell (col and row index)
      w00, w10, w01, w11: floats
          corresponding weights
  """
  cdef double w[4]
  cdef double Lx, Ly, dx, dy
  cdef int Lx0 = iCell % ncols
  cdef int Ly0 = iCell / ncols
  cdef double xllc = 0.
  cdef double yllc = 0.

  # find coordinates in normalized ref (origin (0,0) and cellsize 1)
  Lx = (x - xllc) / csz
  Ly = (y - yllc) / csz

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
  w[0] = (1-dx)*(1-dy)
  # lower right
  w[1] = dx*(1-dy)
  # uper left
  w[2] = (1-dx)*dy
  # and uper right
  w[3] = dx*dy

  return Lx0, Ly0, w[0], w[1], w[2], w[3]


cpdef (int, int, int, double, double, double, double) getCellAndWeights(double x, double y, int ncols, int nrows, double csz, int interpOption):
  cdef int Lx0, Ly0, iCell
  cdef double w[4]
  iCell = getCells(x, y, ncols, nrows, csz)
  Lx0, Ly0, w[0], w[1], w[2], w[3] = getWeights(x, y, iCell, csz, ncols, interpOption)
  return Lx0, Ly0, iCell, w[0], w[1], w[2], w[3]


cpdef (double, double, int, int, int, double, double, double, double) normalProjectionIteratrive(
  double xOld, double yOld, double zOld, double[:,:] ZDEM, double[:,:] nxArray, double[:,:] nyArray,
  double[:,:] nzArray, double csz, int ncols, int nrows, int interpOption,
  int reprojectionIterations, double threshold):
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
      reprojectionIterations: int
          maximum number or iterations
      threshold: double
          stop criterion for reprojection, stops when the distance between two iterations
          is smaller than threshold * csz

    Returns
    -------
    xNew: float
        x coordinate of the projected point
    yNew: float
        y coordinate of the projected point
    iCell: int
        index of the nearest lower left cell
    Lx0: int
        colomn of the nearest lower left cell
    Ly0: int
        row of the nearest lower left cell
    w: float[4]
        corresponding weights
    """
  cdef double zTemp, zn, xNew, yNew, zNew, dist, xPrev, yPrev, zPrev
  cdef double nx, ny, nz
  cdef int Lxi, Lyi, Lxj, Lyj
  cdef int Lx0, Ly0, iCell
  cdef double w[4]
  xNew = xOld
  yNew = yOld
  zNew = zOld
  iCell = getCells(xNew, yNew, ncols, nrows, csz)
  if iCell < 0:
    # if not on the DEM exit with iCell=-1
    return xNew, yNew, iCell, -1, -1, 0, 0, 0, 0
  Lx0, Ly0, w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, iCell, csz, ncols, interpOption)
  dist = csz
  # iterate
  while reprojectionIterations > 0 and dist > threshold*csz:
    reprojectionIterations = reprojectionIterations - 1
    xPrev = xNew
    yPrev = yNew
    zPrev = zNew
    # vertical projection of the point
    zTemp = getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], ZDEM)
    # normal vector at this vertical projection location
    # if we take the normal at the new particle position)
    nx, ny, nz = getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
    nx, ny, nz = normalize(nx, ny, nz)
    # What I think is better for normal projection
    # Normal component of the vector between the initial and projected point
    zn = (xNew-xOld) * nx + (yNew-yOld) * ny + (zTemp-zOld) * nz
    # correct position with the normal part
    xNew = xOld + zn * nx
    yNew = yOld + zn * ny
    zNew = zOld + zn * nz
    # get cell
    iCell = getCells(xNew, yNew, ncols, nrows, csz)
    if iCell < 0:
      # if not on the DEM exit with iCell=-1
      return xNew, yNew, iCell, -1, -1, 0, 0, 0, 0
    Lx0, Ly0, w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, iCell, csz, ncols, interpOption)

    dist = norm(xNew-xPrev, yNew-yPrev, zNew-zPrev)

  return xNew, yNew, iCell, Lx0, Ly0, w[0], w[1], w[2], w[3]


cpdef (double, double, int, int, int, double, double, double, double) samosProjectionIteratrive(
  double xOld, double yOld, double zOld, double[:,:] ZDEM, double[:,:] nxArray, double[:,:] nyArray,
  double[:,:] nzArray, double csz, int ncols, int nrows, int interpOption, int reprojectionIterations):
  """ Find the projection of a point on a mesh (comes from samos)

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
      reprojectionIterations: int
          maximum number or iterations

    Returns
    -------
    xNew: float
        x coordinate of the projected point
    yNew: float
        y coordinate of the projected point
    iCell: int
        index of the nearest lower left cell
    Lx0: int
        colomn of the nearest lower left cell
    Ly0: int
        row of the nearest lower left cell
    w: float[4]
        corresponding weights
    """
  cdef double zTemp, zn, xNew, yNew, zNew
  cdef double nx, ny, nz
  cdef int Lxi, Lyi, Lxj, Lyj
  cdef int Lx0, Ly0, iCell
  cdef double w[4]
  xNew = xOld
  yNew = yOld
  zNew = zOld
  iCell = getCells(xNew, yNew, ncols, nrows, csz)
  if iCell < 0:
    # if not on the DEM exit with iCell=-1
    return xNew, yNew, iCell, -1, -1, 0, 0, 0, 0
  Lx0, Ly0, w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, iCell, csz, ncols, interpOption)
  # get the cell location of the point
  Lxi = Lx0
  Lyi = Ly0
  # iterate
  while reprojectionIterations > 0:
    reprojectionIterations = reprojectionIterations - 1
    # vertical projection of the point
    zTemp = getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], ZDEM)
    # normal vector at this vertical projection location
    # if we take the normal at the new particle position)
    # nx, ny, nz = getVector(Lxy[0], Lxy[1], Lxy[2], Lxy[3], w[0], w[1], w[2], w[3], Nx, Ny, Nz)
    # What Peter does: get the normal of the cell center (but still not exactly giving the same results)
    nx, ny, nz = getVector(Lx0, Ly0, 0.25, 0.25, 0.25, 0.25, nxArray, nyArray, nzArray)
    nx, ny, nz = normalize(nx, ny, nz)
    zn = (xNew-xNew) * nx + (yNew-yNew) * ny + (zTemp-zNew) * nz
    # correct position with the normal part
    xNew = xNew + zn * nx
    yNew = yNew + zn * ny
    zNew = zNew + zn * nz
    # get cell
    iCell = getCells(xNew, yNew, ncols, nrows, csz)
    if iCell < 0:
      # if not on the DEM exit with iCell=-1
      return xNew, yNew, iCell, -1, -1, 0, 0, 0, 0
    Lx0, Ly0, w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, iCell, csz, ncols, interpOption)

    Lxj = Lx0
    Lyj = Ly0
    # are we in the same cell?
    if Lxi==Lxj and Lyi==Lyj:
      return xNew, yNew, iCell, Lx0, Ly0, w[0], w[1], w[2], w[3]
    Lxj = Lxi
    Lyj = Lyi
  return xNew, yNew, iCell, Lx0, Ly0, w[0], w[1], w[2], w[3]


cpdef (double, double, double, int, int, int, double, double, double, double) distConservProjectionIteratrive(
  double xPrev, double yPrev, double zPrev, double[:,:] ZDEM, double[:,:] nxArray, double[:,:] nyArray,
  double[:,:] nzArray, double xOld, double yOld, double zOld, double csz, int ncols, int nrows, int interpOption,
  int reprojectionIterations, double threshold):
  """ Find the projection of a point on a mesh conserving the distance
  with the previous time step position

  Iterative method to find the projection of a point on a surface defined by its mesh

  Parameters
  ----------
      xPrev: float
          x coordinate of the point at the previous time step
      yPrev: float
          y coordinate of the point at the previous time step
      zPrev: float
          z coordinate of the point at the previous time step
      ZDEM: 2D array
          z component of the DEM field at the grid nodes
      Nx: 2D array
          x component of the normal vector field at the grid nodes
      Ny: 2D array
          y component of the normal vector field at the grid nodes
      Nz: 2D array
          z component of the normal vector field at the grid nodes
      xOld: float
          x coordinate of the point to project
      yOld: float
          y coordinate of the point to project
      zOld: float
          z coordinate of the point to project
      csz: float
          cellsize of the grid
      interpOption: int
          -0: nearest neighbour interpolation
          -1: equal weights interpolation
          -2: bilinear interpolation
      reprojectionIterations: int
          maximum number or iterations
      threshold: double
          stop criterion for reprojection, stops when the error on the distance after projection
          is smaller than threshold * (dist + csz)

    Returns
    -------
    xNew: float
        x coordinate of the projected point
    yNew: float
        y coordinate of the projected point
    zNew: float
        z coordinate of the projected point
    iCell: int
        index of the nearest lower left cell
    Lx0: int
        colomn of the nearest lower left cell
    Ly0: int
        row of the nearest lower left cell
    w: float[4]
        corresponding weights
    """
  cdef double zTemp, zn, dist, distn, xNew, yNew, zNew
  cdef double nx, ny, nz
  cdef int Lx0, Ly0, iCell
  cdef double w[4]
  xNew = xOld
  yNew = yOld
  zNew = zOld
  # distance traveled by the point within the time step
  dist = norm(xNew-xPrev, yNew-yPrev, zNew-zPrev)
  # vertical projection of the point on the DEM
  iCell = getCells(xNew, yNew, ncols, nrows, csz)
  if iCell < 0:
    # if not on the DEM exit with iCell=-1
    return xNew, yNew, zNew, iCell, -1, -1, 0, 0, 0, 0
  Lx0, Ly0, w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, iCell, csz, ncols, interpOption)

  zTemp = getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], ZDEM)
  # measure distance between the projected point and the previous time step position
  distn = norm(xNew-xPrev, yNew-yPrev, zTemp-zPrev)
  # while iterate and error too big (aiming to conserve the distance dist during the reprojection)
  while reprojectionIterations > 0 and abs(distn-dist) > threshold*(dist + csz):
    reprojectionIterations = reprojectionIterations - 1
    # first step: orthogonal reprojection
    # get normal vector
    nx, ny, nz = getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
    nx, ny, nz = normalize(nx, ny, nz)
    zn = (xNew-xNew) * nx + (yNew-yNew) * ny + (zTemp-zNew) * nz
    # correct position with the normal part
    xNew = xNew + zn * nx
    yNew = yNew + zn * ny
    zNew = zNew + zn * nz

    # second step: conserve the distance dist (correction ov the position)
    # measure distance between this new point and the previous time step position
    distn = norm(xNew-xPrev, yNew-yPrev, zNew-zPrev)
    # adjust position on the Xprev, Xnew line
    xNew = xPrev + (xNew-xPrev) * dist / distn
    yNew = yPrev + (yNew-yPrev) * dist / distn
    zNew = zPrev + (zNew-zPrev) * dist / distn

    # third step
    # vertical projection of the point on the DEM
    iCell = getCells(xNew, yNew, ncols, nrows, csz)
    if iCell < 0:
      # if not on the DEM exit with iCell=-1
      return xNew, yNew, zNew, iCell, -1, -1, 0, 0, 0, 0
    Lx0, Ly0, w[0], w[1], w[2], w[3] = getWeights(xNew, yNew, iCell, csz, ncols, interpOption)
    zTemp = getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], ZDEM)
    # measure distance between this new point and the previous time step position
    distn = norm(xNew-xPrev, yNew-yPrev, zTemp-zPrev)

  return xNew, yNew, zTemp, iCell, Lx0, Ly0, w[0], w[1], w[2], w[3]


cpdef double[:] projOnRaster(double[:] xArray, double[:] yArray, double[:, :] vArray, double csz, int ncols,
                 int nrows, int interpOption):
  """ Interpolate from raster data to points

  Parameter
  ---------
    xArray: 1D float array
      x coordinates of the new points location
    yArray: 1D float array
      y coordinates of the new points location
    vArray: 2D float array
      raster values
    csz: float
      raster cell size
    ncols: int
      raster number of columns
    nrows: int
      raster number of rows
    interpOption: int
      interpolation option
      (0: nearest neighbour interpolation, 1: equal weights interpolation, 2: bilinear interpolation)

  Returns
  ---------
    v: 1D float array
      raster values interpolated at points (xArray, yArray) location
  """
  cdef int N = xArray.shape[0]
  cdef double x, y
  cdef int Lx0, Ly0, iCell
  cdef double w[4]
  cdef int k
  cdef double[:] v = np.zeros(N)
  for k in range(N):
    x = xArray[k]
    y = yArray[k]

    Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = getCellAndWeights(x, y, ncols, nrows, csz, interpOption)

    v[k] = getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], vArray)

  return v


cpdef double getScalar(int Lx0, int Ly0, double w0, double w1, double w2, double w3, double[:, :] V):
  """ Interpolate vector field from grid to single point location

  Originaly created to get the normal vector at location (x,y) given the
  normal vector field on the grid. Grid has its origin in (0,0).
  Can be used to interpolate any vector field.
  Interpolation using a bilinear interpolation

  Parameters
  ----------
    Lx0: int
        colomn of the nearest lower left cell
    Ly0: int
        row of the nearest lower left cell
    w: float[4]
        corresponding weights
    V: 2D numpy array
        scalar field at the grid nodes

  Returns
  -------
      v: float
          interpolated scalar at position (x, y)
  """
  cdef double v = (V[Ly0, Lx0]*w0 +
                   V[Ly0, Lx0+1]*w1 +
                   V[Ly0+1, Lx0]*w2 +
                   V[Ly0+1, Lx0+1]*w3)


  return v


cpdef (double, double, double) getVector(
  int Lx0, int Ly0, double w0, double w1, double w2, double w3,
  double[:, :] Nx, double[:, :] Ny, double[:, :] Nz):
  """ Interpolate vector field from grid to single point location

  Originaly created to get the normal vector at location (x,y) given the
  normal vector field on the grid. Grid has its origin in (0,0).
  Can be used to interpolate any vector field.

  Parameters
  ----------
      Lx0: int
          colomn of the nearest lower left cell
      Ly0: int
          row of the nearest lower left cell
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
  cdef double nx = (Nx[Ly0, Lx0]*w0 +
                   Nx[Ly0, Lx0+1]*w1 +
                   Nx[Ly0+1, Lx0]*w2 +
                   Nx[Ly0+1, Lx0+1]*w3)
  cdef double ny = (Ny[Ly0, Lx0]*w0 +
                   Ny[Ly0, Lx0+1]*w1 +
                   Ny[Ly0+1, Lx0]*w2 +
                   Ny[Ly0+1, Lx0+1]*w3)
  cdef double nz = (Nz[Ly0, Lx0]*w0 +
                   Nz[Ly0, Lx0+1]*w1 +
                   Nz[Ly0+1, Lx0]*w2 +
                   Nz[Ly0+1, Lx0+1]*w3)
  return nx, ny, nz


cpdef double SamosATfric(double rho, double tau0, double Rs0, double mu, double kappa, double B, double R,
                         double v, double p, double h):
  cdef double Rs = rho * v * v / (p + 0.001)
  cdef double div = h / R
  if div < 1.0:
    div = 1.0
  div = math.log(div) / kappa + B
  cdef double tau = tau0 + p * mu * (1.0 + Rs0 / (Rs0 + Rs)) + rho * v * v / (div * div)
  return tau


def computeIniMovement(cfg, particles, dem, dT, fields):
  """ add artifical viscosity effect on velocity

      Parameters
      ------------
      cfg: configparser
          configuration for DFA simulation
      particles : dict
          particles dictionary at t
      dem : dict
          dictionary with dem information
      dT : float
          time step
      fields: dict
        fields dictionary

      Returns
      --------
      particles: dict
        updated particle dictionary
      force: dict
        force dictionary

  """

  cdef int interpOption = cfg.getint('interpOption')
  cdef double rho = cfg.getfloat('rho')
  cdef double subgridMixingFactor = cfg.getfloat('subgridMixingFactorIni')
  cdef double gravAcc = cfg.getfloat('gravAcc')
  cdef double dt = dT
  cdef int nPart = particles['nPart']
  cdef double csz = dem['header']['cellsize']
  cdef int nrows = dem['header']['nrows']
  cdef int ncols = dem['header']['ncols']
  cdef double[:, :] ZDEM = dem['rasterData']
  cdef double[:, :] nxArray = dem['Nx']
  cdef double[:, :] nyArray = dem['Ny']
  cdef double[:, :] nzArray = dem['Nz']
  cdef double[:, :] areaRatser = dem['areaRaster']
  cdef double[:] mass = particles['m']
  cdef double[:] hArray = particles['h']
  cdef double[:] xArray = particles['x']
  cdef double[:] yArray = particles['y']
  cdef double[:] zArray = particles['z']
  cdef double[:] uxArray = particles['ux']
  cdef double[:] uyArray = particles['uy']
  cdef double[:] uzArray = particles['uz']
  cdef int[:] indXDEM = particles['indXDEM']
  cdef int[:] indYDEM = particles['indYDEM']
  cdef double[:, :] VX = fields['Vx']
  cdef double[:, :] VY = fields['Vy']
  cdef double[:, :] VZ = fields['Vz']
  cdef double[:] forceX = np.zeros(nPart, dtype=np.float64)
  cdef double[:] forceY = np.zeros(nPart, dtype=np.float64)
  cdef double[:] forceZ = np.zeros(nPart, dtype=np.float64)
  cdef double[:] forceFrict = np.zeros(nPart, dtype=np.float64)
  cdef double[:] dM = np.zeros(nPart, dtype=np.float64)
  cdef double[:] curvAcc = np.zeros(nPart, dtype=np.float64)
  cdef double[:] gEff = np.zeros(nPart, dtype=np.float64)

  cdef int indCellX, indCellY
  cdef double areaPart, uMag, m, h
  cdef double vMeanx, vMeany, vMeanz, vMeanNorm, dvX, dvY, dvZ
  cdef double x, y, z, xEnd, yEnd, zEnd, ux, uy, uz, uxDir, uyDir, uzDir
  cdef double nx, ny, nz, nxEnd, nyEnd, nzEnd, nxAvg, nyAvg, nzAvg
  # variables for interpolation
  cdef int Lx0, Ly0, LxEnd0, LyEnd0, iCell, iCellEnd
  cdef double w[4]
  cdef double wEnd[4]
  cdef int k
  force = {}

  for k in range(nPart):
      m = mass[k]
      x = xArray[k]
      y = yArray[k]
      z = zArray[k]
      h = hArray[k]
      ux = uxArray[k]
      uy = uyArray[k]
      uz = uzArray[k]
      indCellX = indXDEM[k]
      indCellY = indYDEM[k]

      # deduce area
      areaPart = m / (h * rho)
      # get cell and weights
      Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = getCellAndWeights(x, y, ncols, nrows, csz, interpOption)

      # get normal at the particle location
      nx, ny, nz = getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
      nx, ny, nz = normalize(nx, ny, nz)

      # add artificial viscosity
      ux, uy, uz = addArtificialViscosity(m, h, dt, rho, ux, uy, uz, subgridMixingFactor, Lx0, Ly0,
                                          w[0], w[1], w[2], w[3], VX, VY, VZ, nx, ny, nz)


      # update velocity
      uxArray[k] = ux
      uyArray[k] = uy
      uzArray[k] = uz
      gEff[k] = gravAcc * nz

  # save results
  force['dM'] = np.asarray(dM)
  force['forceX'] = np.asarray(forceX)
  force['forceY'] = np.asarray(forceY)
  force['forceZ'] = np.asarray(forceZ)
  force['forceFrict'] = np.asarray(forceFrict)
  particles['gEff'] = np.asarray(gEff)
  particles['curvAcc'] = np.asarray(curvAcc)
  particles['ux'] = np.asarray(uxArray)
  particles['uy'] = np.asarray(uyArray)
  particles['uz'] = np.asarray(uzArray)
  particles['m'] = np.asarray(mass)

  return particles, force
