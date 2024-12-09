#!python
# cython: boundscheck=False, wraparound=False, cdivision=True
""" manage Dams in DFA simulation
"""

import logging
import math
import pathlib
import numpy as np
import scipy as sp
import scipy.interpolate
import copy
import cython

# Local imports
import avaframe.in3Utils.geoTrans as gT
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.com1DFA.DFAtools as DFAtls
cimport avaframe.com1DFA.DFAToolsCython as DFAtlsC


# create local logger
log = logging.getLogger(__name__)


cpdef (int, double, double, double, double, double, double, double, double, double, double) getWallInteraction(
                                                                          double xOld, double yOld, double zOld,
                                                                          double xNew, double yNew, double zNew,
                                                                          double uxNew, double uyNew, double uzNew,
                                                                          int nDamPoints, double[:] xFootArray, double[:] yFootArray, double[:] zFootArray,
                                                                          double[:] xCrownArray, double[:] yCrownArray, double[:] zCrownArray,
                                                                          double[:] xTangentArray, double[:] yTangentArray, double[:] zTangentArray,
                                                                          int ncols, int nrows, double csz, int interpOption, double restitutionCoefficient,
                                                                          int nIterDam, double[:,:] nxArray, double[:,:] nyArray, double[:,:] nzArray,
                                                                          double[:,:] ZDEM, double[:,:] FT):
  """ Check if the particle trajectory intersects the dam lines and compute intersection coordinates

  the particle trajectory is given by the start and end points (in 3D)
  the dam line is given by x, y, z arrays of points (in 3D)
  the intersection is a ratio (r value between 0 and 1) of where the lines intersect
  xIntersection = (1.0-r)*xF1 + r*xF2
  yIntersection = (1.0-r)*yF1 + r*yF2

  Parameters
  ----------
  xOld: float
    x coordinate of the old particle position
  yOld: float
    y coordinate of the old particle position
  zOld: float
    z coordinate of the old particle position
  xNew: float
    x coordinate of the new particle position
  yNew: float
    y coordinate of the new particle position
  zNew: float
    z coordinate of the new particle position
  uxNew: float
    x component of the new particle velocity
  uyNew: float
    y component of the new particle velocity
  uzNew: float
    z component of the new particle velocity
  nDamPoints: int
    number of points in the dam line (length of the xFootArray... arrays)
  xFootArray: 1D array
    x coordinates of the dam foot line points
  yFootArray: 1D array
    y coordinates of the dam foot line points
  zFootArray: 1D array
    z coordinates of the dam foot line points
  xCrownArray: 1D array
    x coordinates of the dam crown line points
  yCrownArray: 1D array
    y coordinates of the dam crown line points
  zCrownArray: 1D array
    z coordinates of the dam crown line points
  xTangentArray: 1D array
    x comonent of the dam tangent vector (tangent to the foot line)
  yTangentArray: 1D array
    y comonent of the dam tangent vector (tangent to the foot line)
  zTangentArray: 1D array
    z comonent of the dam tangent vector (tangent to the foot line)
  ncols: int
    number of columns
  nrows: int
    number of rows
  csz: float
    cellsize of the raster
  interpOption: int
    -0: nearest neighbour interpolation
    -1: equal weights interpolation
    -2: bilinear interpolation
  restitutionCoefficient: float
    value between 0 and 1, 0, for a complete dissipation of the normal energy, 1 for a bounce with no dissipation
  nIterDam: int
    max number of interactuion points with the dam
  nxArray : 2D numpy array
    x component of the normal vector of the DEM
  nyArray : 2D numpy array
    y component of the normal vector of the DEM
  nzArray : 2D numpy array
    z component of the normal vector of the DEM
  ZDEM: 2D array
    z component of the DEM raster
  FT: 2D array
    flow thickness raster

  Returns
  -------
  foundIntersection: int
    1 if there is an interaction with the dam, 0 otherwise
  xNew: float
    x coordinate of the new particle position (after dam interaction)
  yNew: float
    y coordinate of the new particle position (after dam interaction)
  zNew: float
    z coordinate of the new particle position (after dam interaction)
  uxNew: float
    x component of the new particle velocity (after dam interaction)
  uyNew: float
    y component of the new particle velocity (after dam interaction)
  uzNew: float
    z component of the new particle velocity (after dam interaction)
  txWall: float
    x component of the tangent vector to the dam at the intersection point
  tyWall: float
    y component of the tangent vector to the dam at the intersection point
  tzWall: float
    z component of the tangent vector to the dam at the intersection point
  dEm: float
    scalar product between the gravity accleration and the dam face tangent vector
  """
  cdef int Lx0, Ly0, LxNew0, LyNew0, iCell, iCellNew, section, sectionNew
  cdef double w[4]
  cdef double wNew[4]
  cdef double nxNew, nyNew, nzNew, uMag, xNewTemp, yNewTemp
  cdef double xFoot, yFoot, zFoot
  cdef double xCrown, yCrown, zCrown
  cdef double nxWall, nyWall, nzWall
  cdef double txWall, tyWall, tzWall
  cdef double normalComponent
  cdef double dissEm = 0
  cdef double velMagMin = 0.000001
  cdef int iterate = nIterDam
  cdef int foundIntersection, foundIntersectionNew
  sectionNew = -1
  # wall interactions
  foundIntersection, section, xFoot, yFoot, zFoot, xCrown, yCrown, zCrown, txWall, tyWall, tzWall = getIntersection(xOld, yOld,
      xNew, yNew, xFootArray, yFootArray, zFootArray, xCrownArray, yCrownArray, zCrownArray, xTangentArray, yTangentArray, zTangentArray, nDamPoints)
  foundIntersectionNew = foundIntersection
  while foundIntersectionNew and iterate>0:
    iterate = iterate - 1
    # get cell and weights of intersection point
    Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = DFAtlsC.getCellAndWeights(xFoot, yFoot, ncols, nrows, csz, interpOption)
    # if(iCell < 0) continue; TODO: do we need to check for this?
    # get intersection foot point z coordinate
    zFoot = DFAtlsC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], ZDEM)
    # get flow thickness at foot point (measured along the surface normal)
    hFoot =  DFAtlsC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], FT)
    # compute vertical flow height from thickness (measured vertically)
    nx, ny, nz = DFAtlsC.getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
    # get average normal between old and new position
    nx, ny, nz = DFAtlsC.normalize(nx, ny, nz)
    hFootVertical = hFoot / nz  # ToDo: hFoot / (nz+0.01) in Peter's code, but nz can never be 0 right?
    # compute wall normal considering filling of the dam
    # update foot z coordinate. ToDo this is very artificial
    zFootFilled = zFoot + 0.5*hFootVertical
    # compute normal vector
    if zCrown-zFootFilled <= 0:
      log.debug('flow depth exceeds dam height')
    # If the snow fills the dam which means zCrown>zFootFilled, we compute the normal the same way
    nxWall, nyWall, nzWall = DFAtlsC.crossProd(xCrown-xFoot, yCrown-yFoot, zCrown-zFootFilled, txWall, tyWall, tzWall)
    # TODO: carefull, if zCrown-zFootFilled = 0 and the slope of the dam is 90° we have a vector of lenght 0...
    # normalizing is impossible
    nxWall, nyWall, nzWall = DFAtlsC.normalize(nxWall, nyWall, nzWall)

    # compute normal component of the trajectory from the foot to the xNew with no dam
    normalComponent = DFAtlsC.scalProd(nxWall, nyWall, nzWall, xNew-xFoot, yNew-yFoot, zNew-zFoot)
    # update position (reflection + dissipation)
    # Samos
    # xNew = xOld - (1.0 + restitutionCoefficient) * normalComponent * nxWall
    # yNew = yOld - (1.0 + restitutionCoefficient) * normalComponent * nyWall
    # zNew = zOld - (1.0 + restitutionCoefficient) * normalComponent * nzWall
    # same as samos but from new position
    # xNew = xNew - (1.0 + restitutionCoefficient) * normalComponent * nxWall
    # yNew = yNew - (1.0 + restitutionCoefficient) * normalComponent * nyWall
    # zNew = zNew - (1.0 + restitutionCoefficient) * normalComponent * nzWall
    # dissiparion in velocity direction and not in wall normal direction
    xNew = (1.0 - restitutionCoefficient) * xFoot + restitutionCoefficient * xNew - 2 * restitutionCoefficient * normalComponent * nxWall
    yNew = (1.0 - restitutionCoefficient) * yFoot + restitutionCoefficient * yNew - 2 * restitutionCoefficient * normalComponent * nyWall
    zNew = (1.0 - restitutionCoefficient) * zFoot + restitutionCoefficient * zNew - 2 * restitutionCoefficient * normalComponent * nzWall

    # update velocity (reflection + dissipation)
    normalComponent = DFAtlsC.scalProd(nxWall, nyWall, nzWall, uxNew, uyNew, uzNew)
    # samos
    # uxNew = uxNew - (1.0 + restitutionCoefficient) * normalComponent * nxWall
    # uyNew = uyNew - (1.0 + restitutionCoefficient) * normalComponent * nyWall
    # uzNew = uzNew - (1.0 + restitutionCoefficient) * normalComponent * nzWall
    # dissiparion in velocity direction and not in wall normal direction
    uxNew = (uxNew - 2 * normalComponent * nxWall) * restitutionCoefficient
    uyNew = (uyNew - 2 * normalComponent * nyWall) * restitutionCoefficient
    uzNew = (uzNew - 2 * normalComponent * nzWall) * restitutionCoefficient

    dissEm = -(zCrown-zFoot)

    if iterate>0:
      # reproject position
      # a simple vertical reprojection is done.
      LxNew0, LyNew0, iCellNew, wNew[0], wNew[1], wNew[2], wNew[3] = DFAtlsC.getCellAndWeights(xNew, yNew, ncols, nrows, csz, interpOption)
      if iCellNew >= 0:
        zNew = DFAtlsC.getScalar(LxNew0, LyNew0, wNew[0], wNew[1], wNew[2], wNew[3], ZDEM)

      # reproject velocity
      # get normal at the new particle location
      nxNew, nyNew, nzNew = DFAtlsC.getVector(LxNew0, LyNew0, wNew[0], wNew[1], wNew[2], wNew[3], nxArray, nyArray, nzArray)
      # normalize new normal
      nxNew, nyNew, nzNew = DFAtlsC.normalize(nxNew, nyNew, nzNew)
      uxNew, uyNew, uzNew, uMag = DFAtlsC.reprojectVelocity(uxNew, uyNew, uzNew, nxNew, nyNew, nzNew, velMagMin, 0)

      # We need to make sure we do not cross the dam again and bounce another time!!!
      # check if we have another intersection after reflexion
      # change the foot to make sure the intersection point is on the left part of the dam
      xFoot = xFoot + 0.0001 * nxWall
      yFoot = yFoot + 0.0001 * nyWall
      foundIntersectionNew, sectionNew, xFoot, yFoot, zFoot, xCrown, yCrown, zCrown, txWall, tyWall, tzWall = getIntersection(xFoot, yFoot,
          xNew, yNew, xFootArray, yFootArray, zFootArray, xCrownArray, yCrownArray, zCrownArray, xTangentArray, yTangentArray, zTangentArray, nDamPoints)

      if foundIntersectionNew and section==sectionNew:
        # crossing the dam on the same section is allowed
        foundIntersectionNew = 0
  return foundIntersection, xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dissEm


cpdef (int, int, double, double, double, double, double, double, double, double, double) getIntersection(double xOld, double yOld,
                                                                                            double xNew, double yNew,
                                                                                            double[:] xFoot,
                                                                                            double[:] yFoot,
                                                                                            double[:] zFoot,
                                                                                            double[:] xCrown,
                                                                                            double[:] yCrown,
                                                                                            double[:] zCrown,
                                                                                            double[:] xTangent,
                                                                                            double[:] yTangent,
                                                                                            double[:] zTangent,
                                                                                            int nDamPoints):
  """ Check if the particle trajectory intersects the dam lines and compute intersection coefficient r

  the particle trajectory is given by the start and end points (in 3D)
  the dam line is given by x, y, z arrays of points (in 3D)
  the intersection is a ratio (r value between 0 and 1) of where the lines intersect
  xIntersection = (1.0-r)*xF1 + r*xF2
  yIntersection = (1.0-r)*yF1 + r*yF2

  Parameters
  ----------
  xOld: float
    x coordinate of the old particle position
  yOld: float
    y coordinate of the old particle position
  xNew: float
    x coordinate of the new particle position
  yNew: float
    y coordinate of the new particle position
  xFoot: 1D array
    x coordinates of the dam foot line points
  yFoot: 1D array
    y coordinates of the dam foot line points
  zFoot: 1D array
    z coordinates of the dam foot line points
  xCrown: 1D array
    x coordinates of the dam crown line points
  yCrown: 1D array
    y coordinates of the dam crown line points
  zCrown: 1D array
    z coordinates of the dam crown line points
  xTangent: 1D array
    x comonent of the dam tangent vector (tangent to the foot line)
  yTangent: 1D array
    y comonent of the dam tangent vector (tangent to the foot line)
  zTangent: 1D array
    z comonent of the dam tangent vector (tangent to the foot line)
  nDamPoints: int
    number of points in the dam line (length of the xFoot... arrays)

  Returns
  -------
  intersection: int
    1 if the lines intersect, 0 otherwise
  section: int
      interaction section of the dam
  xF: float
    x coordinate of the foot intersection point
  yF: float
    y coordinate of the foot intersection point
  zF: float
    z coordinate of the foot intersection point
  xC: float
    x coordinate of the crown point corresponding the intersection point
  yC: float
    y coordinate of the crown point corresponding the intersection point
  zC: float
    z coordinate of the crown point corresponding the intersection point
  xT: float
    x component of the tangent vector to the dam at the intersection point
  yT: float
    y component of the tangent vector to the dam at the intersection point
  zT: float
    z component of the tangent vector to the dam at the intersection point
  """
  cdef int i, intersection
  cdef double xF1, yF1, zF1, xF2, yF2, zF2
  cdef double xF, yF, zF
  cdef double xC, yC, zC, xC1, yC1, zC1, xC2, yC2, zC2
  cdef double xT, yT, zT

  for i in range(nDamPoints-1):
    # get end points of the considered wall section
    xF1 = xFoot[i]
    yF1 = yFoot[i]
    zF1 = zFoot[i]
    xF2 = xFoot[i+1]
    yF2 = yFoot[i+1]
    zF2 = zFoot[i+1]
    # does the particle trajectory intersect with the crown line of the wall
    intersection, r = linesIntersect(xOld, yOld, xNew, yNew, xF1, yF1, xF2, yF2)
    # if yes compute coordinates and tangent at intersection
    if intersection:
      # get crown points of wall segment
      xC1 = xCrown[i]
      xC2 = xCrown[i+1]
      yC1 = yCrown[i]
      yC2 = yCrown[i+1]
      zC1 = zCrown[i]
      zC2 = zCrown[i+1]
      # get tangent vectors of wall segment
      # which is the same as the tangent vector at the intersection
      xT = xTangent[i]
      yT = yTangent[i]
      zT = zTangent[i]
      # compute intersection
      xF = (1.0-r)*xF1 + r*xF2
      yF = (1.0-r)*yF1 + r*yF2
      zF = (1.0-r)*zF1 + r*zF2
      # get crown at intersecion
      xC = (1.0-r)*xC1 + r*xC2
      yC = (1.0-r)*yC1 + r*yC2
      zC = (1.0-r)*zC1 + r*zC2
      return intersection, i, xF, yF, zF, xC, yC, zC, xT, yT, zT
  return intersection, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0




cpdef (int, double) linesIntersect(double xOld, double yOld, double xNew, double yNew,
                                   double xF1, double yF1, double xF2, double yF2):
  """ Check if two lines intersect and compute intersection

  the lines are given by the start and end points (in 2D)
  the intersection is a ratio (r value between 0 and 1) of where the lines intersect
  xIntersection = (1.0-r)*xF1 + r*xF2
  yIntersection = (1.0-r)*yF1 + r*yF2

  Parameters
  ----------
  xOld: float
    x coordinate of the start point of the first line
  yOld: float
    y coordinate of the start point of the first line
  xNew: float
    x coordinate of the end point of the first line
  yNew: float
    y coordinate of the end point of the first line
  xF1: float
    x coordinate of the start point of the second line
  yF1: float
    y coordinate of the start point of the second line
  xF2: float
    x coordinate of the end point of the second line
  yF2: float
    y coordinate of the end point of the second line
  Returns
  -------
  intersection: int
    1 if the lines intersect, 0 otherwise
  r: float
    intersection ration between 0 and 1 of where the lines intersect (on the second line)
  """
  cdef double ax, ay, bx, bY, cx, cy
  cdef double det, u, v, r
  ax = xF2-xF1
  ay = yF2-yF1
  bx = xOld-xNew
  bY = yOld-yNew
  cx = xOld-xF1
  cy = yOld-yF1
  det = ax*bY - ay*bx
  u = cx*bY - cy*bx
  v = ax*cy - ay*cx
  if(det == 0.0):
    return 0, 0
  if(det < 0.0):
    if(u > 0.0):
      return 0, 0
    if(v > 0.0):
      return 0, 0
    if(u < det):
      return 0, 0
    if(v < det):
      return 0, 0
  else:
    if(u < 0.0):
      return 0, 0
    if(v < 0.0):
      return 0, 0
    if(u > det):
      return 0, 0
    if(v > det):
      return 0, 0
  r = u / det
  return 1, r


def initializeWallLines(cfg, dem, wallLineDict, savePath=''):
  """Initialize dam dictionary

  Parameters:
  -----------
  cfg: configparser
    configuration with slope, height and restitutionCoefficient of the dam
  dem: dict
    dem dictionary
  wallLineDict: dict
    dam dictionary with (x,y) coordinates of the centerline
  savePath: pathlib path
    save dam foot line to this path (if savePath is not '')
  Returns
  -------
  wallLineDict: dict
    dam dictionary updated with z coordinate, foot line, crown, tangent....
  """
  cdef double w[4]
  cdef int Lx0, Ly0
  if wallLineDict is not None:
    log.debug('Initializing dam line from: %s' % str(wallLineDict['fileName'][0]))
    log.debug('Dam line feature: %s' % str(wallLineDict['Name'][0]))
    # get z coordinate of the dam polyline
    wallLineDict['x'] = wallLineDict['x'] - dem['originalHeader']['xllcenter']
    wallLineDict['y'] = wallLineDict['y'] - dem['originalHeader']['yllcenter']
    # the z coordinate corresponds to the crown
    wallLineDict['zCrown'] = copy.deepcopy(wallLineDict['z'])
    # get the z of the centerline by projection on the topography
    wallLineDict, _ = gT.projectOnRaster(dem, wallLineDict, interp='bilinear')
    #ToDo: maybe we need to resample!
    # wallLineDict, _ = gT.prepareLine(dem, wallLineDict, distance=dem['header']['cellsize'], Point=None)
    nDamPoints = np.size(wallLineDict['x'])
    wallLineDict['nPoints'] = nDamPoints
    wallLineDict['restitutionCoefficient'] = cfg.getfloat('restitutionCoefficient')
    wallLineDict['nIterDam'] = cfg.getfloat('nIterDam')
    try:
      log.debug('Using dam slope from shape file: %s °' % wallLineDict['slope'])
      wallLineDict['slope'] = np.ones(nDamPoints) * np.radians(wallLineDict['slope'])
    except TypeError:
      message = 'Provide a valid slope value for the dam (\'slope\' attribute in the dam line shape file)'
      log.error(message)
      raise TypeError(message)

    # compute wall tangent vector
    tangentsX = np.zeros(nDamPoints-1)
    tangentsY = np.zeros(nDamPoints-1)
    tangentsZ = np.zeros(nDamPoints-1)
    tangentsTempX = np.zeros(nDamPoints)
    tangentsTempY = np.zeros(nDamPoints)
    tangentsTempZ = np.zeros(nDamPoints)
    # tangent between i and i+1
    for i in range(nDamPoints-1):
      tx = wallLineDict['x'][i+1] - wallLineDict['x'][i]
      ty = wallLineDict['y'][i+1] - wallLineDict['y'][i]
      tz = wallLineDict['z'][i+1] - wallLineDict['z'][i]
      tx, ty, tz = DFAtlsC.normalize(tx, ty, tz)
      # add it to i
      tangentsX[i] = tx
      tangentsY[i] = ty
      tangentsZ[i] = tz
      # add it to i and i+1
      tangentsTempX[i] = tangentsTempX[i] + tx
      tangentsTempY[i] = tangentsTempY[i] + ty
      tangentsTempZ[i] = tangentsTempZ[i] + tz
      tangentsTempX[i+1] = tangentsTempX[i+1] + tx
      tangentsTempY[i+1] = tangentsTempY[i+1] + ty
      tangentsTempZ[i+1] = tangentsTempZ[i+1] + tz
    # normalize
    tangentsTempX, tangentsTempY, tangentsTempZ = DFAtls.normalize(tangentsTempX, tangentsTempY, tangentsTempZ)
    # add it to the wallLineDict
    wallLineDict['xTangent'] = tangentsX
    wallLineDict['yTangent'] = tangentsY
    wallLineDict['zTangent'] = tangentsZ
    # get the foot line (on the left side of the dam) Find the intersection between the dam and the botom surface
    footLineX = np.zeros(nDamPoints)
    footLineY = np.zeros(nDamPoints)
    footLineZ = np.zeros(nDamPoints)
    crownX = np.zeros(nDamPoints)
    crownY = np.zeros(nDamPoints)
    height = np.zeros(nDamPoints)

    for i in range(nDamPoints):
      x = wallLineDict['x'][i]
      y = wallLineDict['y'][i]
      z = wallLineDict['z'][i]
      h = wallLineDict['zCrown'][i] - wallLineDict['z'][i]
      height[i] = h
      slope = wallLineDict['slope'][i]
      # get cell and weights of point
      Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = DFAtlsC.getCellAndWeights(x, y, dem['header']['ncols'],
                                                                          dem['header']['nrows'], dem['header']['cellsize'],
                                                                          2)
      # get the normal vector to the dem surface at the points location
      surfaceNormalX = DFAtlsC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], dem['Nx'])
      surfaceNormalY = DFAtlsC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], dem['Ny'])
      surfaceNormalZ = DFAtlsC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], dem['Nz'])
      # compute crown points
      crownX[i] = x
      crownY[i] = y
      # compute the normal to the dam in 2D ("top view")
      # (0, 0, 1) and (tangentsX[i], tangentsY[i], tangentsZ[i]) but they are not ortogonal, d is not of norm 1
      dx, dy, dz = DFAtls.crossProd(0, 0, 1, tangentsTempX[i], tangentsTempY[i], tangentsTempZ[i])
      d = DFAtlsC.norm(dx, dy, dz)
      # add the z component to get the tangent vector to the sloped wall
      dz = - np.tan(slope) * d
      # get the intersection between the dam side slope and the bottom surface
      r = -h*surfaceNormalZ / DFAtls.scalProd(dx, dy, dz, surfaceNormalX, surfaceNormalY, surfaceNormalZ)
      # compute foot points
      footLineX[i] = x + r * dx
      footLineY[i] = y + r * dy
      # get cell and weights of foot line
      Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = DFAtlsC.getCellAndWeights(footLineX[i], footLineY[i], dem['header']['ncols'],
                                                                          dem['header']['nrows'], dem['header']['cellsize'],
                                                                          2)
      # Samos computes the z as if it was an inclined plane
      # footLineZ[i] = z + h + r * dz
      # Reproject to get the z coord of the footLine
      footLineZ[i] = DFAtlsC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], dem['rasterData'])

    # save foot and crown in the dict
    wallLineDict['x'] = footLineX
    wallLineDict['y'] = footLineY
    wallLineDict['z'] = footLineZ
    wallLineDict['xCrown'] = crownX
    wallLineDict['yCrown'] = crownY
    wallLineDict['height'] = height

    # locate cells around the foot line (then we will only activate the dam effect for particles in the surroudings
    # of the dam)
    wallLineDict = gT.getCellsAlongLine(dem['header'], wallLineDict, addBuffer=True)
    wallLineDict['dam'] = 1
    # FSO: turned off due to being unsafe for parallel computation
    # if savePath != '':
    #   fileName = shpConv.writeLine2SHPfile(wallLineDict, 'dam foot line', savePath, header=dem['originalHeader'])
  else:
    # crete a dummy dict (needed so that cython runs)
    wallLineDict = {'dam': 0, 'cellsCrossed': np.zeros((dem['header']['ncols']*dem['header']['nrows'])).astype(
        np.int64)}
    for key in ['x', 'y', 'z', 'xCrown', 'yCrown', 'zCrown', 'xTangent', 'yTangent', 'zTangent']:
      wallLineDict[key] = np.ones((1))*1.0
    for key in ['nPoints', 'height', 'slope', 'restitutionCoefficient', 'nIterDam']:
      wallLineDict[key] = 0

  return wallLineDict
