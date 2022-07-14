""" manage Dams in DFA simulation
"""

import logging
import math
import pathlib
import numpy as np
import scipy as sp
import scipy.interpolate
import copy

# Local imports
import avaframe.in3Utils.geoTrans as gT
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.com1DFA.DFAfunctionsCython as DFAfunC


# create local logger
log = logging.getLogger(__name__)


cpdef (double, double, double, double, double, double, double, double, double, double) getWallInteraction(
                                                                          double xOld, double yOld, double zOld,
                                                                          double xNew, double yNew, double zNew,
                                                                          double uxNew, double uyNew, double uzNew,
                                                                          int nDamPoints, double[:] xFootArray, double[:] yFootArray, double[:] zFootArray,
                                                                          double[:] xCrownArray, double[:] yCrownArray, double[:] zCrownArray,
                                                                          double[:] xTangentArray, double[:] yTangentArray, double[:] zTangentArray,
                                                                          int ncols, int nrows, double csz, int interpOption, double restitutionCoefficient,
                                                                          double[:,:] nxArray, double[:,:] nyArray, double[:,:] nzArray,
                                                                          double[:,:] ZDEM, double[:,:] FT):
  """ Check if the particle trajectory intersects the dam lines and compute intersection

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
  cdef int Lx0, Ly0, LxNew0, LyNew0, iCell, iCellNew
  cdef double w[4]
  cdef double wNew[4]
  cdef double xFoot, yFoot, zFoot
  cdef double xCrown, yCrown, zCrown
  cdef double nxWall, nyWall, nzWall
  cdef double txWall, tyWall, tzWall
  cdef double normalComponent, dEm
  # wall interactions
  foundIntersection, xFoot, yFoot, zFoot, xCrown, yCrown, zCrown, txWall, tyWall, tzWall = getIntersection(xOld, yOld,
      xNew, yNew, xFootArray, yFootArray, zFootArray, xCrownArray, yCrownArray, zCrownArray, xTangentArray, yTangentArray, zTangentArray, nDamPoints)
  if foundIntersection:
    # get cell and weights of intersection point
    Lx0, Ly0, iCell, w[0], w[1], w[2], w[3] = DFAfunC.getCellAndWeights(xFoot, yFoot, ncols, nrows, csz, interpOption)
    # if(iCell < 0) continue; TODO: do we need to check for this?
    # get intersection foot point z coordinate
    zFoot = DFAfunC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], ZDEM)
    # get flow thickness at foot point (measured along the surface normal)
    hFoot =  DFAfunC.getScalar(Lx0, Ly0, w[0], w[1], w[2], w[3], FT)
    # compute vertical flow height from thickness (measured verticaly)
    nx, ny, nz = DFAfunC.getVector(Lx0, Ly0, w[0], w[1], w[2], w[3], nxArray, nyArray, nzArray)
    # get average normal between old and new position
    nx, ny, nz = DFAfunC.normalize(nx, ny, nz)
    hFootVertical = hFoot / nz  # hFoot / (nz+0.01) in Peter's code, but nz can never be 0 right?
    # compute wall normal considering filling of the dam
    # update foot z coordinate
    zFootFilled = zFoot + 0.5*hFootVertical
    # compute normal vector
    # TODO: what happens if snow fills the dam? which means zCrown>zFootFilled
    nxWall, nyWall, nzWall = DFAfunC.crossProd(xCrown-xFoot, yCrown-yFoot, zCrown-zFootFilled, txWall, tyWall, tzWall)
    # TODO: carefull, if zCrown-zFootFilled = and the slope of the dam is 90° we have a vector of lenght 0...
    # normalizing is impossible
    # TODO: if the angle is 90°, snow cant go through even if zFootFilled>zCrown...
    # TODO: in general, I think there is a broblem if zFootFilled>zCrown...
    nxWall, nyWall, nzWall = DFAfunC.normalize(nxWall, nyWall, nzWall)

    # if there is an interaction with the dam
    normalComponent = DFAfunC.scalProd(nxWall, nyWall, nzWall, xNew-xFoot, yNew-yFoot, zNew-zFoot)
    # update position (reflection + dissipation)
    # ToDo: why take xold????!!! I would use xNew
    # xNew = xOld - (1.0 + restitutionCoefficient) * normalComponent * nxWall
    # yNew = yOld - (1.0 + restitutionCoefficient) * normalComponent * nyWall
    # zNew = zOld - (1.0 + restitutionCoefficient) * normalComponent * nzWall
    xNew = xNew - (1.0 + restitutionCoefficient) * normalComponent * nxWall
    yNew = yNew - (1.0 + restitutionCoefficient) * normalComponent * nyWall
    zNew = zNew - (1.0 + restitutionCoefficient) * normalComponent * nzWall
    # update velocity (reflection + dissipation)
    normalComponent = DFAfunC.scalProd(nxWall, nyWall, nzWall, uxNew, uyNew, uzNew)
    uxNew = uxNew - (1.0 + restitutionCoefficient) * normalComponent * nxWall
    uyNew = uyNew - (1.0 + restitutionCoefficient) * normalComponent * nyWall
    uzNew = uzNew - (1.0 + restitutionCoefficient) * normalComponent * nzWall
    # ToDo: we need to make sure we do not cross the dam again and bounce another time!!!
    dEm = DFAfunC.scalProd(0, 0, -1, xCrown-xFoot, yCrown-yFoot, zCrown-zFoot)

  return xNew, yNew, zNew, uxNew, uyNew, uzNew, txWall, tyWall, tzWall, dEm


cpdef (int, double, double, double, double, double, double, double, double, double) getIntersection(double xOld, double yOld,
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
  """ Check if the particle trajectory intersects the dam lines and compute intersection

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
  cdef int i
  cdef double xF1, yF1, xF2, yF2
  cdef double xF, yF
  cdef double xC, yC, zC, xC1, yC1, zC1, xC2, yC2, zC2
  cdef double xT, yT, zT, xT1, yT1, zT1, xT2, yT2, zT2

  for i in range(nDamPoints-1):
    # get end points of the considered wall section
    xF1 = xFoot[i]
    yF1 = yFoot[i]
    zF1 = zFoot[i]
    xF2 = xFoot[i+1]
    yF2 = yFoot[i+1]
    zF2 = zFoot[i+1]
    # does the particle trajectory intersect with the foot line of the wall
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
      xT1 = xTangent[i]
      xT2 = xTangent[i+1]
      yT1 = yTangent[i]
      yT2 = yTangent[i+1]
      zT1 = zTangent[i]
      zT2 = zTangent[i+1]
      # compute intersection
      xF = (1.0-r)*xF1 + r*xF2
      yF = (1.0-r)*yF1 + r*yF2
      zF = (1.0-r)*zF1 + r*zF2
      # get crown at intersecion
      xC = (1.0-r)*xC1 + r*xC2
      yC = (1.0-r)*yC1 + r*yC2
      zC = (1.0-r)*zC1 + r*zC2
      # get tangent vector at intersection
      xT = (1.0)*xT1 + r*xT2
      yT = (1.0)*yT1 + r*yT2
      zT = (1.0)*zT1 + r*zT2
      # normalize tangent vector
      xT, yT, zT = DFAfunC.normalize(xT, yT, zT)
      return intersection, xF, yF, zF, xC, yC, zC, xT, yT, zT
  return intersection, 0, 0, 0, 0, 0, 0, 0, 0, 0




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
  if wallLineDict is not None:
    # get z coordinate of the dam polyline
    wallLineDict['x'] = wallLineDict['x'] - dem['originalHeader']['xllcenter']
    wallLineDict['y'] = wallLineDict['y'] - dem['originalHeader']['yllcenter']
    wallLineDict, _ = gT.projectOnRaster(dem, wallLineDict, interp='bilinear')
    #ToDo: maybe we need to ressample!
    # wallLineDict, _ = gT.prepareLine(dem, wallLineDict, distance=dem['header']['cellsize'], Point=None)
    nDamPoints = np.size(wallLineDict['x'])
    wallLineDict['nPoints'] = nDamPoints
    wallLineDict['restitutionCoefficient'] = cfg.getfloat('restitutionCoefficient')
    wallLineDict['height'] = np.ones(nDamPoints) * cfg.getfloat('damHeight')
    wallLineDict['slope'] = np.ones(nDamPoints) * np.radians(cfg.getfloat('damSlope'))
    # compute wall tangent vector
    tangentsX = np.zeros(nDamPoints)
    tangentsY = np.zeros(nDamPoints)
    tangentsZ = np.zeros(nDamPoints)
    # tangent between i and i+1
    for i in range(nDamPoints-1):
      tx = wallLineDict['x'][i+1] - wallLineDict['x'][i]
      ty = wallLineDict['y'][i+1] - wallLineDict['y'][i]
      tz = wallLineDict['z'][i+1] - wallLineDict['z'][i]
      tx, ty, tz = DFAfunC.normalize(tx, ty, tz)
      # add it to i and i+1
      tangentsX[i] = tangentsX[i] + tx
      tangentsY[i] = tangentsY[i] + ty
      tangentsZ[i] = tangentsZ[i] + tz
      tangentsX[i+1] = tangentsX[i+1] + tx
      tangentsY[i+1] = tangentsY[i+1] + ty
      tangentsZ[i+1] = tangentsZ[i+1] + tz
    # normalize
    tangentsX, tangentsY, tangentsZ = DFAtls.normalize(tangentsX, tangentsY, tangentsZ)
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
    crownZ = np.zeros(nDamPoints)
    # get the normal vector to the dem surface at the polyline points location
    surfaceNormalX, _ = gT.projectOnRaster(dem, wallLineDict, interp='bilinear', what='Nx', where='nx')
    surfaceNormalX = surfaceNormalX['nx']
    surfaceNormalY, _ = gT.projectOnRaster(dem, wallLineDict, interp='bilinear', what='Ny', where='ny')
    surfaceNormalY = surfaceNormalY['ny']
    surfaceNormalZ, _ = gT.projectOnRaster(dem, wallLineDict, interp='bilinear', what='Nz', where='nz')
    surfaceNormalZ = surfaceNormalZ['nz']
    # compute crown points
    for i in range(nDamPoints):
      x = wallLineDict['x'][i]
      y = wallLineDict['y'][i]
      z = wallLineDict['z'][i]
      h = wallLineDict['height'][i]
      crownX[i] = x
      crownY[i] = y
      crownZ[i] = z + h
      # compute the normal to the dam in 2D ("top view")
      # (0, 0, 1) and (tangentsX[i], tangentsY[i], tangentsZ[i]) but they are not ortogonal, d is not of norm 1
      dx, dy, dz = DFAtls.crossProd(0, 0, 1, tangentsX[i], tangentsY[i], tangentsZ[i])
      d = DFAfunC.norm(dx, dy, dz)
      # add the z component to get the tangent vector to the sloped wall
      dz = - np.tan(wallLineDict['slope'][i]) * d
      dx, dy, dz = DFAfunC.normalize(dx, dy, dz)
      # get the intersection between the dam side slope and the bottom surface
      r = -h*surfaceNormalZ[i] / DFAtls.scalProd(dx, dy, dz, surfaceNormalX[i], surfaceNormalY[i], surfaceNormalZ[i])
      # compute foot points
      footLineX[i] = x + r * dx
      footLineY[i] = y + r * dy
      # ToDo: should we reproject?
      footLineZ[i] = z + h + r * dz

    # save foot and crown in the dict
    wallLineDict['x'] = footLineX
    wallLineDict['y'] = footLineY
    wallLineDict['z'] = footLineZ
    wallLineDict['xCrown'] = crownX
    wallLineDict['yCrown'] = crownY
    wallLineDict['zCrown'] = crownZ

    # locate cells around the foot line (then we will only activate the dam effect for particles in the surroudings
    # of the dam)
    wallLineDict = gT.getCellsAlongLine(dem['header'], wallLineDict, addBuffer=True)
    wallLineDict['flagDam'] = 1
    if savePath != '':
      fileName = shpConv.writeLine2SHPfile(wallLineDict, 'dam foot line', savePath, header=dem['originalHeader'])
  else:
    # crete a dummy dict (needed so that cython runs)
    wallLineDict = {'flagDam': 0, 'cellsCrossed': np.zeros((dem['header']['ncols']*dem['header']['nrows'])).astype(int)}
    for key in ['x', 'y', 'z', 'xCrown', 'yCrown', 'zCrown', 'xTangent', 'yTangent', 'zTangent']:
      wallLineDict[key] = np.ones((1))*1.0
    for key in ['nPoints', 'height', 'slope', 'restitutionCoefficient']:
      wallLineDict[key] = 0

  return wallLineDict
