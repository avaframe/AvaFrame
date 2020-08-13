"""
Opperations and transformations of rasters and lines
"""

import numpy as np
import logging


# create local logger
log = logging.getLogger(__name__)


def projectOnRaster(dem, Points):
    """ Projects the points Points on Raster using a bilinear interpolation
    and returns the z coord
    Input :
    Points: list of points (x,y) 2 rows as many columns as Points
    Output:
    PointsZ: list of points (x,y,z) 3 rows as many columns as Points

    TODO: test
    """
    header = dem['header']
    rasterdata = dem['rasterdata']
    xllcorner = header.xllcorner
    yllcorner = header.yllcorner
    cellsize = header.cellsize
    xcoor = Points['x']
    ycoor = Points['y']
    zcoor = np.array([])
    for i in range(np.shape(xcoor)[0]):
        Lx = (xcoor[i] - xllcorner) / cellsize
        Ly = (ycoor[i] - yllcorner) / cellsize
        Lx0 = int(np.floor(Lx))
        Ly0 = int(np.floor(Ly))
        Lx1 = int(np.floor(Lx)) + 1
        Ly1 = int(np.floor(Ly)) + 1
        dx = Lx - Lx0
        dy = Ly - Ly0
        f11 = rasterdata[Ly0][Lx0]
        f12 = rasterdata[Ly1][Lx0]
        f21 = rasterdata[Ly0][Lx1]
        f22 = rasterdata[Ly1][Lx1]
        # using bilinear interpolation on the cell
        value = f11*(1-dx)*(1-dy) + f21*dx*(1-dy) + f12*(1-dx)*dy + f22*dx*dy
        zcoor = np.append(zcoor, value)

    Points['z'] = zcoor
    return Points


def prepareLine(dem, avapath, splitPoint, distance):
    """ 1- Resample the avapath line with a max intervall of distance=10m
    between points (projected distance on the horizontal plane).
    2- Make avalanch profile out of the path (affect a z value using the dem)
    3- Get projected split point on the profil (closest point)

    TODO: test
    """
    xcoor = avapath['x']
    ycoor = avapath['y']
    xcoornew = np.array([xcoor[0]])
    ycoornew = np.array([ycoor[0]])
    s = np.array([0])  # curvilinear coordinate allong the path
    # loop on the points of the avapath
    for i in range(np.shape(xcoor)[0] - 1):
        Vx = xcoor[i + 1] - xcoor[i]
        Vy = ycoor[i + 1] - ycoor[i]
        D = np.sqrt(Vx**2 + Vy**2)
        nd = int(np.floor(D / distance) + 1)
        # Resample each segment
        S0 = s[-1]
        for j in range(1, nd + 1):
            xn = j / (nd) * Vx + xcoor[i]
            yn = j / (nd) * Vy + ycoor[i]
            xcoornew = np.append(xcoornew, xn)
            ycoornew = np.append(ycoornew, yn)
            s = np.append(s, S0 + D * j / nd)

    # test = np.transpose(np.array([[header.xllcorner,header.yllcorner],
    # [header.xllcorner+header.cellsize*(header.ncols-1),header.yllcorner],
    # [header.xllcorner,header.yllcorner+header.cellsize*(header.nrows-1)],
    # [header.xllcorner+header.cellsize*(header.ncols-1),header.yllcorner+
    # header.cellsize*(header.nrows-1)]]))
    # Test = ProjectOnRaster(header,rasterdata,test)
    ResampAvaPath = avapath
    ResampAvaPath['x'] = xcoornew
    ResampAvaPath['y'] = ycoornew
    ResampAvaPath = projectOnRaster(dem, ResampAvaPath)
    ResampAvaPath['s'] = s
    AvaProfile = ResampAvaPath
    # find split point by computing the distance to the line
    if splitPoint:
        projSplitPoint, splitPoint = findSplitPoint(AvaProfile, splitPoint)
    else:
        projSplitPoint = None

    return AvaProfile, projSplitPoint, splitPoint


def findSplitPoint(AvaProfile, splitPoint):
    """ find split point by computing the distance between splitpoints and the line
        selects the closest splitpoint (if several were given in input)
    """
    xcoor = AvaProfile['x']
    ycoor = AvaProfile['y']
    Dist = np.empty((0))
    IndSplit = np.empty((0))
    for i in range(len(splitPoint['x'])):
        dist = np.sqrt((xcoor - splitPoint['x'][i])**2 +
                       (ycoor - splitPoint['y'][i])**2)
        indSplit = np.argmin(dist)
        IndSplit = np.append(IndSplit, indSplit)
        Dist = np.append(Dist, dist[indSplit])

    ind = np.argmin(Dist)
    indSplit = int(IndSplit[ind])
    projSplitPoint = {}
    projSplitPoint['x'] = AvaProfile['x'][indSplit]
    projSplitPoint['y'] = AvaProfile['y'][indSplit]
    projSplitPoint['z'] = AvaProfile['z'][indSplit]
    projSplitPoint['s'] = AvaProfile['s'][indSplit]
    projSplitPoint['indSplit'] = indSplit
    return projSplitPoint, splitPoint

def checkProfile(AvaProfile, projSplitPoint):
    """ check that the avalanche profiles goes from top to bottom """
    if projSplitPoint:
        indSplit = projSplitPoint['indSplit']
    if AvaProfile['z'][-1] > AvaProfile['z'][0]:
        log.info('Profile reversed')
        L = AvaProfile['s'][-1]
        AvaProfile['x'] = np.flip(AvaProfile['x'])
        AvaProfile['y'] = np.flip(AvaProfile['y'])
        AvaProfile['z'] = np.flip(AvaProfile['z'])
        AvaProfile['s'] = L - np.flip(AvaProfile['s'])
        if projSplitPoint:
            indSplit = len(AvaProfile['x']) - indSplit
            projSplitPoint['indSplit'] = indSplit
        else:
            projSplitPoint = None

    return projSplitPoint, AvaProfile
