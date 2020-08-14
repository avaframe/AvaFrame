"""
Opperations and transformations of rasters and lines
"""

import math
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
    rasterdata = dem['rasterData']
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

    ResampAvaPath = avapath
    ResampAvaPath['x'] = xcoornew
    ResampAvaPath['y'] = ycoornew
    ResampAvaPath = projectOnRaster(dem, ResampAvaPath)
    ResampAvaPath['s'] = s
    AvaProfile = ResampAvaPath
    # find split point by computing the distance to the line
    if splitPoint:
        projSplitPoint = findSplitPoint(AvaProfile, splitPoint)
    else:
        projSplitPoint = None

    return AvaProfile, projSplitPoint


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
    return projSplitPoint


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


def bresenham(x0, y0, x1, y1, cs):
    """
    RASTERIZE - bresenham algorithmus - JT 2011

    input: x0, y0, x1, y1,cellsize
    output: array of x y coodinates of cells hit inbetween

    C IMPLEMENTIERUNG von http://de.wikipedia.org/wiki/Bresenham-Algorithmus
    void line(int x0, int y0, int x1, int y1)
     {
       int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
       int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
       int err = dx+dy, e2; /* error value e_xy */

       for(;;){  /* loop */
         setPixel(x0,y0);
         if (x0==x1 && y0==y1) break;
         e2 = 2*err;
         if (e2 > dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
         if (e2 < dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
       }
     }
    """
    # normalize Cellsize cs to 1
    x0 = round(x0/cs)
    x1 = round(x1/cs)
    y0 = round(y0/cs)
    y1 = round(y1/cs)

    dx = abs(x1-x0)
    dy = abs(y1-y0)
    sx = np.sign(x1-x0)  # step in x direction
    sy = np.sign(y1-y0)  # step in y direction
    err = dx-dy

    z = []
    while True:
        z.append([x0*cs, y0*cs])
        if x0 == x1 and y0 == y1:  # if no step exists we are already there
            break
        e2 = 2*err
        if (e2 > -dy):
            err -= dy
            x0 += sx
        if (e2 < dx):
            err += dx
            y0 += sy

    return z


def path2domain(x, y, w, csz):
    """
    path2domain
    Usage:
        [xp, yp, ] = path2domain(x, y, w, csz)
       Input:
           x, y:   Polyline Coordinates
           w:      Domain width
           csz:    cell size
       Output:
           xp, yp: Arrays determining a path of width w along a polyline

       Uwe Schlifkowitz/ BFW, June 2011
    """
#    Difference between x- bzw. y-Coordinates of Polyline
#    first and last  Vertex: Difference between this and the next
#    other vertices: Difference between previous and next
    dx = np.array((x[1]-x[0]))
    dy = np.array((y[1]-y[0]))
    for i in range(2, len(x)):
        dx = np.append(dx, (x[i]-x[i-2])/2.)
        dy = np.append(dy, (y[i]-y[i-2])/2.)

    dx = np.append(dx, x[len(x)-1]-x[len(x)-2])
    dy = np.append(dy, y[len(x)-1]-y[len(x)-2])

#    Direction of normal vector of difference,
#    a.k.a. bisecting line of angle
    d = np.arctan2(dy, dx) + math.pi/2

#    x- and y-Coordinates (left and right) of path edges,
#    total width w
#    x-KOO[left right]
    OX = np.array((x + w * np.cos(d),
                   x + w * np.cos(d + math.pi)))
#    y-KOO[left right]
    OY = np.array((y + w * np.sin(d),
                   y + w * np.sin(d + math.pi)))

#    AK 2013
#    x- and y-Coordinates (left and right) of path edges,
#    total width w + shift for area/rastersize
    # x-KOO[[left],[right]]
    OOX = np.array((x + (w+csz/2) * np.cos(d),
                    x + (w+csz/2) * np.cos(d + math.pi)))
    # y-KOO[[left],[right]]
    OOY = np.array((y + (w+csz/2) * np.sin(d),
                    y + (w+csz/2) * np.sin(d + math.pi)))

#    x-KOO[[left right], ... ,[left right]]
    OOXX = np.zeros((len(OOX[0])*2, len(OOX)))
#    y-KOO[[left right], ... ,[left right]]
    OOYY = np.zeros((len(OOY[0])*2, len(OOY)))
#        vorwärts
    OOXX[0:-1:2, 0] = OOX[0] + csz/2. * np.cos(d + 1./2*math.pi)
    OOXX[0:-1:2, 1] = OOX[1] + csz/2. * np.cos(d + 1./2*math.pi)
    OOYY[0:-1:2, 0] = OOY[0] + csz/2. * np.sin(d + 1./2*math.pi)
    OOYY[0:-1:2, 1] = OOY[1] + csz/2. * np.sin(d + 1./2*math.pi)
#        rückwärts
    OOXX[1::2, 0] = OOX[0] + csz/2. * np.cos(d + 3./2*math.pi)
    OOXX[1::2, 1] = OOX[1] + csz/2. * np.cos(d + 3./2*math.pi)
    OOYY[1::2, 0] = OOY[0] + csz/2. * np.sin(d + 3./2*math.pi)
    OOYY[1::2, 1] = OOY[1] + csz/2. * np.sin(d + 3./2*math.pi)

    return OX, OY, OOXX, OOYY


def poly2mask_simple(ydep, xdep, ncols, nrows):
    """
    poly2mask_simple
    Create a mask from a polyline
    Usage:
        mask = poly2mask_simple(ydep, xdep, ncols, nrows)
       Input:
           ydep, xdep:      Polyline Coordinates
           ncols, nrows:    Raster size
       Output:
           mask:            Raster of the polyline mask

    """
    mask = np.zeros((nrows, ncols))
    xyframe = bresenham(xdep[0], ydep[0], xdep[1], ydep[1], 1)
    xyframe = np.delete(xyframe, -1, 0)
    xyframe = np.transpose(xyframe)
    for i in range(1, len(xdep)-1):
        xyline = bresenham(xdep[i], ydep[i], xdep[i+1], ydep[i+1], 1)
        # last point is first point of the next line
        xyline = np.delete(xyline, -1, 0)
        xyline = np.transpose(xyline)
        xyframe = np.hstack((xyframe, xyline))

    xyline = bresenham(xdep[-1], ydep[-1], xdep[0], ydep[0], 1)
    xyline = np.delete(xyline, -1, 0)
    xyline = np.transpose(xyline)
    xyframe = np.hstack((xyframe, xyline))
    for i in range(0, len(xyframe[0, :])):
        mask[xyframe[0, i], xyframe[1, i]] = 1

    # filling the inside of the polygon with ones
    i = xyframe[0]
    j = xyframe[1]
    mv, nv = np.meshgrid(np.linspace(0, nrows-1, nrows),
                         np.linspace(0, ncols-1, ncols))  # create index space
    mask = inpolygon(mv, nv, i, j)
    mask = np.transpose(mask)
    return mask


def inpolygon(X, Y, xv, yv):
    """
    inpolygon
    For a polygon defined by vertex points (xv, yv),
    returns a np array of size X with ones if the points (X, Y)
    are inside (or on the boundary) of the polygon;
    Otherwise, returns zeros.
    Usage:
        mask = inpolygon(X, Y, xv, yv)
       Input:
           X, Y:      Set of points to check
           xv, yv:    polygon vertex points
       Output:
           mask:      np array of zeros and ones

    Octave Implementation [IN, ON] = inpolygon (X, Y, xv, yv)
    """
    npol = len(xv)
    lx = np.shape(X)[0]
    ly = np.shape(Y)[1]
    IN = np.zeros(np.shape(X))
    j = npol-1
    for i in range(npol-1):
        delta_xv = xv[j] - xv[i]
        delta_yv = yv[j] - yv[i]
        # distance = [distance from (X,Y) to edge] * length(edge)
        distance = delta_xv*(Y-yv[i]) - (X-xv[i])*delta_yv
        # is Y between the y-values of edge i,j
        # AND (X,Y) on the left of the edge ?
        for ii in range(lx):
            for jj in range(ly):
                if (((yv[i] <= Y[ii][jj] and Y[ii][jj] < yv[j]) or (yv[j] <= Y[ii][jj] and Y[ii][jj] < yv[i])) and 0 < distance[ii][jj]*delta_yv):
                    if IN[ii][jj] == 0:
                        IN[ii][jj] = 1
                    else:
                        IN[ii][jj] = 0
        j = i
    for i in range(npol-1):
        IN[yv[i]][xv[i]] = 1

    return IN
