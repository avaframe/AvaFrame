""" Opperations and transformations of rasters and lines
"""

import logging
import math
import pathlib
import numpy as np
import scipy as sp
import scipy.interpolate
import copy

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU


# create local logger
log = logging.getLogger(__name__)


def projectOnRaster(dem, Points, interp='bilinear'):
    """Projects Points on raster
    using a bilinear or nearest interpolation and returns the z coord (no for loop)

    Parameters
    -------------
    dem: dict
        dem dictionary
    Points: dict
        Points dictionary (x,y)
    interp: str
        interpolation option, between nearest or bilinear
    Returns
    -------
    Points: dict
        Points dictionary with z coordinate added or updated

    """
    header = dem['header']
    rasterdata = dem['rasterData']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    cellsize = header['cellsize']
    xcoor = Points['x']
    ycoor = Points['y']

    zcoor, ioob = projectOnGrid(xcoor, ycoor, rasterdata, csz=cellsize, xllc=xllc, yllc=yllc, interp=interp)
    Points['z'] = zcoor
    return Points, ioob


def projectOnGrid(x, y, Z, csz=1, xllc=0, yllc=0, interp='bilinear'):
    """Projects Z onto points (x,y)
    using a bilinear or nearest interpolation and returns the z coord

    Parameters
    -------------
    x: array
        x coord of the points to project
    y: array
        y coord of the points to project
    Z : 2D numpy array
        raster data
    csz: float
        cellsize corresponding to the raster data
    xllc: float
        x coord of the lower left center of the raster
    yllc: float
        y coord of the lower left center of the raster
    interp: str
        interpolation option, between nearest or bilinear
    Returns
    -------
    z : 2D numpy array
        projected data on the raster data
    ioob: int
        number of out of bounds indexes
    """
    nrow, ncol = np.shape(Z)
    # initialize outputs
    z = np.ones((np.shape(x)))*np.NaN
    dx = np.ones((np.shape(x)))*np.NaN
    dy = np.ones((np.shape(x)))*np.NaN
    f11 = np.ones((np.shape(x)))*np.NaN
    f12 = np.ones((np.shape(x)))*np.NaN
    f21 = np.ones((np.shape(x)))*np.NaN
    f22 = np.ones((np.shape(x)))*np.NaN

    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    Lxx = (x - xllc) / csz
    Lyy = (y - yllc) / csz
    Lx = copy.deepcopy(Lxx)
    Ly = copy.deepcopy(Lyy)

    # find out of bound indexes
    if interp == 'nearest':
        Lx[np.where((Lxx <= -0.5))] = np.NaN
        Ly[np.where((Lxx <= -0.5))] = np.NaN
        Lx[np.where(Lxx >= (ncol-0.5))] = np.NaN
        Ly[np.where(Lxx >= (ncol-0.5))] = np.NaN
        Lx[np.where(Lyy <= -0.5)] = np.NaN
        Ly[np.where(Lyy <= -0.5)] = np.NaN
        Lx[np.where(Lyy >= (nrow-0.5))] = np.NaN
        Ly[np.where(Lyy >= (nrow-0.5))] = np.NaN
    elif interp == 'bilinear':
        Lx[np.where((Lxx < 0))] = np.NaN
        Ly[np.where((Lxx < 0))] = np.NaN
        Lx[np.where(Lxx >= (ncol-1))] = np.NaN
        Ly[np.where(Lxx >= (ncol-1))] = np.NaN
        Lx[np.where(Lyy < 0)] = np.NaN
        Ly[np.where(Lyy < 0)] = np.NaN
        Lx[np.where(Lyy >= (nrow-1))] = np.NaN
        Ly[np.where(Lyy >= (nrow-1))] = np.NaN

    # find index of index of not nan value
    mask = ~np.isnan(Lx+Ly)
    maskInd = np.argwhere(~np.isnan(Lx+Ly))[:, 0]
    itot = len(Lx)
    iinb = len(maskInd)
    ioob = itot - iinb

    # find coordinates of the 4 nearest cornes on the raster
    Lx0 = np.floor(Lx).astype('int')
    Ly0 = np.floor(Ly).astype('int')
    Lx1 = Lx0 + 1
    Ly1 = Ly0 + 1
    # prepare for bilinear interpolation(do not take out of bound into account)
    if interp == 'nearest':
        dx[mask] = np.round(Lx[mask])
        dy[mask] = np.round(Ly[mask])
        z[mask] = Z[dy[mask].astype('int'), dx[mask].astype('int')]
    elif interp == 'bilinear':
        dx[mask] = Lx[mask] - Lx0[mask]
        dy[mask] = Ly[mask] - Ly0[mask]
        f11[mask] = Z[Ly0[mask], Lx0[mask]]
        f12[mask] = Z[Ly1[mask], Lx0[mask]]
        f21[mask] = Z[Ly0[mask], Lx1[mask]]
        f22[mask] = Z[Ly1[mask], Lx1[mask]]
        # using bilinear interpolation on the cell
        z = f11*(1-dx)*(1-dy) + f21*dx*(1-dy) + f12*(1-dx)*dy + f22*dx*dy

    return z, ioob


def resizeData(raster, rasterRef):
    """
    Reproject raster on a grid of shape rasterRef

    Parameters
    ----------
    raster : dict
        raster dictionary
    rasterRef : dict
        reference raster dictionary

    Returns
    -------
    data : 2D numpy array
        reprojected data
    dataRef : 2D numpy array
        reference data
    """
    if IOf.isEqualASCheader(raster['header'], rasterRef['header']):
        return raster['rasterData'], rasterRef['rasterData']
    else:
        headerRef = rasterRef['header']
        ncols = headerRef['ncols']
        nrows = headerRef['nrows']
        csz = headerRef['cellsize']
        xllc = headerRef['xllcenter']
        yllc = headerRef['yllcenter']
        xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
        ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
        X, Y = np.meshgrid(xgrid, ygrid)
        Points = {'x': X, 'y': Y}
        Points, _ = projectOnRaster(raster, Points, interp='bilinear')
        raster['rasterData'] = Points['z']
        return raster['rasterData'], rasterRef['rasterData']


def getMeshXY(rasterDict, cellSizeNew=None):
    """ Get x and y vectors for mesh with given number of rows and columns, lowerleftcenter and cellSize

        Parameters
        -----------
        rasterDict: dict
            'rasterData' : 2D numpy array of data
            'header': class with info on ncols, nrows, csz, xllcenter, yllcenter, noDataValue

        Returns
        --------
        x, y: 1D numpy arrays
            vector of x and y values for mesh center coordinates
        xExtent, yExtent: float
            extent of mesh from first cell center to last cell center in x/y

    """

    header = rasterDict['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    ncols = header['ncols']
    nrows = header['nrows']
    xExtent = (ncols-1) * header['cellsize']
    yExtent = (nrows-1) * header['cellsize']
    x = np.linspace(0, xExtent, ncols) + xllc
    y = np.linspace(0, yExtent, nrows) + yllc

    if cellSizeNew != None:

        nColsNew = int(xExtent/cellSizeNew+1)
        nRowsNew = int(yExtent/cellSizeNew+1)
        xNew = np.linspace(0, (nColsNew-1)*cellSizeNew, nColsNew) + xllc
        yNew = np.linspace(0, (nRowsNew-1)*cellSizeNew, nRowsNew) + yllc
        diffExtentX = xExtent - (nColsNew-1)*cellSizeNew
        diffExtentY = yExtent - (nRowsNew-1)*cellSizeNew
        return x, y, xNew, yNew, diffExtentX, diffExtentY

    else:
        return x, y, xExtent, yExtent


def remeshData(rasterFile, cellSize):
    """ compute raster data on a new mesh with cellSize using scipy RectBivariateSpline

        the new mesh is as big or smaller as the original mesh

    Parameters
    ----------
    rasterFile : str
        path to raster file
    cellSize : float
        mesh size of new mesh

    Returns
    -------
    data : dict
        remeshed data dict with data as numpy array and header info

    """

    # load data
    raster = IOf.readRaster(rasterFile, noDataToNan=True)
    header = raster['header']

    # fetch shape info and get new mesh info
    x, y, xNew, yNew, diffExtentX, diffExtentY = getMeshXY(raster, cellSizeNew=cellSize)
    data = raster['rasterData']

    log.info('Remeshed data extent difference x: %f and y %f' % (diffExtentX, diffExtentY))

    # use scipy interpolate to compute data on points of new mesh and save to raster dict
    rasterNew = sp.interpolate.RectBivariateSpline(y, x, data)(yNew, xNew, grid=True)
    raster['rasterData'] = rasterNew

    # set new header
    header['ncols'] = len(xNew)
    header['nrows'] = len(yNew)
    header['cellsize'] = cellSize

    return raster


def remeshDEM(demFile, cfgSim):
    """ change DEM cell size by reprojecting on a new grid - first check if remeshed DEM available

    the new DEM is as big or smaller as the original DEM and saved to Inputs/DEMremshed as remeshedDEMcellSize

    Interpolation is based on griddata with a cubic method. Here would be the place
    to change the order of the interpolation or to switch to another interpolation method.

    Parameters
    ----------
    demFile: str or pathlib path
        path to DEM in Inputs/
    cfgSim : configParser
        meshCellSizeThreshold : threshold under which no remeshing is done
        meshCellSize : desired cell size

    Returns
    -------
    pathDem : str
        path of DEM with desired cell size relative to Inputs/

    """

    # first check if remeshed DEM is available
    pathDem, DEMFound, allDEMNames = searchRemeshedDEM(demFile.stem, cfgSim)
    if DEMFound:
        return pathDem

    #-------- if no remeshed DEM found - remesh

    # fetch info on dem file
    dem = IOf.readRaster(demFile)
    headerDEM = dem['header']

    # fetch info on desired meshCellSize
    meshCellSize = float(cfgSim['GENERAL']['meshCellSize'])
    meshCellSizeThreshold = float( cfgSim['GENERAL']['meshCellSizeThreshold'])

    # read dem header info
    cszDEM = headerDEM['cellsize']

    # start remesh
    log.info('Remeshing the input DEM (of cell size %.4g m) to a cell size of %.4g m' % (cszDEM, meshCellSize))
    x, y, xNew, yNew, diffExtentX, diffExtentY = getMeshXY(dem, cellSizeNew=meshCellSize)
    xGrid, yGrid = np.meshgrid(x, y)
    xGrid = xGrid.flatten()
    yGrid = yGrid.flatten()
    z = dem['rasterData']
    zCopy = np.copy(z)
    zCopy = zCopy.flatten()
    mask = np.where(~np.isnan(zCopy))
    xGrid = xGrid[mask]
    yGrid = yGrid[mask]
    z = zCopy[mask]

    # create header of remeshed DEM
    headerRemeshed = {}
    headerRemeshed['cellsize'] = meshCellSize
    headerRemeshed['ncols'] = len(xNew)
    headerRemeshed['nrows'] = len(yNew)
    headerRemeshed['xllcenter'] = headerDEM['xllcenter']
    headerRemeshed['yllcenter'] = headerDEM['yllcenter']
    headerRemeshed['noDataValue'] = headerDEM['noDataValue']

    # write new DEM dictionary
    remeshedDEM = {'header': headerRemeshed}
    xNewGrid, yNewGrid = np.meshgrid(xNew, yNew)
    zNew = sp.interpolate.griddata((xGrid, yGrid), z, (xNewGrid, yNewGrid), method='cubic',
                                   fill_value=headerDEM['noDataValue'])
    log.info('Remeshed data extent difference x: %f and y %f' % (diffExtentX, diffExtentY))
    remeshedDEM['rasterData'] = zNew

    # save remeshed DEM
    pathToDem = pathlib.Path(cfgSim['GENERAL']['avalancheDir'], 'Inputs', 'DEMremeshed')
    fU.makeADir(pathToDem)
    outFile = pathToDem / ('%s_remeshedDEM%.2f.asc' % (demFile.stem, remeshedDEM['header']['cellsize']))
    if outFile.name in allDEMNames:
        message = 'Name for saving remeshedDEM already used: %s' % outFile.name
        log.error(message)
        raise FileExistsError(message)

    IOf.writeResultToAsc(remeshedDEM['header'], remeshedDEM['rasterData'], outFile, flip=True)
    log.info('Saved remeshed DEM to %s' % outFile)
    pathDem = str(pathlib.Path('DEMremeshed', outFile.name))

    return pathDem


def searchRemeshedDEM(demName, cfgSim):
    """ search if remeshed DEM with correct name and cell size already available

        Parameters
        -----------
        demName: str
            name of DEM file in Inputs/
        cfgSim: configparser object
            configuration settings: avaDir, meshCellSize, meshCellSizeThreshold

        Returns
        --------
        remshedDEM: dict
            dictionary of remeshed DEM if not found empty dict
        DEMFound: bool
            flag if dem is found
        allDEMNames: list
            list of all names of dems found in Inputs/DEMremeshed
    """

    # path to remeshed DEM folder
    pathToDems = pathlib.Path(cfgSim['GENERAL']['avalancheDir'], 'Inputs', 'DEMremeshed')
    DEMFound = False
    pathDem = ''
    allDEMNames = []

    # fetch info on desired meshCellSize
    meshCellSize = float(cfgSim['GENERAL']['meshCellSize'])
    meshCellSizeThreshold = float( cfgSim['GENERAL']['meshCellSizeThreshold'])

    # check if DEM is available
    if pathToDems.is_dir():
        # look for dems and check if cellSize within tolerance and origin matches
        demFiles = list(pathToDems.glob('*.asc'))
        allDEMNames = [d.name for d in demFiles]
        for demF in demFiles:
            headerDEM = IOf.readASCheader(demF)
            if abs(meshCellSize - headerDEM['cellsize']) < meshCellSizeThreshold and demName in demF.stem:
                log.info('Remeshed DEM found: %s cellSize: %.5f' % (demF.name, headerDEM['cellsize']))
                DEMFound = True
                pathDem = str(pathlib.Path('DEMremeshed', demF.name))
                continue
            else:
                log.debug('Remeshed dem found %s with cellSize %.2f - not used' %
                    (demF, headerDEM['cellsize']))

    else:
        log.debug('Directory %s does not exist' % pathToDems)

    return pathDem, DEMFound, allDEMNames


def prepareLine(dem, avapath, distance=10, Point=None):
    """Resample and project line on dem
    1- Resample the avapath line with a max intervall of distance=10m
    between points (projected distance on the horizontal plane).
    2- Make avalanche profile out of the path (affect a z value using the dem)
    3- Get projection of points on the profil (closest point)

    Parameters
    -----------
    dem: dict
        dem dictionary
    avapath: dict
        line dictionary
    distance: float
        resampling distance
    Point: dict
        a point dictionary (optional, can contain several point)

    Returns
    -------
    AvaProfile: dict
        the resampled avapath with the z coordinate
    projPoint: dict
        point dictionary projected on the profile (if several points
        were give in input, only the closest point to the profile
        is projected)
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
    ResampAvaPath, _ = projectOnRaster(dem, ResampAvaPath)
    ResampAvaPath['s'] = s
    AvaProfile = ResampAvaPath
    # find split point by computing the distance to the line
    if Point:
        projPoint = findSplitPoint(AvaProfile, Point)
    else:
        projPoint = None

    return AvaProfile, projPoint


def findPointOnDEM(dem, vDirX, vDirY, vDirZ, zHighest, xFirst, yFirst, zFirst):
    """ find point on dem given a direction and a z value to reach

    Parameters
    -----------
    dem: dict
        dem dict
    vDirX, vDirY, vDirZ: floats
        x, y and z components of the direction in which to extend
    zHighest: float
        z value to reach
    xFirst, yFirst, zFirst: floats
        x, y and z coordinates of the starting point

    Returns
    --------
    xExtTop, yExtTop, zExtTop:floats
        x, y and z coordinates of the point found
    """
    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    zRaster = dem['rasterData']
    gamma = (zHighest - zFirst) / vDirZ * np.linspace(0.25, 2, 100)
    xArray = xFirst + gamma * vDirX
    yArray = yFirst + gamma * vDirY
    zArray, _ = projectOnGrid(xArray, yArray, zRaster, csz=csz, xllc=xllc, yllc=yllc, interp='bilinear')
    idx = np.nanargmin(np.abs(zArray - np.array([zHighest])))
    xExtTop = np.array([xFirst + gamma[idx] * vDirX])
    yExtTop = np.array([yFirst + gamma[idx] * vDirY])
    zExtTop = np.array([zArray[idx]])
    return xExtTop, yExtTop, zExtTop


def findSplitPoint(AvaProfile, Points):
    """ Finds the closest point in Points to the AvaProfile and returns
    its projection on AvaProfile.

    Parameters
    -----------
    AvaProfile: dict
        line dictionary with x and y coordinates
    Point: dict
        a point dictionary

    Returns
    -------
    projPoint: dict
        point dictionary projected on the profile (if several points
        were give in input, only the closest point to the profile
        is projected)
    """
    xcoor = AvaProfile['x']
    ycoor = AvaProfile['y']
    Dist = np.empty((0))
    IndSplit = np.empty((0))
    for i in range(len(Points['x'])):
        dist = np.sqrt((xcoor - Points['x'][i]) ** 2 + (ycoor - Points['y'][i])**2)
        indSplit = np.argmin(dist)
        IndSplit = np.append(IndSplit, indSplit)
        Dist = np.append(Dist, dist[indSplit])

    ind = np.argmin(Dist)
    indSplit = int(IndSplit[ind])
    projPoint = {}
    projPoint['x'] = AvaProfile['x'][indSplit]
    projPoint['y'] = AvaProfile['y'][indSplit]
    projPoint['z'] = AvaProfile['z'][indSplit]
    projPoint['s'] = AvaProfile['s'][indSplit]
    projPoint['indSplit'] = indSplit
    return projPoint


def checkProfile(AvaProfile, projSplitPoint=None):
    """ check that the avalanche profiles goes from top to bottom
    flip it if not and adjust the splitpoint in consequence

    Parameters
    -----------
    AvaProfile: dict
        line dictionary with x and y coordinates
    projSplitPoint: dict
        a point dictionary already projected on the AvaProfile

    Returns
    -------
    AvaProfile: dict
        avaprofile, fliped if needed
    projSplitPoint: dict
        point dictionary
    """
    if projSplitPoint:
        indSplit = projSplitPoint['indSplit']
    if AvaProfile['z'][-1] > AvaProfile['z'][0]:
        log.info('Profile reversed')
        AvaProfile['x'] = np.flip(AvaProfile['x'])
        AvaProfile['y'] = np.flip(AvaProfile['y'])
        AvaProfile['z'] = np.flip(AvaProfile['z'])
        try:
            L = AvaProfile['s'][-1]
            AvaProfile['s'] = L - np.flip(AvaProfile['s'])
        except KeyError:
            pass

        if projSplitPoint:
            indSplit = len(AvaProfile['x']) - indSplit - 1
            projSplitPoint['indSplit'] = indSplit
            AvaProfile['indSplit'] = indSplit
        else:
            projSplitPoint = None
            AvaProfile['indSplit'] = None

    return projSplitPoint, AvaProfile


def findAngleProfile(tmp, ds, dsMin):
    """
    Find the beta point: first point under the beta value given in
    prepareAngleProfile. Make sure that at least dsMin meters behind the point
    are also under the beta value otherwise keep searching

    Parameters
    ----------
    tmp: 1D numpy array
        index array of point in profile with slope bellow the given beta angle
        and bellow the splitPoint
    ds: 1D numpy array
        distance between points discribed in tmp
    dsMin: float
        threshold distance [m] for looking for the beta point (at least dsMin meters below
        beta degres)

    Returns
    -------
    idsAnglePoint: int
        index of beta point
    """
    noBetaFoundMessage = 'No Beta point found. Check your pathAB.shp and splitPoint.shp.'
    i = 0
    condition = True
    if np.size(tmp) == 0:
        raise IndexError(noBetaFoundMessage)
    while (i <= np.size(tmp) and condition):
        ind = tmp[i]
        j = 0
        dist = 0
        while dist < dsMin:
            try:
                condition = condition and (tmp[i+j+1] == ind+j+1)
                dist = dist + ds[i + j]
            except IndexError:
                raise IndexError(noBetaFoundMessage)
            if not condition:
                i = i + j + 1
                break
            j = j + 1
        if condition:
            idsAnglePoint = ind
            break
        condition = True
    return idsAnglePoint


def prepareAngleProfile(beta, AvaProfile):
    """Prepare inputs for findAngleProfile function
    Read profile (s, z), compute the slope Angle
    look for points for which the slope is under the given Beta value and
    that are located downstream of the splitPoint

    Parameters
    ----------
    beta: float
        beta angle in degrees
    AvaProfile: dict
        profile dictionary, s, z and a split point(optional)
    Returns
    -------
    angle: 1D numpy array
        profile angle
    tmp: 1D numpy array
        index array of point in profile with slope bellow the given beta angle
        and bellow the splitPoint
    ds: 1D numpy array
        distance between points discribed in tmp
    """

    s = AvaProfile['s']
    z = AvaProfile['z']
    try:
        indSplit = AvaProfile['indSplit']
        CuSplit = s[indSplit]
    except KeyError:
        log.warning('No split Point given!')
        CuSplit = 0
    ds = np.abs(s - np.roll(s, 1))
    dz = np.roll(z, 1) - z
    ds[0] = 0.0
    dz[0] = 0.0
    angle = np.rad2deg(np.arctan2(dz, ds))
    # get all values where Angle < 10 but >0
    # get index of first occurance and go one back to get previous value
    # (i.e. last value above 10 deg)
    # tmp = x[(angle < 10.0) & (angle > 0.0) & (x > 450)]
    tmp = np.where((angle <= beta) & (s > CuSplit))
    tmp = np.asarray(tmp).flatten()
    ds = ds[tmp]
    return angle, tmp, ds


def isCounterClockWise(path):
    """ Determines if a polygon path is mostly clockwise or counter clockwise

    https://stackoverflow.com/a/45986805/15887086

    Parameters
    ----------
    path: matplotlib.path
        polygon path
    Returns
    -------
    isCounterCloc1: int
        1 if the path is counter clockwise, 0 otherwise
    """
    v = path.vertices-path.vertices[0, :]
    a = np.arctan2(v[1:, 1], v[1:, 0])
    isCounterClock = (a[1:] >= a[:-1]).astype(int).mean() >= 0.5
    return isCounterClock


def path2domain(xyPath, rasterTransfo):
    """Creates a domain (irregular raster) along a path,
    given the path xyPath, a domain width and a raster cellsize

    Parameters:
    -------------
    xyPath: dict
        line dictionary with coordinates x and y
    rasterTransfo: dict
        rasterTransfo['w']: float
            Domain width
        rasterTransfo['cellSizeSL']: float
            cellsize expected for the new raster

    Returns:
    ---------
    rasterTransfo: dict
        rasterTransfo updated with xp, yp Arrays determining a path of width w along a line

        rasterTransfo['DBXl']:
            x coord of the left boundary
        rasterTransfo['DBXr']:
            x coord of the right boundary
        rasterTransfo['DBYl']:
            y coord of the left boundary
        rasterTransfo['DBYr']:
            y coord of the right boundary

    [Fischer2013] Fischer, Jan-Thomas. (2013).
    A novel approach to evaluate and compare computational snow avalanche
    simulation.
    Natural Hazards and Earth System Sciences.
    13. 1655-. 10.5194/nhess-13-1655-2013.
    Uwe Schlifkowitz/ BFW, June 2011
    """
    csz = rasterTransfo['cellSizeSL']
    x = xyPath['x']
    y = xyPath['y']
    # compute the non dimensional width
    w = rasterTransfo['domainWidth']/2/csz
    # remove scaling due to cellsize
    x = x/csz
    y = y/csz

    # Difference between x- bzw. y-Coordinates of Polyline
    # first and last  Vertex: Difference between this and the next
    # other vertices: Difference between previous and next
    dx = np.array((x[1]-x[0]))
    dy = np.array((y[1]-y[0]))
    n = len(x)
    for i in range(2, n):
        dx = np.append(dx, (x[i]-x[i-2])/2.)
        dy = np.append(dy, (y[i]-y[i-2])/2.)

    dx = np.append(dx, x[-1]-x[-2])
    dy = np.append(dy, y[-1]-y[-2])

    # Direction of normal vector of difference,
    # a.k.a. bisecting line of angle
    d = np.arctan2(dy, dx) + math.pi/2

    # x- and y-Coordinates (left and right) of path edges,
    # total width w
    # x-KOO[left right]
    DBXl = np.array((x + w * np.cos(d)))
    DBXr = np.array((x + w * np.cos(d + math.pi)))
    # y-KOO[left right]
    DBYl = np.array((y + w * np.sin(d)))
    DBYr = np.array((y + w * np.sin(d + math.pi)))

    rasterTransfo['DBXl'] = DBXl
    rasterTransfo['DBXr'] = DBXr
    rasterTransfo['DBYl'] = DBYl
    rasterTransfo['DBYr'] = DBYr

    return rasterTransfo


def areaPoly(X, Y):
    """Gauss's area formula to calculate polygon area

    Parameters
    ----------
    X: 1D numpy array
        x coord of the vertices
    Y: 1D numpy array
        y coord of the vertices
    (Without repeating the first vertex!!!)
    Returns
    -------
    area: float
        Area of the polygon
    """

    X = np.append(X, X[0])
    Y = np.append(Y, Y[0])
    area = 0
    for i in range(np.size(X)-1):
        area = area + (X[i]*Y[i+1]-Y[i]*X[i+1])/2
    return area


def checkOverlap(toCheckRaster, refRaster, nameToCheck, nameRef, crop=False):
    """Check if two rasters overlap

    Parameters
    ----------
    toCheckRaster : 2D numpy array
        Raster to check
    refRaster : 2D numpy array
        refference Raster
    nameToCheck: str
        name of raster that might overlap
    nameRef: str
        name of reference raster
    crop : boolean
        if True, remove overlaping part and send a warning
    Returns
    -------
    toCheckRaster: 2D numpy array
        if crop is True, return toCheckRaster without the overlaping part and send a
        warning if needed
        if crop is False, return error if Rasters overlap otherwise return toCheckRaster
    """
    mask = (toCheckRaster > 0) & (refRaster > 0)
    if mask.any():
        if crop:
            toCheckRaster[mask] = 0
            message = '%s area feature overlapping with %s area - removing the overlapping part' % (nameToCheck, nameRef)
            log.warning(message)
        else:
            message = '%s area features overlapping with %s area - this is not allowed' % (nameToCheck, nameRef)
            log.error(message)
            raise AssertionError(message)

    return toCheckRaster
