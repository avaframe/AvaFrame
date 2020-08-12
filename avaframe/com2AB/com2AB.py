#!/usr/bin/env python
# coding: utf-8
""" Main file for module com2AB - Alpha Beta
"""

import os
import glob
import pickle
import logging
import numpy as np
import matplotlib.pyplot as plt
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)
debugPlot = False

# Local imports


def setEqParameters(smallAva, customParam):
    """Set alpha beta equation parameters to
    - standard (default)
    - small avalanche
    - custom
    TODO: test
    """

    eqParameters = {}

    if smallAva is True:
        log.info('Using small Avalanche Setup')
        eqParameters['k1'] = 0.933
        eqParameters['k2'] = 0.0
        eqParameters['k3'] = 0.0088
        eqParameters['k4'] = -5.02
        eqParameters['SD'] = 2.36

        ParameterSet = "Small avalanches"

    elif customParam:
        log.info('Using custom Avalanche Setup')
        eqParameters['k1'] = customParam['k1']
        eqParameters['k2'] = customParam['k2']
        eqParameters['k3'] = customParam['k3']
        eqParameters['k4'] = customParam['k4']
        eqParameters['SD'] = customParam['SD']

        ParameterSet = "Custom"

    else:
        log.info('Using standard Avalanche Setup')
        eqParameters['k1'] = 1.05
        eqParameters['k2'] = -3130.0
        eqParameters['k3'] = 0.0
        eqParameters['k4'] = -2.38
        eqParameters['SD'] = 1.25

        ParameterSet = "Standard"
    eqParameters['ParameterSet'] = ParameterSet
    return eqParameters


def com2ABMain(DGM, Avapath, SplitPoint, saveOutPath, cfgsetup):
    """ Loops on the given Avapath and runs com2AB to compute AlpahBeta model
    Inputs : DEM header and rater (as np array),
            Avapath and Split points .shp file,
            optional output save path,
            avalanche type,
            reamplind lenght for the Avapath
    Outputs : writes raw results to saveOutPath
    """
    smallAva = cfgsetup.getboolean('smallAva')
    customParam = cfgsetup.getboolean('customParam')
    # customParam = str(customParam or None)
    distance = float(cfgsetup['distance'])

    NameAva = Avapath['Name']
    StartAva = Avapath['Start']
    LengthAva = Avapath['Length']

    for i in range(len(NameAva)):
        name = NameAva[i]
        start = StartAva[i]
        end = start + LengthAva[i]
        avapath = {}
        avapath['x'] = Avapath['x'][int(start):int(end)]
        avapath['y'] = Avapath['y'][int(start):int(end)]
        avapath['Name'] = name
        com2AB(DGM, avapath, SplitPoint, saveOutPath,
               smallAva, customParam, distance)


def com2AB(DGM, avapath, splitPoint, OutPath,
           smallAva, customParam, distance):
    """ Computes the AlphaBeta model given an input raster (of the DEM),
    an avalanche path and split points
    Inputs : DEM header and rater (as np array),
            single avapath as np array,
            Split points as np array,
            output save path,
            avalanche type,
            resamplind lenght for the Avapath
    Outputs : writes raw results to OutPath
    """
    name = avapath['Name']
    abVersion = '4.1'
    log.info('Running Version: %s ', abVersion)
    log.info('Running Alpha Beta model on: %s ', name)
    eqParams = setEqParameters(smallAva, customParam)

    # TODO: make rest work with dict

    # read inputs, ressample ava path
    # make pofile and project split point on path
    AvaProfile, SplitPoint, splitPoint = prepareLine(
        DGM, avapath, splitPoint, distance)

    # Sanity check if first element of AvaProfile[3,:]
    # (i.e z component) is highest:
    # if not, flip all arrays
    splitPoint, AvaProfile = checkProfile(splitPoint, AvaProfile)

    AvaProfile['indSplit'] = splitPoint['indSplit']  # index of split point

    eqOut = calcAB(AvaProfile, eqParams)
    savename = name + '_com2AB_eqparam.pickle'
    save_file = os.path.join(OutPath, savename)
    with open(save_file, 'wb') as handle:
        pickle.dump(eqParams, handle, protocol=pickle.HIGHEST_PROTOCOL)
    savename = name + '_com2AB_eqout.pickle'
    save_file = os.path.join(OutPath, savename)
    with open(save_file, 'wb') as handle:
        pickle.dump(eqOut, handle, protocol=pickle.HIGHEST_PROTOCOL)


def projectOnRaster(DGM, Points):
    """ Projects the points Points on Raster using a bilinear interpolation
    and returns the z coord
    Input :
    Points: list of points (x,y) 2 rows as many columns as Points
    Output:
    PointsZ: list of points (x,y,z) 3 rows as many columns as Points

    TODO: test
    """
    header = DGM['header']
    rasterdata = DGM['rasterdata']
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


def prepareLine(DGM, avapath, splitPoint, distance):
    """ 1- Resample the avapath line with a max intervall of distance=10m
    between points (projected distance on the horizontal plane).
    2- Make avalanch profile out of the path (affect a z value using the DEM)
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
    ResampAvaPath = projectOnRaster(DGM, ResampAvaPath)
    ResampAvaPath['s'] = s
    AvaProfile = ResampAvaPath
    # find split point by computing the distance to the line
    SplitPoint, splitPoint = findSplitPoint(AvaProfile, splitPoint)

    return AvaProfile, SplitPoint, splitPoint


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
    SplitPoint = {}
    SplitPoint['x'] = AvaProfile['x'][indSplit]
    SplitPoint['y'] = AvaProfile['y'][indSplit]
    SplitPoint['z'] = AvaProfile['z'][indSplit]
    SplitPoint['s'] = AvaProfile['s'][indSplit]
    splitPoint['indSplit'] = indSplit
    return SplitPoint, splitPoint



def readABinputs(cfgINPATH):
    cfgAva = cfgINPATH['AvaLancheRoot']
    cfgPath = {}

    ProfileLayer = glob.glob(cfgAva + '/Inputs/LINES/*.shp')
    cfgPath['ProfileLayer'] = ''.join(ProfileLayer)

    DGMSource = glob.glob(cfgAva + '/Inputs/*.asc')
    try:
        assert len(DGMSource)==1, 'There should be only one and only one DEM .asc file in ' + cfgAva + '/Inputs/'
    except AssertionError as e:
        raise

    cfgPath['DGMSource'] = ''.join(DGMSource)

    SplitPointSource = glob.glob(cfgAva + '/Inputs/POINTS/*.shp')
    cfgPath['SplitPointSource'] = ''.join(SplitPointSource)

    saveOutPath = os.path.join(cfgAva, 'Outputs/')
    if not os.path.exists(saveOutPath):
        # log.info('Creating output folder %s', saveOutPath)
        os.makedirs(saveOutPath)
    cfgPath['saveOutPath'] = saveOutPath

    DefaultName = str(cfgAva).split('/')[-1]
    cfgPath['DefaultName'] = DefaultName

    return cfgPath



def readRaster(fname):
    """ Read raster file (.asc)"""

    log.info('Reading DEM : %s', fname)
    header = IOf.readASCheader(fname)
    rasterdata = IOf.readASCdata2numpyArray(fname, header)
    rasterdata[rasterdata == header.noDataValue] = np.NaN
    DGM = {}
    DGM['header'] = header
    DGM['rasterdata'] = np.flipud(rasterdata)
    return DGM


def readAvaPath(fname, defname, header):
    """ Read avalanche path from  .shp"""

    log.info('Reading avalanche path : %s ', fname)
    Avapath = shpConv.SHP2Array(fname, defname)
    coordx = Avapath['x']
    coordy = Avapath['y']
    if len(coordx) < 2:
        raise ValueError('Ava path file should contain at least 2 columns')
    for i in range(len(coordx)):
        Lx = int(np.floor((coordx[i] - header.xllcorner) /
                          header.cellsize))
        Ly = int(np.floor((coordy[i] - header.yllcorner) /
                          header.cellsize))
        if (Ly < 0 or Ly > header.nrows or Lx < 0 or Lx > header.ncols):
            raise ValueError('Avalanche path exceeds DEM extend')
    return Avapath


def readSplitPoint(fname, header):
    """ Read split point path from .shp"""

    log.info('Reading split point : %s ', fname)
    defname = 'SHP'
    SplitPoint = shpConv.SHP2Array(fname, defname)
    splitPointx = SplitPoint['x']
    splitPointy = SplitPoint['y']
    if len(splitPointx) < 2:
        raise ValueError('split point file should contain at least 2 columns')
    for i in range(len(splitPointx)):
        Lx = int(np.floor((splitPointx[i] - header.xllcorner) /
                          header.cellsize))
        Ly = int(np.floor((splitPointy[i] - header.yllcorner) /
                          header.cellsize))
        if (Ly < 0 or Ly > header.nrows or Lx < 0 or Lx > header.ncols):
            raise ValueError('Split point is not on the DEM')
    return SplitPoint


def checkProfile(SplitPoint, AvaProfile):
    """ check that the avalanche profiles goes from top to bottom """
    indSplit = SplitPoint['indSplit']
    if AvaProfile['z'][-1] > AvaProfile['z'][0]:
        log.info('Profile reversed')
        L = AvaProfile['s'][-1]
        AvaProfile['x'] = np.flip(AvaProfile['x'])
        AvaProfile['y'] = np.flip(AvaProfile['y'])
        AvaProfile['z'] = np.flip(AvaProfile['z'])
        AvaProfile['s'] = L - np.flip(AvaProfile['s'])
        indSplit = len(AvaProfile['x']) - indSplit

    return SplitPoint, AvaProfile


def find_10Point(tmp, delta_ind):
    """ find the beta point: first point under 10°
     (make sure that the delta_ind next indexes are also under 10°)
     otherwise keep looking
     """
    i = 0
    while True:
        ind = tmp[0][i]
        condition = True
        for j in range(delta_ind):
            condition = condition and (tmp[0][i+j+1] == ind+j+1)
            if not condition:
                i = i + j + 1
                break
        if condition:
            ids_10Point = ind - 1
            break
    return ids_10Point


def calcAB(AvaProfile, eqParameters):
    """
    Calculate Alpha Beta for data in eqInput according to chosen eqParameters
    """
    log.info("Calculating alpha beta")
    k1 = eqParameters['k1']
    k2 = eqParameters['k2']
    k3 = eqParameters['k3']
    k4 = eqParameters['k4']
    SD = eqParameters['SD']

    s = AvaProfile['s']
    z = AvaProfile['z']
    distance = s[1] - s[0]
    delta_ind = max(int(np.floor(30/distance)), 1)
    indSplit = AvaProfile['indSplit']
    ds = np.abs(s - np.roll(s, 1))
    dz = np.abs(z - np.roll(z, 1))
    ds[0] = 0.0
    dz[0] = 0.0
    angle = np.rad2deg(np.arctan2(dz, ds))
    CuSplit = s[indSplit]
    # TODO SPLIT POINT READING
    # get all values where Angle < 10 but >0
    # get index of first occurance and go one back to get previous value
    # (i.e. last value above 10 deg)
    # tmp = x[(angle < 10.0) & (angle > 0.0) & (x > 450)]

    tmp = np.where((angle < 10.0) & (angle > 0.0) & (s > CuSplit))

    # find the beta point: first point under 10°
    # (make sure that the 30 next meters are also under 10°)
    ids_10Point = find_10Point(tmp, delta_ind)
    if debugPlot:
        plt.figure(figsize=(10, 6))
        plt.plot(s, angle)
        plt.plot(s[ids_10Point], angle[ids_10Point], 'or')
        plt.axhline(y=10, color='0.8',
                    linewidth=1, linestyle='-.', label='10^\circ line')
        plt.show()

    # Do a quadtratic fit and get the polynom for 2nd derivative later
    zQuad = np.polyfit(s, z, 2)
    poly = np.poly1d(zQuad)
    # Get H0: max - min for parabola
    H0 = max(poly(s)) - min(poly(s))
    # get beta
    dz_beta = z[0] - z[ids_10Point]
    beta = np.rad2deg(np.arctan2(dz_beta, s[ids_10Point]))
    # get Alpha
    alpha = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4

    # get Alpha standard deviations
    SDs = [SD, -1*SD, -2*SD]
    alphaSD = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4 + SDs

    AvaProfile['CuSplit'] = CuSplit
    AvaProfile['ids_10Point'] = ids_10Point
    AvaProfile['poly'] = poly
    AvaProfile['beta'] = beta
    AvaProfile['alpha'] = alpha
    AvaProfile['SDs'] = SDs
    AvaProfile['alphaSD'] = alphaSD
    return AvaProfile
