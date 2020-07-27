#!/usr/bin/env python
# coding: utf-8
""" Main file for module com2AB - Alpha Beta
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
colors = ["#393955", "#8A8A9B", "#E9E940"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

# Local imports
import avaframe
import avaframe.com2AB.com2ABCfg as conf
import avaframe.com2AB.IO_functionality as IOf

def setEqParameters(smallAva=False,customParam=[]):
    """Set alpha beta equation parameters to
    - standard (default)
    - small avalanche
    - custom
    TODO: test
    """

    eqParameters = []

    if smallAva is True:
        print('[com2AB] Using small Avalanche Setup')
        eqParameters['k1'] = 0.933
        eqParameters['k2'] = 0.0
        eqParameters['k3'] = 0.0088
        eqParameters['k4'] = -5.02
        eqParameters['SD'] = 2.36

        ParameterSet = "Small avalanches"

    elif customParam:
        print('[com2AB] Using standard Avalanche Setup')
        eqParameters['k1'] = customParam['k1']
        eqParameters['k2'] = customParam['k2']
        eqParameters['k3'] = customParam['k3']
        eqParameters['k4'] = customParam['k4']
        eqParameters['SD'] = customParam['SD']

        ParameterSet = "Custom"

    else:
        print('[com2AB] Using standard Avalanche Setup')
        eqParameters['k1'] = 1.05
        eqParameters['k2'] = -3130.0
        eqParameters['k3'] = 0.0
        eqParameters['k4'] = -2.38
        eqParameters['SD'] = 1.25

        ParameterSet = "Standard"

    return eqParameters, ParameterSet


def com2ABMain(header, rasterdata, avapath, splitPoint, saveOutPath='./',
               smallAva=False):
    """ Computes the AlphaBeta model given an input raster (of the DEM),
    an avalanche path and a split point
    """

    abVersion = '4.1'
    print('Running Version: ', abVersion)

    k1, k2, k3, k4, ParameterSet = setEqParameters(smallAva)

    # read inputs, ressample ava path
    # make pofile and project split point on path
    AvaProfile, SplitPoint, indSplit = PrepareLine(
        header, rasterdata, avapath, splitPoint, distance=10)

    # Sanity check if first element of AvaProfile[3,:]
    # (i.e z component) is highest:
    # if not, flip all arrays
    # TODO: break out to testable function
    if AvaProfile[2, -1] > AvaProfile[2, 0]:
        print('[com2AB] Profile reversed')
        L = AvaProfile[3, -1]
        AvaProfile = np.fliplr(AvaProfile)
        AvaProfile[3, :] = L - AvaProfile[3, :]
        indSplit = np.shape(AvaProfile)[1]-indSplit

    # END TODO


    # TODO: breackout to testable function
    # TODO: more descriptiv: what is s? x,y,z are clear
    s = AvaProfile[3, :]
    x = AvaProfile[0, :]
    y = AvaProfile[1, :]
    z = AvaProfile[2, :]
    ds = np.abs(s - np.roll(s, 1))
    dz = np.abs(z - np.roll(z, 1))
    ds[0] = 0.0
    dz[0] = 0.0
    angle = np.rad2deg(np.arctan2(dz, ds))
    CuSplit = AvaProfile[3, indSplit]
    # TODO SPLIT POINT READING
    # get all values where Angle < 10 but >0
    # get index of first occurance and go one back to get previous value
    # (i.e. last value above 10 deg)
    # tmp = x[(angle < 10.0) & (angle > 0.0) & (x > 450)]
    tmp = np.where((angle < 10.0) & (angle > 0.0) & (s > CuSplit))
    ids_10Point = tmp[0][0] - 1

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

    # TODO END: breackout to testable function


    # TODO: the rest is plotting + analytical stuff
    # I (FSO) suggest we just write it out and move these functions to IO module?

    # Line down to alpha
    f = z[0] + np.tan(np.deg2rad(-alpha)) * s
    fplus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[0])) * s
    fminus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[1])) * s
    fminus2SD = z[0] + np.tan(np.deg2rad(-alphaSD[2])) * s

    # First it calculates f - g and the corresponding signs
    # using np.sign. Applying np.diff reveals all
    # the positions, where the sign changes (e.g. the lines cross).
    ids_alpha = np.argwhere(np.diff(np.sign(f - z))).flatten()
    ids_alphaP1SD = np.argwhere(np.diff(np.sign(fplus1SD - z))).flatten()
    ids_alphaM1SD = np.argwhere(np.diff(np.sign(fminus1SD - z))).flatten()
    ids_alphaM2SD = np.argwhere(np.diff(np.sign(fminus2SD - z))).flatten()

    # Only get the first index past the splitpoint
    try:
        ids_alpha = ids_alpha[s[ids_alpha] > CuSplit][0]
    except print('Alpha out of profile'):
        ids_alpha = None

    try:
        ids_alphaP1SD = ids_alphaP1SD[s[ids_alphaP1SD] > CuSplit][0]
    except print('+1 SD above beta point'):
        ids_alphaP1SD = None

    try:
        ids_alphaM1SD = ids_alphaM1SD[s[ids_alphaM1SD] > CuSplit][0]
    except print('-1 SD out of profile'):
        ids_alphaM1SD = None

    try:
        ids_alphaM2SD = ids_alphaM2SD[s[ids_alphaM2SD] > CuSplit][0]
    except print('-2 SD out of profile'):
        ids_alphaM2SD = None

    # Plot raster and path
    fig1, ax1 = plt.subplots()
    cmap = mpl.cm.Greys
    cmap.set_bad(color='white')
    im1 = plt.imshow(rasterdata, cmap, origin='lower')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cb1 = fig1.colorbar(im1, cax=cax)
    path1 = ax1.plot((x-header.xllcorner)/header.cellsize,
                     (y-header.yllcorner)/header.cellsize)
    ax1.plot((AvaProfile[0]-header.xllcorner)/header.cellsize,
             (AvaProfile[1]-header.yllcorner)/header.cellsize, 'k')
    ax1.plot((splitPoint[0]-header.xllcorner)/header.cellsize,
             (splitPoint[1]-header.yllcorner)/header.cellsize, '.',
             color='0.6', label='Split point')
    ax1.plot((SplitPoint[0]-header.xllcorner)/header.cellsize,
             (SplitPoint[1]-header.yllcorner)/header.cellsize, '.',
             color='0.3', label='Projection of Split Point on ava path')
    plt.show()

    # Plot the whole profile with beta, alpha ... points and lines
    plotSaveResults(beta, alpha, x, y, s, z, f, poly, indSplit, ids_10Point,
                    ids_alpha, ids_alphaM1SD, ids_alphaM2SD, abVersion,
                    ParameterSet, saveOutPath)

    print('Parameter Set %s\n' % ParameterSet)
    print('Alpha point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_alpha], y[ids_alpha], z[ids_alpha], s[ids_alpha], alpha))
    print('Beta point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_10Point], y[ids_10Point], z[ids_10Point], s[ids_10Point], beta))
    print('alphaM1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_alphaM1SD], y[ids_alphaM1SD], z[ids_alphaM1SD], s[ids_alphaM1SD], alphaSD[1]))
    print('alphaM2SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_alphaM2SD], y[ids_alphaM2SD], z[ids_alphaM2SD], s[ids_alphaM2SD], alphaSD[2]))
    print('alphaP1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle  in [°] : %.2f\n' % (
        x[ids_alphaP1SD], y[ids_alphaP1SD], z[ids_alphaP1SD], s[ids_alphaP1SD], alphaSD[0]))
    WriteResults(beta, alpha, x, y, s, z, ids_10Point, ids_alpha, ids_alphaM1SD,
                 ids_alphaM2SD, ids_alphaP1SD, alphaSD, abVersion, ParameterSet, saveOutPath)


def plotSaveResults(beta, alpha, x, y, s, z, f, poly, indSplit,
                    ids_10Point, ids_alpha, ids_alphaM1SD,
                    ids_alphaM2SD, abVersion, ParameterSet, saveOutPath):
    """TODO: - doc
    - split into save and plot"""

    # Plot the whole profile with beta, alpha ... points and lines
    plt.close("all")
    fig = plt.figure(4, figsize=(10, 6))
    titleText = 'Profile'
    plt.title(titleText)

    xlabelText = 'Distance [m]\nBeta: '+str(round(beta, 1)) + '$^\circ$' + \
        '  Alpha: '+str(round(alpha, 1)) + '$^\circ$'
    plt.xlabel(xlabelText, multialignment='center')

    plt.ylabel('Height [m]')

    plt.plot(s, z, '-', label='DEM')
    plt.plot(s, poly(s), ':', label='QuadFit')
    plt.axvline(x=s[indSplit], color='0.7',
                linewidth=1, linestyle='--', label='Split point')
    plt.axvline(x=s[ids_10Point], color='0.8',
                linewidth=1, linestyle='-.', label='Beta')

    plt.plot(s, f, '-', label='AlphaLine')
    if ids_alpha:
        plt.plot(s[ids_alphaM1SD], z[ids_alphaM1SD], 'x', markersize=8,
                 label='Alpha - 1SD')
    if ids_alphaM2SD:
        plt.plot(s[ids_alphaM2SD], z[ids_alphaM2SD], 'x', markersize=8,
                 label='Alpha - 2SD')

    ax = plt.gca()

    versionText = datetime.datetime.now().strftime("%d.%m.%y") + \
        '; ' + 'AlphaBeta ' + abVersion + ' ' + ParameterSet
    plt.text(0, 0, versionText, fontsize=8, verticalalignment='bottom',
             horizontalalignment='left', transform=ax.transAxes,
             color='0.5')
    # plt.text(-0.2, 0, 'matplotlib -2', \
    #          verticalalignment='center', transform=ax.transAxes)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(linestyle=':', color='0.9')
    plt.legend(frameon=False)
    plt.draw()
    # plt.show()
    print('[com2AB] Saving profile figure to:', saveOutPath)
    save_file = os.path.join(saveOutPath, 'AlphaBeta_.pdf')
    plt.savefig(save_file)
    plt.close(fig)
    plt.close("all")

def WriteResults(beta, alpha, x, y, s, z, ids_10Point, ids_alpha,
                 ids_alphaM1SD, ids_alphaM2SD, ids_alphaP1SD,
                 alphaSD, abVersion, ParameterSet, saveOutPath):
    """ Write com2AB results to file """

    FileName_ext = saveOutPath + 'results_python.txt'
    with open(FileName_ext, 'w') as outfile:
        outfile.write('Parameter Set %s\n' % ParameterSet)
        outfile.write('Alpha point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alpha], y[ids_alpha], z[ids_alpha], s[ids_alpha], alpha))
        outfile.write('Beta point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_10Point], y[ids_10Point], z[ids_10Point], s[ids_10Point], beta))
        outfile.write('alphaM1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaM1SD], y[ids_alphaM1SD], z[ids_alphaM1SD], s[ids_alphaM1SD], alphaSD[1]))
        outfile.write('alphaM2SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaM2SD], y[ids_alphaM2SD], z[ids_alphaM2SD], s[ids_alphaM2SD], alphaSD[2]))
        outfile.write('alphaP1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaP1SD], y[ids_alphaP1SD], z[ids_alphaP1SD], s[ids_alphaP1SD], alphaSD[0]))


def ProjectOnRaster(header, rasterdata, Points):
    """ Projects the points Points on Raster and returns the z coord
    Input :
    Points: list of points (x,y) 2 rows as many columns as Points
    Output:
    PointsZ: list of points (x,y,z) 3 rows as many columns as Points

    TODO: test
    """
    xllcorner = header.xllcorner
    yllcorner = header.yllcorner
    cellsize = header.cellsize
    ncol = header.ncols
    nrow = header.nrows
    xcoor = Points[0]
    ycoor = Points[1]
    zcoor = np.array([])
    for i in range(np.shape(xcoor)[0]):
        Lx = int(np.round((xcoor[i] - xllcorner) / cellsize))
        Ly = int(np.round((ycoor[i] - yllcorner) / cellsize))
        zcoor = np.append(zcoor, rasterdata[Ly][Lx])
    PointsZ = np.vstack((Points, zcoor))
    return (PointsZ)

def PrepareLine(header, rasterdata, avapath, splitPoint, distance=10):
    """ 1- Resample the avapath line with a max intervall of distance=10m
    between points (projected distance on the horizontal plane).
    2- Make avalanch profile out of the path (affect a z value using the DEM)
    3- Get projected split point on the profil (closest point)

    TODO: test
    """

    xcoor = avapath[0]
    ycoor = avapath[1]
    xcoornew = np.array([xcoor[0]])
    ycoornew = np.array([ycoor[0]])
    s = np.array([0])  # curvilinear coordinate
    # loop on the points of the avapath
    for i in range(np.shape(xcoor)[0]-1):
        Vx = xcoor[i+1] - xcoor[i]
        Vy = ycoor[i+1] - ycoor[i]
        D = np.sqrt(Vx**2 + Vy**2)
        nd = int(np.round(D / distance) + 1)
        # Resample each segment
        S0 = s[-1]
        for j in range(1, nd):
            xn = j / (nd - 1) * Vx + xcoor[i]
            yn = j / (nd - 1) * Vy + ycoor[i]
            xcoornew = np.append(xcoornew, xn)
            ycoornew = np.append(ycoornew, yn)
            s = np.append(s, S0 + D * j / nd)

    # test = np.transpose(np.array([[header.xllcorner,header.yllcorner],
    # [header.xllcorner+header.cellsize*(header.ncols-1),header.yllcorner],
    # [header.xllcorner,header.yllcorner+header.cellsize*(header.nrows-1)],
    # [header.xllcorner+header.cellsize*(header.ncols-1),header.yllcorner+
    # header.cellsize*(header.nrows-1)]]))
    # Test = ProjectOnRaster(header,rasterdata,test)
    ResampAvaPath = np.vstack((xcoornew, ycoornew))
    AvaProfile = ProjectOnRaster(header, rasterdata, ResampAvaPath)
    AvaProfile = np.vstack((AvaProfile, s))

    # find split point by computing the distance to the line
    dist = np.sqrt((xcoornew - splitPoint[0])**2 + (ycoornew - splitPoint[1])**2)
    indSplit = np.argmin(dist)
    SplitPoint = AvaProfile[:, indSplit]
    SplitPoint = np.append(SplitPoint, s[indSplit])

    return AvaProfile, SplitPoint, indSplit


def ReadRaster(fname):
    """ Read raster file from raster filename (.asc)"""

    print('[com2AB] Reading DEM :', fname)
    header = IOf.readASCheader(fname)
    rasterdata = IOf.readASCdata2numpyArray(fname, header)
    rasterdata[rasterdata == header.noDataValue] = np.NaN
    return [header, np.flipud(rasterdata)]


def ReadAvaPath(fname, header):
    """ Read avalanche path file from .txt or .xyz"""

    print('[com2AB] Reading avalanche path :', fname)
    avapath = IOf.readpathdata2numpyArray(fname)
    if np.shape(avapath)[0] < 2:
        raise ValueError('Ava path file should contain at least 2 columns')
    for i in range(np.shape(avapath)[1]):
        Lx = int(np.round((avapath[0, i] - header.xllcorner) / header.cellsize))
        Ly = int(np.round((avapath[1, i] - header.yllcorner) / header.cellsize))
        if (Ly < 0 or Ly > header.nrows or Lx < 0 or Lx > header.ncols):
            raise ValueError('Avalanche path exceeds DEM extend')
    return avapath


def ReadSplitPoint(fname, header):
    """ Read split point path file from .txt or .xyz"""

    print('[com2AB] Reading split point :', fname)
    splitPoint = IOf.readpathdata2numpyArray(fname)
    if np.shape(splitPoint)[0] < 2:
        raise ValueError('split point file should contain at least 2 columns')
    Lx = int(np.round((splitPoint[0] - header.xllcorner) / header.cellsize))
    Ly = int(np.round((splitPoint[1] - header.yllcorner) / header.cellsize))
    if (Ly < 0 or Ly > header.nrows or Lx < 0 or Lx > header.ncols):
        raise ValueError('Split point is not on the DEM')
    return splitPoint

