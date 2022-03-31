""" functions to plot range-time diagrams """

import numpy as np
from matplotlib import pyplot as plt
import logging
import pathlib
from cmcrameri import cm
from matplotlib.colors import LightSource

# Local imports

import avaframe.out3Plot.plotUtils as pU
import avaframe.ana5Utils.distanceTimeAnalysis as dtAna
import avaframe.in3Utils.geoTrans as gT

log = logging.getLogger(__name__)


def plotRangeTime(mtiInfo, cfgRangeTime):
    """ plot range-time diagram with avalanche front and colorcoded average values of result parameter

        Parameters
        -----------
        mtiInfo: dict
            dictionary with average values for range to reference point (mti), timeStep list,
            list with distance to reference point of avalanche front (rangeList)
        cfgRangeTime: configparser object
            configuration settings for range time diagram - here used avalancheDir, rangeTimeResType,
            simHash

    """

    # fetch required input info
    mti =  mtiInfo['mti']
    rangeGates = mtiInfo['rangeGates']
    timeList = mtiInfo['timeList']
    rangeList = mtiInfo['rangeList']
    rangeTimeResType = cfgRangeTime['rangeTimeResType']
    maxVel, rangeVel, timeVel = dtAna.approachVelocity(mtiInfo, cfgRangeTime.getfloat('minVelTimeStep'))

    # in case time steps are not ordered - the colormesh x and y need to be ordered
    timeIndex = np.argsort(np.array(timeList))
    timeListNew = np.array(timeList)[timeIndex]
    mti = mti[:, timeIndex]

    # create plot
    fig = plt.figure(figsize=(pU.figW, pU.figH))
    ax = fig.add_subplot(1, 1, 1)
    plt.title(mtiInfo['plotTitle'])
    pc = plt.pcolormesh(timeListNew, rangeGates, mti, cmap=pU.cmapRangeTime)
    plt.plot(timeList, rangeList , '.', color='black', markersize=4,
        label='avalanche front')
    plt.xlabel('Time [s]')
    plt.ylabel('Distance to %s [m]' % mtiInfo['referencePointName'])

    # add colorbar and infobox
    unit = pU.cfgPlotUtils['unit' + rangeTimeResType]
    cName = 'avg. ' + pU.cfgPlotUtils['name' + rangeTimeResType]
    pU.addColorBar(pc, ax, None, unit, title=cName)
    pU.putAvaNameOnPlot(ax, cfgRangeTime['avalancheDir'])
    rangeTimeVelocityLegend(timeListNew, rangeGates, ax, maxVel, width=0.25, height=0.25, lw=0.5, textsize=7)
    # add max velocity location
    ax.plot(timeVel, rangeVel, 'r*', label='max velocity location')

    # add info on avalanche front in legend
    plt.legend(facecolor='grey', framealpha=0.2, loc='lower right', fontsize=8)

    # if tt-diagram add beta point info
    if mtiInfo['type'] == 'thalwegTime':
        # invert y axis as ava flow starts from minus distance to beta point
        ax.invert_yaxis()
        ax.axhline(y=0.0, color='gray', linestyle='--', linewidth=1, alpha=0.5,
            label='beta point: %.1f°' % mtiInfo['betaPointAngle'])

    # set path for saving figure
    outDir = pathlib.Path(cfgRangeTime['avalancheDir'], 'Outputs', 'ana5Utils')
    outFileName = mtiInfo['type'] + '_' + rangeTimeResType + '_' + cfgRangeTime['simHash']
    plotPath = pU.saveAndOrPlot({'pathResult': outDir}, outFileName, fig)


def radarFieldOfViewPlot(radarFov, radarRange, cfgRangeTime, rangeGates, dem):
    """ Create radar field of view plot

        Parameters
        -----------
        radarFov: numpy array
            list with radar location and end point of field of view, x and y coors
        radarRange: masked array
            masked array of DEM with radar field of view - showing distance to radar
        cfgRangeTime: configparser object
            configuration settings section - here used avalancheDir, simHash,
                aperture angle [degree]
        rangeGates: numpy array
            range gates of radar field of view
        dem: dict
            dictionary with dem header and data
    """

    # get input parameters
    aperture = cfgRangeTime['GENERAL'].getfloat('aperture')

    # fetch header info - required for creating coordinate grid
    xllc = dem['header']['xllcenter']
    yllc = dem['header']['yllcenter']
    ncols = dem['header']['ncols']
    nrows = dem['header']['nrows']
    cellSize = dem['header']['cellsize']

    # Set coordinate grid with given origin
    X, Y = gT.makeCoordinateGrid(xllc, yllc, cellSize, ncols, nrows)

    # load required input parameters for contour plot
    gateContours = cfgRangeTime['PLOTS'].getint('gateContours')

    # get field of view with radar location and aperture angle
    xR = gT.rotate(radarFov, aperture)
    xL = gT.rotate(radarFov, -aperture)
    radarFovArrow = np.diff(radarFov)*0.2

    # create plot of radar range values and field of view
    fig = plt.figure(figsize=(pU.figW, pU.figH))
    ax1 = plt.subplot(111)
    cmapRadar = pU.cmapRadarFOV
    cmapRadar.set_bad('w', alpha=0.)
    pc = plt.pcolormesh(X, Y, radarRange, cmap=cmapRadar)
    ax1.plot(radarFov[0], radarFov[1], marker='o', color='k', zorder=1)
    ax1.contour(X, Y, radarRange, rangeGates[::gateContours], colors=['brown'], alpha=0.5,
        linewidths=0.5)
    ax1.plot(xR[0], xR[1], marker='o', color='k')
    ax1.plot(xL[0], xL[1], marker='o', color='k')
    ax1.arrow(x=radarFov[0][0], y=radarFov[1][0], dx=radarFovArrow[0][0], dy=radarFovArrow[1][0],
        fc='r', ec='r', width=10, zorder=2)
    ax1.set_title('Radar range')
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    plt.xticks(rotation=45)
    # add colorbar and infoboxes
    unit = 'm'
    cName = 'range'
    pU.addColorBar(pc, ax1, None, unit, title=cName)
    pU.putAvaNameOnPlot(ax1, cfgRangeTime['GENERAL']['avalancheDir'])
    pU.putInfoBox(ax1, '- range gates [%d]' % gateContours, location='lowerRight', color='brown')

    # create plot of DEM and radar field of view
    fig2 = plt.figure(figsize=(pU.figW, pU.figH))
    # add second panel with DEM and radar location
    ax2 = plt.subplot(111)
    #pc0 = plt.pcolormesh(X, Y, dem['rasterData'], cmap=cm.grayC)
    ls = LightSource(azdeg=315, altdeg=45)
    ax2.imshow(ls.hillshade(dem['rasterData'], vert_exag=10, dx=dem['rasterData'].shape[1],
        dy=dem['rasterData'].shape[0]), cmap='gray', extent=[np.amin(X), np.amax(X), np.amin(Y), np.amax(Y)],
        origin='lower')
    CS =  ax2.contour(X, Y, dem['rasterData'], colors=['brown'], alpha=0.75,
        linewidths=0.5)
    # add radar field of view
    ax2.plot(radarFov[0], radarFov[1], marker='o', color='k', zorder=1)
    ax2.plot(xR[0], xR[1], marker='o', color='k')
    ax2.plot(xL[0], xL[1], marker='o', color='k')
    ax2.set_title('DEM with radar location')
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    plt.xticks(rotation=45)
    # add infoboxes
    pU.putAvaNameOnPlot(ax2, cfgRangeTime['GENERAL']['avalancheDir'])
    pU.putInfoBox(ax2, '- elevation [m]', location='lowerRight', color='brown')
    ax2.clabel(CS, CS.levels, inline=True, fontsize=10)

    # set path for saving figure
    outDir = pathlib.Path(cfgRangeTime['GENERAL']['avalancheDir'], 'Outputs', 'ana5Utils')
    outFileName = 'radarFieldOfViewPlot_%s' % cfgRangeTime['GENERAL']['simHash']
    plotPath = pU.saveAndOrPlot({'pathResult': outDir}, outFileName, fig)
    outFileName = 'radarFieldOfViewPlotDEM_%s' % cfgRangeTime['GENERAL']['simHash']
    plotPath = pU.saveAndOrPlot({'pathResult': outDir}, outFileName, fig2)


def plotRangeRaster(slRaster, rasterTransfo, cfgRangeTime, rangeRaster, cLower):
    """ create plot of distance from start of runout area point extracted avalanche front and result field

        Parameters
        -----------
        slRaster: numpy array
            transformed result field
        rasterTransfo: dict
            info on coordinate transformation
        cfgRangeTime: configparser object
            configuration settings of range time diagram
        rangeRaster: np array
            distance to runout area point in s coordinate
        cLower: int
            index of avalanche front
    """

    # fetch unit for result parameter
    resType = cfgRangeTime['rangeTimeResType']
    unitRes = pU.cfgPlotUtils['unit' + resType]

    # fetch s. l coordinates
    l = rasterTransfo['l']
    s = rasterTransfo['s']
    indStartOfRunout = rasterTransfo['indStartOfRunout']

    # create figure
    fig = plt.figure(figsize=(pU.figW*1.5, pU.figH))
    ax1 = plt.subplot(111)
    plt.title('Distance to reference point and %s field' % pU.cfgPlotUtils['name' + resType])
    # create colormap for range field
    cmapRange, _, ticksRange, normRange = pU.makeColorMap(cm.grayC.reversed(), np.nanmin(
        rangeRaster), np.nanmax(rangeRaster), continuous=pU.contCmap)
    # add imshow plot of range
    ref0, im1 = pU.NonUnifIm(ax1, l, s, rangeRaster, 'l [m]', 's [m]',
                            extent=[l.min(), l.max(), s.min(), s.max()], cmap=cmapRange)
    # add horizontal line indicating location of start of runout area and avalanche front
    ax1.axhline(y=s[indStartOfRunout], color='w', linestyle='--',
                label='start of run-out area point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])
    ax1.axhline(y=s[cLower], color='y', linestyle='--',
                label='avalanche front based on %s limit %s' % (resType, cfgRangeTime['thresholdResult']))

    # add result field masked with threshold
    bmaskRes = slRaster < cfgRangeTime.getfloat('thresholdResult')
    maskResType = np.ma.array(slRaster, mask=bmaskRes)
    # create colormap
    cmapRes, _, ticksRes, normRes = pU.makeColorMap(pU.colorMaps[resType], np.nanmin(
        maskResType), np.nanmax(maskResType), continuous=True)
    cmapRes.set_bad(alpha=0.0)
    ref2, im2 = pU.NonUnifIm(ax1, l, s, maskResType, 'l [m]', 's [m]',
                            extent=[l.min(), l.max(), s.min(), s.max()], cmap=cmapRes)

    # add legend and colorbar
    ax1.legend(loc='center', facecolor='black', framealpha=0.2)
    pU.addColorBar(im1, ax1, ticksRange, 'm', pad=0.005)
    pU.addColorBar(im2, ax1, ticksRes, unitRes)

    # create plot name and location for saving and save/plot
    outDir = pathlib.Path(cfgRangeTime['avalancheDir'], 'Outputs', 'ana5Utils')
    pathDict = {'pathResult': outDir}
    outFileName = 'distanceToStartOfRunoutArea_%s_%s' % (resType, cfgRangeTime['simHash'])
    pU.saveAndOrPlot(pathDict, outFileName, fig)


def plotMaskForMTI(cfgRangeTime, bmaskRange, bmaskAvaRadar, bmaskAvaRadarRangeslice, mtiInfo):
    """ plot masks used to compute averages for range-time diagram with radar field of view

        Parameters
        -----------
        cfgRangeTime: configparser settings
            configuration settings for range time diagram, avalancheDir
        bmaskRange: numpy mask
            mask for radar range slice for range gate and rgWidth
        bmaskAvaRadar: numpy mask
            mask for result parameter below threshold and out of radar field of view
        bmaskAvaRadarRangeslice: numpy mask
            full mask combination
    """

    # get extent of masks
    headerInfo = mtiInfo['demOriginal']['header']
    x0 = headerInfo['xllcenter']
    y0 = headerInfo['yllcenter']
    x1 = headerInfo['xllcenter'] + headerInfo['cellsize'] * bmaskRange.shape[1]
    y1 = headerInfo['yllcenter'] + headerInfo['cellsize'] * bmaskRange.shape[0]
    # dimensions
    nx = bmaskRange.shape[1]
    ny = bmaskRange.shape[0]

    # create figure with three panels
    fig = plt.figure(figsize=(pU.figW*3, pU.figH))
    # add radar range gate mask
    ax1 = fig.add_subplot(131)
    ax1.set_title('Radar range gate mask')
    im1 = ax1.imshow(bmaskRange, extent=[x0, x1, y0, y1], origin='lower', aspect=nx/ny)
    # add masked avalanche result field with threshold and radar field of view
    ax2 = fig.add_subplot(132)
    ax2.set_title('RadarFOV and resultThreshold mask')
    im2 = ax2.imshow(bmaskAvaRadar, extent=[x0, x1, y0, y1], origin='lower', aspect=nx/ny)
    # add combined mask
    ax3 = fig.add_subplot(133)
    ax3.set_title('Combined mask')
    im3 = ax3.imshow(bmaskAvaRadarRangeslice, extent=[x0, x1, y0, y1], origin='lower', aspect=nx/ny)

    # create plot name and location for saving and save/plot
    outDir = pathlib.Path(cfgRangeTime['avalancheDir'], 'Outputs', 'ana5Utils')
    pathDict = {'pathResult': outDir}
    outFileName = 'rangeTimeMasks_%s' % cfgRangeTime['simHash']
    pU.saveAndOrPlot(pathDict, outFileName, fig)


def rangeTimeVelocityLegend(timeSteps, rangeGates, ax, maxVel, width=0.25, height=0.25, lw=0.5, textsize=7):
    """ set legend in range time diagram for velocity in terms of steepness of front position

        Parameters
        -----------
        timeSteps: list
            time steps - needs to be ordered
        rangeGates: list
            list of range gates - needs to be ordered
        ax: matplotlib axis object
            axes where the legend shall be accomodated to
        width: float
            width of legend in fraction of full plot width
        height: float
            height of legend in fraction of full plot height
        lw: float
            line width
        textsize: float
            size of text

    """

    # compute slope of max velocity line for transformed axes
    xMax = timeSteps[-1]
    yMax = np.amax(rangeGates)# - np.amin(rangeGates)
    xinterval = xMax - timeSteps[0]
    yinterval = np.amax(rangeGates) - np.amin(rangeGates)
    kPlot = (yinterval /xinterval)

    # fetch points for drawing velocity legend
    point = getVelocityPoints(1, 1, width, height, 1, 1, maxVel, kPlot)

    # add title
    ax.text(point[1][0], point[1][1], 'Approach \n velocity [m/s]', size=6,
        transform=ax.transAxes, va='top', ha='right', fontweight='bold',
        bbox=dict(facecolor='white', alpha=0.5))

    # plot boundary of velocity box
    ax.plot([point[0][0], point[1][0]], [point[0][1], point[1][1]], color='k', linewidth=lw, transform=ax.transAxes)
    ax.plot([point[1][0], point[5][0]], [point[1][1], point[5][1]], color='k', linewidth=lw, transform=ax.transAxes)
    ax.plot([point[5][0], point[10][0]], [point[5][1], point[10][1]], color='k', linewidth=lw, transform=ax.transAxes)
    ax.plot([point[10][0], point[0][0]], [point[10][1], point[0][1]], color='k', linewidth=lw, transform=ax.transAxes)

    # these are the diagonal lines referring to velocity
    for j in range(2, len(point) - 1):
        x1 = [point[0][0], point[j][0]]
        y1 = [point[0][1], point[j][1]]
        if j == 9: # max velocity make red
            ax.plot(x1, y1, color='r', linewidth=1.0, transform=ax.transAxes)
        else:
            ax.plot(x1, y1, color='k', linewidth=0.5, transform=ax.transAxes)

    # compute velocity values
    dataPoint = getVelocityPoints(xMax, yMax, width, height, xinterval, yinterval, maxVel, kPlot)
    velocity = []
    for i in range(2, len(dataPoint)-1):
        deltaS = dataPoint[0][1] - dataPoint[i][1]
        deltaT = dataPoint[i][0] - dataPoint[0][0]
        velocity.append(deltaS/deltaT)

    # add velocity labels
    for j in range(2, len(point)-1):
        x1 = [point[0][0], point[j][0]]
        y1 = [point[0][1], point[j][1]]

        if j < 5: # text along vertical axes
            ax.text(x1[1]+0.01, y1[1]+0.002, '%.1f' % velocity[j-2], size=5, fontweight='bold',
                color='k', transform=ax.transAxes, gid=j, va='center',ha='left')
        elif j == 9: # text for max velocity
            if velocity[j-2] < velocity[3]:
                xLoc = x1[1]-0.04
                yLoc =  y1[1]-0.004
            else:
                xLoc = x1[1]
                yLoc =  y1[1]+0.004
            ax.text(xLoc, yLoc, 'max %.1f' % velocity[j-2], size=5, fontweight='bold',
                color='r', transform=ax.transAxes, gid=j, va='center',ha='left',
                bbox=dict(facecolor='white', alpha=0.5))
            if abs(velocity[j-2] - maxVel) > 1.e7:
                log.warning('velocity legend wrong')

        else: # text along horizontal axes
            ax.text(x1[1], y1[1]-0.01, '%.1f' % velocity[j-2], size=5, fontweight='bold',
                color='k', transform=ax.transAxes, gid=j, va='top', ha='center')


def getVelocityPoints(xMax, yMax, width, height, xinterval, yinterval, maxVel, kPlot):
    """ get points for legend creation (box and sloping lines)

        Parameters
        -----------
        xMax: float
            max x extent of data (timeSteps)
        yMax: float
            max y extent of data (rangeGates)

        Returns
        --------
        points: list
            list of point coordinates
    """

    points = []
    points.append([xMax - width*xinterval, yMax])
    points.append([xMax, yMax])
    points.append([xMax, yMax - height/4.*yinterval])
    points.append([xMax, yMax - height/2.*yinterval])
    points.append([xMax, yMax - height/4.*3.*yinterval])
    points.append([xMax, yMax - height*yinterval])
    points.append([xMax - width/4*xinterval, yMax - height*yinterval])
    points.append([xMax - width/2*xinterval, yMax - height*yinterval])
    points.append([xMax - width/4*3*xinterval, yMax - height*yinterval])

    if maxVel > kPlot:
        points.append([(xMax - width*xinterval) + width*xinterval/maxVel*kPlot, yMax-height*yinterval])
    else:
        points.append([xMax, yMax - height*yinterval*maxVel/kPlot])

    points.append([xMax - width*xinterval, yMax - height*yinterval])

    return points
