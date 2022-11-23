import numpy as np
import pathlib
import matplotlib.pyplot as plt
import logging
from matplotlib.animation import FuncAnimation, PillowWriter

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outQuickPlot as oQ
import avaframe.out3Plot.outAIMEC as oA
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.Tools import PlotTools

cfgMain = cfgUtils.getGeneralConfig()
cfgFlags = cfgMain['FLAGS']
# create local logger
log = logging.getLogger(__name__)


def plotTrackParticle(outDirData, particlesList, trackedPartProp, cfg, dem):
    """ Plot time series of tracked particles

    Parameters
    ----------
    outDirData: str
        path to output directory
    particlesList: list
        list or particles dictionaries
    trackedPartProp: dict
        dictionary with time series of the wanted properties for tracked
        particles
    cfg : dict
        configuration read from ini file
    dem: dict
        dem dictionary with normal information

    """
    cfgTrackPart = cfg['TRACKPARTICLES']
    radius = cfgTrackPart.getfloat('radius')
    centerList = cfgTrackPart['centerTrackPartPoint']
    centerList = centerList.split('|')
    center = {'x': np.array([float(centerList[0])]),
              'y': np.array([float(centerList[1])])}
    center, _ = geoTrans.projectOnRaster(dem, center, interp='bilinear')
    time = trackedPartProp['time']

    # do some ploting
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))

    # ax1 : particle trajectory
    #fig.suptitle('Tracked particles')
    ax1 = plt.subplot(121)
    #ax1 = addDem2Plot(ax1, dem, what='slope')

    simu_number = 0
    plotDict = PlotTools.demParticles(simu_number)
    # uncomment this to set the under value for discrete cmap transparent
    # cmap.set_under(alpha=0)
    xllcenter = plotDict['raster']['header']['xllcenter']
    yllcenter = plotDict['raster']['header']['yllcenter']
    rowsMinPlot = plotDict['rowsMin']*plotDict['cellSize']  + yllcenter
    rowsMaxPlot = (plotDict['rowsMax']+1)*plotDict['cellSize']  + yllcenter
    colsMinPlot = plotDict['colsMin']*plotDict['cellSize']  + xllcenter
    colsMaxPlot = (plotDict['colsMax']+1)*plotDict['cellSize']  + xllcenter

    #extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]
    extent = [0, colsMaxPlot-colsMinPlot, 0, rowsMaxPlot-rowsMinPlot]

    # add DEM hillshade with contour lines
    ls, CS = pU.addHillShadeContours(ax1, plotDict['demConstrained'], plotDict['cellSize'], extent)

    # add the particles
    ax1.plot(trackedPartProp['x']+xllcenter-colsMinPlot, trackedPartProp['y']+yllcenter-rowsMinPlot)
    #circle1 = plt.Circle((center['x'], center['y']), radius, color='r')
    #ax1.add_patch(circle1)

    ax1.set_xlabel("x [m]", fontsize=22)
    ax1.set_ylabel("y [m]", fontsize=22)
    ax1.tick_params(axis='both', labelsize=15)
    #ax1.set_title('Trajectory')

    # ax2 = plt.subplot(222)
    # ax2.plot(time, trackedPartProp['m'])
    # ax2.set_xlabel('t [s]')
    # ax2.set_ylabel('m [kg]')
    # ax2.set_title('Mass')

    ax2 = plt.subplot(122)
    velocity = DFAtls.norm(trackedPartProp['ux'], trackedPartProp['uy'],
                           trackedPartProp['uz'])
    ax2.plot(time, velocity)

    ax2.set_xlabel("Time [s]\n\n", fontsize=22)
    ax2.set_ylabel("Velocity [m/s]", fontsize=22)
    ax2.tick_params(axis='both', labelsize=15)
    # add AvaNode legend
    #ax2.legend(loc='lower right', ncol=1, fancybox=False, shadow=False)
    fig.legend(loc='lower center', ncol=3, fancybox=True, shadow=True, fontsize=14)

    # ax4 = plt.subplot(224)
    # ax4.plot(time, trackedPartProp['h'])
    # ax4.set_xlabel('t [s]')
    # ax4.set_ylabel('h [m]')
    # ax4.set_title('Flow thickness')

    pathDict = {}
    pathDict['pathResult'] = outDirData
    outFileName = 'trackedParticles'
    pU.saveAndOrPlot(pathDict, outFileName, fig)

    if cfgFlags.getboolean('showPlot'):
        fig2 = plt.figure()
        ax1 = plt.subplot(111)
        for count in range(len(particlesList)):
            particles = particlesList[count]
            ax1 = updateTrackPart(particles, ax1, dem)
            ax1 = addDem2Plot(ax1, dem, what='slope')
            plt.pause(0.1)
        plt.show()

        # ani = FuncAnimation(fig2, update, round(len(Particles)),
        #                     fargs=(Particles, xllc, yllc, ax1, XX, YY, dem))
        # # plt.show()
        #
        # writer = PillowWriter(fps=4)
        # # ani.save("MalSecRel.gif", writer=writer)
        # ani.save("testTrackAlr1.gif", writer=writer)


def plotTrackParticleAcceleration(outDirData,trackedPartProp, cfg):
    """ Plot time series of tracked particles
    Parameters
    ----------
    outDirData: str
        path to output directory
    trackedPartProp: dict
        dictionary with time series of the wanted properties for tracked
        particles
    cfg : dict
        configuration read from ini file

    """


    # do some ploting
    fig = plt.figure(figsize=(pU.figW, pU.figH))
    fig.suptitle('Tracked particles acceleration')
    ax1 = plt.subplot(111)

    ax1.plot(trackedPartProp['time'], trackedPartProp['uAcc'])
    ax1.set_xlabel('time step [s]')
    ax1.set_ylabel('acceleration [ms-2]')

    pathDict = {}
    pathDict['pathResult'] = outDirData
    outFileName = 'trackedParticles_acceleration'
    pU.saveAndOrPlot(pathDict, outFileName, fig)


def updateTrackPart(particles, ax, dem):
    """Update axes with particles (tracked particles are highlighted in red)
    """

    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']

    X = particles['x'] + xllc
    Y = particles['y'] + yllc

    ax.clear()
    ax.set_title('t=%.2f s' % particles['t'])
    variable = particles['trackedParticles']
    # set range and steps of colormap
    cc = np.where(variable == 1, True, False)
    ax.scatter(X, Y, c='b', cmap=None, marker='.')
    ax.scatter(X[cc], Y[cc], c='r', cmap=None, marker='.', s=5)
    return ax


def addParticles2Plot(particles, ax, dem, whatS='m', whatC='h', colBarResType=''):
    """Update axes with particles

    Parameters
    ----------
    particles: dict
        particles dictionary
    ax: matplotlib ax object
    dem: dict
        dem dictionary with normal information
    whatS: str
        which particle property should be used for the markersize
    whatC: str
        which particle property should be used for the marker color

    """

    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']

    X = particles['x'] + xllc
    Y = particles['y'] + yllc
    cmap = pU.cmapT
    if colBarResType != '':
        unit = pU.cfgPlotUtils['unit%s' % colBarResType]
    else:
        unit = ''
    variableC = particles[whatC]
    variableS = particles[whatS]
    minMax = np.nanmax(variableS)-np.nanmin(variableS)
    if minMax>0:
        variableS = ((variableS-np.nanmin(variableS))/(np.nanmax(variableS)-np.nanmin(variableS)) + 1)*pU.ms
    else:
        variableS = pU.ms
    cmap, _, ticks, norm = pU.makeColorMap(cmap, np.amin(variableC), np.amax(variableC), continuous=pU.contCmap)
    # set range and steps of colormap
    sc = ax.scatter(X, Y, c=variableC, s=variableS, cmap=cmap, marker='.', zorder=15)
    cb = pU.addColorBar(sc, ax, ticks, unit)
    return ax, cb


def addDem2Plot(ax, dem, what='slope', extent=''):
    """ Add dem to the background of a plot

    Parameters
    ----------
    ax: matplotlib ax object
    dem: dict
        dem dictionary with normal information
    what: str
        what information about the dem will be plotted?
        slope: use the dem slope (computed from the normals) to color the plot
        z : use the elevation to color the plot
    """
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    xArray = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    yArray = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    cmap = pU.cmapGreys
    cmap.set_bad(color='white')

    if what == 'slope':
        value = dem['Nz'] / DFAtls.norm(dem['Nx'], dem['Ny'], dem['Nz'])
    elif what == 'z':
        value = dem['rasterData']
    elif what == 'hillshade':
        ls = pU.LightSource(azdeg=pU.azimuthDegree, altdeg=pU.elevationDegree)
        value = dem['rasterData']
        value = ls.hillshade(value, vert_exag=pU.vertExag, dx=value.shape[1],
                             dy=value.shape[0])
    else:
        value = dem['rasterData']
    if extent == '':
        extent = [xArray.min(), xArray.max(), yArray.min(), yArray.max()]

    ref0, im = pU.NonUnifIm(ax, xArray, yArray, value, 'x [m]', 'y [m]',
                            # extent=[2400, 2700, YY.min(), YY.max()],
                            extent=extent,
                            cmap=cmap, norm=None, zorder=0)
    CS = ax.contour(xArray, yArray, dem['rasterData'], levels=pU.hillshadeContLevs, colors='k', linewidths=0.5,
                    zorder=3)

    # add labels
    ax.clabel(CS, CS.levels, inline=True, fontsize=8, zorder=3)

    # add info box with indication of label meaning
    pU.putInfoBox(ax, '- elevation [m]', location='upperLeft', color='gray',
                  hAlignment='left', alphaF=1.0)
    return ax


def plotParticles(particlesList, cfg, dem):
    """ Plot particles on dem

    Parameters
    ----------
    particlesList: list
        list or particles dictionaries
    cfg : dict
        configuration read from ini file
    dem: dict
        dem dictionary with normal information

    """
    if cfgFlags.getboolean('showPlot'):
        for count in range(len(particlesList)):
            fig2 = plt.figure()
            ax1 = plt.subplot(111)
            ax1.clear()
            particles = particlesList[count]
            ax1.set_title('t=%.2f s' % particles['t'])
            ax1 = addParticles2Plot(particles, ax1, dem, whatS='h', whatC='m')
            ax1 = addDem2Plot(ax1, dem, what='slope')
            plt.show()


def addResult2Plot(ax, header, rasterData, resType, colorbar=True, contour=False):
    """ Add raster data to a plot

    Parameters
    ----------
    ax: matplotlib ax object
    header: dict
        raster header dictionary
    rasterData: 2D numpy array
        data to plot
    resType: str
        what kinf of result is it? ppr, pft...
    colorbar: bool
        If true add the colorbar
    """
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    rowsMin, rowsMax, colsMin, colsMax, rasterData = pU.constrainPlotsToData(rasterData, csz, extentOption=False, constrainedData=True)
    xArray = np.linspace(xllc+colsMin*csz, xllc+colsMax*csz, colsMax-colsMin+1)
    yArray = np.linspace(yllc+rowsMin*csz, yllc+rowsMax*csz, rowsMax-rowsMin+1)
    extent = [xArray.min(), xArray.max(), yArray.min(), yArray.max()]
    unit = pU.cfgPlotUtils['unit%s' % resType]
    contourLevels = pU.cfgPlotUtils['contourLevels%s' % resType]
    contourLevels = fU.splitIniValueToArraySteps(contourLevels)
    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.amin(rasterData), np.amax(rasterData),
                                           continuous=pU.contCmap)
    cmap.set_bad(alpha=0)
    rasterData = np.ma.masked_where(rasterData == 0, rasterData)

    ref0, im = pU.NonUnifIm(ax, xArray, yArray, rasterData, 'x [m]', 'y [m]',
                            extent=extent,
                            cmap=cmap, norm=norm, zorder=9)
    if colorbar:
        cb = pU.addColorBar(im, ax, ticks, unit)
    else:
        cb = None
    if contour:
        CS = ax.contour(xArray, yArray, rasterData, levels=contourLevels, colors='k')
        ax.clabel(CS, inline=1, fontsize=8, zorder=10)
    else:
        CS = None
    return ax, extent, cb, CS


def createContourPlot(reportDictList, avalancheDir, simDF):
    """ create a contour line plot of all simulations of current run

        Parameters
        -----------
        reportDictList: list
            list of com1DFA dictionary with info on contour dicts for each sim
        avalancheDir: str or pathlib path
            path to avalanche directory
        simDF: pandas DataFrame
            dataframe with one row per simulation performed and its parameter settings

        Returns
        -------
        reportDictList: list
            updated reportDictList - deleted contours dict
    """

    modName = 'com1DFA'
    contourD = {}

    # fetch coordinates of contour line for each sim in reportDict and create contourD
    for cont in reportDictList:
        contourD[list(cont['contours'].keys())[0]] = cont['contours'][list(cont['contours'].keys())[0]]
        del cont['contours']

    # create contours directory in Ouputs to save plot
    contourPlotDir = pathlib.Path(avalancheDir, 'Outputs', modName, 'contours')
    fU.makeADir(contourPlotDir)

    # create pathDict and fetch simName of first sim found
    pathDict = {'pathResult': contourPlotDir, 'plotScenario': pathlib.Path(avalancheDir).stem,
        'avaDir': avalancheDir}
    simName1 = simDF['simName'].iloc[0]
    # check if consistent settings throughout all sims
    if len(np.unique(simDF['contourResType'])) != 1 or len(np.unique(simDF['thresholdValue'])) != 1:
        log.warning('Contou result type or thresholdValue are not identical for all sims performed - so cannot create contour plot')
        return reportDictList, False

    # generate plot
    oQ.plotContours(contourD,
                    contourD[simName1]['contourResType'],
                    contourD[simName1]['thresholdValue'],
                    pathDict)

    log.info('Saved contour plot of %d sims to %s' % (len(simDF), contourPlotDir))
    log.info('Sim names are: %s' % simDF['simName'].to_list())

    return reportDictList, True


def fetchContCoors(demHeader, flowF, cfgVisu, simName):
    """ fetch coordinates of contour line

        Parameters
        ------------
        demHeader: dict
            dictionary of dem nrows, ncols, cellsize
        flowF: np.ndarray
            field data used to compute contour line
        cfgVisu: configparser object
            configuration settings for visualisation, here used:
            contourResType, thresholdValue
        simName: str
            simName

        Rertuns
        --------
        contDictXY: dict
            dictionary with simName and subDict with x, y coordinates of contour line
            contourResType and thresholdValue
    """
    # create coordinate grid first
    xGrid, yGrid, _, _ = geoTrans.makeCoordGridFromHeader(demHeader)

    # fetch contour line
    x, y = pU.fetchContourCoords(xGrid, yGrid, flowF,cfgVisu.getfloat('thresholdValue'))

    # setup dict
    contDictXY = {simName: {'x': x, 'y': y, 'contourResType': cfgVisu['contourResType'],
        'thresholdValue': cfgVisu.getfloat('thresholdValue')}}

    return contDictXY
