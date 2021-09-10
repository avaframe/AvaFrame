import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.plotUtils as pU


cfgMain = cfgUtils.getGeneralConfig()
cfgFlags = cfgMain['FLAGS']


def plotTrackParticle(outDirData, Particles, trackedPartProp, cfg, demOri, dem):
    """ Plot time series of tracked partcles"""
    cfgTrackPart = cfg['TRACKPARTICLES']
    radius = cfgTrackPart.getfloat('radius')
    centerList = cfgTrackPart['centerTrackPartPoint']
    centerList = centerList.split('|')
    center = {'x': np.array([float(centerList[0])]),
              'y': np.array([float(centerList[1])])}
    center, _ = geoTrans.projectOnRaster(dem, center, interp='bilinear')
    time = trackedPartProp['time']

    # do some ploting
    # ToDo: put this in a plotting folder (here just to demonstrate how this works)
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))
    # fig.suptitle('This is a somewhat long figure title')
    ax1 = plt.subplot(221)
    ax1 = addDem2Plot(ax1, dem, what='slope')
    circle1 = plt.Circle((center['x'], center['y']), radius, color='r')
    ax1.plot(trackedPartProp['x'], trackedPartProp['y'])
    ax1.add_patch(circle1)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.set_title('Tracked particles trajectory')

    ax2 = plt.subplot(222)
    ax2.plot(time, trackedPartProp['m'])
    ax2.set_xlabel('t [s]')
    ax2.set_ylabel('m [kg]')
    ax2.set_title('Tracked particles mass')

    ax3 = plt.subplot(223)
    velocity = DFAtls.norm(trackedPartProp['ux'], trackedPartProp['uy'],
                           trackedPartProp['uz'])
    ax3.plot(time, velocity)
    ax3.set_xlabel('t [s]')
    ax3.set_ylabel('v [m/s]')
    ax3.set_title('Tracked particles velocity')

    ax4 = plt.subplot(224)
    ax4.plot(time, trackedPartProp['h'])
    ax4.set_xlabel('t [s]')
    ax4.set_ylabel('h [m]')
    ax4.set_title('Tracked particles flow depth')

    pathDict = {}
    pathDict['pathResult'] = outDirData
    outFileName = 'trackedParticles'
    pU.saveAndOrPlot(pathDict, outFileName, fig)

    if cfgFlags.getboolean('showPlot'):
        fig2 = plt.figure()
        ax1 = plt.subplot(111)
        for count in range(len(Particles)):
            update(count, Particles, ax1, dem)
        plt.show()

        # ani = FuncAnimation(fig2, update, round(len(Particles)),
        #                     fargs=(Particles, xllc, yllc, ax1, XX, YY, dem))
        # # plt.show()
        #
        # writer = PillowWriter(fps=4)
        # # ani.save("MalSecRel.gif", writer=writer)
        # ani.save("testTrackAlr1.gif", writer=writer)


def update(count, Particles, ax, dem):
    particles = Particles[count]

    header = dem['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']

    X = particles['x'] + xllc
    Y = particles['y'] + yllc

    ax.clear()
    ax.set_title('t=%.2f s' % particles['t'])
    variable = particles['trackedParticles']
    ax = addDem2Plot(ax, dem, what='slope')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapPres, np.amin(variable),
                                           np.amax(variable), continuous=True)
    # set range and steps of colormap
    cc = np.where(variable == 1, True, False)
    ax.scatter(X, Y, c='b', cmap=None, marker='.')
    ax.scatter(X[cc], Y[cc], c='r', cmap=None, marker='.', s=5)

    plt.pause(0.1)


def addDem2Plot(ax, dem, what='slope'):
    """ Add dem to the background of a plot"""
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
    else:
        value = dem['rasterData']

    ref0, im = pU.NonUnifIm(ax, xArray, yArray, value, 'x [m]', 'y [m]',
                            # extent=[2400, 2700, YY.min(), YY.max()],
                            extent=[xArray.min(), xArray.max(),
                                    yArray.min(), yArray.max()],
                            cmap=cmap, norm=None)
    ax.contour(xArray, yArray, dem['rasterData'], levels=10, colors='k')
    return ax
