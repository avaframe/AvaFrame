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


def plotTrackParticle(outDirData, particlesList, trackedPartProp, cfg, dem):
    """ Plot time series of tracked partcles
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
    fig.suptitle('Tracked particles')
    ax1 = plt.subplot(221)
    ax1 = addDem2Plot(ax1, dem, what='slope')
    circle1 = plt.Circle((center['x'], center['y']), radius, color='r')
    ax1.plot(trackedPartProp['x'], trackedPartProp['y'])
    ax1.add_patch(circle1)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.set_title('Trajectory')

    ax2 = plt.subplot(222)
    ax2.plot(time, trackedPartProp['m'])
    ax2.set_xlabel('t [s]')
    ax2.set_ylabel('m [kg]')
    ax2.set_title('Mass')

    ax3 = plt.subplot(223)
    velocity = DFAtls.norm(trackedPartProp['ux'], trackedPartProp['uy'],
                           trackedPartProp['uz'])
    ax3.plot(time, velocity)
    ax3.set_xlabel('t [s]')
    ax3.set_ylabel('v [m/s]')
    ax3.set_title('Velocity')

    ax4 = plt.subplot(224)
    ax4.plot(time, trackedPartProp['h'])
    ax4.set_xlabel('t [s]')
    ax4.set_ylabel('h [m]')
    ax4.set_title('Flow depth')

    pathDict = {}
    pathDict['pathResult'] = outDirData
    outFileName = 'trackedParticles'
    pU.saveAndOrPlot(pathDict, outFileName, fig)

    if cfgFlags.getboolean('showPlot'):
        fig2 = plt.figure()
        ax1 = plt.subplot(111)
        for count in range(len(particlesList)):
            particles = particlesList[count]
            ax1 = update(particles, ax1, dem)
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


def update(particles, ax, dem):
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


def addDem2Plot(ax, dem, what='slope'):
    """ Add dem to the background of a plot

    Parameters
    ----------
    ax: mathplotlix ax object
    dem: dict
        dem dictionary with normal information
    what: str
        what information abour the dem will be plotted?
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
    else:
        value = dem['rasterData']

    ref0, im = pU.NonUnifIm(ax, xArray, yArray, value, 'x [m]', 'y [m]',
                            # extent=[2400, 2700, YY.min(), YY.max()],
                            extent=[xArray.min(), xArray.max(),
                                    yArray.min(), yArray.max()],
                            cmap=cmap, norm=None)
    ax.contour(xArray, yArray, dem['rasterData'], levels=10, colors='k')
    return ax
