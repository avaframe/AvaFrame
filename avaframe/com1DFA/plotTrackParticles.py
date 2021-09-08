import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.plotUtils as pU


def plotTrackParticle(Particles, trackedPartProp, cfgTrackPart, demOri, dem):
    radius = cfgTrackPart.getfloat('radius')
    centerList = cfgTrackPart['center']
    centerList = centerList.split('|')
    center = {'x': np.array([float(centerList[0])]),
              'y': np.array([float(centerList[1])])}
    center, _ = geoTrans.projectOnRaster(dem, center, interp='bilinear')

    header = demOri['header']
    ncols = header['ncols']
    nrows = header['nrows']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    XX = PointsX[0, :]
    YY = PointsY[:, 0]

    # do some ploting
    # ToDo: put this in a plotting folder (here just to demonstrate how this works)
    fig = plt.figure()
    ax = plt.subplot(131)
    ax.plot(trackedPartProp['time'], DFAtls.norm(
                              trackedPartProp['ux'], trackedPartProp['uy'], trackedPartProp['uz']))

    ax.set_xlabel('t [s]')
    ax.set_ylabel('v [m/s]')
    ax1 = plt.subplot(132)
    ax1.plot(trackedPartProp['time'], trackedPartProp['m'])

    ax1.set_xlabel('t [s]')
    ax1.set_ylabel('m [kg]')

    ax2 = plt.subplot(133)
    circle1 = plt.Circle((center['x'], center['y']), radius, color='r')
    ax2.plot(trackedPartProp['x'], trackedPartProp['y'])
    ax2.add_patch(circle1)
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    plt.show()

    fig2 = plt.figure()
    ax1 = plt.subplot(111)
    for count in range(len(Particles)):
        update(count, Particles, xllc, yllc, ax1, XX, YY, dem)
    plt.show()

    # ani = FuncAnimation(fig2, update, round(len(Particles)),
    #                     fargs=(Particles, xllc, yllc, ax1, XX, YY, dem))
    # # plt.show()
    #
    # writer = PillowWriter(fps=4)
    # # ani.save("MalSecRel.gif", writer=writer)
    # ani.save("testTrackAlr1.gif", writer=writer)


def update(count, Particles, xllc, yllc, ax1, XX, YY, dem):
    particles = Particles[count]
    X = particles['x'] + xllc
    Y = particles['y'] + yllc

    ax1.clear()
    ax1.set_title('t=%.2f s' % particles['t'])
    variable = particles['trackedParticles']
    cmap = pU.cmapGreys
    cmap.set_bad(color='white')
    ref0, im = pU.NonUnifIm(ax1, XX, YY, dem['Nz']/DFAtls.norm(dem['Nx'], dem['Ny'], dem['Nz']), 'x [m]', 'y [m]',
                            # extent=[2400, 2700, YY.min(), YY.max()],
                            extent=[XX.min(), XX.max(),
                                    YY.min(), YY.max()],
                            cmap=cmap, norm=None)
    Cp1 = ax1.contour(XX, YY, dem['rasterData'], levels=10, colors='k')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapPres, np.amin(variable),
                                           np.amax(variable), continuous=True)
    # set range and steps of colormap
    cc = np.where(variable == 1, True, False)
    sc = ax1.scatter(X, Y, c='b', cmap=None, marker='.')
    sc = ax1.scatter(X[cc], Y[cc],
                     c='r', cmap=None, marker='.', s=5)

    plt.pause(0.1)
