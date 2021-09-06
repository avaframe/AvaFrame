import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in1Data import getInput as gI
import avaframe.in2Trans.ascUtils as IOf
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.plotUtils as pU


def runTrackParticle():
    # +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
    # log file name; leave empty to use default runLog.log
    logName = 'runTrackParticles'
    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']

    modCfg = cfgUtils.getDefaultModuleConfig(com1DFA)

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT Track Particles')
    log.info('Current avalanche: %s', avalancheDir)

    # get dem
    inputSimFiles = gI.getInputDataCom1DFA(avalancheDir, modCfg['FLAGS'])
    demOri = IOf.readRaster(inputSimFiles['demFile'], noDataToNan=True)
    demOri, dem = com1DFA.initializeMesh(modCfg['GENERAL'], demOri, 1)
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

    # get particles
    inDirPart = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'particles')
    Particles, TimeStepInfo = com1DFA.readPartFromPickle(inDirPart)

    # location of particles to track
    # ToDo: put this in an ini file
    radius = 5
    # center = {'x': np.array([2800])-xllc, 'y': np.array([-4000])-yllc}
    center = {'x': np.array([255865])-xllc, 'y': np.array([381383])-yllc}
    center, _ = geoTrans.projectOnRaster(dem, center)
    # properties to track
    # ToDo: put this in an ini file
    properties = ['x', 'y', 'z', 'ux', 'uy', 'uz']

    # get particles and related properties of particles to track
    Particles, trackedPartProp = com1DFA.mainTrackParticles(
        Particles, TimeStepInfo, center, radius, properties)

    # do some ploting
    # ToDo: put this in a plotting folder (here just to demonstrate how this works)
    fig = plt.figure()
    ax = plt.subplot(121)
    ax.plot(trackedPartProp['time'], DFAtls.norm(
        trackedPartProp['ux'], trackedPartProp['uy'], trackedPartProp['uz']))

    ax.set_xlabel('t [s]')
    ax.set_ylabel('v [m/s]')

    ax2 = plt.subplot(122)
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


if __name__ == "__main__":
    runTrackParticle()
