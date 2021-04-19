"""
    trackedPath
    get information from Outputs/com1DFAPy/particles/particlesxxxx.xxxx.p
    create shape file as Inputs/LINES/pathAB for com2AB
"""

import shapefile
import numpy as np
import matplotlib.pyplot as plt

# Local imports
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.ascUtils as IOf
import avaframe.com1DFAPy.DFAtools as DFAtls
from avaframe.out3Plot.plotUtils import *
from avaframe.out3Plot.makePalette import *


def writeLine2SHPfile(part, lineName, fileName):
    """ copied from
    https://pypi.org/project/pyshp/#writing-shapefiles
    section: Adding a LineString shape
    """

    w = shapefile.Writer(fileName)
    w.field('name', 'C')
    w.line([part])
    w.record(lineName)
    w.close()


def trackParticle():

    # inDir = '/home/marie/ava0/AvaFrame/avaframe/data/avaParabola/Outputs/com1DFAPy/particles'
    inDir = '/home/marie/ava0/AvaFrame/avaframe/data/avaAlr/Outputs/com1DFAPy/particles'
    # inDEM = '/home/marie/ava0/AvaFrame/avaframe/data/avaParabola/Inputs/DEM_PF_Topo.asc'
    inDEM = '/home/marie/ava0/AvaFrame/avaframe/data/avaAlr/Inputs/avaAlr.asc'

    # inDir = '/home/matthias/Documents/github/AvaFrame/avaframe/data/avaParabola/Outputs/com1DFAPy/particles'
    # inDEM = '/home/matthias/Documents/github/AvaFrame/avaframe/data/avaParabola/Inputs/DEM_PF_Topo.asc'

    inDir = '/home/matthias/Documents/github/AvaFrame/avaframe/data/avaAlr0/Outputs/com1DFAPy/particles'
    inDEM = '/home/matthias/Documents/github/AvaFrame/avaframe/data/avaAlr0/Inputs/avaAlr.asc'

    # inDir = '/home/matthiastonnel/Documents/github/AvaFrame/avaframe/data/avaInclinedPlane/Outputs/com1DFAPy/particles'
    # inDEM = '/home/matthiastonnel/Documents/github/AvaFrame/avaframe/data/avaInclinedPlane/Inputs/DEM_IP_Topo.asc'
    header = IOf.readASCheader(inDEM)
    dem = IOf.readRaster(inDEM)
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    XX = PointsX[0, :]
    YY = PointsY[:, 0]
    ZZ = dem['rasterData']
    Particles, TimeStepInfo = com1DFA.readPartFromPickle(inDir, flagAvaDir=False)

    # listt = []
    # nTimeSteps = len(TimeStepInfo)
    # save pairs of x,y values
    # (this will be a 2d array with 2 cols and as many rows as time steps)
    xPath = np.empty((0, 1))
    yPath = np.empty((0, 1))
    zPath = np.empty((0, 1))
    sPath = np.empty((0, 1))
    lPath = np.empty((0, 1))
    V2Path = np.empty((0, 1))
    EkinPath = np.empty((0, 1))
    EpotPath = np.empty((0, 1))
    count = 0
    fig2 = plt.figure(figsize=(2*figW, figH))
    ax = plt.subplot(111)
    for t in TimeStepInfo:
        particles = Particles[count]
        # force = Force[count]
        m = particles['m']
        X = particles['x'] + xllc
        Y = particles['y'] + yllc
        Z = particles['z']
        ux = particles['ux']
        uy = particles['uy']
        uz = particles['uz']
        u = DFAtls.norm(ux, uy, uz)
        U2 = u*u
        Npart = particles['Npart']
        S = particles['s']
        L = particles['l']
        kineticEne = 0.5*m*u*u
        kineticEneSum = np.sum(kineticEne)
        # option to take not one particle but a mass-averaged path
        if kineticEneSum <= 100:
            pond = np.ones(np.shape(kineticEne))
            pondSum = Npart
        else:
            pond = kineticEne
            pondSum = kineticEneSum

        # pond = np.zeros(np.shape(kineticEne))
        # pond[13000] = 1
        # pondSum = 1
        pond = m
        pondSum = np.sum(m)
        # pond = np.ones(np.shape(kineticEne))
        # pondSum = Npart

        zcoE = np.nansum(pond*Z)/pondSum
        # calculate XY locations of computation step
        xcoE = np.sum(pond*X)/pondSum
        ycoE = np.sum(pond*Y)/pondSum
        # calculate global distance of coE coordinates
        scoE = np.sum(pond*S)/pondSum
        lcoE = np.sum(pond*L)/pondSum
        v2coE = np.sum(pond*U2)/pondSum

        # append x, y and update distance
        xPath = np.append(xPath, xcoE)
        yPath = np.append(yPath, ycoE)
        zPath = np.append(zPath, zcoE)
        sPath = np.append(sPath, scoE)
        lPath = np.append(lPath, lcoE)
        V2Path = np.append(V2Path, v2coE)

        # update energy
        EkinPath = np.append(EkinPath, kineticEneSum)
        # EpotPath = np.append(EpotPath, EpotSumCoE)
        count = count + 1

        ax.clear()
        ax.set_title('t=%.2f s' % particles['t'])
        variable = particles['h']
        Cp1 = ax.contour(XX, YY, ZZ, levels=10, colors='k')
        cmap, _, _, norm, ticks = makeColorMap(
            cmapDepth, np.amin(variable), np.amax(variable), continuous=True)
        # set range and steps of colormap
        cc = variable
        sc = ax.scatter(X, Y, c=cc, cmap=cmap, marker='.')
        plt.pause(0.01)


    print(np.min(m), np.mean(m), np.max(m), np.std(m))
    mu = 0.55
    g = 9.81
    avapath = {}
    avapath['x'] = np.array([xPath[0], xPath[-1]])
    avapath['y'] = np.array([yPath[0], yPath[-1]])
    header = IOf.readASCheader(inDEM)
    dem = IOf.readRaster(inDEM)
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    XX = PointsX[0, :]
    YY = PointsY[:, 0]
    ZZ = dem['rasterData']
    AvaProfile, projPoint = geoTrans.prepareLine(dem, avapath, distance=10, Point=None)

    fig = plt.figure(figsize=(2*figW, figH))
    ax1 = plt.subplot(211)
    cmap = cmapPlasma
    cmap.set_under(color='w')
    ax1.plot(xPath, yPath, 'k--', label='avapath')
    sc = ax1.scatter(X, Y, c=cc, cmap=cmap, marker='.')
    Cp1 = ax.contour(XX, YY, ZZ, levels=10, colors='k')
    # ax1.axis('equal')
    # ax1.set_xlim([x.min(), x.max()])
    # ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel(r'$x\;[m]$')
    ax1.set_ylabel(r'$y\;[m]$')
    ax1.legend()

    ax2 = plt.subplot(212)

    ax2.plot(sPath, zPath, 'k-', label='Avalanche profile')
    # ax2.plot(AvaProfile['s'], AvaProfile['z'], 'r-', label='Avalanche profile')
    Zene = zPath + V2Path/(2*g)
    f = zPath[5] - mu * sPath + (V2Path[5]/(2*g) + mu * sPath[5])
    ax2.plot(sPath, f, '-', color='b', label='AlphaLine')
    scat = ax2.scatter(sPath, Zene, marker='s', cmap=cmap, s=2*ms, c= EkinPath, label='Total energy height')
    cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('Kinetic Energy [J]')
    ax2.axvline(x=sPath[-1], color='k',
    linewidth=1, linestyle='-.', label='Run out point')
    ax2.set_xlabel('s [m]', fontsize=fs)
    ax2.set_ylabel('Altitude [m]', fontsize=fs)
    # ax1.set_title('title subplot 1')
    ax2.legend()
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    trackParticle()
