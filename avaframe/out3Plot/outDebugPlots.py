import numpy as np
import pathlib
import copy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import avaframe.com1DFA.DFAtools as DFAtls


import avaframe.out3Plot.plotUtils as pU


def plotPartIni(particles, dem):
    header = dem['header']
    x = np.arange(header['ncols']) * header['cellsize']
    y = np.arange(header['nrows']) * header['cellsize']
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    ref0, im = pU.NonUnifIm(ax, x, y, dem['areaRaster'], 'x [m]', 'y [m]',
                            extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, norm=None)

    ax.plot(particles['x'], particles['y'], 'or', linestyle='None')
    pU.addColorBar(im, ax, None, 'm²')
    plt.show()


def plotAreaDebug(dem, avapath, Raster):
    ncols = dem['header']['ncols']
    nrows = dem['header']['nrows']
    cellsize = dem['header']['cellsize']
    x = np.arange(ncols) * cellsize
    y = np.arange(nrows) * cellsize
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Release area')
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    ref0, im = pU.NonUnifIm(ax, x, y, Raster, 'x [m]', 'y [m]',
                            extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, norm=None)
    ax.plot(avapath['x'] * cellsize, avapath['y'] * cellsize, 'r', label='release polyline')
    plt.legend()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax)
    plt.show()


def plotRemovePart(xCoord0, yCoord0, header, X, Y, Mask, mask):
    x = np.arange(header['ncols']) * header['cellsize']
    y = np.arange(header['nrows']) * header['cellsize']
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Release area')
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    ref0, im = pU.NonUnifIm(ax, x, y, Mask, 'x [m]', 'y [m]',
                            extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, norm=None)
    ax.plot(xCoord0 * header['cellsize'], yCoord0 * header['cellsize'], 'r', label='release polyline')
    ax.plot(X[mask] * header['cellsize'], Y[mask] * header['cellsize'], '.b')
    plt.legend()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax)
    plt.show()


def plotPartAfterRemove(points, xCoord0, yCoord0, mask):
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    ax.set_title('Release area')
    ax.plot(xCoord0, yCoord0, 'r', label='release polyline')
    ax.plot(points['x'], points['y'], '.b')
    ax.plot(points['x'][mask], points['y'][mask], '.g')
    plt.legend()
    plt.show()


def analysisPlots(particlesList, fieldsList, cfg, demOri, dem, outDir):
    """ create analysis plots during simulation run """

    cfgGen = cfg['GENERAL']
    partRef = particlesList[0]
    Z0 = partRef['z'][0]
    rho = cfgGen.getfloat('rho')
    gravAcc = cfgGen.getfloat('gravAcc')
    mu = cfgGen.getfloat('mu')
    repeat = True
    while repeat:
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        T = np.array([0])
        Z = np.array([0])
        U = np.array([0])
        S = np.array([0])
        for part, field in zip(particlesList, fieldsList):
            T = np.append(T, part['t'])
            S = np.append(S, part['s'][0])
            Z = np.append(Z, part['z'][0])
            U = np.append(U, DFAtls.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
            fig, ax = plotPosition(
                fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True)
            fig.savefig(pathlib.Path(outDir, 'particlest%f.%s' % (part['t'], pU.outputFormat)))

        fig, ax = plotPosition(
                fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True, last=True)
        fig.savefig(pathlib.Path(outDir, 'particlesFinal.%s' % (pU.outputFormat)))
        value = input("[y] to repeat:\n")
        if value != 'y':
            repeat = False

    fieldEnd = fieldsList[-1]
    partEnd = particlesList[-1]
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig3, ax3 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig1, ax1 = plotPosition(
        fig1, ax1, partEnd, demOri, fieldEnd['FD'], pU.cmapPres, 'm', plotPart=False)
    fig2, ax2 = plotPosition(
        fig2, ax2, partEnd, demOri, fieldEnd['FV'], pU.cmapPres, 'm/s', plotPart=False)
    fig3, ax3 = plotPosition(
        fig3, ax3, partEnd, demOri, fieldEnd['P']/1000, pU.cmapPres, 'kPa', plotPart=False)
    plt.show()


def plotPosition(fig, ax, particles, dem, data, Cmap, unit, plotPart=False, last=False):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    X = PointsX[0, :]
    Y = PointsY[:, 0]
    Z = dem['rasterData']
    x = particles['x'] + xllc
    y = particles['y'] + yllc
    xx = np.arange(ncols) * csz + xllc
    yy = np.arange(nrows) * csz + yllc
    try:
        # Get the images on an axis
        cb = ax.images[-1].colorbar
        if cb:
            cb.remove()
    except IndexError:
        pass

    ax.clear()
    ax.set_title('t=%.2f s' % particles['t'])
    cmap, _, ticks, norm = pU.makeColorMap(Cmap, np.nanmin(data), np.nanmax(data), continuous=True)
    cmap.set_under(color='w')
    ref0, im = pU.NonUnifIm(ax, xx, yy, data, 'x [m]', 'y [m]',
                         extent=[x.min(), x.max(), y.min(), y.max()],
                         cmap=cmap, norm=norm)

    Cp1 = ax.contour(X, Y, Z, levels=10, colors='k')
    pU.addColorBar(im, ax, ticks, unit)
    if plotPart:
        # ax.plot(x, y, '.b', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 1], y[NPPC == 1], '.c', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 4], y[NPPC == 4], '.b', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 9], y[NPPC == 9], '.r', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 16], y[NPPC == 16], '.m', linestyle='None', markersize=1)
        # load variation colormap
        variable = particles['h']
        cmap, _, ticks, norm = pU.makeColorMap(pU.cmapDepth, np.nanmin(data), np.amax(variable), continuous=True)
        # set range and steps of colormap
        cc = variable
        sc = ax.scatter(x, y, c=cc, cmap=cmap, marker='.')

        if last:
            pU.addColorBar(sc, ax, ticks, 'm', 'Flow Depth')

    plt.pause(0.1)
    return fig, ax


def plotContours(fig, ax, particles, dem, data, Cmap, unit, last=False):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    X = PointsX[0, :]
    Y = PointsY[:, 0]
    Z = dem['rasterData']
    try:
        # Get the images on an axis
        cb = ax.images[-1].colorbar
        if cb:
            cb.remove()
    except IndexError:
        pass

    ax.clear()
    ax.set_title('t=%.2f s' % particles['t'])
    cmap, _, ticks, norm = pU.makeColorMap(Cmap, np.nanmin(data), np.nanmax(data), continuous=True)
    cmap.set_under(color='w')

    CS = ax.contour(X, Y, data, levels=8, origin='lower', cmap=cmap,
                    linewidths=2)
    lev = CS.levels

    if last:
        # pU.addColorBar(im, ax, ticks, unit, 'Flow Depth')
        CB = fig.colorbar(CS)
        ax.clabel(CS, inline=1, fontsize=8)

    plt.pause(0.1)
    return fig, ax, cmap, lev