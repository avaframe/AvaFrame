import pathlib
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from avaframe.in3Utils import geoTrans
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outCom1DFA as outCom1DFA

def hybridProfilePlot(avalancheDir, resultsHybrid):
    """Update profile plot with result of curent iteration"""
    fig = plt.figure(figsize=(3*pU.figW, 2*pU.figH))
    ax = plt.subplot(111)
    nIter = len(resultsHybrid.keys())
    i = 0
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, 0, nIter, continuous=pU.contCmap)
    for key, dict in resultsHybrid.items():
        avaProfileMassExt = dict['path']
        avaProfileMassExt = geoTrans.computeS(avaProfileMassExt)
        alpha = dict['alpha']
        sBetaPoint = dict['sBetaPoint']
        col = cmap(norm(i))
        # Plot the whole profile with beta, alpha ... points and lines
        ax.plot(avaProfileMassExt['s'], avaProfileMassExt['z'], linestyle='-', color=col, label='profile (iteration %d)' % i)
        ax.axvline(x=sBetaPoint, color=col, linestyle=':', linewidth=1, label='Beta point (iteration %d)' % i)
        s = avaProfileMassExt['s'][[0, -1]]
        z = avaProfileMassExt['z'][0] - s*np.tan(np.deg2rad(alpha))
        ax.plot(s, z, color=col, linestyle='--', label='AlphaLine (iteration %d)' % i)
        i = i+1

    titleText = r'Profiles extracted from the DFA simulations with corresponding $\alpha-\beta$ model results'
    ax.set_title(titleText)
    ax.set_xlabel('projectd length s [m]')
    ax.set_ylabel('Height [m]')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(linestyle=':', color='0.9')
    ax.legend(frameon=False)
    title = ('com3HybProfPlot')
    l = ax.legend(loc='lower left')
    l.set_zorder(40)
    pU.putAvaNameOnPlot(ax, avalancheDir)
    path = pathlib.Path(avalancheDir, 'Outputs', 'com3Hybrid')
    pU.saveAndOrPlot({'pathResult': path}, title, fig)


def hybridPathPlot(avalancheDir, dem, resultsHybrid, fields, particles, muArray):
    """Update path plot with result of curent iteration"""
    headerOri = dem['originalHeader']
    xllcOri = headerOri['xllcenter']
    yllcOri = headerOri['yllcenter']
    fig = plt.figure(figsize=(3*pU.figW, 2*pU.figH))
    ax = plt.subplot(111)
    nIter = len(resultsHybrid.keys())
    i = 0
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, 0, nIter, continuous=pU.contCmap)
    for key, dict in resultsHybrid.items():
        mu = muArray[i]
        avaProfileMassExt = dict['path']
        xAB = dict['xAB']
        yAB = dict['yAB']
        col = cmap(norm(i))
        ax.plot(avaProfileMassExt['x'], avaProfileMassExt['y'], color=col,
                label='Center of mass path iteration %d ($\mu$ = %.2f )' % (i, mu), zorder = 20)
        ax.plot(xAB - xllcOri, yAB - yllcOri, 'X', color=col, markersize=8,
                label=r'com2AB $\alpha$ point iteration %d' % i, zorder = 20)
        i = i+1
    titleText = 'Avalanche path for each iteration with peak travel angle field \n and particles flow thickness in final time step'
    ax.set_title(titleText)
    ax.set_ylabel('x [m]')
    ax.set_ylabel('y [m]')
    ax, extent, cb, CS = outCom1DFA.addResult2Plot(ax, dem['header'], fields['pta'], 'pta')
    cb.ax.set_ylabel(pU.cfgPlotUtils['namepta'])
    ax = outCom1DFA.addDem2Plot(ax, dem, what='slope', extent=extent)
    ax, cb2 = outCom1DFA.addParticles2Plot(particles, ax, dem, whatS='h')
    cb2.ax.set_ylabel('particle ' + pU.cfgPlotUtils['nameFT'])
    ax.set_ylim(extent[2:])
    ax.set_xlim(extent[:2])
    title = ('com3HybRasterPlot')
    l = ax.legend(loc='lower left')
    l.set_zorder(40)
    pU.putAvaNameOnPlot(ax, avalancheDir)
    path = pathlib.Path(avalancheDir, 'Outputs', 'com3Hybrid')
    pU.saveAndOrPlot({'pathResult': path}, title, fig)
