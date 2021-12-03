import pathlib
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from avaframe.in3Utils import geoTrans
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outCom1DFA as outCom1DFA


def initializeFigures():
    figPath, axPath = plt.subplots(figsize=(1.5*pU.figW, 1*pU.figH))
    # figPath = plt.figure(figsize=(1.5*pU.figW, 1*pU.figH))
    # axPath = plt.subplot(111)
    figProf = plt.figure(figsize=(1.5*pU.figW, 1*pU.figH))
    axProf = plt.subplot(111)
    figDict = {'figPath': figPath, 'axPath': axPath, 'figProf': figProf, 'axProf': axProf}
    return figDict


def updateProfilePlot(resAB, name, figDict, iteration):
    """ Plot and save results depending on flags options"""
    fig = figDict['figProf']
    ax = figDict['axProf']
    s = resAB[name]['s']
    z = resAB[name]['z']
    ids10Point = resAB[name]['ids10Point']
    beta = resAB[name]['beta']
    alpha = resAB[name]['alpha']
    f = resAB[name]['f']
    # Plot the whole profile with beta, alpha ... points and lines
    ax.plot(s, z, '-', label='profile (iteration %d)' % iteration)
    ax.axvline(x=s[ids10Point], color='0.8', linewidth=1, linestyle='-.', label='Beta point (iteration %d)' % iteration)
    ax.plot(s, f, '-', label='AlphaLine (iteration %d)' % iteration)
    figDict['figProf'] = fig
    figDict['axProf'] = ax
    return figDict


def updatePathPlot(demOri, avaProfile, resAB, name, figDict, iteration):
    """ Plot and save results depending on flags options"""
    fig = figDict['figPath']
    ax = figDict['axPath']
    headerOri = demOri['header']
    xllcOri = headerOri['xllcenter']
    yllcOri = headerOri['yllcenter']

    alpha = resAB[name]['alpha']
    ids_alpha = resAB[name]['ids_alpha']
    xAB = resAB[name]['x'][ids_alpha]
    yAB = resAB[name]['y'][ids_alpha]

    ax.plot(xAB - xllcOri, yAB - yllcOri, 'X', markersize=8, label=r'com2AB $\alpha$ point iteration %d' % iteration, zorder = 20)
    ax.plot(avaProfile['x'] - xllcOri, avaProfile['y'] - yllcOri, '-', label='Center of mass path iteration %d' % iteration, zorder = 20)
    figDict['figPath'] = fig
    figDict['axPath'] = ax
    return figDict


def finalizePathPlot(avalancheDir, figDict, resAnalysis, indSim, dem, demOri, particles, fields):
    headerOri = demOri['header']
    xllcOri = headerOri['xllcenter']
    yllcOri = headerOri['yllcenter']
    xAIMEC = resAnalysis['runout'][1][indSim]
    yAIMEC = resAnalysis['runout'][2][indSim]

    fig = figDict['figPath']
    ax = figDict['axPath']
    titleText = 'Path plot'
    ax.set_title(titleText)
    ax.set_ylabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.plot(xAIMEC - xllcOri, yAIMEC - yllcOri, 'X', markersize=8, color='0.7', label='Runout point com1DFA (AIMEC pfd=0m)')
    ax, extent = outCom1DFA.addResult2Plot(ax, dem, fields['pta'], 'pta')
    ax = outCom1DFA.addDem2Plot(ax, dem, what='slope', extent=extent)
    ax = outCom1DFA.addParticles2Plot(particles, ax, dem, whatS='h')
    ax.set_ylim(extent[2:])
    ax.set_xlim(extent[:2])
    title = ('ana5HybRasterPlot')
    l = ax.legend(loc='lower left')
    l.set_zorder(40)
    path = pathlib.Path(avalancheDir, 'Outputs', 'ana5Hybrid')
    pU.saveAndOrPlot({'pathResult': path}, title, fig)


def finalizeProfilePlot(avalancheDir, figDict, resAnalysis, indSim):
    sAIMEC = resAnalysis['runout'][0][indSim]

    fig = figDict['figProf']
    ax = figDict['axProf']
    plt.axvline(x=sAIMEC, color='0.7', linewidth=1, linestyle='--', label='RunOut distance AIMEC')
    titleText = 'Profile plot'
    ax.set_title(titleText)
    ax.set_ylabel('projectd length s [m]')
    ax.set_ylabel('Height [m]')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(linestyle=':', color='0.9')
    ax.legend(frameon=False)
    title = ('ana5HybProfPlot')
    l = ax.legend(loc='lower left')
    l.set_zorder(40)
    path = pathlib.Path(avalancheDir, 'Outputs', 'ana5Hybrid')
    pU.saveAndOrPlot({'pathResult': path}, title, fig)


def plotHybridRes(avalancheDir, cfg, resAB, name, simID, demOri, avaProfileMassNew):

    alpha = resAB[name]['alpha']

    V2Path = avaProfileMassNew['v2']
    EkinPath = avaProfileMassNew['ekin']

    g = cfg['GENERAL'].getfloat('gravAcc')

    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    cmap = pU.cmapPlasma
    cmap.set_under(color='w')
    projectedZ = 'yes'
    if projectedZ == 'yes':
        avaProfileMassNew, _ = geoTrans.projectOnRaster(demOri, avaProfileMassNew, interp='bilinear')
    Zene = avaProfileMassNew['z'] + V2Path/(2*g)
    # Colorbar: kinietische Energie [J]
    scat = ax2.scatter(avaProfileMassNew['sCor'], Zene, marker='s', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Gesamtenergie(s_mod)')
    scat = ax2.scatter(avaProfileMassNew['s'], Zene, marker='o', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Gesamtenergie(s_real)')
    cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('kinetische Energie [J]')

    ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'], 'b-', label='Z_av(s_real)')
    # ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'] + 2*avaProfileMassNew['zstd'], 'b:')
    # ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'] - 2*avaProfileMassNew['zstd'], 'b:')
    # ax2.plot(avaProfileMassNew['s'] + 2*avaProfileMassNew['sstd'], avaProfileMassNew['z'], 'b--')
    # ax2.plot(avaProfileMassNew['s'] - 2*avaProfileMassNew['sstd'], avaProfileMassNew['z'], 'b--')
    ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'], 'k-', label='Z_av(s_mod)')
    # ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'] + 2*avaProfileMassNew['zstd'], 'k:')
    # ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'] - 2*avaProfileMassNew['zstd'], 'k:')
    # ax2.plot(avaProfileMassNew['sCor'] + 2*avaProfileMassNew['sstd'], avaProfileMassNew['z'], 'k--')
    # ax2.plot(avaProfileMassNew['sCor'] - 2*avaProfileMassNew['sstd'], avaProfileMassNew['z'], 'k--')
    GK = avaProfileMassNew['sCor'][-1] * np.tan(alpha*np.pi/180)
    z_enda = avaProfileMassNew['z'][0] - GK
    s_geomL = [avaProfileMassNew['sCor'][0], avaProfileMassNew['sCor'][-1]]
    z_geomL = [avaProfileMassNew['z'][0], z_enda]
    ax2.plot(s_geomL, z_geomL, 'r-', linewidth=0.3, label='alpha line from Z_av func s_mod')
    GK = avaProfileMassNew['s'][-1] * np.tan(alpha*np.pi/180)
    z_enda = avaProfileMassNew['z'][0] - GK
    s_geomL = [avaProfileMassNew['s'][0], avaProfileMassNew['s'][-1]]
    z_geomL = [avaProfileMassNew['z'][0], z_enda]
    ax2.plot(s_geomL, z_geomL, 'g-', linewidth=0.3, label='alpha line from Z_true func s_mod')

    ax2.set_xlabel('s [m]', fontsize=pU.fs)
    ax2.set_ylabel('z [m]', fontsize=pU.fs)
    ax2.legend(loc='lower left')

    # set titel of output png
    title = ('ana5HybEnergyProfilePlot')
    path = pathlib.Path(avalancheDir, 'Outputs', 'ana5Hybrid')
    pU.saveAndOrPlot({'pathResult': path}, title, fig2)
