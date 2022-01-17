import pathlib
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from avaframe.in3Utils import geoTrans
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outCom1DFA as outCom1DFA


def initializeFigures():
    """ Initialize the two figure from com3Hybrid main loop"""
    figPath, axPath = plt.subplots(figsize=(1.5*pU.figW, 1*pU.figH))
    # figPath = plt.figure(figsize=(1.5*pU.figW, 1*pU.figH))
    # axPath = plt.subplot(111)
    figProf = plt.figure(figsize=(1.5*pU.figW, 1*pU.figH))
    axProf = plt.subplot(111)
    figDict = {'figPath': figPath, 'axPath': axPath, 'figProf': figProf, 'axProf': axProf}
    return figDict


def updateProfilePlot(resAB, name, figDict, iteration):
    """Update profile plot with result of curent iteration"""
    fig = figDict['figProf']
    ax = figDict['axProf']
    s = resAB[name]['s']
    z = resAB[name]['z']
    ids10Point = resAB[name]['ids10Point']
    beta = resAB[name]['beta']
    alpha = resAB[name]['alpha']
    f = resAB[name]['f']
    # Plot the whole profile with beta, alpha ... points and lines
    ax.plot(s, z, linestyle=pU.ls[iteration], color='k', label='profile (iteration %d)' % iteration)
    ax.axvline(x=s[ids10Point], color='0.8', linewidth=1, linestyle=pU.ls[iteration], label='Beta point (iteration %d)' % iteration)
    ax.plot(s, f, linestyle=pU.ls[iteration], color='b', label='AlphaLine (iteration %d)' % iteration)
    figDict['figProf'] = fig
    figDict['axProf'] = ax
    return figDict


def updatePathPlot(demOri, avaProfile, resAB, name, figDict, iteration):
    """Update path plot with result of curent iteration"""
    fig = figDict['figPath']
    ax = figDict['axPath']
    headerOri = demOri['header']
    xllcOri = headerOri['xllcenter']
    yllcOri = headerOri['yllcenter']

    alpha = resAB[name]['alpha']
    ids_alpha = resAB[name]['ids_alpha']
    xAB = resAB[name]['x'][ids_alpha]
    yAB = resAB[name]['y'][ids_alpha]

    ax.plot(xAB - xllcOri, yAB - yllcOri, 'X', color=str(0.8-(iteration+1)/4), markersize=8, label=r'com2AB $\alpha$ point iteration %d' % iteration, zorder = 20)
    ax.plot(avaProfile['x'] - xllcOri, avaProfile['y'] - yllcOri, color=str(0.8-(iteration+1)/4), linestyle=pU.ls[iteration], label='Center of mass path iteration %d' % iteration, zorder = 20)
    figDict['figPath'] = fig
    figDict['axPath'] = ax
    return figDict


def finalizePathPlot(avalancheDir, figDict, resAnalysisDF, refSimulation, dem, demOri, particles, fields):
    """ Add title ... and save path plot"""
    headerOri = demOri['header']
    xllcOri = headerOri['xllcenter']
    yllcOri = headerOri['yllcenter']
    xAIMEC = resAnalysisDF.loc[refSimulation, 'xRunout']
    yAIMEC = resAnalysisDF.loc[refSimulation, 'yRunout']

    fig = figDict['figPath']
    ax = figDict['axPath']
    titleText = 'Avalanche Path with peak travel angle field \n and particles flow depth in final time step'
    ax.set_title(titleText)
    ax.set_ylabel('x [m]')
    ax.set_ylabel('y [m]')
    ax, extent = outCom1DFA.addResult2Plot(ax, dem, fields['pta'], 'pta')
    ax = outCom1DFA.addDem2Plot(ax, dem, what='slope', extent=extent)
    ax = outCom1DFA.addParticles2Plot(particles, ax, dem, whatS='h')
    ax.plot(xAIMEC - xllcOri, yAIMEC - yllcOri, 'P', markersize=8, color='r',
            label='Runout point com1DFA (AIMEC pfd=0m)', zorder=40)
    ax.set_ylim(extent[2:])
    ax.set_xlim(extent[:2])
    title = ('com3HybRasterPlot')
    l = ax.legend(loc='lower left')
    l.set_zorder(40)
    pU.putAvaNameOnPlot(ax, avalancheDir)
    path = pathlib.Path(avalancheDir, 'Outputs', 'com3Hybrid')
    pU.saveAndOrPlot({'pathResult': path}, title, fig)


def finalizeProfilePlot(avalancheDir, figDict, resAnalysisDF, refSimulation):
    """Add title ... and save profile plot"""
    sAIMEC = resAnalysisDF.loc[refSimulation, 'sRunout']

    fig = figDict['figProf']
    ax = figDict['axProf']
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


def plotEnergyProfile(avalancheDir, cfg, resAB, name, simID, demOri, avaProfileMassNew):
    """ Make energy line profile plot"""
    alpha = resAB[name]['alpha']

    V2Path = avaProfileMassNew['v2']
    EkinPath = avaProfileMassNew['ekin']

    g = cfg['GENERAL'].getfloat('gravAcc')

    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    titleText = r"""Center of mass profile extracted from the DFA simulation
    function of the real and modified s with corresponding
    $\alpha$ line and energy height points"""
    ax2.set_title(titleText)
    cmap = pU.cmapPlasma
    cmap.set_under(color='w')
    projectedZ = 'no'
    if projectedZ == 'yes':
        avaProfileMassNew, _ = geoTrans.projectOnRaster(demOri, avaProfileMassNew, interp='bilinear')
    Zene = avaProfileMassNew['z'] + V2Path/(2*g)
    # Colorbar: kinietische Energie [J]
    scat = ax2.scatter(avaProfileMassNew['sCor'], Zene, marker='s', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Energy altitude(s_mod)')
    scat = ax2.scatter(avaProfileMassNew['s'], Zene, marker='o', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Energy altitude(s_real)')
    cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('kinetic energy [J]')

    ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'], 'b-', label='Z_av(s_real)')
    ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'], 'k-', label='Z_av(s_mod)')
    GK = 1.05*avaProfileMassNew['s'][-1] * np.tan(alpha*np.pi/180)
    z_end = avaProfileMassNew['z'][0] - GK
    s_geomL = [avaProfileMassNew['s'][0], 1.05*avaProfileMassNew['s'][-1]]
    z_geomL = [avaProfileMassNew['z'][0], z_end]
    ax2.plot(s_geomL, z_geomL, 'r-', label='alpha line %.2f °' % alpha)

    ax2.set_xlabel('s [m]', fontsize=pU.fs)
    ax2.set_ylabel('z [m]', fontsize=pU.fs)
    ax2.legend(loc='lower left')
    pU.putAvaNameOnPlot(ax2, avalancheDir)
    plt.show()

    # set titel of output png
    title = ('com3HybEnergyProfilePlot')
    path = pathlib.Path(avalancheDir, 'Outputs', 'com3Hybrid')
    pU.saveAndOrPlot({'pathResult': path}, title, fig2)
