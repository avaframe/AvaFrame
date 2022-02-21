"""
Energy line test

This module runs a DFA simulation extracts the center of mass path
and compares it to he analytic geometric/alpha line solution
"""
import pathlib
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import interpolate

# Local imports
# import config and init tools
from avaframe.in3Utils import cfgUtils
from avaframe.in1Data import getInput

# import computation modules
from avaframe.com1DFA import com1DFA, particleTools

# import analysis tools

# import plotting tools
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outCom1DFA as outCom1DFA

# create local logger
log = logging.getLogger(__name__)


def mainEnergyLineTest(cfgMain):
    """This is the core function of the com3Hybrid module

    Here the DFA and AB models are run one after the other to first produce (from com1DFA)
    the avalanche path, then get the friction angle corresponding to the topography (from com2AB)
    and finally run the DFA simultaion one last time with the correct friction angle and get a 3D
    output of the avalanche
    """
    avalancheDir = cfgMain['MAIN']['avalancheDir']
    demOri = getInput.readDEM(avalancheDir)
    # get comDFA configuration path for the energy line test
    energyLineTestCfgFile = pathlib.Path('ana1Tests', 'ana1EnergyLineTest_com1DFACfg.ini')
    # Run dense flow with coulomb friction
    energyLineTestCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=energyLineTestCfgFile, modInfo=True)
    dem, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=energyLineTestCfgFile)
    particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName='', flagAvaDir=True,
                                                                   comModule='com1DFA')
    # postprocess to extract path and energy line
    _, avaProfileMass, _ = particleTools.getCom1DFAPath(particlesList, dem)

    simName = simDF.index[0]
    fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['pfv'], simName=simName, flagAvaDir=True,
                                                           comModule='com1DFA')

    fig1 = plt.figure(figsize=(pU.figW*3, pU.figH*2))
    ax1 = plt.subplot2grid((2, 2), (1, 0))
    ax1, extent, cbar0, cs1 = outCom1DFA.addResult2Plot(ax1, dem, fieldsList[-1]['pfv'], 'pfv')
    cbar0.ax.set_ylabel('peak flow velocity')
    ax1 = outCom1DFA.addDem2Plot(ax1, dem, what='slope', extent=extent)
    ax1, cbar1 = outCom1DFA.addParticles2Plot(particlesList[-1], ax1, dem, whatS='h')
    cbar1.ax.set_ylabel('particle flow thickness')
    ax1.plot(avaProfileMass['x'], avaProfileMass['y'], '-k.', zorder = 20, label='Center of mass path')
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(fieldsList[-1]['pfv'], 5, extentOption=True,
                                                                 constrainedData=False, buffer='')
    ax1.set_ylim([rowsMin, rowsMax])
    ax1.set_xlim([colsMin, colsMax])
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    title = ('com3HybRasterPlot')
    l = ax1.legend(loc='upper right')
    l.set_zorder(40)
    pU.putAvaNameOnPlot(ax1, avalancheDir)

    V2Path = avaProfileMass['v2']
    EkinPath = avaProfileMass['ekin']

    g = energyLineTestCfg['GENERAL'].getfloat('gravAcc')
    mu = energyLineTestCfg['GENERAL'].getfloat('mu')
    alphaRad = np.arctan(mu)
    alphaDeg = np.rad2deg(alphaRad)
    deltaS = avaProfileMass['s'][-1] - avaProfileMass['s'][0]
    deltaZ = avaProfileMass['z'][-1] - avaProfileMass['z'][0]
    runOutAngleRad = np.arctan(np.abs(deltaZ/deltaS))
    runOutAngleDeg = np.rad2deg(runOutAngleRad)
    p1 = np.polyfit(avaProfileMass['s'][-3:], avaProfileMass['z'][-3:], 1)

    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps['pfv'], np.amin(V2Path/(2*g)), np.amax(V2Path/(2*g)),
                                           continuous=pU.contCmap)
    cmap.set_under(color='w')

    # compute mass average velocity elevation
    zEne = avaProfileMass['z'] + V2Path/(2*g)
    # Create the alpha line
    GK = 1.05*avaProfileMass['s'][-1] * mu
    z_end = avaProfileMass['z'][0] - GK
    zEneTarget = avaProfileMass['z'][0] - avaProfileMass['s'] * mu
    s_geomL = [avaProfileMass['s'][0], 1.05*avaProfileMass['s'][-1]]
    z_geomL = [avaProfileMass['z'][0], z_end]

    # find intersection between alpha line and profile
    sExt = np.append(avaProfileMass['s'], 1.05*avaProfileMass['s'][-1])
    zExt = np.append(avaProfileMass['z'], avaProfileMass['z'][-1] + p1[0]*avaProfileMass['s'][-1]*0.05)
    alphaLine = avaProfileMass['z'][0] - sExt * mu
    idx = np.argwhere(np.diff(np.sign(zExt - alphaLine))).flatten()
    idx = idx[-1]
    s0 = sExt[idx]
    s1 = sExt[idx+1]
    zP0 = zExt[idx]
    zP1 = zExt[idx+1]
    zA0 = alphaLine[idx]
    zA1 = alphaLine[idx+1]
    sIntersection = s0 + (s1-s0)*(zA0-zP0)/((zP1-zP0)-(zA1-zA0))
    zIntersection = zP0 + (sIntersection-s0) * (zP1-zP0) / (s1-s0)

    # compute errors
    # rmse of the energy height
    rmseVelocityElevation = rmse(zEne, zEneTarget)
    # error on s runout
    runOutSError = avaProfileMass['s'][-1] - sIntersection
    # error on z runout
    runOutZError = avaProfileMass['z'][-1] - zIntersection
    # error on angle runout
    runOutAngleError = runOutAngleDeg - alphaDeg

    ax2 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    # plot mass averaged center of mass
    ax2.plot(avaProfileMass['s'], avaProfileMass['z'], '-k.', label='Center of mass altitude')
    # ax2.plot(avaProfileMass['sCor'], avaProfileMass['z'], '--k.', label='Center of mass altitude (corrected s)')
    # extend this curve towards the bottom using a linear regression on the last x points
    ax2.plot(avaProfileMass['s'][-1]*np.array([1, 1.05]), avaProfileMass['z'][-1] + p1[0]*avaProfileMass['s'][-1]*np.array([0, 0.05]),
             ':k', label='Center of mass altitude extrapolation')

    # add center of mass velocity points and run-out line
    ax2.plot(avaProfileMass['s'][[0, -1]], zEne[[0, -1]], '-r', label='com1dfa energy line (%.2f°)' % runOutAngleDeg)
    scat = ax2.scatter(avaProfileMass['s'], zEne, marker='o', cmap=cmap, s=2*pU.ms, c=V2Path/(2*g),
                       label='Center of mass velocity altitude')
    cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('Center of mass velocity altitude [m]')

    # add alpha line
    ax2.plot(s_geomL, z_geomL, 'b-', label=r'$\alpha$ line (%.2f°)' % alphaDeg)

    ax2.set_xlabel('s [m]')
    ax2.set_ylabel('z [m]')
    ax2.legend(loc='upper right')
    ax2.set_title('Energy line test')

    ax3 = plt.subplot2grid((2, 2), (1, 1))
    # same plot as on ax2 but zoomed in
    # plot mass averaged center of mass
    ax3.plot(avaProfileMass['s'], avaProfileMass['z'], '-k.', label='Center of mass altitude')
    # extend this curve towards the bottom using a linear gegression on the last x points
    ax3.plot(avaProfileMass['s'][-1]*np.array([1, 1.05]), avaProfileMass['z'][-1] + p1[0]*avaProfileMass['s'][-1]*np.array([0, 0.05]),
             ':k', label='Center of mass altitude extrapolation')

    # add center of mass velocity points and run-out line
    ax3.plot(avaProfileMass['s'][[0, -1]], zEne[[0, -1]], '-r', label='com1dfa energy line (%.2f°)' % runOutAngleDeg)
    scat = ax3.scatter(avaProfileMass['s'], zEne, marker='o', cmap=cmap, s=2*pU.ms, c=V2Path/(2*g),
                       label='Center of mass velocity altitude extrapolation')
    cbar3 = ax3.figure.colorbar(scat, ax=ax3, use_gridspec=True)
    cbar3.ax.set_ylabel('Center of mass velocity altitude [m]')

    # add alpha line
    ax3.plot(s_geomL, z_geomL, 'b-', label=r'$\alpha$ line (%.2f°)' % alphaDeg)

    # add horizontal line at the final mass averaged position
    errorS = abs(sIntersection - avaProfileMass['s'][-1])
    ax3.vlines(x=avaProfileMass['s'][-1], ymin=avaProfileMass['z'][-1]-mu*2*errorS, ymax=avaProfileMass['z'][-1], color='r', linestyle='--')
    ax3.vlines(x=sIntersection, color='b', ymin=avaProfileMass['z'][-1]-mu*2*errorS, ymax=zIntersection, linestyle='--')

    ax3.set_xlabel('s [m]')
    ax3.set_ylabel('z [m]')
    # ax3.legend(loc='upper right')
    ax3.set_ylim([avaProfileMass['z'][-1]-mu*2*errorS, avaProfileMass['z'][-1]+mu*2*errorS])
    ax3.set_xlim([avaProfileMass['s'][-1]-2*errorS, avaProfileMass['s'][-1]+2*errorS])
    ax3.get_yaxis().set_ticks([])
    ax3.tick_params(axis='x', which='major', labelsize=8, rotation=45)
    ax3.set_xticks([avaProfileMass['s'][-1], sIntersection])
    text = ('Run-out s diff : %.2f m \nRun-out z diff : %.2f m \nRun-out angle diff : %.4f° \nvelocity height rsme : %.2f m \n(energy line - alpha line)' %
            (runOutSError, runOutZError, runOutAngleError, rmseVelocityElevation))
    text_box = AnchoredText(text, frameon=False, loc=1, pad=0.5)
    plt.setp(text_box.patch, facecolor='white', alpha=0.5)
    ax3.add_artist(text_box)

    plt.show()


def rmse(predictions, targets):
    """ Compute the root mean square error between two numpy 1D arrays"""
    return np.sqrt(((predictions - targets) ** 2).mean())
