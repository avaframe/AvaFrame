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

    # read peak field
    simName = simDF.index[0]
    fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['pfv'], simName=simName, flagAvaDir=True,
                                                           comModule='com1DFA')

    generateEnergyPlot(avalancheDir, energyLineTestCfg, avaProfileMass, dem, particlesList, fieldsList, simName)


def generateEnergyPlot(avalancheDir, energyLineTestCfg, avaProfileMass, dem, particlesList, fieldsList, simName):
    """ Make energy test analysis and plot results

    Parameters
    -----------
    avalancheDir: pathlib
        avalanche directory pathlib path
    energyLineTestCfg: configParser
        com1DFA energy test config
    avaProfileMass: dict
        particle mass average properties
    dem: dict
        com1DFA simulation dictionnary
    particlesList: list
        particles dictionary list
    fieldsList: list
        field dictionary list
    simName: str
        sime name
    """
    # read pysical parameters from DFA configuration
    g = energyLineTestCfg['GENERAL'].getfloat('gravAcc')
    mu = energyLineTestCfg['GENERAL'].getfloat('mu')
    alphaRad = np.arctan(mu)
    alphaDeg = np.rad2deg(alphaRad)
    runOutAngleRad, runOutAngleDeg = getRunOutAngle(avaProfileMass)
    slopeExt, sIntersection, zIntersection = getAlphaProfileIntersection(avaProfileMass, mu)
    zEne, V2Path, s_geomL, z_geomL, errorEnergyTest = getEnergyInfo(avaProfileMass, g, mu, sIntersection, zIntersection,
                                                                    runOutAngleDeg, alphaDeg)
    z0 = avaProfileMass['z'][0]
    # create figures and plots
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))
    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps['pfv'], np.amin(V2Path/(2*g)), np.amax(V2Path/(2*g)),
                                           continuous=pU.contCmap)
    cmap.set_under(color='w')
    # make the bird view plot
    ax1 = plt.subplot2grid((2, 2), (1, 0))
    ax1, extent, cbar0, cs1 = outCom1DFA.addResult2Plot(ax1, dem, fieldsList[-1]['pfv'], 'pfv')
    cbar0.ax.set_ylabel('peak flow velocity')
    ax1 = outCom1DFA.addDem2Plot(ax1, dem, what='slope', extent=extent)
    ax1, cbar1 = outCom1DFA.addParticles2Plot(particlesList[-1], ax1, dem, whatS='h')
    cbar1.ax.set_ylabel('final particle flow thickness')
    ax1.plot(avaProfileMass['x'], avaProfileMass['y'], '-k.', zorder = 20, label='Center of mass path')
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(fieldsList[-1]['pfv'], 5, extentOption=True,
                                                                 constrainedData=False, buffer='')
    ax1.set_ylim([rowsMin, rowsMax])
    ax1.set_xlim([colsMin, colsMax])
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    l1 = ax1.legend(loc='upper right')
    l1.set_zorder(40)
    pU.putAvaNameOnPlot(ax1, avalancheDir)

    # make profile plot, zoomed out
    ax2 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    # plot mass averaged center of mass
    ax2.plot(avaProfileMass['s'], avaProfileMass['z'], '-k.', label='Center of mass altitude')
    # ax2.plot(avaProfileMass['sCor'], avaProfileMass['z'], '--k.', label='Center of mass altitude (corrected s)')
    # extend this curve towards the bottom using a linear regression on the last x points
    ax2.plot(avaProfileMass['s'][-1]*np.array([1, 1.05]), avaProfileMass['z'][-1] +
             slopeExt*avaProfileMass['s'][-1]*np.array([0, 0.05]), ':k', label='Center of mass altitude extrapolation')

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

    # make profile plot, zoomed in
    ax3 = plt.subplot2grid((2, 2), (1, 1))
    # plot mass averaged center of mass
    ax3.plot(avaProfileMass['s'], avaProfileMass['z']-z0, '-k.', label='Center of mass altitude')
    # ax3.plot(avaProfileMass['sCor'], avaProfileMass['z'], '--k.', label='Center of mass altitude (corrected s)')
    # extend this curve towards the bottom using a linear gegression on the last x points
    ax3.plot(avaProfileMass['s'][-1]*np.array([1, 1.05]), avaProfileMass['z'][-1] +
             slopeExt*avaProfileMass['s'][-1]*np.array([0, 0.05])-z0, ':k',
             label='Center of mass altitude extrapolation')

    # add center of mass velocity points and run-out line
    ax3.plot(avaProfileMass['s'][[0, -1]], zEne[[0, -1]]-z0, '-r',
             label='com1dfa energy line (%.2f°)' % runOutAngleDeg)
    scat = ax3.scatter(avaProfileMass['s'], zEne-z0, marker='o', cmap=cmap, s=2*pU.ms, c=V2Path/(2*g),
                       label='Center of mass velocity altitude extrapolation')
    cbar3 = ax3.figure.colorbar(scat, ax=ax3, use_gridspec=True)
    cbar3.ax.set_ylabel('Center of mass velocity altitude [m]')

    # add alpha line
    ax3.plot(s_geomL, z_geomL-z0, 'b-', label=r'$\alpha$ line (%.2f°)' % alphaDeg)

    # add horizontal line at the final mass averaged position
    # compute plot limits
    runOutSError = errorEnergyTest['runOutSError']
    runOutZError = errorEnergyTest['runOutZError']
    runOutAngleError = errorEnergyTest['runOutAngleError']
    rmseVelocityElevation = errorEnergyTest['rmseVelocityElevation']
    errorS = abs(runOutSError)
    sMin = min(avaProfileMass['s'][-1], sIntersection) - errorS
    sMax = max(avaProfileMass['s'][-1], sIntersection) + errorS
    zMin = avaProfileMass['z'][-1] + slopeExt*(sMax - avaProfileMass['s'][-1])-z0
    zMax = avaProfileMass['z'][0] - sMin*np.tan(min(runOutAngleRad, alphaRad))-z0
    ax3.vlines(x=avaProfileMass['s'][-1], ymin=zMin, ymax=avaProfileMass['z'][-1]-z0,
               color='r', linestyle='--')
    ax3.vlines(x=sIntersection, color='b', ymin=zMin, ymax=zIntersection-z0, linestyle='--')

    ax3.hlines(y=avaProfileMass['z'][-1]-z0, xmin=sMin, xmax=avaProfileMass['s'][-1],
               color='r', linestyle='--')
    ax3.hlines(y=zIntersection-z0, color='b', xmin=sMin, xmax=sIntersection, linestyle='--')

    ax3.set_xlabel('s [m]')
    ax3.set_ylabel(r'$\Delta z$ [m]')
    ax3.yaxis.set_label_coords(0, 0.9)
    # ax3.legend(loc='upper right')
    ax3.set_xlim([sMin, sMax])
    ax3.set_ylim([zMin, zMax])
    ax3.tick_params(axis='x', which='major', labelsize=8, rotation=45)
    ax3.tick_params(axis='y', which='major', labelsize=8, rotation=45)
    ax3.set_xticks([avaProfileMass['s'][-1], sIntersection])
    ax3.set_yticks([avaProfileMass['z'][-1]-z0, zIntersection-z0])
    text = ('Run-out s diff : %.2f m \nRun-out z diff : %.2f m \nRun-out angle diff : %.4f° \nvelocity height rmse : %.2f m \n(energy line - alpha line)' %
            (runOutSError, runOutZError, runOutAngleError, rmseVelocityElevation))
    text_box = AnchoredText(text, frameon=False, loc=1, pad=0.5)
    plt.setp(text_box.patch, facecolor='white', alpha=0.5)
    ax3.add_artist(text_box)
    outFileName = '_'.join([simName, 'EnergyTest'])
    outDir = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests')
    pU.saveAndOrPlot({'pathResult': outDir}, outFileName, fig)


def getRunOutAngle(avaProfileMass):
    """Compute the Center of Mass runout angle

    Parameters
    -----------
    avaProfileMass: dict
        particle mass average properties

    Returns
    --------
    runOutAngleRad: float
        Center of Mass runout angle in radians
    runOutAngleDeg: float
        Center of Mass runout angle in degrees
    """
    # compute horizontal and vertical traveled distance
    deltaS = avaProfileMass['s'][-1] - avaProfileMass['s'][0]
    deltaZ = avaProfileMass['z'][-1] - avaProfileMass['z'][0]
    # extract simulation run-out
    runOutAngleRad = np.arctan(np.abs(deltaZ/deltaS))
    runOutAngleDeg = np.rad2deg(runOutAngleRad)
    return runOutAngleRad, runOutAngleDeg


def getAlphaProfileIntersection(avaProfileMass, mu):
    """Extend the  profile path and compute the intersection
    between the theoretical energy line and the path profile

    Parameters
    -----------
    avaProfileMass: dict
        particle mass average properties
    mu: float
        friction coefficient

    Returns
    --------
    slopeExt: float
        slope of the extrapolation line
    sIntersection: float
        s coord of the intersection betwee the line of slope mu and the
        mass average path profile
    zIntersection: float
        z coord of the intersection betwee the line of slope mu and the
        mass average path profile
    """
    # get slope of the last profile foints to extend the profile
    p1 = np.polyfit(avaProfileMass['s'][-3:], avaProfileMass['z'][-3:], 1)
    slopeExt = p1[0]
    # extend profile
    sExt = np.append(avaProfileMass['s'], 1.05*avaProfileMass['s'][-1])
    zExt = np.append(avaProfileMass['z'], avaProfileMass['z'][-1] + slopeExt*avaProfileMass['s'][-1]*0.05)
    # find intersection between alpha line and profile
    alphaLine = avaProfileMass['z'][0] - sExt * mu
    # find the intersection segment
    idx = np.argwhere(np.diff(np.sign(zExt - alphaLine))).flatten()
    idx = idx[-1]
    # find the exact intersection point
    s0 = sExt[idx]
    s1 = sExt[idx+1]
    zP0 = zExt[idx]
    zP1 = zExt[idx+1]
    zA0 = alphaLine[idx]
    zA1 = alphaLine[idx+1]
    sIntersection = s0 + (s1-s0)*(zA0-zP0)/((zP1-zP0)-(zA1-zA0))
    zIntersection = zP0 + (sIntersection-s0) * (zP1-zP0) / (s1-s0)
    return slopeExt, sIntersection, zIntersection


def getEnergyInfo(avaProfileMass, g, mu, sIntersection, zIntersection, runOutAngleDeg, alphaDeg):
    """Compute energy dots and errors

    Parameters
    -----------
    avaProfileMass: dict
        particle mass average properties
    g: float
        gravity
    mu: float
        friction coefficient
    sIntersection: float
        s coord of the intersection betwee the line of slope mu and the
        mass average path profile
    zIntersection: float
        z coord of the intersection betwee the line of slope mu and the
        mass average path profile
    runOutAngleRad: float
        Center of Mass runout angle in radians
    runOutAngleDeg: float
        Center of Mass runout angle in degrees

    Returns
    --------
    zEne: numpy 1D array
        energy height of the particle averaged time steps
    V2Path: numpy 1D array
        kinetic energy of the particle averaged time steps
    s_geomL: 2 element list
        s coord (start and end) of the run out angle line
    z_geomL: 2 element list
        z coord (start and end) of the run out angle line
    errorEnergyTest: dict
        rmseVelocityElevation, runOutSError, runOutZError, runOutAngleError
        between theoretical solution and simulation result
    """
    # read mass average quantities
    V2Path = avaProfileMass['v2']
    # compute mass average velocity elevation
    # extract energy altitude
    zEne = avaProfileMass['z'] + V2Path/(2*g)

    # Create the alpha line
    GK = sIntersection * mu
    z_end = avaProfileMass['z'][0] - GK
    zEneTarget = avaProfileMass['z'][0] - avaProfileMass['s'] * mu
    s_geomL = [avaProfileMass['s'][0], sIntersection]
    z_geomL = [avaProfileMass['z'][0], z_end]
    # compute errors
    # rmse of the energy height
    rmseVelocityElevation = rmse(zEne, zEneTarget)
    # error on s runout
    runOutSError = avaProfileMass['s'][-1] - sIntersection
    # error on z runout
    runOutZError = avaProfileMass['z'][-1] - zIntersection
    # error on angle runout
    runOutAngleError = runOutAngleDeg - alphaDeg
    errorEnergyTest = {'rmseVelocityElevation': rmseVelocityElevation, 'runOutSError': runOutSError,
                       'runOutZError': runOutZError, 'runOutAngleError': runOutAngleError}
    return zEne, V2Path, s_geomL, z_geomL, errorEnergyTest


def rmse(predictions, targets):
    """ Compute the root mean square error between two numpy 1D arrays

    Parameters
    -----------
    predictions: 1D numpy array
        prediction array
    targets: 1D numpy array
        targets array

    Returns
    --------
    rmse: float
        rmse value
    """
    return np.sqrt(((predictions - targets) ** 2).mean())
