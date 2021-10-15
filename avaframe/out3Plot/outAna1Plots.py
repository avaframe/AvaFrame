# imports
import math
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import logging

# local imports
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outQuickPlot as outQuickPlot
import avaframe.out3Plot.outDebugPlots as outDebugPlots

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def showSaveTimeSteps(cfgMain, cfgSimi, particlesList, fieldsList, solSimi, Tsave, header, outDirTest, simHash, simDFrow):
    """ Generate plots of the comparison of DFA solution and simiSol
    """

    relTh = simDFrow['relTh']
    gravAcc = simDFrow['gravAcc']

    # user interaction?
    if cfgSimi.getboolean('interaction'):
        value = input("give time step to plot (float in s):\n")
    else:
        value = cfgSimi.getfloat('tSave')

    try:
        value = float(value)
    except ValueError:
        value = 'n'
    while isinstance(value, float):

        # determine index for time step
        ind_t = min(np.searchsorted(Tsave, value), min(len(Tsave)-1, len(fieldsList)-1))
        ind_time = np.absolute(solSimi['Time'] - Tsave[ind_t]).argmin()
        # get similartiy solution h, u at reuired time step
        simiDict = simiSolTest.getSimiSolParameters(solSimi, header, ind_time, cfgSimi, relTh, gravAcc)

        # get DFA simulation
        # load fields and particles of required time step described by ind_t
        fields = fieldsList[ind_t]
        particles = particlesList[ind_t]
        for axis in ['xaxis', 'yaxis']:
            # get particle parameters
            comSol = simiSolTest.prepareParticlesFieldscom1DFA(fields, particles, header, simiDict, axis)
            comSol['outDirTest'] = outDirTest
            comSol['showPlot'] = cfgMain['FLAGS'].getboolean('showPlot')
            comSol['Tsave'] = Tsave[ind_t]
            comSol['dt'] = simDFrow['dt']
            comSol['deltaTh'] = simDFrow['deltaTh']
            comSol['sphKernelRadius'] = simDFrow['sphKernelRadius']

            # make plot
            plotProfilesSimiSol(ind_time, simHash, comSol, simiDict, solSimi, axis)

        # create flow depth raster difference plot
        cellSize = header['cellsize']
        hSimi = simiDict['hSimi']
        hNumerical = fields['FD']
        dataDict = {'name2': 'solAnalytic', 'name1': simHash, 'data2': hSimi, 'data1': hNumerical, 'cellSize': cellSize,
                    'suffix': 'FD', 'compareType': 'compToRef', 'simName': simHash}
        avaName = 'FDComparison' + simHash
        outQuickPlot.generatePlot(dataDict, avaName, outDirTest / 'pics', cfgMain, {'plots': [], 'difference': [], 'stats': [],
                                  'differenceZoom': []}, crossProfile=False)

        # create flow momentum per surface unit raster difference plot
        vhSimi = simiDict['vSimi'] * hSimi
        vhNumerical = fields['FV'] * hNumerical
        dataDict = {'name2': 'solAnalytic', 'name1': simHash, 'data2': vhSimi, 'data1': vhNumerical, 'cellSize': cellSize,
                    'suffix': 'FDV', 'compareType': 'compToRef', 'simName': simHash}
        avaName = 'FDVComparison' + simHash
        outQuickPlot.generatePlot(dataDict, avaName, outDirTest / 'pics', cfgMain, {'plots': [], 'difference': [], 'stats': [],
                                  'differenceZoom': []}, crossProfile=False)

        # # option for user interaction
        if cfgSimi.getboolean('interaction'):
            value = input("give time step to plot (float in s):\n")
            try:
                value = float(value)
            except ValueError:
                value = 'n'
        else:
            value = 'n'


def plotProfilesSimiSol(ind_time, outputName, comSol, simiDict, solSimi, axis):
    """ Plot flow depth and velocity for similarity solution and simulation results

        Parameters
        -----------
        ind_time: int
            time index for simiSol
        outputName: str
            outputName
        comSol: dict
            dictionary of simulation results and info (particles, fields, indices, time step)
        simiDict: dict
            dictionary with similiarty solution for h, u, and xCenter at required time step
        solSimi: dict
            dictionary with similiarty solution
        axis: str

    """

    # com1DFA results
    fields = comSol['fields']
    xArrayFields = comSol['xArrayFields']
    yArrayFields = comSol['yArrayFields']
    # particle properties
    x = comSol['x']
    y = comSol['y']
    h = comSol['h']
    v = comSol['v']
    vx = comSol['vx']
    vy = comSol['vy']
    vz = comSol['vz']
    outDirTest = comSol['outDirTest']
    indFinal = comSol['indFinal']
    Tsave = comSol['Tsave']
    dt = comSol['dt']
    deltaTh = comSol['deltaTh']
    sphKernelRadius = comSol['sphKernelRadius']

    # similarity solution results
    vSimi = simiDict['vSimi']
    vxSimi = simiDict['vxSimi']
    vySimi = simiDict['vySimi']
    vzSimi = simiDict['vzSimi']
    hSimi = simiDict['hSimi']
    xCenter = simiDict['xCenter']
    Time = solSimi['Time']

    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, pU.figH))
    ax2 = ax1.twinx()

    if axis == 'xaxis':
        ax1.axvline(x=xCenter, linestyle=':')
        FD = fields['FD'][indFinal, :]
        hSimi = hSimi[indFinal, :]
        # DFA simulation
        ax1.plot(xArrayFields, FD, 'k', label='Field flow depth')
        ax2.plot(xArrayFields, FD*fields['FV'][indFinal, :], 'g', label='Field flow velocity')
        ax2.plot(xArrayFields, FD*fields['Vx'][indFinal, :], 'm', label='Field x velocity')
        ax2.plot(xArrayFields, FD*fields['Vy'][indFinal, :], 'b', label='Field y velocity')
        ax2.plot(xArrayFields, FD*fields['Vz'][indFinal, :], 'c', label='Field z velocity')
        ax1.plot(x, h, '.k', linestyle='None', label='Part flow depth')
        ax2.plot(x, h*v, '.g', linestyle='None', label='Part flow velocity')
        # similarity solution
        ax1.plot(xArrayFields, hSimi, '--k', label='SimiSol flow depth')
        ax2.plot(xArrayFields, hSimi*vSimi[indFinal, :], '--g', label='SimiSol flow velocity')
        ax2.plot(xArrayFields, hSimi*vxSimi[indFinal, :], '--m', label='SimiSol x velocity')
        ax2.plot(xArrayFields, hSimi*vySimi[indFinal, :], '--b', label='SimiSol y velocity')
        ax2.plot(xArrayFields, hSimi*vzSimi[indFinal, :], '--c', label='SimiSol z velocity')
        ax1.set_title('Profile along flow at t=%.2f (com1DFA), %.2f s (simiSol) (csz = %s m, dt = %s s, deltaTh = %s m)'
                      % (Tsave, Time[ind_time], sphKernelRadius, dt, deltaTh))
        ax1.set_xlabel('x in [m]')
        indStart = min(first_nonzero(hSimi, 0), first_nonzero(FD, 0)) - 2
        indEnd = max(last_nonzero(hSimi, 0), last_nonzero(FD, 0)) + 2
        ax1.set_xlim([xArrayFields[indStart], xArrayFields[indEnd]])

    elif axis == 'yaxis':
        FD = fields['FD'][:, indFinal]
        hSimi = hSimi[:, indFinal]
        # DFA simulation
        ax1.plot(yArrayFields, FD, 'k', label='Field flow depth')
        ax2.plot(yArrayFields, FD*fields['FV'][:, indFinal], 'g', label='Field flow velocity')
        ax2.plot(yArrayFields, FD*fields['Vx'][:, indFinal], 'm', label='Field x velocity')
        ax2.plot(yArrayFields, FD*fields['Vy'][:, indFinal], 'b', label='Field y velocity')
        ax2.plot(yArrayFields, FD*fields['Vz'][:, indFinal], 'c', label='Field z velocity')
        ax1.plot(y, h, '.k', linestyle='None', label='Part flow depth')
        ax2.plot(y, h*v, '.g', linestyle='None', label='Part flow velocity')
        # similarity solution
        ax1.plot(yArrayFields, hSimi, '--k', label='SimiSol flow depth')
        ax2.plot(yArrayFields, hSimi*vSimi[:, indFinal], '--g', label='SimiSol flow velocity')
        ax2.plot(yArrayFields, hSimi*vxSimi[:, indFinal], '--m', label='SimiSol x velocity')
        ax2.plot(yArrayFields, hSimi*vySimi[:, indFinal], '--b', label='SimiSol y velocity')
        ax2.plot(yArrayFields, hSimi*vzSimi[:, indFinal], '--c', label='SimiSol z velocity')
        ax1.set_title('Profile across flow at t=%.2f (com1DFA), %.2f s (simiSol) (csz = %s m, dt = %s s, deltaTh = %s m)'
                      % (Tsave, Time[ind_time], sphKernelRadius, dt, deltaTh))
        ax1.set_xlabel('y in [m]')
        indStart = min(first_nonzero(hSimi, 0), first_nonzero(FD, 0)) - 2
        indEnd = max(last_nonzero(hSimi, 0), last_nonzero(FD, 0)) + 2
        ax1.set_xlim([yArrayFields[indStart], yArrayFields[indEnd]])

    ax1.set_ylabel('flow depth [m]')
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    color = 'tab:green'
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.grid(color='tab:green', linestyle='-', linewidth=0.25, alpha=0.5)
    ax2.set_ylabel('flow velocity [ms-1]', color=color)
    ax2.legend(loc='upper right')
    ax1.legend(loc='upper left')

    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'profile_' + str(outputName) + '_%sCutSol_T%.2f.' % (axis, Tsave) + pU.outputFormat, fig1)


def plotErrorTime(time, hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray, outDirTest, outputName,
                  simDFrow, relativ):
    """plot error between given com1DFA sol and analytic sol
    function of time
    """
    dt = simDFrow['dt']
    deltaTh = simDFrow['deltaTh']
    sphKernelRadius = simDFrow['sphKernelRadius']
    title = (' between similarity solution and com1DFA \n(csz = %s m, dt = %s s, deltaTh = %s m)'
                     % (sphKernelRadius, dt, deltaTh))
    title = getTitleError(relativ, title)
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    ax2 = ax1.twinx()
    ax1.plot(time, hErrorL2Array, 'k-', label='Flow depth L2 error')
    ax1.plot(time, hErrorLMaxArray, 'k--', label='Flow depth LMax error')
    ax1.set_title(title)
    ax1.set_xlabel('time in [s]')
    ax1.set_ylabel(getTitleError(relativ, ' on flow depth'))
    ax1.legend(loc='upper left')
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

    color = 'tab:green'
    ax2.plot(time, vhErrorL2Array, 'g-', label='Momentum L2 error')
    ax2.plot(time, vhErrorLMaxArray, 'g--', label='Momentum LMax error')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel(getTitleError(relativ, ' on momentum'), color=color)
    ax2.legend(loc='lower right')
    ax2.grid(color='tab:green', linestyle='-', linewidth=0.25, alpha=0.5)

    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'Error_Time_' + str(outputName), fig1)


def plotError(simDF, outDirTest, cfgSimi):
    """plot error between all com1DFA sol and analytic sol
    function of nParts for all dt
    """
    relativ = cfgSimi.getboolean('relativError')
    title = getTitleError(relativ, ' between similarity solution and com1DFA')
    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, 2*pU.figH))
    ax2 = ax1.twinx()
    ax1 = sns.pointplot(x='Npart', y='hErrorL2', hue='dt', data=simDF, ax=ax1, markers=['o', 's', 'd', 'v', '^', '<', '>', '.', '+', '*'], palette=['k'])
    ax2 = sns.pointplot(x='Npart', y='vhErrorL2', hue='dt', data=simDF, ax=ax2, markers=['o', 's', 'd', 'v', '^', '<', '>', '.', '+', '*'], palette=['g'])
    ax1 = sns.pointplot(x='Npart', y='hErrorLMax', hue='dt', data=simDF, ax=ax1, linestyles='--', markers=['o', 's', 'd', 'v', '^', '<', '>', '.', '+', '*'], palette=['k'])
    ax2 = sns.pointplot(x='Npart', y='vhErrorLMax', hue='dt', data=simDF, ax=ax2, linestyles='--', markers=['o', 's', 'd', 'v', '^', '<', '>', '.', '+', '*'], palette=['g'])
    ax1.set_title(title)
    ax1.set_xlabel('Number of particles')

    ax1.set_ylabel(getTitleError(relativ, ' on flow depth'))
    color = 'tab:green'
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel(getTitleError(relativ, ' on momentum'), color=color)
    ax2.legend(loc='lower right')
    ax1.legend(loc='upper left')
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'Error', fig1)


def plotErrorLog(simDF, outDirTest, cfgSimi):
    """plot error between all com1DFA sol and analytic sol
    function of nParts for all dt
    """
    sphKernelRadiusList = simDF["sphKernelRadius"].unique()
    dt = simDF["dt"].unique()[0]
    tSave = cfgSimi.getfloat('tSave')
    relativ = cfgSimi.getboolean('relativError')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, min(simDF["sphKernelRadius"])*0.25, max(simDF["sphKernelRadius"])*2, continuous=pU.contCmap)
    cmap = 'viridis'
    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, 2*pU.figH))
    ax2 = ax1.twinx()
    scatter = ax1.scatter(simDF["Npart"], simDF["hErrorL2"], c=simDF["sphKernelRadius"], s=simDF["dt"]*200, cmap=cmap, marker='o', alpha=1, edgecolors='k')
    scatte2 = ax2.scatter(simDF["Npart"], simDF["vhErrorL2"], c=simDF["sphKernelRadius"], s=simDF["dt"]*200, cmap=cmap, marker='s', alpha=0.8, edgecolors='k')
    for sphKernelRadius in sphKernelRadiusList:
        simDFNew = simDF[(simDF['sphKernelRadius'] == sphKernelRadius) & (simDF['dt'] == dt)]
        Npart = simDFNew["Npart"]
        hErrorL2 = simDFNew["hErrorL2"]
        vErrorL2 = simDFNew["vhErrorL2"]
        p = np.polyfit(np.log(simDFNew["Npart"]), np.log(hErrorL2), deg=1)
        p1H = p[0]
        p0H = np.exp(p[1])
        p = np.polyfit(np.log(simDFNew["Npart"]), np.log(vErrorL2), deg=1)
        p1U = p[0]
        p0U = np.exp(p[1])
        ax1.plot(Npart, p0H*Npart**p1H, 'r')
        ax2.plot(Npart, p0U*Npart**p1U, 'g')
        log.info('power law fit sphKernelRadius = %.2f m: hErrorL2 = %.1f * Npart^{%.2f}' % (sphKernelRadius, p0H, p1H))
        log.info('power law fit sphKernelRadius = %.2f m: vhErrorL2 = %.1f * Npart^{%.2f}' % (sphKernelRadius, p0U, p1U))
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_title('Convergence of DFA simulation for the similarity solution test at t = %.2fs' % tSave)
    ax1.set_xlabel('number of particles')
    ax1.set_ylabel(getTitleError(relativ, r' L2 on flow depth ($\bullet$)'))
    ax2.set_ylabel(getTitleError(relativ, r' L2 on momentum ($\blacksquare$)'))
    legend1 = ax1.legend(*scatter.legend_elements(), loc="lower left", title="sphKernelRadius")
    ax1.add_artist(legend1)

    # produce a legend with a cross section of sizes from the scatter
    kw = dict(prop="sizes", color=scatter.cmap(0.7),
          func=lambda s: s/200)
    legend2 = ax1.legend(*scatter.legend_elements(**kw), loc="upper right", title="dt")
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    ax1.grid(color='grey', which='minor', linestyle='--', linewidth=0.25, alpha=0.5)
    b1, t1 = ax1.get_ylim()
    b2, t2 = ax2.get_ylim()
    ax1.set_ylim([min(b1, b2), max(t1, t2)])
    ax2.set_ylim([min(b1, b2), max(t1, t2)])
    plt.show()
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'ErrorLog%ds' % int(tSave), fig1)


def plotContoursSimiSol(particlesList, fieldsList, solSimi, relDict, cfgSimi, Hini, outDirTest):
    """ Make a contour plot of flow depth for analytical solution and simulation result """

    # load parameters
    bedFrictionAngleDeg = cfgSimi.getfloat('bedFrictionAngle')
    planeinclinationAngleDeg = cfgSimi.getfloat('planeinclinationAngle')
    L_x = cfgSimi.getfloat('L_x')
    L_y = cfgSimi.getfloat('L_y')
    X1 = relDict['X1']
    Y1 = relDict['Y1']
    X = relDict['X']
    Y = relDict['Y']
    dem = relDict['dem']
    demOri = relDict['demOri']
    dem['header']['xllcenter'] = demOri['header']['xllcenter']
    dem['header']['yllcenter'] = demOri['header']['yllcenter']

    # Set parameters
    Pi = math.pi
    zeta = planeinclinationAngleDeg * Pi /180       # plane inclination
    delta = bedFrictionAngleDeg * Pi /180           # basal angle of friction
    # A-C
    A = np.sin(zeta)
    C = np.cos(zeta) * np.tan(delta)
    AminusC = A - C
    # make plot
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    for part, field in zip(particlesList, fieldsList):
        t = part['t']
        ind_time = np.searchsorted(solSimi['Time'], t)
        hSimi = simiSolTest.computeH(solSimi, X1, Y1, ind_time, L_y, L_x, Hini, AminusC)
        hSimi = np.where(hSimi <= 0, 0, hSimi)
        fig, ax, cmap, lev = outDebugPlots.plotContours(
            fig, ax, part, dem, field['FD'], pU.cmapDepth, 'm')
        CS = ax.contour(X, Y, hSimi, levels=lev, origin='lower', cmap=cmap,
                        linewidths=2, linestyles='dashed')
        plt.pause(1)
        fig.savefig(pathlib.Path(outDirTest, 'ContourSimiSol%f.%s' % (t, pU.outputFormat)))

    fig, ax, cmap, lev = outDebugPlots.plotContours(
        fig, ax, part, dem, field['FD'], pU.cmapDepth, 'm', last=True)
    CS = ax.contour(X, Y, hSimi, levels=lev, origin='lower', cmap=cmap,
                    linewidths=2, linestyles='dashed')
    ax.clabel(CS, inline=1, fontsize=8)
    fig.savefig(pathlib.Path(outDirTest, 'ContourSimiSolFinal.%s' % (pU.outputFormat)))


def getTitleError(relativ, ending=''):
    """Get error plot title (relativ error or not?)"""
    if relativ:
        return 'Relativ error' + ending
    else:
        return 'Error between' + ending


def last_nonzero(arr, axis, invalid_val=-1):
    """Get index of last non zero value
        Parameters
        -----------
        arr: numpy array
            data array
        axis: int
            axis along which you want to get the index

        Returns
        --------
        index of last non zero in axis direction

    """
    mask = arr != 0
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)


def first_nonzero(arr, axis, invalid_val=-1):
    """Get index of first non zero value
        Parameters
        -----------
        arr: numpy array
            data array
        axis: int
            axis along which you want to get the index

        Returns
        --------
        index of first non zero in axis direction
    """
    mask = arr != 0
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)
