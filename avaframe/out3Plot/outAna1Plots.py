# imports
import math
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns

# local imports
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outQuickPlot as outQuickPlot
import avaframe.out3Plot.outDebugPlots as outDebugPlots


def showSaveTimeSteps(cfgMain, cfgSimi, particlesList, fieldsList, solSimi, Tsave, header, outDirTest, simHash, simDFrow):
    """ Generate plots of the comparison of DFA solution and simiSol
    """

    relTh = simDFrow['relTh']
    gravAcc = simDFrow['gravAcc']

    # user interaction?
    if cfgSimi.getboolean('flagInteraction'):
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

            # make plot
            plotProfilesSimiSol(ind_time, simHash, comSol, simiDict, solSimi, axis)

        cellSize = header['cellsize']
        hSimi = simiDict['hSimi']
        hNumerical = fields['FD']
        dataDict = {'name2': 'solAnalytic', 'name1': simHash, 'data2': hSimi, 'data1': hNumerical, 'cellSize': cellSize,
                    'suffix': 'FD', 'compareType': 'compToRef', 'simName': simHash}
        avaName = 'FDComparison' + simHash
        outQuickPlot.generatePlot(dataDict, avaName, outDirTest / 'pics', cfgMain, {'plots': [], 'difference': [], 'stats': [],
                                  'differenceZoom': []}, crossProfile=False)

        vSimi = DFAtls.norm(simiDict['vxSimi'], simiDict['vySimi'], simiDict['vzSimi'])
        vNumerical = DFAtls.norm(fields['Vx'], fields['Vy'], fields['Vz'])
        dataDict = {'name2': 'solAnalytic', 'name1': simHash, 'data2': vSimi, 'data1': vNumerical, 'cellSize': cellSize,
                    'suffix': 'FV', 'compareType': 'compToRef', 'simName': simHash}
        avaName = 'FVComparison' + simHash
        outQuickPlot.generatePlot(dataDict, avaName, outDirTest / 'pics', cfgMain, {'plots': [], 'difference': [], 'stats': [],
                                  'differenceZoom': []}, crossProfile=False)

        # # option for user interaction
        if cfgSimi.getboolean('flagInteraction'):
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

    # similarity solution results
    vSimi = simiDict['vSimi']
    vxSimi = simiDict['vxSimi']
    vySimi = simiDict['vySimi']
    vzSimi = simiDict['vzSimi']
    hSimi = simiDict['hSimi']
    xCenter = simiDict['xCenter']

    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, pU.figH))
    ax2 = ax1.twinx()

    if axis == 'xaxis':
        ax1.axvline(x=xCenter, linestyle=':')
        # DFA simulation
        ax1.plot(xArrayFields, fields['FD'][indFinal, :], 'k', label='Field flow depth')
        ax2.plot(xArrayFields, fields['FV'][indFinal, :], 'g', label='Field flow velocity')
        ax2.plot(xArrayFields, fields['Vx'][indFinal, :], 'm', label='Field x velocity')
        ax2.plot(xArrayFields, fields['Vy'][indFinal, :], 'b', label='Field y velocity')
        ax2.plot(xArrayFields, fields['Vz'][indFinal, :], 'c', label='Field z velocity')
        ax1.plot(x, h, '.k', linestyle='None', label='Part flow depth')
        ax2.plot(x, v, '.g', linestyle='None', label='Part flow velocity')
        # similarity solution
        ax1.plot(xArrayFields, hSimi[indFinal, :], '--k', label='SimiSol flow depth')
        ax2.plot(xArrayFields, vSimi[indFinal, :], '--g', label='SimiSol flow velocity')
        ax2.plot(xArrayFields, vxSimi[indFinal, :], '--m', label='SimiSol x velocity')
        ax2.plot(xArrayFields, vySimi[indFinal, :], '--b', label='SimiSol y velocity')
        ax2.plot(xArrayFields, vzSimi[indFinal, :], '--c', label='SimiSol z velocity')
        ax1.set_title('Profile along flow at t=%.2f (com1DFA), %.2f s (simiSol)' % (Tsave, solSimi['Time'][ind_time]))
        ax1.set_xlabel('x in [m]')
        indStart = min(first_nonzero(hSimi[indFinal, :], 0), first_nonzero(fields['FD'][indFinal, :], 0)) - 2
        indEnd = max(last_nonzero(hSimi[indFinal, :], 0), last_nonzero(fields['FD'][indFinal, :], 0)) + 2
        ax1.set_xlim([xArrayFields[indStart], xArrayFields[indEnd]])

    elif axis == 'yaxis':
        # DFA simulation
        ax1.plot(yArrayFields, fields['FD'][:, indFinal], 'k', label='Field flow depth')
        ax2.plot(yArrayFields, fields['FV'][:, indFinal], 'g', label='Field flow velocity')
        ax2.plot(yArrayFields, fields['Vx'][:, indFinal], 'm', label='Field x velocity')
        ax2.plot(yArrayFields, fields['Vy'][:, indFinal], 'b', label='Field y velocity')
        ax2.plot(yArrayFields, fields['Vz'][:, indFinal], 'c', label='Field z velocity')
        ax1.plot(y, h, '.k', linestyle='None', label='Part flow depth')
        ax2.plot(y, v, '.g', linestyle='None', label='Part flow velocity')
        # similarity solution
        ax1.plot(yArrayFields, hSimi[:, indFinal], '--k', label='SimiSol flow depth')
        ax2.plot(yArrayFields, vSimi[:, indFinal], '--g', label='SimiSol flow velocity')
        ax2.plot(yArrayFields, vxSimi[:, indFinal], '--m', label='SimiSol x velocity')
        ax2.plot(yArrayFields, vySimi[:, indFinal], '--b', label='SimiSol y velocity')
        ax2.plot(yArrayFields, vzSimi[:, indFinal], '--c', label='SimiSol z velocity')
        ax1.set_title('Profile across flow at t=%.2f (com1DFA), %.2f s (simiSol)' % (Tsave, solSimi['Time'][ind_time]))
        ax1.set_xlabel('y in [m]')
        indStart = min(first_nonzero(hSimi[:, indFinal], 0), first_nonzero(fields['FD'][:, indFinal], 0)) - 2
        indEnd = max(last_nonzero(hSimi[:, indFinal], 0), last_nonzero(fields['FD'][:, indFinal], 0)) + 2
        ax1.set_xlim([yArrayFields[indStart], yArrayFields[indEnd]])

    ax1.set_ylabel('flow depth [m]')
    color = 'tab:green'
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel('flow velocity [ms-1]', color=color)
    ax2.legend(loc='upper right')
    ax1.legend(loc='upper left')

    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'profile_' + outputName + '_%sCutSol_T%.2f.' % (axis, Tsave) + pU.outputFormat, fig1)


def plotErrorTime(time, hErrorL2Array, hErrorLMaxArray, vErrorL2Array, vErrorLMaxArray, outDirTest, outputName):
    """plot error between given com1DFA sol and analytic sol
    function of time
    """
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    ax2 = ax1.twinx()
    ax1.plot(time, hErrorL2Array, 'k-', label='Flow depth L2 error')
    ax1.plot(time, hErrorLMaxArray, 'k--', label='Flow depth LMax error')
    ax1.set_title('Error between similarity solution and com1DFA')
    ax1.set_xlabel('time in [s]')
    ax1.set_ylabel('error on flow depth')
    ax1.legend(loc='upper left')

    color = 'tab:green'
    ax2.plot(time, vErrorL2Array, 'g-', label='Velocity L2 error')
    ax2.plot(time, vErrorLMaxArray, 'g--', label='Velocity LMax error')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel('error on velocity', color=color)
    ax2.legend(loc='lower right')

    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'Error_Time_' + outputName, fig1)


def plotError(simDF, outDirTest):
    """plot error between all com1DFA sol and analytic sol
    function of nParts for all dt
    """
    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, 2*pU.figH))
    ax2 = ax1.twinx()
    ax1 = sns.pointplot(x='Npart', y='hErrorL2', hue='dt', data=simDF, ax=ax1, markers=['o', 's', 'd', 'v', '^', '<', '>'], palette=['k'])
    ax2 = sns.pointplot(x='Npart', y='vErrorL2', hue='dt', data=simDF, ax=ax2, markers=['o', 's', 'd', 'v', '^', '<', '>'], palette=['g'])
    ax1 = sns.pointplot(x='Npart', y='hErrorLMax', hue='dt', data=simDF, ax=ax1, linestyles='--', markers=['o', 's', 'd', 'v', '^', '<', '>'], palette=['k'])
    ax2 = sns.pointplot(x='Npart', y='vErrorLMax', hue='dt', data=simDF, ax=ax2, linestyles='--', markers=['o', 's', 'd', 'v', '^', '<', '>'], palette=['g'])
    ax1.set_title('Error between similarity solution and com1DFA')
    ax1.set_xlabel('Number of particles')

    ax1.set_ylabel('error on flow depth')
    color = 'tab:green'
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel('error on velocity', color=color)
    ax2.legend(loc='lower right')
    ax1.legend(loc='upper left')
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'Error', fig1)


def plotErrorLog(simDF, outDirTest):
    # fig1, ax1 = plt.subplots(figsize=(2*pU.figW, 2*pU.figH))
    # ax2 = ax1.twinx()
    # g = sns.relplot(data=simDF, x="Npart", y="hErrorL2", hue="sphKernelRadius", style="dt",  # hue="dt",
    #                 markers=['o', 's'], palette=['k', 'b', 'r'])
    # g.set(xscale="log", yscale="log")
    # g.ax.xaxis.grid(True, "minor", linewidth=.25)
    # g.ax.yaxis.grid(True, "minor", linewidth=.25)
    # # g.despine(left=True, bottom=True)

    ax = simDF.plot.scatter(x="Npart", y="hErrorL2", c="sphKernelRadius", cmap="viridis", s=simDF["dt"] * 200)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'ErrorLog', g)


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
