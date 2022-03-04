# imports
import math
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import pandas as pd
import logging
import copy

# local imports
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.outQuickPlot as outQuickPlot
import avaframe.out3Plot.outDebugPlots as outDebugPlots

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def showSaveTimeStepsSimiSol(cfgMain, cfgSimi, fieldsList, solSimi, timeList, header, outDirTest, simHash, simDFrow):
    """ Generate plots of the comparison of DFA solution and simiSol
    Parameters
    -----------
    cfgMain: configparser
        main cfg
    cfgSimi: configparser
        simiSol cfg
    fieldsList: list
        list of fields Dictionaries
    solSimi: dict
        dictionary with similiarty solution
    timeList: list
        time coresponding to fields in fieldsList
    header: dict
        field header dictionary
    outDirTest: str or pathlib
        output directory
    simHash: str
        hash of the curent simulation
    simDFrow: dataFrame
        row of the simDF dataFrame coresponding to simHash
    """

    relTh = cfgSimi.getfloat('relTh')
    gravAcc = simDFrow['gravAcc']

    # user interaction?
    if cfgSimi.getboolean('interaction'):
        value = input("give time step to plot (float in s):\n")
    else:
        value = cfgSimi['tSave']
        if value == '':
            value = timeList[-1]

    try:
        value = float(value)
    except ValueError:
        value = 'n'
    while isinstance(value, float):

        # determine index for time step
        ind_t = min(np.searchsorted(timeList, value), min(len(timeList)-1, len(fieldsList)-1))
        ind_time = np.absolute(solSimi['Time'] - timeList[ind_t]).argmin()
        # get similartiy solution h, u at reuired time step
        simiDict = simiSolTest.getSimiSolParameters(solSimi, header, ind_time, cfgSimi, relTh, gravAcc)

        # get DFA simulation
        # load fields and particles of required time step described by ind_t
        fields = fieldsList[ind_t]
        for axis in ['xaxis', 'yaxis']:
            # get particle parameters
            xCenter = simiDict['xCenter']
            # get info on DEM extent
            ncols = header['ncols']
            nrows = header['nrows']
            xllc = header['xllcenter']
            yllc = header['yllcenter']
            csz = header['cellsize']

            xArrayFields = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
            yArrayFields = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
            if axis == 'xaxis':
                indFinal = int(nrows * 0.5) -1
            elif axis == 'yaxis':
                indFinal = int(np.round((xCenter - xllc)/csz) + 1)

            comSol = {'indFinal': indFinal, 'xArrayFields': xArrayFields, 'yArrayFields': yArrayFields, 'fields': fields}
            comSol['outDirTest'] = outDirTest
            comSol['showPlot'] = cfgMain['FLAGS'].getboolean('showPlot')
            comSol['Tsave'] = timeList[ind_t]

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


def plotProfilesSimiSol(ind_time, outputName, comSol, simiDict, solSimi, axis, particles=False):
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
    if particles:
        x = comSol['x']
        y = comSol['y']
        h = comSol['h']
        v = comSol['v']
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
    Time = solSimi['Time']

    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, pU.figH))
    ax2 = ax1.twinx()

    if axis == 'xaxis':
        ax1.axvline(x=xCenter, linestyle=':')
        FD = fields['FD'][indFinal, :]
        hSimi = hSimi[indFinal, :]
        # DFA simulation
        ax1.plot(xArrayFields, FD, 'k', label='Field flow depth')
        ax2.plot(xArrayFields, FD*fields['FV'][indFinal, :], 'g', label=getLabel('Field', '', dir='', vert=True))
        ax2.plot(xArrayFields, FD*fields['Vx'][indFinal, :], 'm', label=getLabel('Field', '', dir='x'))
        ax2.plot(xArrayFields, FD*fields['Vy'][indFinal, :], 'b', label=getLabel('Field', '', dir='y'))
        ax2.plot(xArrayFields, FD*fields['Vz'][indFinal, :], 'c', label=getLabel('Field', '', dir='z'))
        if particles:
            ax1.plot(x, h, '.k', linestyle='None', label='Part flow depth')
            ax2.plot(x, h*v, '.g', linestyle='None', label=getLabel('Part', '', dir='', vert=True))
        # similarity solution
        ax1.plot(xArrayFields, hSimi, '--k', label='SimiSol flow depth')
        ax2.plot(xArrayFields, hSimi*vSimi[indFinal, :], '--g', label=getLabel('SimiSol', '', dir='', vert=True))
        ax2.plot(xArrayFields, hSimi*vxSimi[indFinal, :], '--m', label=getLabel('SimiSol', '', dir='x'))
        ax2.plot(xArrayFields, hSimi*vySimi[indFinal, :], '--b', label=getLabel('SimiSol', '', dir='y'))
        ax2.plot(xArrayFields, hSimi*vzSimi[indFinal, :], '--c', label=getLabel('SimiSol', '', dir='z'))
        ax1.set_title('Profile along flow at t=%.2f (com1DFA), %.2f s (simiSol) (simulation %s)'
                      % (Tsave, Time[ind_time], outputName))
        ax1.set_xlabel('x in [m]')
        indStart = min(first_nonzero(hSimi, 0), first_nonzero(FD, 0)) - 2
        indEnd = max(last_nonzero(hSimi, 0), last_nonzero(FD, 0)) + 2
        ax1.set_xlim([xArrayFields[indStart], xArrayFields[indEnd]])

    elif axis == 'yaxis':
        FD = fields['FD'][:, indFinal]
        hSimi = hSimi[:, indFinal]
        # DFA simulation
        ax1.plot(yArrayFields, FD, 'k', label='Field flow depth')
        ax2.plot(yArrayFields, FD*fields['FV'][:, indFinal], 'g', label=getLabel('Field', '', dir='', vert=True))
        ax2.plot(yArrayFields, FD*fields['Vx'][:, indFinal], 'm', label=getLabel('Field', '', dir='x'))
        ax2.plot(yArrayFields, FD*fields['Vy'][:, indFinal], 'b', label=getLabel('Field', '', dir='y'))
        ax2.plot(yArrayFields, FD*fields['Vz'][:, indFinal], 'c', label=getLabel('Field', '', dir='z'))
        if particles:
            ax1.plot(y, h, '.k', linestyle='None', label='Part flow depth')
            ax2.plot(y, h*v, '.g', linestyle='None', label=getLabel('Part', '', dir='', vert=True))
        # similarity solution
        ax1.plot(yArrayFields, hSimi, '--k', label='SimiSol flow depth')
        ax2.plot(yArrayFields, hSimi*vSimi[:, indFinal], '--g', label=getLabel('SimiSol', '', dir='', vert=True))
        ax2.plot(yArrayFields, hSimi*vxSimi[:, indFinal], '--m', label=getLabel('SimiSol', '', dir='x'))
        ax2.plot(yArrayFields, hSimi*vySimi[:, indFinal], '--b', label=getLabel('SimiSol', '', dir='y'))
        ax2.plot(yArrayFields, hSimi*vzSimi[:, indFinal], '--c', label=getLabel('SimiSol', '', dir='z'))
        ax1.set_title('Profile across flow at t=%.2f (com1DFA), %.2f s (simiSol) (simulation %s)'
                      % (Tsave, Time[ind_time], outputName))
        ax1.set_xlabel('y in [m]')
        indStart = min(first_nonzero(hSimi, 0), first_nonzero(FD, 0)) - 2
        indEnd = max(last_nonzero(hSimi, 0), last_nonzero(FD, 0)) + 2
        ax1.set_xlim([yArrayFields[indStart], yArrayFields[indEnd]])

    ax1.set_ylabel('flow depth [m]')
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    color = 'tab:green'
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.grid(color='tab:green', linestyle='-', linewidth=0.25, alpha=0.5)
    ax2.set_ylabel(getLabel('', '', dir='') + r'$[ms-1]$', color=color)
    ax2.legend(loc='upper right')
    ax1.legend(loc='upper left')

    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'profile_' + str(outputName) + '_%sCutSol_T%.2f.' % (axis, Tsave) + pU.outputFormat, fig1)


def plotErrorTime(time, hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray, outDirTest, outputName,
                  simDFrow, relativ):
    """plot error between a given com1DFA sol and the analytic sol
    function of time
    Parameters
    -----------
    time: 1D numpy array
        time array
    hErrorL2Array: 1D numpy array
        flow thickness L2 error array
    hErrorLMaxArray: 1D numpy array
        flow thickness LMax error array
    vhErrorL2Array: 1D numpy array
        flow momentum L2 error array
    vhErrorLMaxArray: 1D numpy array
        flow momentum LMax error array
    outDirTest: pathlib path
        path to outpute folder
    outputName: str
        oupute file name
    simDFrow: dataFrame
        simDF row corresponding to the curent simulation analyzed
    relativ: str
    """
    title = (' between similarity solution and com1DFA \n(simulation %s)'
                     % (outputName))
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
    ax2.plot(time, vhErrorL2Array, 'g-', label=getLabel('L2 error', '', dir=''))
    ax2.plot(time, vhErrorLMaxArray, 'g--', label=getLabel('LMax error', '', dir=''))
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel(getTitleError(relativ, getLabel(' on', '', dir='')), color=color)
    ax2.legend(loc='lower right')
    ax2.grid(color='tab:green', linestyle='-', linewidth=0.25, alpha=0.5)

    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'Error_Time_' + str(outputName), fig1)


def plotErrorConvergence(simDF, outDirTest, cfgSimi, xField, yField, coloredBy, sizedBy, logScale=False, fit=False):
    """plot error between all com1DFA sol and analytic sol
    function of whatever you want

    Parameters
    -----------
    simDF: dataFrame
        the simulation data with the postprocessing results
    outDirTest: str or pathlib
        output directory
    cfgSimi: configparser
        the cfg
    xField: str
        column of the simDF to use for the x axis
    yField: str
        column of the simDF to use for the y axis
    coloredBy: str
        column of the simDF to use for the colors
    sizedBy: str
        column of the simDF to use for the marker size
    logScale: boolean
        If you want a loglog scale
    fit: boolean
        if True add power law regression
    """
    tSave = cfgSimi.getfloat('tSave')
    relativ = cfgSimi.getboolean('relativError')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, min(simDF[coloredBy]), max(simDF[coloredBy]), continuous=pU.contCmap)
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(3*pU.figW, 2*pU.figH))
    # get the sizing function
    sizeList = simDF[sizedBy].unique()
    lenSize = len(sizeList)
    minSize = np.nanmin(sizeList)
    maxSize = np.nanmax(sizeList)
    if lenSize > 1:
        sizeList = (simDF[sizedBy].to_numpy() - minSize) / (maxSize - minSize) * 70 + 10
    else:
        sizeList = np.array([100])
    # make the scatter plot
    scatter1 = ax1.scatter(simDF[xField], simDF[yField], c=simDF[coloredBy], s=sizeList, cmap=cmap, norm=norm,
                          marker=pU.markers[0], alpha=1)#, edgecolors='k')

    # #########################################
    # If you want to add some regression lines
    slopeU = pd.DataFrame()
    slopeH = pd.DataFrame()
    colorValueList = simDF[coloredBy].unique()
    lenColor = len(colorValueList)
    if fit:
        for colorValue in colorValueList:
            simDFNew = simDF[simDF[coloredBy] == colorValue]
            for sizeValue in simDFNew[sizedBy].unique():
                simDFNewNew = simDFNew[simDFNew[sizedBy] == sizeValue]
                xArray = simDFNewNew[xField]
                hErrorL2 = simDFNewNew["hErrorL2"]
                vErrorL2 = simDFNewNew["vhErrorL2"]
                p, rSquaredH, _, _, _ = np.polyfit(np.log(xArray), np.log(hErrorL2), deg=1, full=True)
                p1H = p[0]
                p0H = np.exp(p[1])
                slopeH.loc[sizeValue, colorValue] = p1H
                ax1.plot(xArray, p0H*xArray**p1H, 'r')
                if np.size(rSquaredH) == 0:
                    rSquaredH = np.nan
                log.info('power law fit sphKernelRadius = %.2f m: hErrorL2 = %.1f * Npart^{%.2f}, r=%.2f' % (colorValue, p0H, p1H, rSquaredH))

    if logScale:
        ax1.set_yscale('log')
        ax1.set_xscale('log')

    fig1.suptitle('Convergence of DFA simulation for the similarity solution test at t = %.2fs' % tSave)
    ax1.set_title(getTitleError(relativ, r' L2 on flow thickness'))
    ax1.set_xlabel(xField)
    ax1.set_ylabel(getTitleError(relativ, r' L2 on flow thickness'))
    if lenColor<=10:
        lenColor = None
    legend1 = ax1.legend(*scatter1.legend_elements(num=lenColor), loc="upper center", title=coloredBy)
    ax1.add_artist(legend1)

    # produce a legend with a cross section of sizes from the scatter
    if lenSize<=10:
        lenSize = None
    kw = dict(prop="sizes", color=scatter1.cmap(0.7),
          func=lambda s: (s-10)*(maxSize - minSize)/70 + minSize)
    legend3 = ax1.legend(*scatter1.legend_elements(num=lenSize, **kw), loc="upper right", title=sizedBy)
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    ax1.grid(color='grey', which='minor', linestyle='--', linewidth=0.25, alpha=0.5)
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'ErrorLog%ds' % int(tSave), fig1)

    return fig1, ax1, slopeH


def plotErrorRef(simDF, outDirTest, cfgSimi, xField, yFieldArray, coloredBy, sizedBy, logScale=False):
    """plot error between all com1DFA sol and analytic sol
    function of whatever (ini parameter given in the simDF) you want

    Parameters
    -----------
    simDF: dataFrame
        the simulation data with the postprocessing results
    outDirTest: str or pathlib
        output directory
    cfgSimi: configparser
        the cfg
    xField: str
        column of the simDF to use for the x axis
    yFieldArray: list
        list of max 2 column of the simDF to use for the y axis
    coloredBy: str
        column of the simDF to use for the colors
    sizedBy: str
        column of the simDF to use for the marker size
    logScale: boolean
        If you want a loglog scale
    """
    tSave = cfgSimi.getfloat('tSave')
    relativ = cfgSimi.getboolean('relativError')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, min(simDF[coloredBy]), max(simDF[coloredBy]), continuous=pU.contCmap)
    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, 2*pU.figH))
    ax2 = ax1.twinx()
    # for the scatter plot (new filtering)
    simDFScatter = simDF[simDF['viscOption'] == 1]
    # get the sizing function
    sizeList = simDFScatter[sizedBy].unique()
    minSize = np.nanmin(sizeList)
    maxSize = np.nanmax(sizeList)
    if len(sizeList) > 1:
        sizeList = (simDFScatter[sizedBy].to_numpy() - minSize) / (maxSize - minSize) * 70 + 10
    else:
        sizeList = np.array([100])
    # make the scatter plot
    scatter = ax1.scatter(simDFScatter[xField], simDFScatter[yFieldArray[0]], c=simDFScatter[coloredBy], s=sizeList, cmap=cmap,
                          marker=pU.markers[0], alpha=1)#, edgecolors='k')
    scatter2 = ax2.scatter(simDFScatter[xField], simDFScatter[yFieldArray[1]], c=simDFScatter[coloredBy], s=sizeList, cmap=cmap,
                           marker=pU.markers[1], alpha=1)#, edgecolors='k')

    # ######################################
    # plot the reference as lines here
    # Add the horizontal lines for the atta viscosity
    simDFhline = simDF[simDF['viscOption'] == 2]
    simDFhline = simDFhline.sort_values(by=coloredBy, ascending=True)
    handles1 = []
    # handles2 = []
    for simHash, simDFrow in simDFhline.iterrows():
        hl = ax1.axhline(simDFrow[yFieldArray[0]], color=cmap(norm(simDFrow[coloredBy])), label='Ata, nPart = %d, csz = %.2f' % (simDFrow['nPart'], simDFrow['sphKernelRadius']))
        handles1.append(hl)
        hl = ax2.axhline(simDFrow[yFieldArray[1]], color=cmap(norm(simDFrow[coloredBy])), linestyle='--', label='Ata, nPart = %d, csz = %.2f' % (simDFrow['nPart'], simDFrow['sphKernelRadius']))
        # handles1.append(hl)

    # #########################################
    # Adding legend and titles
    legend = ax1.legend(handles=handles1, loc="upper left")
    ax1.add_artist(legend)
    # legend = ax2.legend(handles=handles2, loc="upper left")
    # ax2.add_artist(legend)
    if logScale:
        ax1.set_yscale('log')
        ax2.set_yscale('log')
        ax1.set_xscale('log')
    ax1.set_title('Convergence of DFA simulation for the similarity solution test at t = %.2fs' % tSave)
    ax1.set_xlabel(xField)
    ax1.set_ylabel(getTitleError(relativ, r' L2 on flow depth ($\bullet$)'))
    ax2.set_ylabel(getTitleError(relativ, getLabel(' L2 on', r'$(\blacksquare)$', dir='')))
    legend1 = ax1.legend(*scatter.legend_elements(), loc="upper center", title=coloredBy)
    ax1.add_artist(legend1)

    # produce a legend with a cross section of sizes from the scatter
    kw = dict(prop="sizes", color=scatter.cmap(0.7),
          func=lambda s: (s-10)*(maxSize - minSize)/70 + minSize)
    legend2 = ax1.legend(*scatter.legend_elements(**kw), loc="upper right", title=sizedBy)
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    ax1.grid(color='grey', which='minor', linestyle='--', linewidth=0.25, alpha=0.5)
    b1, t1 = ax1.get_ylim()
    b2, t2 = ax2.get_ylim()
    ax1.set_ylim([min(b1, b2), max(t1, t2)])
    ax2.set_ylim([min(b1, b2), max(t1, t2)])
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'ErrorRef%ds' % int(tSave), fig1)
    return fig1, ax1, ax2


def plotPresentation(simDF, outDirTest, cfgSimi, xField, yField, coloredBy, sizedBy, logScale=False, fit=False):
    """plot error between all com1DFA sol and analytic sol
    function of whatever you want

    Parameters
    -----------
    simDF: dataFrame
        the simulation data with the postprocessing results
    outDirTest: str or pathlib
        output directory
    cfgSimi: configparser
        the cfg
    xField: str
        column of the simDF to use for the x axis
    yField: str
        column of the simDF to use for the y axis
    coloredBy: str
        column of the simDF to use for the colors
    sizedBy: str
        column of the simDF to use for the marker size
    logScale: boolean
        If you want a loglog scale
    fit: boolean
        if True add power law regression
    """
    tSave = cfgSimi.getfloat('tSave')
    relativ = cfgSimi.getboolean('relativError')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, min(simDF[coloredBy]), max(simDF[coloredBy]), continuous=pU.contCmap)
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(4*pU.figW, 2*pU.figH))
    if logScale:
        ax1.set_yscale('log')
        ax1.set_xscale('log')
    # get the sizing function
    sizeList = simDF[sizedBy].unique()
    lenSize = len(sizeList)
    minSize = np.nanmin(sizeList)
    maxSize = np.nanmax(sizeList)
    if lenSize > 1:
        sizeList = (simDF[sizedBy].to_numpy() - minSize) / (maxSize - minSize) * 70 + 10
    else:
        sizeList = np.array([100])

    # #########################################
    # If you want to add some regression lines
    slopeH = pd.DataFrame()
    colorValueList = -np.sort(-simDF[coloredBy].unique())
    lenColor = len(colorValueList)
    if lenColor<=10:
        lenColor = None
    if lenSize<=10:
        lenSize = None
    count = 0
    simDFOld = ''
    for colorValue in colorValueList:
        simDFNew = simDF[simDF[coloredBy] == colorValue]
        for sizeValue in simDFNew[sizedBy].unique():
            simDFNewNew = simDFNew[simDFNew[sizedBy] == sizeValue]
            xArray = simDFNewNew[xField]
            hErrorL2 = simDFNewNew[yField]
            if count >= 1:
                simDFNewNew = simDFNewNew.append(simDFOld)
            colorList = simDFNewNew[coloredBy].unique()
            lenColor = len(colorList)
            if lenColor<=10:
                lenColor = None
            # make the scatter plot
            scatter1 = ax1.scatter(simDFNewNew[xField], simDFNewNew[yField], c=simDFNewNew[coloredBy], s=sizeList, cmap=cmap, norm=norm,
                                   marker=pU.markers[0], alpha=1)
            fig1.suptitle('Convergence of DFA simulation for the similarity solution test at t = %.2fs' % tSave)
            ax1.set_title(getTitleError(relativ, r' L2 on flow thickness'))
            ax1.set_xlabel(xField)
            ax1.set_ylabel(getTitleError(relativ, r' L2 on flow thickness'))
            legend1 = ax1.legend(*scatter1.legend_elements(num=lenColor), loc="upper center", title=coloredBy)
            ax1.add_artist(legend1)

            # produce a legend with a cross section of sizes from the scatter
            kw = dict(prop="sizes", color=scatter1.cmap(0.7),
                  func=lambda s: (s-10)*(maxSize - minSize)/70 + minSize)
            legend3 = ax1.legend(*scatter1.legend_elements(num=lenSize, **kw), loc="upper right", title=sizedBy)
            ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
            ax1.grid(color='grey', which='minor', linestyle='--', linewidth=0.25, alpha=0.5)
            pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'ErrorPresentation%d' % count, fig1)
            if fit:
                p, rSquaredH, _, _, _ = np.polyfit(np.log(xArray), np.log(hErrorL2), deg=1, full=True)
                p1H = p[0]
                p0H = np.exp(p[1])
                slopeH.loc[sizeValue, colorValue] = p1H
                color = cmap(norm(simDFNewNew[coloredBy][0]))
                ax1.plot(xArray, p0H*xArray**p1H, color=color)
                infoText = '%s = %.2f' % (coloredBy, simDFNewNew[coloredBy][0])
                ax1.text(max(1.05*xArray), p0H*max(xArray)**p1H, infoText, color=color,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.5))

                if np.size(rSquaredH) == 0:
                    rSquaredH = np.nan
                log.info('power law fit sphKernelRadius = %.2f m: hErrorL2 = %.1f * Npart^{%.2f}, r=%.2f' % (colorValue, p0H, p1H, rSquaredH))

                pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'ErrorPresentation%dFit' % count, fig1)

            count = count + 1
            simDFOld = copy.deepcopy(simDFNewNew)
            legend1.remove()
    return fig1, ax1


def plotTimeCPULog(simDF, outDirTest, cfgSimi, xField, coloredBy, sizedBy, logScale=False):
    """plot computation time function of nParts
    function of whatever (ini parameter given in the simDF) you want
    Parameters
    -----------
    simDF: dataFrame
        the simulation data with the postprocessing results
    outDirTest: str or pathlib
        output directory
    cfgSimi: configparser
        the cfg
    xField: str
        column of the simDF to use for the x axis
    coloredBy: str
        column of the simDF to use for the colors
    sizedBy: str
        column of the simDF to use for the marker size
    logScale: boolean
        If you want a loglog scale
    """
    colorList = simDF[coloredBy].unique()
    tSave = cfgSimi.getfloat('tSave')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, min(simDF[coloredBy]), max(simDF[coloredBy]), continuous=pU.contCmap)
    # get the sizing function
    sizeList = simDF[sizedBy].unique()
    minSize = np.nanmin(sizeList)
    maxSize = np.nanmax(sizeList)
    if len(sizeList) > 1:
        sizeList = (simDF[sizedBy].to_numpy() - minSize) / (maxSize - minSize) * 70 + 10
    else:
        sizeList = np.array([100])
    fig1, ax1 = plt.subplots(figsize=(2*pU.figW, 2*pU.figH))
    nameList = ['timeLoop', 'timeForce', 'timeForceSPH', 'timePos', 'timeNeigh', 'timeField']
    for count, name in enumerate(nameList):
        scatter = ax1.scatter(simDF[xField], simDF[name], c=simDF[coloredBy], s=sizeList, cmap=cmap,
                              marker=pU.markers[count], alpha=1, edgecolors='k')

    slopeTsph = pd.DataFrame()
    for colorValue in colorList:
        simDFNew = simDF[simDF[coloredBy] == colorValue]
        sizeValue = max(simDFNew[sizedBy].unique())
        for sizeValue in simDFNew[sizedBy].unique():
            simDFNewNew = simDFNew[simDFNew[sizedBy] == sizeValue]
            nPart = simDFNewNew[xField]
            timeLoop = simDFNewNew["timeLoop"]
            timeForce = simDFNewNew["timeForce"]
            timeForceSPH = simDFNewNew["timeForceSPH"]
            timePos = simDFNewNew["timePos"]
            timeNeigh = simDFNewNew["timeNeigh"]
            timeField = simDFNewNew["timeField"]
            p = np.polyfit(np.log(simDFNewNew[xField]), np.log(timeLoop), deg=1)
            p11 = p[0]
            p01 = np.exp(p[1])
            p = np.polyfit(np.log(simDFNewNew[xField]), np.log(timeForce), deg=1)
            p12 = p[0]
            p02 = np.exp(p[1])
            p = np.polyfit(np.log(simDFNewNew[xField]), np.log(timeForceSPH), deg=1)
            p13 = p[0]
            p03 = np.exp(p[1])
            p = np.polyfit(np.log(simDFNewNew[xField]), np.log(timePos), deg=1)
            p14 = p[0]
            p04 = np.exp(p[1])
            p = np.polyfit(np.log(simDFNewNew[xField]), np.log(timeNeigh), deg=1)
            p15 = p[0]
            p05 = np.exp(p[1])
            p = np.polyfit(np.log(simDFNewNew[xField]), np.log(timeField), deg=1)
            p16 = p[0]
            p06 = np.exp(p[1])
            slopeTsph.loc[sizeValue, colorValue] = p13
            handles1 = []
            hl = ax1.plot(nPart, p01*nPart**p11, 'k', label='timeLoop')
            handles1.append(hl[0])
            hl = ax1.plot(nPart, p02*nPart**p12, 'g', label='timeForce')
            handles1.append(hl[0])
            hl = ax1.plot(nPart, p03*nPart**p13, 'r', label='timeForceSPH')
            handles1.append(hl[0])
            hl = ax1.plot(nPart, p04*nPart**p14, 'b', label='timePos')
            handles1.append(hl[0])
            hl = ax1.plot(nPart, p05*nPart**p15, 'm', label='timeNeigh')
            handles1.append(hl[0])
            hl = ax1.plot(nPart, p06*nPart**p16, 'c', label='timeField')
            handles1.append(hl[0])
            log.info('power law fit sphKernelRadius = %.2f m: timeLoop = %.1f * nPart^{%.2f}' % (colorValue, p01, p11))
            log.info('power law fit sphKernelRadius = %.2f m: timeForce = %.1f * nPart^{%.2f}' % (colorValue, p02, p12))
            log.info('power law fit sphKernelRadius = %.2f m: timeForceSPH = %.1f * nPart^{%.2f}' % (colorValue, p03, p13))
            log.info('power law fit sphKernelRadius = %.2f m: timePos = %.1f * nPart^{%.2f}' % (colorValue, p04, p14))
            log.info('power law fit sphKernelRadius = %.2f m: timeNeigh = %.1f * nPart^{%.2f}' % (colorValue, p05, p15))
            log.info('power law fit sphKernelRadius = %.2f m: timeField = %.1f * nPart^{%.2f}' % (colorValue, p06, p16))

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_title('CPU time')
    ax1.set_xlabel('number of particles')
    ax1.set_ylabel('Time [s]')
    # Adding legend and titles
    legend = ax1.legend(handles=handles1, loc="upper left")
    ax1.add_artist(legend)
    legend1 = ax1.legend(*scatter.legend_elements(), loc="lower left", title=coloredBy)
    ax1.add_artist(legend1)

    # produce a legend with a cross section of sizes from the scatter
    kw = dict(prop="sizes", color=scatter.cmap(0.7),
          func=lambda s: (s-10)*(maxSize - minSize)/70 + minSize)
    legend2 = ax1.legend(*scatter.legend_elements(**kw), loc="upper right", title=sizedBy)
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    ax1.grid(color='grey', which='minor', linestyle='--', linewidth=0.25, alpha=0.5)
    plt.show()
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'timeCPU%ds' % int(tSave), fig1)
    fig2, ax = plt.subplots(figsize=(2*pU.figW, 2*pU.figH))
    for sizeValue, simDFrow in slopeTsph.iterrows():
        ax.scatter(slopeTsph.columns, simDFrow, label=sizeValue)
    plt.legend()
    plt.show()


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
        return 'Error' + ending


def getLabel(start, end, dir='', vert=True):
    """Get error plot title (relativ error or not?)"""
    if dir:
        return start + r' $h \bar{u}_' + dir + r'$ ' + end
    else:
        if vert:
            return start + r' $\vert h \mathbf{\bar{u}} \vert $' + end
        else:
            return start + r' $h \mathbf{\bar{u}}$ ' + end


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
