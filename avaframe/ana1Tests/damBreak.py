"""

Simple python script to reproduce analytic solution for a Riemann problem,
following the derivations in Faccanoni and Mangeney (2012), Test 2, Case 1.2.
but scaled up in size.

Here the instantanous release of fluid from rest is described using incompressible,
depth-avaeraged mass and momentum conservation equations and a Coulomb-tpye friction law.


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import logging

# local imports
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.out3Plot.plotUtils as pU
import avaframe.in2Trans.ascUtils as IOf
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.com1DFA.DFAfunctionsCython as DFAfunC
import avaframe.com1DFA.particleTools as particleTools
from avaframe.com1DFA import DFAtools
import avaframe.ana1Tests.analysisTools as anaTools
import avaframe.out3Plot.outAna1Plots as outAna1Plots
import avaframe.out3Plot.outCom1DFA as outCom1DFA
from avaframe.in1Data import getInput as gI

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def _plotVariable(var, x, xMiddle, dtInd, dtStep, label):

    fig = plt.figure(figsize=(pU.figW, pU.figH))
    plt.title('Dry-Bed')
    plt.plot(x, var[:,0], 'k--', label='t_init')
    plt.plot(x, var[:,dtInd], label='t = %.1fs' % dtStep)
    plt.vlines(xMiddle[dtInd], ymin=0, ymax=max(var[:,0]), color='r', linestyles='--')
    plt.xlabel('x-coordinate [m]')
    plt.ylabel(label)
    plt.legend()

    return fig


def plotResults(t, x, xMiddle, h, u, tSave, cfg, outDirTest):
    """ Create plots of the analytical solution for the given settings,
        including an animation

    """

    # index of time steps
    dtInd = np.searchsorted(t, tSave)

    fig = _plotVariable(h, x, xMiddle, dtInd, tSave, 'Flow depth [m]')
    outputName = 'damBreakFlowDepth'
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outputName, fig)

    fig = _plotVariable(u, x, xMiddle, dtInd, tSave, 'Flow velocity [ms-1]')
    outputName = 'damBreakFlowVelocity'
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outputName, fig)


def _plotMultVariables(x, y, nx_loc, dtAnalysis, data1, data2, xR, dataR, tR, label, unit):
    """ generate plots """

    fig, ax = plt.subplots(nrows=1, sharex=True, figsize=(pU.figW*3, pU.figH*2))
    ax.plot(x, y, 'grey', linestyle='--')
    ax.plot(x, data1[nx_loc, :], 'k--', label='init')
    ax.plot(x, data2[nx_loc, :], 'b', label='com1DFAPy')
    ax.plot(xR, dataR[:,tR], 'r-', label='analyt')
    ax.set_xlabel('Along track [ncols]')
    ax.set_ylabel('%s [%s]' % (label, unit))
    ax.set_xlim([-200, 200])
    plt.legend()
    ax.set_title('%s at time step %.02f s' % (label, dtAnalysis))

    return fig, ax


def plotComparison(cfg, simName, fields0, fieldsT, header, solDam, tSave, outDirTest):
    """ Generate plots that compare the simulation results to the analytical solution

        Parameters
        -----------

        dataComSol: dataFrame
            dataframe of simulation results (including name, file path, result type, etc.)
        hL: float
            initial release thickness
        xR: numpy array
            x extent of domain
        hR: numpy array
            flow depth of analytical solution
        uR: numpy array
            flow velocity of analytical solution
        dtAnalysis: float
            time step of analaysis
        outDirTest: pathli path
            path to output directory

    """
    phi = cfg.getfloat('phi')
    # Load data
    dataIniFD = fields0['FD']
    dataAnaFD = fieldsT['FD']
    dataIniVx = fields0['Vx']
    dataIniVy = fields0['Vy']
    dataIniVz = fields0['Vz']
    dataAnaVx = fieldsT['Vx']
    dataAnaVy = fieldsT['Vy']
    dataAnaVz = fieldsT['Vz']
    # project velocity on inclined plane
    dataIniV = DFAtools.scalProd(dataIniVx, dataIniVy, dataIniVz, np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))
    dataAnaV = DFAtools.scalProd(dataAnaVx, dataAnaVy, dataAnaVz, np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))

    # Location of Profiles
    cellSize = header['cellsize']
    ny = dataAnaFD.shape[0]
    nx = dataAnaFD.shape[1]
    xllc = header['xllcenter']
    nx_loc = int(ny *0.5)

    # set x Vector
    x = np.arange(xllc, xllc + nx*cellSize, cellSize)
    y = np.zeros(len(x))
    y[x<0] = solDam['h0']
    y[x>=0] = 0.0
    y[x<-120] = 0.0

    # setup index for time of analyitcal solution
    indTime = np.searchsorted(solDam['tAna'], tSave)
    t = solDam['tAna'][indTime]

    fig, ax1 = _plotMultVariables(x, y, nx_loc, tSave, dataIniFD, dataAnaFD, solDam['xAna'], solDam['hAna'],
                             indTime, 'Flow depth', 'm')
    outputName = 'CompareDamBreakH%s_%.02f' % (simName, tSave)
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outputName, fig)

    y = np.zeros(len(x))
    fig, ax2 = _plotMultVariables(x, y, nx_loc, tSave, dataIniFD*dataIniV, dataAnaFD*dataAnaV, solDam['xAna'],
                             solDam['hAna']*solDam['uAna'], indTime, 'Flow momentum in slope directrion', 'm/s')
    outputName = 'CompareDamBreakHVel%s_%.02f' % (simName, tSave)
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outputName, fig)

    return ax1, ax2

def damBreakSol(avaDir, cfg, cfgC, outDirTest):
    """ Compute flow depth for dam break for granular flow over a dry rough sloping bed with the Savage Hutter model

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        cfg: configParser object
            main avaframe configuration - here showPlot flag used
        cfgC: configParser object
            configuration setting for avalanche simulation including DAMBREAK section

        Returns
        -------
        hL: float
            release thickness
        h: numpy array
            flow depth
        u: numpy array
            flow velocity
        x: numpy array
            extent of domain in the horizontal plane coordinate system
        xMiddle: numpy array
            middle of the material in x dir in the horizontal plane coordinate system
    """

    # Set Parameters
    # Coordinate system chosen in the direction of the inclined plane
    g = cfgC['GENERAL'].getfloat('gravAcc')       # acceleration due to gravity [ms-2]
    phiDeg = cfgC['DAMBREAK'].getfloat('phi')       # acceleration due to gravity [ms-2]
    deltaDeg = cfgC['DAMBREAK'].getfloat('delta')       # acceleration due to gravity [ms-2]
    phi = np.radians(phiDeg)                          # slope angle [°]
    delta = np.radians(deltaDeg)                        # bed friction angle [°]
    gz = g * np.cos(phi)                          # projection of g perpendicular to the inclined plane
    m0 = gz * (np.tan(phi) - np.tan(delta))       # constant x-acceleration resulting from gravity and friction force
    hL = cfgC['GENERAL'].getfloat('relTh')        # initial height [m] in Riemann problem in state 1 (x<0), hR (x>0)=0
    cL = np.sqrt(gz * hL)
    hR = cfgC['GENERAL'].getfloat('relTh')        # initial height [m] in Riemann problem in state 1 (x<0), hR (x>0)=0
    cR = np.sqrt(gz * hR)
    x0R = -120 / np.cos(phi)                           # wave celeritiy
    tSave = cfgC['DAMBREAK'].getfloat('tSave')
    # Define time [0-1] seconds and space [-2,2] meters domains multiplied times 100
    t = np.linspace(0, 30, 3000)
    x = np.linspace(-200 / np.cos(phi) , 200 / np.cos(phi) , 800)
    nt = len(t)
    nx = len(x)
    # Initialise flow depth solution and velocity
    h = np.zeros((nx, nt))
    u = np.zeros((nx, nt))
    xMiddle = np.zeros(nt)

    # Compute exact solution for case: 'dry bed' - including three different states
    for m in range(nt):
        for k in range(nx):
            cond1R = x0R + ((m0*t[m]) / 2.0 - 2*cR) * t[m]
            cond2R = x0R + ((m0*t[m]) / 2.0 + cR) * t[m]
            cond1 = ((m0*t[m]) / 2.0 - cL) * t[m]
            cond2 = (2.0 *cL + ((m0*t[m]) / 2.0)) * t[m]
                # if x[k] <= cond1R:
                #     u[k,m] = 0
                #     h[k,m] = 0
                # elif cond1R < x[k] <= cond2R:
                #     u[k,m] = (2./3.) * (-cR + ((x[k]-x0R) / t[m]) + m0 * t[m])
                #     h[k,m] = ((2.* cR + ((x[k]-x0R) / t[m]) - ((m0 * t[m]) / 2.))**2) / (9. * gz)
            if x[k] <= cond1:
                u[k,m] = m0 * t[m]
                h[k,m] = hL
            elif cond1 < x[k] <= cond2:
                u[k,m] = (2./3.) * (cL + (x[k] / t[m]) + m0 * t[m])
                h[k,m] = ((2.* cL - (x[k] / t[m]) + ((m0 * t[m]) / 2.))**2) / (9. * gz)
            elif x[k] > cond2:
                u[k,m] = 0
                h[k,m] = 0.0
        xMiddle[m] = (2*cond2R + cond1)/3
    #-----------------------------Plot results --------------
    # Reproduce figure 6, case 1.2.1 - Test 2
    plotResults(t, x, xMiddle, h, u, tSave, cfg, outDirTest)

    x = x * np.cos(phi)  # projected on the horizontal plane
    xMiddle = xMiddle * np.cos(phi)  # projected on the horizontal plane
    solDam = {'tAna': t, 'h0': hL, 'xAna': x, 'hAna': h, 'uAna': u, 'xMidAna': xMiddle}
    return solDam


def postProcessDamBreak(avalancheDir, cfgMain, cfgDam, simDF, solDam, outDirTest):
    """ loop on all DFA simulations and compare then to the anlytic solution

    Parameters
    -----------
    avalancheDir: str or pathlib path
        avalanche directory
    cfgMain: confiparser
        avaframeCfg configuration
    cfgDam: dict
        configuration for dambreak sol
    simDF: pandas dataFrame
        configuration DF
    solDam: dict
        dambreak analytical solution
    outDirTest: pathlib path
        path to output directory

    Returns
    --------
    simDF: pandas dataFrame
        configuration DF appended with the analysis results
    """
    # loop on all the simulations and make the comparison to reference
    for simHash, simDFrow in simDF.iterrows():
        simName = simDFrow['simName']
        # fetch the simulation results
        # particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName=simHash, flagAvaDir=True, comModule='com1DFA')
        fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['FD', 'FV', 'Vx', 'Vy', 'Vz'],
                                                            simName=simName, flagAvaDir=True, comModule='com1DFA')
        # analyze and compare results
        if cfgDam['DAMBREAK']['tSave'] == '':
            indTime = -1
            tSave = timeList[-1]
        else:
            tSave = cfgDam['DAMBREAK'].getfloat('tSave')
            indTime = min(np.searchsorted(timeList, tSave), min(len(timeList)-1, len(fieldsList)-1))
            # fieldsList = [fieldsList[indTime]]
            # timeList = [timeList[indTime]]
            # indTime = 0

        # computeAndPlotGradient(avalancheDir, particlesList, timeList, solDam, cfgDam, outDirTest, simHash, simDFrow)
        hEL2Array, hELMaxArray, vhEL2Array, vhELMaxArray = analyzeResults(avalancheDir, fieldsList, timeList, solDam, fieldHeader,
                                                                          cfgDam, outDirTest, simHash, simDFrow)
        # add result of error analysis
        # save results in the simDF
        simDF.loc[simHash, 'hErrorL2'] = hEL2Array[indTime]
        simDF.loc[simHash, 'vhErrorL2'] = vhEL2Array[indTime]
        simDF.loc[simHash, 'hErrorLMax'] = hELMaxArray[indTime]
        simDF.loc[simHash, 'vhErrorLMax'] = vhELMaxArray[indTime]
        # +++++++++POSTPROCESS++++++++++++++++++++++++
        # -------------------------------
    name = 'results' + str(round(tSave)) + '.p'
    simDF.to_pickle(outDirTest / name)

    return simDF


def analyzeResults(avalancheDir, fieldsList, timeList, solDam, fieldHeader, cfgDam, outDirTest, simHash, simDFrow):
    """Compare analytical and com1DFA results
        Parameters
        -----------
        fieldsList: list
            list of fields dictionaries
        solSimi: dictionary
            similarity solution:
                time: time array (without dimention)
                Time: time array (with dimention)
                g_sol: g array
                g_p_sol: first derivativ of g array
                f_sol: f array
                f_p_sol: first derivativ of f array
        fieldHeader: dict
            header dictionary with info about the extend and cell size
        cfgSimi: dict
            simisol confguration section
        outDirTest: pathlib path
            path output directory (where to save the figures)
        simHash: str
            com1DFA simulation id

        Returns
        --------
        hErrorL2Array: numpy array
            L2 error on flow depth for saved time steps
        hErrorLMaxArray: numpy array
            LMax error on flow depth for saved time steps
        vErrorL2Array: numpy array
            L2 error on flow velocity for saved time steps
        vErrorLMaxArray: numpy array
            LMax error on flow velocity for saved time steps

    """
    phi = cfgDam['DAMBREAK'].getfloat('phi')
    cellSize = fieldHeader['cellsize']
    xMiddle = solDam['xMidAna']
    tAna = solDam['tAna']
    xAna = solDam['xAna']
    hAna = solDam['hAna']
    uAna = solDam['uAna']
    hErrorL2Array = np.zeros((len(fieldsList)))
    vhErrorL2Array = np.zeros((len(fieldsList)))
    hErrorLMaxArray = np.zeros((len(fieldsList)))
    vhErrorLMaxArray = np.zeros((len(fieldsList)))
    count = 0
    # run the comparison routine for each saved time step
    for t, field in zip(timeList, fieldsList):
        # get similartiy solution h, u at required time step
        indTime = np.searchsorted(tAna, t)
        xM = xMiddle[indTime]
        # get extend where the sol should be compared
        xDamPlus, nColMid, nColMax, nRowMin, nRowMax = getDamExtend(fieldHeader, xM, cfgDam['DAMBREAK'])
        # get analytical solution (in same format as the simulation results)
        hDamPlus = np.interp(xDamPlus, xAna, hAna[:, indTime])
        uDamPlus = np.interp(xDamPlus, xAna, uAna[:, indTime])
        # get numirical sol on good interval
        hNumPlus = field['FD'][nRowMin:nRowMax, nColMid:nColMax]
        dataNumV = DFAtools.scalProd(field['Vx'], field['Vy'], field['Vz'], np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))
        uNumPlus = dataNumV[nRowMin:nRowMax, nColMid:nColMax]

        # make analytical results 2D
        nrows, ncols = np.shape(hNumPlus)
        O = np.ones((nrows, ncols))
        hDamPlus = O*hDamPlus
        uDamPlus = O*uDamPlus

        hEL2Plus, hL2RelPlus, hLmaxPlus, hLmaxRelPlus = anaTools.normL2Scal(hDamPlus, hNumPlus, cellSize, np.cos(np.radians(phi)))
        vhL2Plus, vhL2RelPlus, vhLmaxPlus, vhLmaxRelPlus = anaTools.normL2Scal(hDamPlus*uDamPlus, hNumPlus*uNumPlus, cellSize, np.cos(np.radians(phi)))
        if cfgDam['DAMBREAK'].getboolean('relativError'):
            hErrorL2Array[count] = hL2RelPlus
            hErrorLMaxArray[count] = hLmaxRelPlus
            vhErrorL2Array[count] = vhL2RelPlus
            vhErrorLMaxArray[count] = vhLmaxRelPlus
        else:
            hErrorL2Array[count] = hEL2Plus
            hErrorLMaxArray[count] = hLmaxPlus
            vhErrorL2Array[count] = vhL2Plus
            vhErrorLMaxArray[count] = vhLmaxPlus
        title = outAna1Plots.getTitleError(cfgDam['DAMBREAK'].getboolean('relativError'))
        log.info("L2 %s error on the Flow Depth at t=%.2f s is : %.4f" % (title, t, hEL2Plus))
        log.info("L2 %s error on the momentum at t=%.2f s is : %.4f" % (title, t, vhL2Plus))
        # if cfgDam['DAMBREAK'].getboolean('plotIntermediate'):
        #     plotComparison(cfgDam['DAMBREAK'], simHash, fieldsList[0], field, fieldHeader, solDam, t, outDirTest)

        count = count + 1
    if cfgDam['DAMBREAK'].getboolean('plotIntermediate') and len(timeList)>1:
        outAna1Plots.plotErrorTime(timeList, hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray, outDirTest,
                                   simHash, simDFrow, cfgDam['DAMBREAK'].getboolean('relativError'))
        plotDamBreakSummary(avalancheDir, timeList, fieldsList, fieldHeader, solDam, hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray, outDirTest,
                            simHash, simDFrow, cfgDam)

    return hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray


def getDamExtend(fieldHeader, xM, cfgDam):
    cellSize = fieldHeader['cellsize']
    xllc = fieldHeader['xllcenter']
    yllc = fieldHeader['yllcenter']
    nCols = fieldHeader['ncols']
    xStart = cfgDam.getfloat('xStart')
    xEnd = cfgDam.getfloat('xEnd')
    yStart = cfgDam.getfloat('yStart')
    yEnd = cfgDam.getfloat('yEnd')
    xArray = np.arange(nCols)*cellSize + xllc
    nColMid = round((xM-xllc)/cellSize)
    nColMax = round((xEnd-xllc)/cellSize)
    nRowMin = round((yStart-yllc)/cellSize)
    nRowMax = round((yEnd-yllc)/cellSize)

    xDamPlus = xArray[nColMid:nColMax]

    return xDamPlus, nColMid, nColMax, nRowMin, nRowMax


def plotDamBreakSummary(avalancheDir, timeList, fieldsList, fieldHeader, solDam, hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray, outDirTest,
                        simHash, simDFrow, cfgDam):

    # Initialise DEM
    demFile = gI.getDEMPath(avalancheDir)
    demOri = IOf.readRaster(demFile, noDataToNan=True)
    demOri, dem = com1DFA.initializeMesh(cfgDam['GENERAL'], demOri, cfgDam['GENERAL'].getint('methodMeshNormal'))
    dem['header']['xllcenter'] = demOri['header']['xllcenter']
    dem['header']['yllcenter'] = demOri['header']['yllcenter']

    phi = cfgDam['DAMBREAK'].getfloat('phi')
    tSave = cfgDam['DAMBREAK'].getfloat('tSave')
    relativ = cfgDam['DAMBREAK'].getboolean('relativError')
    indT = min(np.searchsorted(timeList, tSave), min(len(timeList)-1, len(fieldsList)-1))
    tSave = timeList[indT]
    # Load data
    fields0 = fieldsList[0]
    fieldsT = fieldsList[indT]
    dataIniFD = fields0['FD']
    dataFD = fieldsT['FD']
    dataIniVx = fields0['Vx']
    dataIniVy = fields0['Vy']
    dataIniVz = fields0['Vz']
    dataVx = fieldsT['Vx']
    dataVy = fieldsT['Vy']
    dataVz = fieldsT['Vz']
    # project velocity on inclined plane
    dataIniV = DFAtools.scalProd(dataIniVx, dataIniVy, dataIniVz, np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))
    dataV = DFAtools.scalProd(dataVx, dataVy, dataVz, np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))

    # Location of Profiles
    cellSize = fieldHeader['cellsize']
    ny = dataFD.shape[0]
    nx = dataFD.shape[1]
    xllc = fieldHeader['xllcenter']
    yllc = fieldHeader['yllcenter']
    nx_loc = int(ny *0.5)

    # set x Vector
    x = np.arange(xllc, xllc + nx*cellSize, cellSize)
    y = np.zeros(len(x))
    y[x<0] = solDam['h0']
    y[x>=0] = 0.0
    y[x<-120] = 0.0

    # setup index for time of analyitcal solution
    indTime = np.searchsorted(solDam['tAna'], tSave)

    # create figures and plots
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))
    fig.suptitle('Simulation %s, t = %.2f s' % (simHash, tSave), fontsize=30)
    # make flow thickness comparison plot
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax1.plot(x, y, 'grey', linestyle='--')
    ax1.plot(x, dataIniFD[nx_loc, :], 'k--', label='init')
    ax1.plot(x, dataFD[nx_loc, :], 'b', label='simulation')
    ax1.plot(solDam['xAna'], solDam['hAna'][:,indTime], 'r-', label='analytic')
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('Flow thickness [m]')
    ax1.set_xlim([-200, 200])
    plt.legend()
    ax1.set_title('Flow thickness')

    # make flow momentum comparison
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax2.plot(x, dataIniV[nx_loc, :], 'k--', label='init')
    ax2.plot(x, dataFD[nx_loc, :]*dataV[nx_loc, :], 'b', label='simulation')
    ax2.plot(solDam['xAna'], solDam['uAna'][:,indTime]*solDam['hAna'][:,indTime], 'r-', label='analytic')
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel(outAna1Plots.getLabel('Flow momentum', '[m²/s]', dir='', vert=True))
    ax2.set_xlim([-200, 200])
    plt.legend()
    ax2.set_title(outAna1Plots.getLabel('Flow momentum', '', dir='', vert=True))

    # make bird view plot
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax3, extent, cbar0, cs1 = outCom1DFA.addResult2Plot(ax3, dem, fieldsList[-1]['FD'], 'FD')
    cbar0.ax.set_ylabel('flow thickness')
    ax3 = outCom1DFA.addDem2Plot(ax3, dem, what='slope', extent=extent)
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(fieldsList[-1]['FD'], dem['header']['cellsize'],
                                                                 extentOption=True, constrainedData=False, buffer='')
    ax3.set_ylim([rowsMin+yllc, rowsMax+yllc])
    ax3.set_xlim([colsMin+xllc, colsMax+xllc])
    ax3.set_xlabel('x [m]')
    ax3.set_ylabel('y [m]')
    # ax3.set_title('Flow thickness')
    pU.putAvaNameOnPlot(ax3, avalancheDir)

    # make error plot
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    title = ' between analytical solution and com1DFA'
    title = outAna1Plots.getTitleError(relativ, title)
    ax5 = ax4.twinx()
    ax4.plot(timeList, hErrorL2Array, 'k-', label='Flow thickness L2 error')
    ax4.plot(timeList, hErrorLMaxArray, 'k--', label='Flow thickness LMax error')
    ax4.set_title(title)
    ax4.set_xlabel('time in [s]')
    ax4.set_ylabel(outAna1Plots.getTitleError(relativ, ' on flow thickness'))
    ax4.legend(loc='upper left')
    ax4.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

    color = 'tab:green'
    ax5.plot(timeList, vhErrorL2Array, 'g-', label=outAna1Plots.getLabel('L2 error', '', dir=''))
    ax5.plot(timeList, vhErrorLMaxArray, 'g--', label=outAna1Plots.getLabel('LMax error', '', dir=''))
    ax5.tick_params(axis='y', labelcolor=color)
    ax5.set_ylabel(outAna1Plots.getTitleError(relativ, outAna1Plots.getLabel(' on', '', dir='')), color=color)
    ax5.legend(loc='lower right')
    ax5.grid(color='tab:green', linestyle='-', linewidth=0.25, alpha=0.5)
    ax4.set_yscale('log')
    ax5.set_yscale('log')

    outFileName = '_'.join([simHash, 'DamBreakTest'])
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outFileName, fig)



def computeAndPlotGradient(avalancheDir, particlesList, timeList, solDam, cfgDam, outDirTest, simHash, simDFrow):

    phi = cfgDam['DAMBREAK'].getfloat('phi')
    # Initialise DEM
    demFile = gI.getDEMPath(avalancheDir)
    demOri = IOf.readRaster(demFile, noDataToNan=True)
    xllc = demOri['header']['xllcenter']
    yllc = demOri['header']['yllcenter']
    demOri, dem = com1DFA.initializeMesh(cfgDam['GENERAL'], demOri, cfgDam['GENERAL'].getint('methodMeshNormal'))
    for t, particles in zip(timeList, particlesList):
        x = particles['x'] + xllc
        y = particles['y'] + yllc
        z = particles['z']
        ux = particles['ux']
        uy = particles['uy']
        uz = particles['uz']
        force2 = {}
        particles, force2 = DFAfunC.computeForceSPHC(cfgDam['GENERAL'], particles, force2, dem,
            cfgDam.getint('GENERAL', 'sphOption'), gradient=1)
        gradNorm = DFAtls.norm(force2['forceSPHX'], force2['forceSPHY'], force2['forceSPHZ'])
        uMag = DFAtls.norm(ux, uy, uz)
        v = DFAtls.scalProd(ux, uy, uz, np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))
        grad = DFAtls.scalProd(force2['forceSPHX'], force2['forceSPHY'], force2['forceSPHZ'], np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))

        fig = plt.figure(figsize=(pU.figW, pU.figH))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        indPartOfInterest = np.where(np.abs(y)<20)
        xPart = x[indPartOfInterest]
        vPart = v[indPartOfInterest]
        gradPart = grad[indPartOfInterest]
        indPartRest = np.where(gradPart<-0.06136465)
        indPartBack = np.where(gradPart<-0.82909272)
        ax1.plot(xPart, vPart, color='b', marker='.', linestyle='None', label='particle flow depth')
        ax2.plot(xPart, gradPart, color='b', marker='o', linestyle='None', label='SPH gradient used')
        ax1.plot(xPart[indPartRest], vPart[indPartRest], color='k', marker='.', linestyle='None', label='particle flow depth')
        ax2.plot(xPart[indPartRest], gradPart[indPartRest], color='k', marker='o', linestyle='None', label='SPH gradient used')
        ax1.plot(xPart[indPartBack], vPart[indPartBack], color='r', marker='.', linestyle='None', label='particle flow depth')
        ax2.plot(xPart[indPartBack], gradPart[indPartBack], color='r', marker='o', linestyle='None', label='SPH gradient used')
        # ax2.plot(r, v, color='b', marker='.', linestyle='None')

        ax1.set_xlabel('r in [m]')
        ax2.set_xlabel('r in [m]')
        ax2.set_title('Gradient of the flow depth')
        ax1.legend()
        ax2.legend()
        plt.show()
