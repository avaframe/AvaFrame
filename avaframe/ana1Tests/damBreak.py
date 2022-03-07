"""

Simple python script to reproduce analytic solution for a Riemann problem,
following the derivations in Faccanoni and Mangeney (2012), Test 2, Case 1.2.
but scaled up in size.

Here the instantanous release of fluid from rest is described using incompressible,
depth-avaeraged mass and momentum conservation equations and a Coulomb-tpye friction law.


"""

import numpy as np
import matplotlib.pyplot as plt
import logging

# local imports
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.out3Plot.plotUtils as pU
from avaframe.com1DFA import DFAtools
import avaframe.ana1Tests.analysisTools as anaTools
import avaframe.out3Plot.outAna1Plots as outAna1Plots

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
        cfgC: configParser object
            configuration setting for avalanche simulation including DAMBREAK section
        solDam: dict
            analytic solution dictionary:
                tAna: 1D numpy array
                    time array
                h0: float
                    release thickness
                hAna: 2D numpy array
                    flow depth (rows for x and columns for time)
                uAna: 2D numpy array
                    flow velocity (rows for x and columns for time)
                xAna: 2D numpy array
                    extent of domain in the horizontal plane coordinate system (rows for x and columns for time)
                xMidAna: 1D numpy array
                    middle of the material in x dir in the horizontal plane coordinate system
                    (used to compute the error)
        fields0: dict
            initial time step field dictionary
        fieldsT: dict
            tSave field dictionary
        tSave: float
            time step of analaysis
        header: dict
            fields header dictionary
        outDirTest: pathli path
            path to output directory

    """
    phi = cfg.getfloat('phi')
    phiRad = np.radians(phi)
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
    dataIniV = DFAtools.scalProd(dataIniVx, dataIniVy, dataIniVz, np.cos(phiRad), 0, -np.sin(phiRad))
    dataAnaV = DFAtools.scalProd(dataAnaVx, dataAnaVy, dataAnaVz, np.cos(phiRad), 0, -np.sin(phiRad))

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

    fig, ax1 = _plotMultVariables(x, y, nx_loc, tSave, dataIniFD, dataAnaFD, solDam['xAna'], solDam['hAna'],
                             indTime, 'Flow depth', 'm')
    outputName = 'CompareDamBreakH%s_%.02f' % (simName, tSave)
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outputName, fig)

    y = np.zeros(len(x))
    fig, ax2 = _plotMultVariables(x, y, nx_loc, tSave, dataIniFD*dataIniV, dataAnaFD*dataAnaV, solDam['xAna'],
                             solDam['hAna']*solDam['uAna'], indTime, 'Flow momentum in slope directrion', 'm²/s')
    outputName = 'CompareDamBreakHVel%s_%.02f' % (simName, tSave)
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outputName, fig)

    y = np.zeros(len(x))
    fig, ax2 = _plotMultVariables(x, y, nx_loc, tSave, dataIniV, dataAnaV, solDam['xAna'],
                             solDam['uAna'], indTime, 'Flow velocity in slope directrion', 'm/s')
    outputName = 'CompareDamBreakVel%s_%.02f' % (simName, tSave)
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, outputName, fig)

    return ax1, ax2

def damBreakSol(avaDir, cfg, cfgC, outDirTest):
    """ Compute analytical flow depth and velocity for dam break test case

        for a granular flow over a dry rough sloping bed with the Savage Hutter model

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        cfg: configParser object
            main avaframe configuration - here showPlot flag used
        cfgC: configParser object
            configuration setting for avalanche simulation including DAMBREAK section
        outDirTest: pathlib path
            path to output directory

        Returns
        -------
        solDam: dict
            analytic solution dictionary:
                tAna: 1D numpy array
                    time array
                h0: float
                    release thickness
                hAna: 2D numpy array
                    flow depth (rows for x and columns for time)
                uAna: 2D numpy array
                    flow velocity (rows for x and columns for time)
                xAna: 2D numpy array
                    extent of domain in the horizontal plane coordinate system (rows for x and columns for time)
                xMidAna: 1D numpy array
                    middle of the material in x dir in the horizontal plane coordinate system
                    (used to compute the error)
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
        cond1R = x0R + ((m0*t[m]) / 2.0 - 2*cR) * t[m]
        cond2R = x0R + ((m0*t[m]) / 2.0 + cR) * t[m]
        cond1 = ((m0*t[m]) / 2.0 - cL) * t[m]
        cond2 = (2.0 *cL + ((m0*t[m]) / 2.0)) * t[m]
        # elif x[k] > cond2:
        u[:, m] = 0
        h[:, m] = 0
        # elif cond1 < x[k] <= cond2:
        if t[m] > 0:
            u[:, m] = np.where(cond2 >= x, (2./3.) * (cL + (x / t[m]) + m0 * t[m]), u[:, m])
            h[:, m] = np.where(cond2 >= x, ((2.* cL - (x / t[m]) + ((m0 * t[m]) / 2.))**2) / (9. * gz), h[:, m])
        # if x[k] <= cond1:
        u[:, m] = np.where(cond1 >= x, m0 * t[m], u[:, m])
        h[:, m] = np.where(cond1 >= x, hL, h[:, m])
        # to get the analytical solution on the back part (not working very well, theory is probably wrong)
        # # elif cond1R < x[k] <= cond2R:
        # if t[m] > 0:
        #     u[:, m] = np.where(cond2R >= x, (2./3.) * (-cR + ((x-x0R) / t[m]) + m0 * t[m]), u[:, m])
        #     h[:, m] = np.where(cond2R >= x, ((2.* cR + ((x-x0R) / t[m]) - ((m0 * t[m]) / 2.))**2) / (9. * gz), h[:, m])
        # # if x[k] <= cond1R:
        # u[:, m] = np.where(cond1R >= x, 0, u[:, m])
        # h[:, m] = np.where(cond1R >= x, 0, h[:, m])


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
    cfgDam: configParser object
        configuration setting for avalanche simulation including DAMBREAK section
    simDF: pandas dataFrame
        configuration DF
    solDam: dict
        analytic solution dictionary
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
            if cfgDam['DAMBREAK'].getboolean('onlyLast'):
                fieldsList = [fieldsList[0], fieldsList[indTime]]
                timeList = [0, timeList[indTime]]
                indTime = 1

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
    """Compare analytical and com1DFA results, compute error

        Parameters
        -----------
        avalancheDir: pathlib path
            path to avalanche directory
        fieldsList: list
            list of fields dictionaries
        timeList: list
            list of time steps
        solDam: dict
            analytic solution dictionary
        fieldHeader: dict
            header dictionary with info about the extend and cell size
        cfgDam: configParser object
            configuration setting for avalanche simulation including DAMBREAK section
        outDirTest: pathlib path
            path output directory (where to save the figures)
        simHash: str
            com1DFA simulation id
        simDFrow: pandas object
            com1DFA simulation rox coresponding to simHash

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
    phiRad = np.radians(phi)
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
        dataNumV = DFAtools.scalProd(field['Vx'], field['Vy'], field['Vz'], np.cos(phiRad), 0, -np.sin(phiRad))
        uNumPlus = dataNumV[nRowMin:nRowMax, nColMid:nColMax]

        # make analytical results 2D
        nrows, ncols = np.shape(hNumPlus)
        O = np.ones((nrows, ncols))
        hDamPlus = O*hDamPlus
        uDamPlus = O*uDamPlus

        hEL2Plus, hL2RelPlus, hLmaxPlus, hLmaxRelPlus = anaTools.normL2Scal(hDamPlus, hNumPlus, cellSize,
                                                                            np.cos(phiRad))
        vhL2Plus, vhL2RelPlus, vhLmaxPlus, vhLmaxRelPlus = anaTools.normL2Scal(hDamPlus*uDamPlus, hNumPlus*uNumPlus,
                                                                               cellSize, np.cos(phiRad))
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
        # Make all individual time step comparison plot
        if cfgDam['DAMBREAK'].getboolean('plotSequence'):
            plotComparison(cfgDam['DAMBREAK'], simHash, fieldsList[0], field, fieldHeader, solDam, t, outDirTest)
        count = count + 1

    if cfgDam['DAMBREAK'].getboolean('plotErrorTime') and len(timeList)>1:
        outAna1Plots.plotErrorTime(timeList, hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray,
                                   outDirTest, simHash, simDFrow, cfgDam['DAMBREAK'].getboolean('relativError'))

    outAna1Plots.plotDamBreakSummary(avalancheDir, timeList, fieldsList, fieldHeader, solDam, hErrorL2Array, hErrorLMaxArray,
                        vhErrorL2Array, vhErrorLMaxArray, outDirTest, simHash, simDFrow, cfgDam)

    return hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray


def getDamExtend(fieldHeader, xM, cfgDam):
    """ Get the extend where the analytic and simulation should be compared

    Parameters
    -----------
    fieldHeader: dict
        header dictionary with info about the extend and cell size
    xM: float
        x coordinate ot the start of the comparizon domain
    cfgDam: configParser object
        configuration setting for the DAMBREAK section

    Returns
    --------
    xDamPlus: 1D numpy array
        x array corresponding to the comparizon domain
    nColMid: int
        index of the column corresponding to xM
    nColMax: int
        index of the column corresponding to xEnd
    nRowMin: int
        index of the row corresponding to yStart
    nRowMax: int
        index of the row corresponding to yEnd
    """
    cellSize = fieldHeader['cellsize']
    xllc = fieldHeader['xllcenter']
    yllc = fieldHeader['yllcenter']
    nCols = fieldHeader['ncols']
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
