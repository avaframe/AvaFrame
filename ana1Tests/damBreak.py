"""
Simple python script to reproduce analytic solution for a Riemann problem,
following the derivations in Faccanoni and Mangeney (2012), Test 2, Case 1.2.
but scaled up in size.

Here the instantanous release of fluid from rest is described using incompressible,
thickness-avaeraged mass and momentum conservation equations and a Coulomb-tpye friction law.
"""

import numpy as np
import logging

# local imports
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.com1DFA import DFAtools
import avaframe.ana1Tests.analysisTools as anaTools
import avaframe.out3Plot.outAna1Plots as outAna1Plots

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def damBreakSol(avaDir, cfg, cfgC, outDirTest):
    """ Compute analytical Flow thickness and velocity for dam break test case

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
                    Flow thickness (rows for x and columns for time)
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
    phiDeg = cfgC['DAMBREAK'].getfloat('phi')       # slope angle [°]
    deltaDeg = cfgC['DAMBREAK'].getfloat('delta')       # friction angle [°]
    phi = np.radians(phiDeg)                          # slope angle [rad]
    delta = np.radians(deltaDeg)                        # bed friction angle [rad]
    gz = g * np.cos(phi)                          # projection of g perpendicular to the inclined plane
    m0 = gz * (np.tan(phi) - np.tan(delta))       # constant x-acceleration resulting from gravity and friction force
    hL = cfgC['GENERAL'].getfloat('relTh')        # initial height [m] in Riemann problem in state 1 (x<0), hR (x>0)=0
    cL = np.sqrt(gz * hL)
    x0R = cfgC['DAMBREAK'].getfloat('xBack') / np.cos(phi)
    tSave = cfgC['DAMBREAK'].getfloat('tSave')
    # Define time [0-1] seconds and space [-2,2] meters domains multiplied times 100
    tAna = np.arange(0, cfgC['DAMBREAK'].getfloat('tEnd'), cfgC['DAMBREAK'].getfloat('dt'))
    x = np.arange(cfgC['DAMBREAK'].getfloat('xStart') / np.cos(phi) , cfgC['DAMBREAK'].getfloat('xEnd') / np.cos(phi),
                  cfgC['DAMBREAK'].getfloat('dx'))
    nt = len(tAna)
    nx = len(x)
    # Initialise Flow thickness solution and velocity
    h = np.zeros((nx, nt))
    u = np.zeros((nx, nt))
    xMiddle = np.zeros(nt)

    # Compute exact solution for case: 'dry bed' - including three different states
    for m in range(nt):
        cond1R = x0R + ((m0*tAna[m]) / 2.0 - 2*cL) * tAna[m]
        cond2R = x0R + ((m0*tAna[m]) / 2.0 + cL) * tAna[m]
        cond1 = ((m0*tAna[m]) / 2.0 - cL) * tAna[m]
        cond2 = (2.0 *cL + ((m0*tAna[m]) / 2.0)) * tAna[m]
        # elif x[k] > cond2:
        u[:, m] = 0
        h[:, m] = 0
        # elif cond1 < x[k] <= cond2:
        if tAna[m] > 0:
            u[:, m] = np.where(cond2 >= x, (2./3.) * (cL + (x / tAna[m]) + m0 * tAna[m]), u[:, m])
            h[:, m] = np.where(cond2 >= x, ((2.* cL - (x / tAna[m]) + ((m0 * tAna[m]) / 2.))**2) / (9. * gz), h[:, m])
        # if x[k] <= cond1:
        u[:, m] = np.where(cond1 >= x, m0 * tAna[m], u[:, m])
        h[:, m] = np.where(cond1 >= x, hL, h[:, m])
        # uncomment this if you also want the rear part of the dam break
        # to get the analytical solution on the back part (not working very well, theory is probably wrong)
        # # elif cond1R < x[k] <= cond2R:
        # if tAna[m] > 0:
        #     u[:, m] = np.where(cond2R >= x, (2./3.) * (-cL + ((x-x0R) / tAna[m]) + m0 * tAna[m]), u[:, m])
        #     h[:, m] = np.where(cond2R >= x, ((2.* cL + ((x-x0R) / tAna[m]) -
        #                        ((m0 * tAna[m]) / 2.))**2) / (9. * gz), h[:, m])
        # # if x[k] <= cond1R:
        # u[:, m] = np.where(cond1R >= x, 0, u[:, m])
        # h[:, m] = np.where(cond1R >= x, 0, h[:, m])


        xMiddle[m] = (2*cond2R + cond1)/3

    #-----------------------------Plot results --------------
    # Reproduce figure 6, case 1.2.1 - Test 2
    outAna1Plots.plotDamAnaResults(tAna, x, xMiddle, h, u, tSave, cfg, outDirTest)

    x = x * np.cos(phi)  # projected on the horizontal plane
    xMiddle = xMiddle * np.cos(phi)  # projected on the horizontal plane
    solDam = {'tAna': tAna, 'h0': hL, 'xAna': x, 'hAna': h, 'uAna': u, 'xMidAna': xMiddle}
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
        fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['FT', 'FV', 'Vx', 'Vy', 'Vz'],
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
        hEL2Array, hELMaxArray, vhEL2Array, vhELMaxArray, t = analyzeResults(avalancheDir, fieldsList, timeList, solDam,
                                                                             fieldHeader, cfgDam, outDirTest, simHash,
                                                                             simDFrow)
        # add result of error analysis
        # save results in the simDF
        simDF.loc[simHash, 'timeError'] = t
        simDF.loc[simHash, 'hErrorL2'] = hEL2Array[indTime]
        simDF.loc[simHash, 'vhErrorL2'] = vhEL2Array[indTime]
        simDF.loc[simHash, 'hErrorLMax'] = hELMaxArray[indTime]
        simDF.loc[simHash, 'vhErrorLMax'] = vhELMaxArray[indTime]
        # +++++++++POSTPROCESS++++++++++++++++++++++++
        # -------------------------------
    name = 'results' + str(round(t)) + '.p'
    simDF.to_pickle(outDirTest / name)

    return simDF


def analyzeResults(avalancheDir, fieldsList, timeList, solDam, fieldHeader, cfg, outDirTest, simHash, simDFrow):
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
        cfg: configParser object
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
            L2 error on Flow thickness for saved time steps
        hErrorLMaxArray: numpy array
            LMax error on Flow thickness for saved time steps
        vErrorL2Array: numpy array
            L2 error on flow velocity for saved time steps
        vErrorLMaxArray: numpy array
            LMax error on flow velocity for saved time steps
        tSave: float
            time corresponding the errors

    """
    cfgDam = cfg['DAMBREAK']
    phi = cfgDam.getfloat('phi')
    phiRad = np.radians(phi)
    cellSize = fieldHeader['cellsize']
    xMiddle = solDam['xMidAna']
    tAna = solDam['tAna']
    xAna = solDam['xAna']
    hAna = solDam['hAna']
    uAna = solDam['uAna']
    # get plots limits
    limits = outAna1Plots.getPlotLimits(cfgDam, fieldsList, fieldHeader)
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
        # get extent where the sol should be compared
        xDamPlus, nColMid, nColMax, nRowMin, nRowMax = getDamExtend(fieldHeader, xM, cfgDam)
        # get analytical solution (in same format as the simulation results)
        hDamPlus = np.interp(xDamPlus, xAna, hAna[:, indTime])
        uDamPlus = np.interp(xDamPlus, xAna, uAna[:, indTime])
        # get numirical sol on good interval
        hNumPlus = field['FT'][nRowMin:nRowMax, nColMid:nColMax]
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
        if cfgDam.getboolean('relativError'):
            hErrorL2Array[count] = hL2RelPlus
            hErrorLMaxArray[count] = hLmaxRelPlus
            vhErrorL2Array[count] = vhL2RelPlus
            vhErrorLMaxArray[count] = vhLmaxRelPlus
        else:
            hErrorL2Array[count] = hEL2Plus
            hErrorLMaxArray[count] = hLmaxPlus
            vhErrorL2Array[count] = vhL2Plus
            vhErrorLMaxArray[count] = vhLmaxPlus
        title = outAna1Plots.getTitleError(cfgDam.getboolean('relativError'))
        log.debug("L2 %s error on the Flow thickness at t=%.2f s is : %.4f" % (title, t, hEL2Plus))
        log.debug("L2 %s error on the momentum at t=%.2f s is : %.4f" % (title, t, vhL2Plus))
        # Make all individual time step comparison plot
        if cfgDam.getboolean('plotSequence'):
            outAna1Plots.plotComparisonDam(cfgDam, simHash, fieldsList[0], field, fieldHeader, solDam, t, limits,
                                           outDirTest)
        count = count + 1

    tSave = cfgDam.getfloat('tSave')
    indT = min(np.searchsorted(timeList, tSave), min(len(timeList)-1, len(fieldsList)-1))
    tSave = timeList[indT]
    if cfgDam.getboolean('plotErrorTime') and len(timeList)>1:
        # Create result plots

        outAna1Plots.plotErrorTime(timeList, hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray,
                                   cfgDam.getboolean('relativError'), tSave, simHash, outDirTest)

    outAna1Plots.plotDamBreakSummary(avalancheDir, timeList, fieldsList, fieldHeader, solDam, hErrorL2Array,
                                     hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray, outDirTest, simDFrow, simHash,
                                     cfg)

    return hErrorL2Array, hErrorLMaxArray, vhErrorL2Array, vhErrorLMaxArray, tSave


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
