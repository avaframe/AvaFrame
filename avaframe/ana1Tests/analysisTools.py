# imports
import numpy as np

# local imports
import avaframe.com1DFA.DFAtools as DFAtls


def normL2Vect(analyticalSol, numericalSol, cellSize, cosAngle):
    """ Compute relativ L2 and Lmax norm of the error between the analytic and numerical solution
    """
    nonZeroIndex = np.where((np.abs(analyticalSol['fx']) > 0) & (np.abs(analyticalSol['fy']) > 0)
                            & (np.abs(analyticalSol['fz']) > 0) & (np.abs(numericalSol['fz']) > 0)
                            & (np.abs(numericalSol['fz']) > 0) & (np.abs(numericalSol['fz']) > 0))
    dvx = analyticalSol['fx'] - numericalSol['fx']
    dvy = analyticalSol['fy'] - numericalSol['fy']
    dvz = analyticalSol['fz'] - numericalSol['fz']
    dv = DFAtls.norm2(dvx, dvy, dvz)
    vAnalyticalSol = DFAtls.norm2(analyticalSol['fx'], analyticalSol['fy'], analyticalSol['fz'])
    localError = dv[nonZeroIndex]
    vAnalyticalSolL2 = cellSize * cellSize / cosAngle * np.nansum(vAnalyticalSol[nonZeroIndex])
    vAnalyticalSolL2 = np.sqrt(vAnalyticalSolL2)
    error2 = cellSize * cellSize / cosAngle * np.nansum(localError)
    error = np.sqrt(error2) / vAnalyticalSolL2
    errorMax = np.sqrt(np.nanmax(np.append(localError, 0))) / np.sqrt(np.nanmax(np.append(vAnalyticalSol[nonZeroIndex], 0)))

    return error, errorMax


def normL2Scal(analyticalSol, numericalSol, cellSize, cosAngle):
    """ Compute relativ L2 and Lmax norm of the error between the analytic and numerical solution
    """
    nonZeroIndex = np.where((analyticalSol > 0) & (numericalSol > 0))
    localError = (analyticalSol[nonZeroIndex] - numericalSol[nonZeroIndex])
    localError = localError * localError
    hAnalyticalSolL2 = cellSize * cellSize / cosAngle * np.nansum(analyticalSol[nonZeroIndex])
    hAnalyticalSolL2 = np.sqrt(hAnalyticalSolL2)
    error2 = cellSize * cellSize / cosAngle * np.nansum(localError)
    error = np.sqrt(error2) / hAnalyticalSolL2
    errorMax = np.sqrt(np.nanmax(localError)) / np.sqrt(np.nanmax(analyticalSol[nonZeroIndex]))

    return error, errorMax
