# imports
import numpy as np

# local imports
import avaframe.com1DFA.DFAtools as DFAtls


def normL2Vect(analyticalSol, numericalSol, cellSize, cosAngle):
    """ Compute L2 and Lmax norm of the error between the analytic and numerical solution
        Parameters
        -----------
        analyticalSol: dictionary
            analytic solution dictionary:
                fx: x component of the vector
                fy: y component of the vector
                fz: z component of the vector
        numericalSol: dictionary
            numericalSol solution dictionary:
                fx: x component of the vector
                fy: y component of the vector
                fz: z component of the vector
        cellSize: float
            grid cell size
        cosAngle: float
            cosine of the slope angle
        Returns
        ---------
        errorL2: float
            L2 error
        errorL2Rel: float
            Relativ L2 error
        errorMax: float
            LMax error
        errorMaxRel: float
            Relativ LMax error
    """
    # only compare location where both solutions are non zero
    nonZeroIndex = np.where((np.abs(analyticalSol['fx']) > 0) & (np.abs(analyticalSol['fy']) > 0)
                            & (np.abs(analyticalSol['fz']) > 0) & (np.abs(numericalSol['fz']) > 0)
                            & (np.abs(numericalSol['fz']) > 0) & (np.abs(numericalSol['fz']) > 0))
    # compute difference in all 3 directions
    dvx = analyticalSol['fx'] - numericalSol['fx']
    dvy = analyticalSol['fy'] - numericalSol['fy']
    dvz = analyticalSol['fz'] - numericalSol['fz']
    # compute the norm2 of the difference
    dv = DFAtls.norm2(dvx, dvy, dvz)
    localError = dv[nonZeroIndex]
    # compute the norm2 of the reference
    analyticalSol2 = DFAtls.norm2(analyticalSol['fx'], analyticalSol['fy'], analyticalSol['fz'])
    analyticalSol2 = analyticalSol2[nonZeroIndex]

    # compute error L2, LMax and relativ errors
    errorL2, errorL2Rel, errorMax, errorMaxRel = computeErrorAndNorm(localError, analyticalSol2, cellSize, cosAngle)

    return errorL2, errorL2Rel, errorMax, errorMaxRel


def normL2Scal(analyticalSol, numericalSol, cellSize, cosAngle):
    """ Compute L2 and Lmax norm of the error between the analytic and numerical solution
        Parameters
        -----------
        analyticalSol: numpy array
            analytic solution array
        numericalSol: numpy array
            numericalSol solution array
        cellSize: float
            grid cell size
        cosAngle: float
            cosine of the slope angle
        Returns
        ---------
        errorL2: float
            L2 error
        errorL2Rel: float
            Relativ L2 error
        errorMax: float
            LMax error
        errorMaxRel: float
            Relativ LMax error
    """
    # only compare location where both solutions are non zero
    nonZeroIndex = np.where((analyticalSol > 0) & (numericalSol > 0))
    # compute the norm2 of the difference
    localError = (analyticalSol[nonZeroIndex] - numericalSol[nonZeroIndex])
    localError = localError * localError
    # compute the norm2 of the reference
    analyticalSol2 = analyticalSol*analyticalSol
    analyticalSol2 = analyticalSol2[nonZeroIndex]

    # compute error L2, LMax and relativ errors
    errorL2, errorL2Rel, errorMax, errorMaxRel = computeErrorAndNorm(localError, analyticalSol2, cellSize, cosAngle)

    return errorL2, errorL2Rel, errorMax, errorMaxRel


def computeErrorAndNorm(localError, analyticalSol2, cellSize, cosAngle):
    """ Compute error between two functions given their norme 2
        Parameters
        -----------
        localError: numpy array
            norm2 of the error function
        analyticalSol2: numpy array
            norm2 of the reference function
        cellSize: float
            grid cell size
        cosAngle: float
            cosine of the slope angle
        Returns
        ---------
        errorL2: float
            L2 error
        errorL2Rel: float
            Relativ L2 error
        errorMax: float
            LMax error
        errorMaxRel: float
            Relativ LMax error
    """
    # compute the L2 norm of the difference
    errorL2 = L2Norm(localError, cellSize, cosAngle)
    # compute the LMax norm of the difference
    errorMax = np.sqrt(np.nanmax(np.append(localError, 0)))
    # compute the L2 norm of the reference
    analyticalSolL2 = L2Norm(analyticalSol2, cellSize, cosAngle)
    # compute the LMax norm of the reference
    vAnalyticalSolLMax = np.sqrt(np.nanmax(np.append(analyticalSol2, 0)))
    # compute relativ error
    if analyticalSolL2 > 0:
        errorL2Rel = errorL2 / analyticalSolL2
    else:
        errorL2Rel = errorL2

    if vAnalyticalSolLMax > 0:
        errorMaxRel = errorMax / vAnalyticalSolLMax
    else:
        errorMaxRel = errorMax

    return errorL2, errorL2Rel, errorMax, errorMaxRel


def L2Norm(norm2Array, cellSize, cosAngle):
    """ Compute L2 norm of an array
        Parameters
        -----------
        norm2Array: numpy array
            norm2 of the function
        cellSize: float
            grid cell size
        cosAngle: float
            cosine of the slope angle
        Returns
        ---------
        normL2: float
            normL2 of the function
    """
    normL2 = cellSize * cellSize / cosAngle * np.nansum(norm2Array)
    normL2 = np.sqrt(normL2)
    return normL2
