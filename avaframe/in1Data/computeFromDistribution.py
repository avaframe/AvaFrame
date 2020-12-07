"""
    Extract sample of values to fit a beta-pert distribution

    This file is part of Avaframe.

"""

# load python modules
import os
import numpy as np
import scipy.special as sc
from scipy import stats
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# load avaframe modules
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.out3Plot import in1DataPlots as iPlot
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def computeParameters(a, b, c):
    """ Compute alpha, beta and mu """

    # computation of paramters
    mu = (a + 4*b + c) / 6.0
    alpha = (4*b + c - 5*a) / (c - a)
    beta = (5*c - a - 4*b) / (c - a)

    return alpha, beta, mu


def computePert(a, b, c, x, alpha, beta):
    """ Compute the CDF and PDF of the Pert distribution using scipy betainc function """

    # Compute pert pdf for testing if the retrieved sample fits the desired distribution
    PDF = (((x - a)**(alpha - 1)) * ((c - x)**(beta-1))) / \
          (sc.beta(alpha, beta) * ((c - a)**(alpha + beta - 1)))

    # compute regularized incomplete beta function for pert distribution using scipy
    z = (x - a) / (c - a)
    CDF = sc.betainc(alpha, beta, z)

    # use scipy.interpolate to create a function that can be used to extract samples
    CDFint = interp1d(CDF, x)

    return PDF, CDF, CDFint


def extractFromCDF(CDF, CDFint, x, cfg):
    """ Extract a sample from the CDF with prescribed steps """

    # extract number of samples using function generated using scipy.interpolate
    # more robust regarding septs vs CDFinterval
    ySampling = np.linspace(0.0, 1.0, int(cfg['sampleSize'])+1)
    sampleVect = CDFint(ySampling)

    return sampleVect


def extractUniform(a, c, x, cfg):
    """ Extract sample of a uniform distriution """

    # load parameters
    sampleSize = int(cfg['sampleSize'])

    # compute interval
    fullInterval = 100000.0 * ((c - a) / sampleSize)
    interval = np.around(fullInterval) * 0.00001

    # Compute sample
    sampleVect = np.zeros(sampleSize+1)
    sampleVect[0] = a
    for m in range(1,sampleSize+1):
        sampleVect[m] = sampleVect[m-1] + interval

    CDF = np.linspace(0, 1, int(cfg['support']))
    CDFInt = interp1d(x, CDF)

    return CDF, CDFInt, sampleVect


def getEmpiricalCDF(sample, CDF, cfg):
    """ Derive empirical CDF using numpy histogram and cumsum """

    binsNo = int(int(cfg['sampleSize']) *0.25)
    hist, binsEd = np.histogram(sample, bins=binsNo)
    CDFEmp = np.cumsum(hist)
    CDFEmpPlot = CDFEmp / CDFEmp[-1]

    return CDFEmpPlot, binsEd[1:]


def getEmpiricalCDFNEW(sample, CDF, cfg):
    """ Derive empirical CDF using sorted sample """

    # sort sample
    sampleSorted = np.sort(sample)
    sampleSize = len(sample)

    ECDF = np.zeros(sampleSize)
    for m in range(sampleSize):
        cumsum = 0
        for l in range(sampleSize):
            if sample[l] <= sampleSorted[m]:
                cumsum = cumsum + 1
        ECDF[m] = cumsum / sampleSize

    return ECDF, sampleSorted
