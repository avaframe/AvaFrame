"""
    Generate sample of values following a beta-pert or a uniform distribution
"""

# load python modules
import numpy as np
import scipy.special as sc
from scipy.interpolate import interp1d
import logging
from scipy.stats import norm


# create local logger
log = logging.getLogger(__name__)


def computeParameters(a, b, c):
    """ Compute alpha, beta and mu """

    # computation of paramters
    if a > b or b > c or a > c:
        message = 'a:%.2f must be smaller than b: %.2f must be smaller than c: %.2f' % (a, b, c)
        log.error(message)
        raise ValueError(message)

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
    if cfg.getboolean('flagMinMax'):
        ySampling = np.linspace(0.0, 1.0, cfg.getint('sampleSize'))
        sampleVect = CDFint(ySampling)
    else:
        ySampling = np.linspace(0.0, 1.0, cfg.getint('sampleSize')+2)
        sampleVect = CDFint(ySampling)
        sampleVect = sampleVect[1:-1]

    return sampleVect


def extractUniform(a, c, x, cfg):
    """ Extract sample of a uniform distriution """

    # load parameters
    sampleSize = cfg.getint('sampleSize')

    if cfg.getboolean('flagMinMax'):
        sampleSize = sampleSize
    else:
        sampleSize = sampleSize + 2

    # compute interval
    fullInterval = 100000.0 * ((c - a) / (sampleSize-1))
    interval = np.around(fullInterval) * 0.00001

    # Compute sample
    sampleVect = np.zeros(sampleSize)
    sampleVect[0] = a
    for m in range(1, sampleSize):
        sampleVect[m] = sampleVect[m-1] + interval

    if cfg.getboolean('flagMinMax'):
        sampleSize = sampleSize
    else:
        sampleVect = sampleVect[1:-1]

    CDF = np.linspace(0, 1, cfg.getint('support'))
    CDFInt = interp1d(CDF, x)

    return CDF, CDFInt, sampleVect


def getEmpiricalCDF(sample):
    """ Derive empirical CDF using numpy histogram and cumsum """

    binsNo = int(len(sample) * 0.25)
    hist, binsEd = np.histogram(sample, bins=binsNo)
    CDFEmp = np.cumsum(hist)
    CDFEmpPlot = CDFEmp / CDFEmp[-1]

    return CDFEmpPlot, binsEd[1:]


def getEmpiricalCDFNEW(sample):
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


def extractNormalDist(cfg):
    """ create a normal distribution from given parameters and draw a sample

        Parameters
        ------------
        cfg: configparser object or dict
            configuration settings for computing distributions, std or ci95, minMaxInterval, flagMinMax

        Returns
        --------
        CDFint: function object
            interpolated CDF function
        sampleVect: numpy array
            sample values
        pdf: numpy array
            pdf values computed for normal distribution and support x
        x: numpy array
            support x
    """

    # load parameters from config
    sampleSize = int(cfg['sampleSize'])
    mean = float(cfg['mean'])
    minMaxInterval = float(cfg['minMaxInterval'])

    # first check if std or ci95 is provided throw error if both
    if cfg['buildType'] == 'ci95':
        std = float(cfg['buildValue']) / 1.96
    elif cfg['buildType'] == 'std':
        std = float(cfg['buildValue'])
    else:
        message = 'buildType: %s not a valid option' % cfg['buildType']
        log.error(message)
        raise AssertionError
    log.info('Compute normal distribution with mean: %.4f and std: %.4f' % (mean, std))

    # compute min and max values of range derived from 99% confidence interval
    # Note: final sample includes these min and max values
    min = norm.interval(alpha=(minMaxInterval/100.), loc=mean, scale=std)[0]
    max = norm.interval(alpha=(minMaxInterval/100.), loc=mean, scale=std)[1]

    # derive normal distribution (pdf and cdf) for range from min to max values
    x = np.linspace(min, max, int(cfg['support']))
    cdf = norm.cdf(x, loc=mean, scale=std)
    pdf = norm.pdf(x, loc=mean, scale=std)

    # create interpolated function of cdf to draw samples from
    CDFint = interp1d(cdf, x)
    # draw sampleSize samples from distribution
    xSample = np.linspace(np.min(cdf), np.max(cdf), sampleSize)
    sampleVect = CDFint(xSample)

    return CDFint, sampleVect, pdf, x
