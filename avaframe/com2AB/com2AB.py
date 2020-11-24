"""
    Main logic for Alpha beta computational module

    This file is part of Avaframe.
"""

import os
import glob
import logging
import numpy as np
import matplotlib.pyplot as plt

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def setEqParameters(smallAva, customParam):
    """Set alpha beta equation parameters to
    - standard (default)
    - small avalanche
    - custom
    """

    eqParameters = {}

    if smallAva is True:
        log.debug('Using small Avalanche Setup')
        eqParameters['k1'] = 0.933
        eqParameters['k2'] = 0.0
        eqParameters['k3'] = 0.0088
        eqParameters['k4'] = -5.02
        eqParameters['SD'] = 2.36

        ParameterSet = "Small avalanches"

    elif customParam:
        log.debug('Using custom Avalanche Setup')
        eqParameters['k1'] = customParam['k1']
        eqParameters['k2'] = customParam['k2']
        eqParameters['k3'] = customParam['k3']
        eqParameters['k4'] = customParam['k4']
        eqParameters['SD'] = customParam['SD']

        ParameterSet = "Custom"

    else:
        log.debug('Using standard Avalanche Setup')
        eqParameters['k1'] = 1.05
        eqParameters['k2'] = -3130.0
        eqParameters['k3'] = 0.0
        eqParameters['k4'] = -2.38
        eqParameters['SD'] = 1.25

        ParameterSet = "Standard"
    eqParameters['ParameterSet'] = ParameterSet
    return eqParameters


def com2ABMain(cfg, avalancheDir):
    """ Run main model
    Loops on the given AvaPath and runs com2AB to compute AlpahBeta model
    Inputs : dem header and rater (as np array),
            AvaPath and Split points .shp file,
            optional output save path,
            avalanche type,
            reamplind lenght for the AvaPath
    Outputs : writes raw results to saveOutPath
    """
    abVersion = '4.1'
    cfgsetup = cfg['ABSETUP']
    smallAva = cfgsetup.getboolean('smallAva')
    customParam = cfgsetup.getboolean('customParam')
    distance = float(cfgsetup['distance'])
    resAB = {}
    # Extract input file locations
    cfgPath = readABinputs(avalancheDir)

    log.info("Running com2ABMain model on DEM \n \t %s \n \t with profile \n \t %s ",
             cfgPath['demSource'], cfgPath['profileLayer'])

    resAB['saveOutPath'] = cfgPath['saveOutPath']
    # Read input data for ALPHABETA
    dem = IOf.readRaster(cfgPath['demSource'])
    resAB['dem'] = dem
    AvaPath = shpConv.readLine(cfgPath['profileLayer'], cfgPath['defaultName'],
                               dem)
    resAB['AvaPath'] = AvaPath
    resAB['splitPoint'] = shpConv.readPoints(cfgPath['splitPointSource'], dem)

    # Read input setup
    eqParams = setEqParameters(smallAva, customParam)
    resAB['eqParams'] = eqParams

    NameAva = AvaPath['Name']
    StartAva = AvaPath['Start']
    LengthAva = AvaPath['Length']

    for i in range(len(NameAva)):
        name = NameAva[i]
        start = StartAva[i]
        end = start + LengthAva[i]
        avapath = {}
        avapath['x'] = AvaPath['x'][int(start):int(end)]
        avapath['y'] = AvaPath['y'][int(start):int(end)]
        avapath['Name'] = name
        log.info('Running Alpha Beta %s on: %s ', abVersion, name)
        resAB = com2AB(avapath, resAB, distance)

    return resAB


def com2AB(avapath, resAB, distance):
    """ Computes the AlphaBeta model given an input raster (of the dem),
    an avalanche path and split points
    Inputs : dem header and rater (as np array),
            single avapath as np array,
            Split points as np array,
            output save path,
            avalanche type,
            resamplind lenght for the Avapath
    Outputs : writes raw results to OutPath
    """
    dem = resAB['dem']
    splitPoint = resAB['splitPoint']
    name = avapath['Name']
    eqParams = resAB['eqParams']

    # read inputs, ressample ava path
    # make pofile and project split point on path
    AvaProfile, projSplitPoint = geoTrans.prepareLine(
        dem, avapath, distance, splitPoint)

    if np.isnan(np.sum(AvaProfile['z'])):
        raise ValueError('The resampled avalanche path exceeds the dem extent.'
                         + 'Try with another path')

    # Sanity check if first element of AvaProfile[3,:]
    # (i.e z component) is highest:
    # if not, flip all arrays
    projSplitPoint, AvaProfile = geoTrans.checkProfile(AvaProfile,
                                                       projSplitPoint)

    AvaProfile['indSplit'] = projSplitPoint['indSplit']  # index of split point

    eqOut = calcAB(AvaProfile, eqParams)
    resAB[name] = eqOut

    return resAB


def readABinputs(cfgAva):

    cfgPath = {}
    # read avalanche paths for AB
    profileLayer = glob.glob(cfgAva + '/Inputs/LINES/*AB*.shp')
    try:
        message = 'There should be exactly one pathAB.shp file containing (multiple) avalanche paths in ' + cfgAva + '/Inputs/LINES/'
        assert len(profileLayer) == 1, message
    except AssertionError:
        raise
    cfgPath['profileLayer'] = ''.join(profileLayer)

    # read DEM
    demSource = glob.glob(cfgAva + '/Inputs/*.asc')
    try:
        assert len(demSource) == 1, 'There should be exactly one topography .asc file in ' + \
            cfgAva + '/Inputs/'
    except AssertionError:
        raise
    cfgPath['demSource'] = ''.join(demSource)

    # read split points
    splitPointSource = glob.glob(cfgAva + '/Inputs/POINTS/*.shp')
    try:
        message = 'There should be exactly one .shp file containing the split points in ' + cfgAva + '/Inputs/POINTS/'
        assert len(splitPointSource) == 1, message
    except AssertionError:
        raise
    cfgPath['splitPointSource'] = ''.join(splitPointSource)

    # make output path
    saveOutPath = os.path.join(cfgAva, 'Outputs/com2AB/')
    if not os.path.exists(saveOutPath):
        # log.info('Creating output folder %s', saveOutPath)
        os.makedirs(saveOutPath)
    cfgPath['saveOutPath'] = saveOutPath

    defaultName = str(cfgAva).split('/')[-1]
    cfgPath['defaultName'] = defaultName

    return cfgPath


def calcAB(AvaProfile, eqParameters):
    """
    Calculate Alpha Beta for data in eqInput according to chosen eqParameters
    """
    log.debug("Calculating alpha beta")
    k1 = eqParameters['k1']
    k2 = eqParameters['k2']
    k3 = eqParameters['k3']
    k4 = eqParameters['k4']
    SD = eqParameters['SD']

    s = AvaProfile['s']
    z = AvaProfile['z']

    # prepare find Beta points
    betaValue = 10
    angle, tmp, deltaInd = geoTrans.prepareAngleProfile(betaValue, AvaProfile)

    # find the beta point: first point under 10°
    # (make sure that the 30 next meters are also under 10°)
    ids10Point = geoTrans.findAngleProfile(tmp, deltaInd)
    if debugPlot:
        plt.figure(figsize=(10, 6))
        plt.plot(s, angle)
        plt.plot(s[ids10Point], angle[ids10Point], 'or')
        plt.axhline(y=10, color='0.8',
                    linewidth=1, linestyle='-.', label='10° line')

        plt.figure(figsize=(10, 6))
        plt.plot(s, z)
        plt.axvline(x=s[ids10Point], color='0.8',
                    linewidth=1, linestyle='-.')
        plt.show()
        plt.close('all')

    # Do a quadtratic fit and get the polynom for 2nd derivative later
    zQuad = np.polyfit(s, z, 2)
    poly = np.poly1d(zQuad)
    # Get H0: max - min for parabola
    H0 = max(poly(s)) - min(poly(s))
    # get beta
    dzBeta = z[0] - z[ids10Point]
    beta = np.rad2deg(np.arctan2(dzBeta, s[ids10Point]))
    # get Alpha
    alpha = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4

    # get Alpha standard deviations
    SDs = [SD, -1*SD, -2*SD]
    alphaSD = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4 + SDs

    AvaProfile['CuSplit'] = s[AvaProfile['indSplit']]
    AvaProfile['ids10Point'] = ids10Point
    AvaProfile['poly'] = poly
    AvaProfile['beta'] = beta
    AvaProfile['alpha'] = alpha
    AvaProfile['SDs'] = SDs
    AvaProfile['alphaSD'] = alphaSD
    return AvaProfile
