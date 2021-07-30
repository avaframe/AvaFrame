"""
    Main module for Alpha beta

"""

import logging
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pathlib

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def setEqParameters(cfg, smallAva):
    """ Set alpha beta equation parameters

    Set alpha beta equation parameters to
    - standard (default)
    - small avalanche (if smallAva==True)
    - custom (if cfgsetup('customParam') is True use the custom values provided in the ini file)

    Parameters
    ----------
    cfg : configParser
        if cfgsetup('customParam') is True, the custom parameters provided
        in the cfg are used in theAlphaBeta equation. (provide all 5 k1, k2, k3, k4 and SD values)
    smallAva : boolean
        True if the small avallanche AlphaBeta equation parameters should be used

    Returns
    -------
    eqParameters : dict
        k1, k2, k3, k4 and SD values to be used in the AlphaBeta equation
    """

    cfgsetup = cfg['ABSETUP']
    eqParameters = {}

    if smallAva is True:
        log.debug('Using small Avalanche Setup')
        eqParameters['k1'] = 0.933
        eqParameters['k2'] = 0.0
        eqParameters['k3'] = 0.0088
        eqParameters['k4'] = -5.02
        eqParameters['SD'] = 2.36

        ParameterSet = "Small avalanches"

    elif cfgsetup.getboolean('customParam'):
        log.debug('Using custom Avalanche Setup:')
        parameterName = ['k1', 'k2', 'k3', 'k4', 'SD']
        for paramName in parameterName:
            if cfg.has_option('ABSETUP', paramName):
                eqParameters[paramName] = cfgsetup.getfloat(paramName)
            else:
                message = 'Custom parameter %s is missing in the configuration file' % paramName
                log.error(message)
                raise KeyError(message)

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
    """ Main AlphaBeta model function

    Loops on the given AvaPaths and runs com2AB to compute AlpahBeta model

    Parameters
    ----------
    cfg : configparser
        configparser with all requiered fields in com2ABCfg.ini
    avalancheDir : str
        path to directory of avalanche to analyze

    Returns
    -------
    resAB : dict
        dictionary with AlphaBeta model results
    """
    abVersion = '4.1'
    cfgsetup = cfg['ABSETUP']
    smallAva = cfgsetup.getboolean('smallAva')
    resAB = {}
    # Extract input file locations
    pathDict = readABinputs(avalancheDir)

    log.info("Running com2ABMain model on DEM \n \t %s \n \t with profile \n \t %s ",
             pathDict['demSource'], pathDict['profileLayer'])

    resAB['saveOutPath'] = pathDict['saveOutPath']
    # Read input data for ALPHABETA
    dem = IOf.readRaster(pathDict['demSource'])
    resAB['dem'] = dem
    AvaPath = shpConv.readLine(pathDict['profileLayer'], pathDict['defaultName'],
                               dem)
    resAB['AvaPath'] = AvaPath
    resAB['splitPoint'] = shpConv.readPoints(pathDict['splitPointSource'], dem)

    # Read input setup
    eqParams = setEqParameters(cfg, smallAva)
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
        resAB = com2ABKern(avapath, resAB, cfgsetup.getfloat('distance'), cfgsetup.getfloat('dsMin'))

        if cfg.getboolean('FLAGS', 'fullOut'):
            # saving results to pickle saveABResults(resAB, name)
            savename = name + '_com2AB_eqparam.pickle'
            save_file = pathlib.Path(pathDict['saveOutPath'], savename)
            pickle.dump(resAB['eqParams'], open(save_file, "wb"))
            log.info('Saving intermediate results to: %s' % (save_file))
            savename = name + '_com2AB_eqout.pickle'
            save_file = pathlib.Path(pathDict['saveOutPath'], savename)
            pickle.dump(resAB[name], open(save_file, "wb"))
            log.info('Saving intermediate results to: %s' % (save_file))

    return resAB


def com2ABKern(avapath, resAB, distance, dsMin):
    """ Compute AlpahBeta model for a given avapath

    Call calcAB to compute the AlphaBeta model given an input raster (of the dem),
    an avalanche path and split points

    Parameters
    ----------
    avapath : dict
        dictionary with the name of the avapath, the x and y coordinates of the
        path
    resAB : dict
        dictionary with AlphaBeta model intput parameters as well as the
        results of already computed avapath.
    distance: float
        line resampling distance
    dsMin: float
        threshold distance [m] for looking for the 10° point

    Returns
    -------
    resAB : dict
        dictionary with AlphaBeta model results
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

    eqOut = calcAB(AvaProfile, eqParams, dsMin)
    resAB[name] = eqOut

    return resAB


def readABinputs(avalancheDir):
    """ Fetch inputs for AlpahBeta model

    Get path to AlphaBeta model inputs (dem raster, avalanche path and split points)

    Parameters
    ----------
    avalancheDir : str
        path to directory of avalanche to analyze

    Returns
    -------
    pathDict : dict
        dictionary with path to AlphaBeta inputs (dem, avaPath, splitPoint)
    """
    pathDict = {}
    avalancheDir = pathlib.Path(avalancheDir)
    # read avalanche paths for AB
    refDir = avalancheDir / 'Inputs' / 'LINES'
    profileLayer = list(refDir.glob('*AB*.shp'))
    try:
        message = 'There should be exactly one pathAB.shp file containing (multiple) avalanche paths in %s /Inputs/LINES/' % avalancheDir
        assert len(profileLayer) == 1, message
    except AssertionError:
        raise
    pathDict['profileLayer'] = profileLayer[0]

    # read DEM
    refDir = avalancheDir / 'Inputs'
    demSource = list(refDir.glob('*.asc'))
    try:
        assert len(demSource) == 1, 'There should be exactly one topography .asc file in %s /Inputs/' % avalancheDir
    except AssertionError:
        raise
    pathDict['demSource'] = demSource[0]

    # read split points
    refDir = avalancheDir / 'Inputs' / 'POINTS'
    splitPointSource = list(refDir.glob('*.shp'))
    try:
        message = 'There should be exactly one .shp file containing the split points in %s /Inputs/POINTS/' %  avalancheDir
        assert len(splitPointSource) == 1, message
    except AssertionError:
        raise
    pathDict['splitPointSource'] = splitPointSource[0]

    # make output path
    saveOutPath = avalancheDir / 'Outputs' / 'com2AB'
    if not saveOutPath.exists():
        # log.info('Creating output folder %s', saveOutPath)
        saveOutPath.mkdir(parents=True, exist_ok=True)
    pathDict['saveOutPath'] = saveOutPath

    defaultName = avalancheDir.stem
    pathDict['defaultName'] = defaultName

    return pathDict


def calcAB(AvaProfile, eqParameters, dsMin):
    """ Kernel function that computes the AlphaBeta model
    for a given AvaProfile and eqParameters

    Parameters
    ----------
    AvaProfile : dict
        dictionary with the name of the avapath, the x, y and z coordinates of
        the path
    eqParameters: dict
        AB parameter dictionary
    dsMin: float
        threshold distance [m] for looking for the 10° point

    Returns
    -------
    AvaProfile : dict
        updated AvaProfile with alpha, beta and other values resulting from the
        AlphaBeta model computation
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
    angle, tmp, ds = geoTrans.prepareAngleProfile(betaValue, AvaProfile)
    # find the beta point: first point under 10°
    # (make sure that the 30 next meters are also under 10°)
    ids10Point = geoTrans.findAngleProfile(tmp, ds, dsMin)
    if debugPlot:
        plt.figure(figsize=(10, 6))
        plt.plot(s, angle, '.k')
        plt.plot(s[ids10Point], angle[ids10Point], 'or')
        plt.axhline(y=10, color='0.8',
                    linewidth=1, linestyle='-.', label='10° line')

        plt.figure(figsize=(10, 6))
        plt.plot(s, z)
        plt.axvline(x=s[ids10Point], color='0.8',
                    linewidth=1, linestyle='-.')
        plt.show()
        plt.close()

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
