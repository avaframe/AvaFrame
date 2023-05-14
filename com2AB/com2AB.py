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
import avaframe.out3Plot.outDebugPlots as debPlot

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

        parameterSet = "Small avalanches"

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

        parameterSet = "Custom"

    else:
        log.debug('Using standard Avalanche Setup')
        eqParameters['k1'] = 1.05
        eqParameters['k2'] = -3130.0
        eqParameters['k3'] = 0.0
        eqParameters['k4'] = -2.38
        eqParameters['SD'] = 1.25

        parameterSet = "Standard"

    eqParameters['parameterSet'] = parameterSet
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
    pathDict : dict
        dictionary with AlphaBeta inputs
    dem: dict
        dem dictionary used to get the avaProfile from the avaPath
    splitPoint: dict
        split point dict
    eqParams: dict
        dict containing the AB model parameters (produced by setEqParameters and depends on the com2ABCfg.ini)
    resAB : dict
        dictionary with AlphaBeta model results
    """
    abVersion = '4.1'
    cfgsetup = cfg['ABSETUP']
    smallAva = cfgsetup.getboolean('smallAva')
    resampleDistance = cfgsetup.getfloat('distance')
    betaThresholdDistance = cfgsetup.getfloat('dsMin')
    resAB = {}
    # Extract input file locations
    pathDict = readABinputs(avalancheDir, path2Line=cfgsetup['path2Line'], path2SplitPoint=cfgsetup['path2SplitPoint'])

    log.info("Running com2ABMain model on DEM \n \t %s \n \t with profile \n \t %s ",
             pathDict['demSource'], pathDict['profileLayer'])
    # Read input data for ALPHABETA
    dem = IOf.readRaster(pathDict['demSource'])
    # read line (may contain multiple lines)
    fullAvaPath = shpConv.readLine(pathDict['profileLayer'], pathDict['defaultName'], dem)
    splitPoint = shpConv.readPoints(pathDict['splitPointSource'], dem)

    # Read input setup
    eqParams = setEqParameters(cfg, smallAva)

    NameAva = fullAvaPath['Name']
    StartAva = fullAvaPath['Start']
    LengthAva = fullAvaPath['Length']
    # loop on each feature in the shape file
    for i in range(len(NameAva)):
        name = NameAva[i]
        start = StartAva[i]
        end = start + LengthAva[i]
        # extract individual line
        avaPath = {'sks': fullAvaPath['sks']}
        avaPath['x'] = fullAvaPath['x'][int(start):int(end)]
        avaPath['y'] = fullAvaPath['y'][int(start):int(end)]
        avaPath['name'] = name
        log.info('Running Alpha Beta %s on: %s ', abVersion, name)
        avaProfile = com2ABKern(avaPath, splitPoint, dem, eqParams, resampleDistance, betaThresholdDistance)
        resAB[name] = avaProfile
        if cfg.getboolean('FLAGS', 'fullOut'):
            # saving results to pickle saveABResults(resAB, name)
            savename = name + '_com2AB_eqparam.pickle'
            save_file = pathlib.Path(pathDict['saveOutPath'], savename)
            pickle.dump(eqParams, open(save_file, "wb"))
            log.info('Saving intermediate results to: %s' % (save_file))
            savename = name + '_com2AB_avaProfile.pickle'
            save_file = pathlib.Path(pathDict['saveOutPath'], savename)
            pickle.dump(avaProfile, open(save_file, "wb"))
            log.info('Saving intermediate results to: %s' % (save_file))

    return pathDict, dem, splitPoint, eqParams, resAB


def com2ABKern(avaPath, splitPoint, dem, eqParams, distance, dsMin):
    """ Compute AlpahBeta model for a given avapath

    Call calcABAngles to compute the AlphaBeta model given an input raster (of the dem),
    an avalanche path and split points

    Parameters
    ----------
    avaPath : dict
        dictionary with the name of the avaPath, the x and y coordinates of the
        path
    splitPoint : dict
        dictionary split points
    dem: dict
        dem dictionary used to get the avaProfile from the avaPath
    eqParams: dict
        dict containing the AB model parameters (produced by setEqParameters and depends on the com2ABCfg.ini)
    distance: float
        line resampling distance
    dsMin: float
        threshold distance [m] when looking for the beta point

    Returns
    -------
    avaProfile : dict
        avaPath dictionary with AlphaBeta model results (path became a profile adding the z and s arrays.
        AB runout angles and distances)
    """
    name = avaPath['name']

    # read inputs, ressample ava path
    # make pofile and project split point on path
    avaProfile, projSplitPoint = geoTrans.prepareLine(dem, avaPath, distance, splitPoint)

    if np.isnan(np.sum(avaProfile['z'])):
        raise ValueError('The resampled avalanche path exceeds the dem extent. Try with another path')

    # Sanity check if first element of avaProfile[3,:]
    # (i.e z component) is highest:
    # if not, flip all arrays
    projSplitPoint, avaProfile = geoTrans.checkProfile(avaProfile, projSplitPoint)

    avaProfile['indSplit'] = projSplitPoint['indSplit']  # index of split point

    # run AB model and get angular results
    avaProfile = calcABAngles(avaProfile, eqParams, dsMin)
    # convert the angular results in distances
    avaProfile = calcABDistances(avaProfile, name)

    return avaProfile


def readABinputs(avalancheDir, path2Line='', path2SplitPoint=''):
    """ Fetch inputs for AlpahBeta model

    Get path to AlphaBeta model inputs (dem raster, avalanche path and split points)

    Parameters
    ----------
    avalancheDir : str
        path to directory of avalanche to analyze
    path2Line : pathlib path
        pathlib path to altrnative line
        (if empty, reading the line from the input directory Inputs/LINES/*AB.shp)
    path2SplitPoint : pathlib path
        pathlib path to altrnative splitPoint
        (if empty, reading the point from the input directory Inputs/LINES/*AB.shp)

    Returns
    -------
    pathDict : dict
        dictionary with path to AlphaBeta inputs (dem, avaPath, splitPoint)
    """
    pathDict = {}
    avalancheDir = pathlib.Path(avalancheDir)
    # read avalanche paths for AB
    if path2Line == '':
        refDir = avalancheDir / 'Inputs' / 'LINES'
        profileLayer = list(refDir.glob('*AB*.shp'))
        try:
            message = ('There should be exactly one pathAB.shp file containing (multiple)'
                       + 'avalanche paths in %s /Inputs/LINES/' % avalancheDir)
            assert len(profileLayer) == 1, message
        except AssertionError:
            log.error(message)
            raise
        pathDict['profileLayer'] = profileLayer[0]
    else:
        path2Line = pathlib.Path(path2Line)
        if not path2Line.is_file():
            message = 'No line called: %s' % (path2Line)
            log.error(message)
            raise FileNotFoundError(message)
        pathDict['profileLayer'] = path2Line

    # read DEM
    refDir = avalancheDir / 'Inputs'
    demSource = list(refDir.glob('*.asc'))
    try:
        assert len(demSource) == 1, 'There should be exactly one topography .asc file in %s /Inputs/' % avalancheDir
    except AssertionError:
        raise
    pathDict['demSource'] = demSource[0]

    # read split points
    if path2SplitPoint == '':
        refDir = avalancheDir / 'Inputs' / 'POINTS'
        splitPointSource = list(refDir.glob('*.shp'))
        try:
            message = 'There should be exactly one .shp file containing the split points in %s /Inputs/POINTS/' %  avalancheDir
            assert len(splitPointSource) == 1, message
        except AssertionError:
            raise
        pathDict['splitPointSource'] = splitPointSource[0]
    else:
        path2SplitPoint = pathlib.Path(path2SplitPoint)
        if not path2SplitPoint.is_file():
            message = 'No line called: %s' % (path2SplitPoint)
            log.error(message)
            raise FileNotFoundError(message)
        pathDict['splitPointSource'] = path2SplitPoint

    # make output path
    saveOutPath = avalancheDir / 'Outputs' / 'com2AB'
    if not saveOutPath.exists():
        # log.info('Creating output folder %s', saveOutPath)
        saveOutPath.mkdir(parents=True, exist_ok=True)
    pathDict['saveOutPath'] = saveOutPath

    defaultName = avalancheDir.stem
    pathDict['defaultName'] = defaultName

    return pathDict


def calcABAngles(avaProfile, eqParameters, dsMin):
    """ Kernel function that computes the AlphaBeta model (angular results)
    for a given avaProfile and eqParameters

    Parameters
    ----------
    avaProfile : dict
        dictionary with the name of the avapath, the x, y and z coordinates of
        the path
    eqParameters: dict
        AB parameter dictionary
    dsMin: float
        threshold distance [m] when looking for the Beta point

    Returns
    -------
    avaProfile : dict
        updated avaProfile with alpha, beta and other values resulting from the
        AlphaBeta model computation
    """
    log.debug("Calculating alpha beta")
    k1 = eqParameters['k1']
    k2 = eqParameters['k2']
    k3 = eqParameters['k3']
    k4 = eqParameters['k4']
    SD = eqParameters['SD']

    s = avaProfile['s']
    z = avaProfile['z']

    # prepare find Beta points
    betaValue = 10
    angle, tmp, ds = geoTrans.prepareAngleProfile(betaValue, avaProfile)
    # find the beta point: first point under the beta angle
    # (make sure that the dsMin next meters are also under te beta angle)
    try:
        indBetaPoint = geoTrans.findAngleProfile(tmp, ds, dsMin)
    except IndexError:
        noBetaFoundMessage = 'No Beta point found. Check your pathAB.shp and splitPoint.shp.'
        raise IndexError(noBetaFoundMessage)
    if debugPlot:
        debPlot.plotSlopeAngle(s, angle, indBetaPoint)
        debPlot.plotProfile(s, z, indBetaPoint)

    # Do a quadtratic fit and get the polycom2ABKernnom for 2nd derivative later
    zQuad = np.polyfit(s, z, 2)
    poly = np.poly1d(zQuad)
    # Get H0: max - min for parabola
    H0 = max(poly(s)) - min(poly(s))
    # get beta
    dzBeta = z[0] - z[indBetaPoint]
    beta = np.rad2deg(np.arctan2(dzBeta, s[indBetaPoint]))
    # get Alpha
    alpha = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4

    # get Alpha standard deviations
    SDs = [SD, -1*SD, -2*SD]
    alphaSD = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4 + SDs

    avaProfile['sSplit'] = s[avaProfile['indSplit']]
    avaProfile['indBetaPoint'] = indBetaPoint
    avaProfile['poly'] = poly
    avaProfile['beta'] = beta
    avaProfile['alpha'] = alpha
    avaProfile['SDs'] = SDs
    avaProfile['alphaSD'] = alphaSD
    return avaProfile


def calcABDistances(avaProfile, name):
    """ Compute runout distances and points from angles computed in calcABAngles

    Parameters
    ----------
    avaProfile : dict
        dictionary with the name of the avapath, the x, y and z coordinates of
        the path
    name: str
        profile name

    Returns
    -------
    avaProfile : dict
        updated avaProfile with s index of alpha, and alphaSD points
    """
    s = avaProfile['s']
    z = avaProfile['z']
    sSplit = avaProfile['sSplit']
    alpha = avaProfile['alpha']
    alphaSD = avaProfile['alphaSD']

    # Line down to alpha
    f = z[0] + np.tan(np.deg2rad(-alpha)) * s
    fplus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[0])) * s
    fminus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[1])) * s
    fminus2SD = z[0] + np.tan(np.deg2rad(-alphaSD[2])) * s

    # First it calculates f - g and the corresponding signs
    # using np.sign. Applying np.diff reveals all
    # the positions, where the sign changes (e.g. the lines cross).
    indAlpha = np.argwhere(np.diff(np.sign(f - z))).flatten()
    indAlphaP1SD = np.argwhere(np.diff(np.sign(fplus1SD - z))).flatten()
    indAlphaM1SD = np.argwhere(np.diff(np.sign(fminus1SD - z))).flatten()
    indAlphaM2SD = np.argwhere(np.diff(np.sign(fminus2SD - z))).flatten()

    # Only get the first index past the splitpoint
    try:
        indAlpha = indAlpha[s[indAlpha] > sSplit][0]
    except IndexError:
        log.warning('Alpha out of profile')
        indAlpha = None

    try:
        indAlphaP1SD = indAlphaP1SD[s[indAlphaP1SD] > sSplit][0]
    except IndexError:
        log.warning('+1 SD above beta point')
        indAlphaP1SD = None

    try:
        indAlphaM1SD = indAlphaM1SD[s[indAlphaM1SD] > sSplit][0]
    except IndexError:
        log.warning('-1 SD out of profile')
        indAlphaM1SD = None

    try:
        indAlphaM2SD = indAlphaM2SD[s[indAlphaM2SD] > sSplit][0]
    except IndexError:
        log.warning('-2 SD out of profile')
        indAlphaM2SD = None

    avaProfile['f'] = f
    avaProfile['indAlpha'] = indAlpha
    avaProfile['indAlphaP1SD'] = indAlphaP1SD
    avaProfile['indAlphaM1SD'] = indAlphaM1SD
    avaProfile['indAlphaM2SD'] = indAlphaM2SD

    return avaProfile
