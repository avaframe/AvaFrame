"""
    Main module for Alpha beta

"""

import logging
import numpy as np
import pickle
import pathlib

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.rasterUtils as IOf
import avaframe.in1Data.getInput as gI
import avaframe.out3Plot.outDebugPlots as debPlot

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def setEqParameters(cfg):
    """ Set alpha beta equation parameters

    Set alpha beta equation parameters to
    - standard (default values from ini file)
    - small avalanche (if smallAva==True)

    Parameters
    ----------
    cfg : configParser
        cfg that's providing all 5 k1, k2, k3, k4 and SD values

    Returns
    -------
    eqParameters : dict
        k1, k2, k3, k4 and SD values to be used in the AlphaBeta equation
    """

    cfgsetup = cfg['ABSETUP']
    eqParameters = {}
    smallAva = cfgsetup.getboolean('smallAva')

    if smallAva is True:
        log.debug('Using small Avalanche Setup')
        eqParameters['k1'] = cfgsetup.getfloat('k1_small')
        eqParameters['k2'] = cfgsetup.getfloat('k2_small')
        eqParameters['k3'] = cfgsetup.getfloat('k3_small')
        eqParameters['k4'] = cfgsetup.getfloat('k4_small')
        eqParameters['SD'] = cfgsetup.getfloat('SD_small')

        parameterSet = "Small avalanches"

    else:
        log.debug('Using configuration file Avalanche Setup')
        eqParameters['k1'] = cfgsetup.getfloat('k1')
        eqParameters['k2'] = cfgsetup.getfloat('k2')
        eqParameters['k3'] = cfgsetup.getfloat('k3')
        eqParameters['k4'] = cfgsetup.getfloat('k4')
        eqParameters['SD'] = cfgsetup.getfloat('SD')

        parameterSet = "Standard"

    eqParameters['parameterSet'] = parameterSet
    return eqParameters


def com2ABMain(cfg, avalancheDir):
    """ Main AlphaBeta model function

    Loops on the given AvaPaths and runs com2AB to compute AlphaBeta model

    Parameters
    ----------
    cfg : configparser
        configparser with all required fields in com2ABCfg.ini
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
    defaultName = 'line' #set generic name in case no 'Name' attribute exist
    fullAvaPath = shpConv.readLine(pathDict['profileLayer'], defaultName, dem)
    splitPoint = shpConv.readPoints(pathDict['splitPointSource'], dem)

    # Read input setup
    eqParams = setEqParameters(cfg)

    NameAva = fullAvaPath['Name']
    StartAva = fullAvaPath['Start']
    LengthAva = fullAvaPath['Length']
    # loop on each feature in the shape file
    for i in range(len(NameAva)):
        name = NameAva[i] + str(i + 1) if NameAva[i] == defaultName else NameAva[i]
        start = StartAva[i]
        end = start + LengthAva[i]
        # extract individual line
        avaPath = {'sks': fullAvaPath['sks'],
                   'x': fullAvaPath['x'][int(start):int(end)],
                   'y': fullAvaPath['y'][int(start):int(end)],
                   'name': name}
        log.info('Running Alpha Beta %s on: %s ', abVersion, name)
        avaProfile = com2ABKern(avaPath, splitPoint, dem, eqParams, resampleDistance, betaThresholdDistance)
        resAB[name] = avaProfile
        if cfg.getboolean('FLAGS', 'fullOut'):
            # saving results to pickle saveABResults(resAB, name)
            savename = name + '_com2AB_eqparam.pickle'
            save_file = pathlib.Path(pathDict['saveOutPath'], savename)
            pickle.dump(eqParams, open(save_file, "wb"))
            log.info('Saving intermediate results to: %s' % save_file)
            savename = name + '_com2AB_avaProfile.pickle'
            save_file = pathlib.Path(pathDict['saveOutPath'], savename)
            pickle.dump(avaProfile, open(save_file, "wb"))
            log.info('Saving intermediate results to: %s' % save_file)

    return pathDict, dem, splitPoint, eqParams, resAB


def com2ABKern(avaPath, splitPoint, dem, eqParams, distance, dsMin):
    """ Compute AlphaBeta model for a given AvaPath

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

    # read inputs, resample ava path
    # make pofile and project split point on path
    avaProfile, projSplitPoint = geoTrans.prepareLineStrict(dem, avaPath, distance, splitPoint)

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
    avaProfile = calcABDistances(avaProfile)

    return avaProfile


def readABinputs(avalancheDir, path2Line='', path2SplitPoint=''):
    """ Fetch inputs for AlphaBeta model

    Get path to AlphaBeta model inputs (dem raster, avalanche path and split points)

    Parameters
    ----------
    avalancheDir : str
        path to directory of avalanche to analyze
    path2Line : pathlib path
        to alternative line
        (if empty, reading the line from the input directory Inputs/LINES/yourNameAB.shp)
    path2SplitPoint : pathlib path
        to alternative splitPoint
        (if empty, reading the point from the input directory Inputs/LINES/yourNameAB.shp)

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
            message = 'No line called: %s' % path2Line
            log.error(message)
            raise FileNotFoundError(message)
        pathDict['profileLayer'] = path2Line

    # read DEM
    pathDict['demSource'] = gI.getDEMPath(avalancheDir)

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
            message = 'No line called: %s' % path2SplitPoint
            log.error(message)
            raise FileNotFoundError(message)
        pathDict['splitPointSource'] = path2SplitPoint

    # make output path
    saveOutPath = avalancheDir / 'Outputs' / 'com2AB'
    if not saveOutPath.exists():
        # log.info('Creating output folder %s', saveOutPath)
        saveOutPath.mkdir(parents=True, exist_ok=True)
    pathDict['saveOutPath'] = saveOutPath

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

    # Do a quadratic fit and get the polycom2ABKernnom for 2nd derivative later
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


def calcABDistances(avaProfile):
    """ Compute runout distances and points from angles computed in calcABAngles

    Parameters
    ----------
    avaProfile : dict
        dictionary with the name of the avapath, the x, y and z coordinates of
        the path

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
