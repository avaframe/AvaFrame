"""
    AIMEC post processing workflow
"""

import logging

# Local imports
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.outAIMEC as outAimec

# create local logger
log = logging.getLogger(__name__)


def AIMEC2Report(pathDict, cfg):
    """ perform AIMEC analysis and generate plots for reports

    Reads the required files location for ana3AIMEC postpocessing
    given a path dictionary to the input files

    Parameters
    ----------
    pathDict : dict
        dictionary with paths to data to analyze
    cfg : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini

    Returns
    -------
    rasterTransfo: dict
        domain transformation information
    newRasters: dict
        raster data expressed in the new coordinates
    resAnalysis: dict
        results of ana3AIMEC analysis
    """

    # Extract input config parameters
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    interpMethod = cfgSetup['interpMethod']

    log.info('Prepare data for post-ptocessing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = aT.makeDomainTransfo(pathDict, cfgSetup)
    ###########################################################################
    # visualisation
    # TODO: needs to be moved somewhere else
    # read reference file
    nRef = pathDict['referenceFile']
    rasterSource = pathDict['ppr'][nRef]

    pressureRaster = IOf.readRaster(rasterSource)
    slRaster = aT.transform(rasterSource, rasterTransfo, interpMethod)
    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = pressureRaster['rasterData']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict, cfgFlags)
    #################################################################

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning pressure data to deskewed raster")
    newRasters['newRasterPPR'] = aT.assignData(pathDict['ppr'],
                                                 rasterTransfo, interpMethod)
    # assign depth data
    log.debug("Assigning depth data to deskewed raster")
    newRasters['newRasterPFD'] = aT.assignData(pathDict['pfd'],
                                              rasterTransfo, interpMethod)
    # assign speed data
    log.debug("Assigning speed data to deskewed raster")
    newRasters['newRasterPFV'] = aT.assignData(pathDict['pfv'],
                                              rasterTransfo, interpMethod)

    # assign dem data
    log.debug("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([pathDict['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.debug('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')

    plotName = outAimec.visuRunoutComp(rasterTransfo, resAnalysis, newRasters, cfgSetup, pathDict, cfgFlags)
    resAnalysis['slCompPlot'] = {'Aimec comparison of mean and max values along path': plotName}
    if cfgFlags.getboolean('flagMass'):
        plotName = outAimec.visuMass(resAnalysis, pathDict, cfgFlags)
        resAnalysis['massAnalysisPlot'] = {'Aimec mass analysis': plotName}

    return rasterTransfo, newRasters, resAnalysis


def mainAIMEC(pathDict, cfg):
    """ Main logic for AIMEC postprocessing

    Reads the required files location for ana3AIMEC postpocessing
    given an avalanche directory
    Make the domain transformation corresponding to the input avalanche path
    Transform 2D fields (pressure, flow depth ...)
    Analyze transformed rasters and mass
    Produce plots and report

    Parameters
    ----------
    pathDict : dict
        dictionary with paths to data to analyze
    cfg : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini

    Returns
    -------
    rasterTransfo: dict
        domain transformation information
    newRasters: dict
        raster data expressed in the new coordinates
    resAnalysis: dict
        results of ana3AIMEC analysis
    """

    # Extract input config parameters
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    interpMethod = cfgSetup['interpMethod']

    log.info('Prepare data for post-ptocessing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = aT.makeDomainTransfo(pathDict, cfgSetup)

    ###########################################################################
    # visualisation
    # TODO: needs to be moved somewhere else
    # read reference file
    nRef = pathDict['referenceFile']
    rasterSource = pathDict[cfgSetup['resType']][nRef]

    raster = IOf.readRaster(rasterSource)
    slRaster = aT.transform(rasterSource, rasterTransfo, interpMethod)
    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = raster['rasterData']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict, cfgFlags)
    #################################################################

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning pressure data to deskewed raster")
    newRasters['newRasterPPR'] = aT.assignData(pathDict['ppr'],
                                                 rasterTransfo, interpMethod)
    # assign depth data
    log.debug("Assigning depth data to deskewed raster")
    newRasters['newRasterPFD'] = aT.assignData(pathDict['pfd'],
                                              rasterTransfo, interpMethod)
    # assign speed data
    if pathDict['pfv']:
        log.debug("Assigning speed data to deskewed raster")
        newRasters['newRasterPFV'] = aT.assignData(pathDict['pfv'],
                                                  rasterTransfo, interpMethod)

    # assign dem data
    log.info("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([pathDict['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.info('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')
    outAimec.visuSimple(rasterTransfo, resAnalysis, newRasters, pathDict, cfgFlags)
    if pathDict['numSim'] == 2:
        outAimec.visuRunoutComp(rasterTransfo, resAnalysis, newRasters, cfgSetup, pathDict, cfgFlags)
        if cfgFlags.getboolean('flagMass'):
            outAimec.visuMass(resAnalysis, pathDict, cfgFlags)
    else:
        outAimec.visuRunoutStat(rasterTransfo, resAnalysis, newRasters, cfgSetup, pathDict, cfgFlags)
    outAimec.resultVisu(cfgSetup, pathDict, cfgFlags, rasterTransfo, resAnalysis)

    # -----------------------------------------------------------
    # write results to file
    # -----------------------------------------------------------
    log.info('Writing results to file')
    flagMass = cfgFlags.getboolean('flagMass')
    outAimec.resultWrite(pathDict, cfgSetup, flagMass, rasterTransfo, resAnalysis)

    return rasterTransfo, newRasters, resAnalysis


def AIMECIndi(pathDict, cfg):
    """ perform AIMEC analysis and generate plots for reports

    Reads the required files location for ana3AIMEC postpocessing
    given a path dictionary to the input files

    Parameters
    ----------
    pathDict : dict
        dictionary with paths to data to analyze
    cfg : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini

    Returns
    -------
    rasterTransfo: dict
        domain transformation information
    newRasters: dict
        raster data expressed in the new coordinates
    resAnalysis: dict
        results of ana3AIMEC analysis
    """

    # Extract input config parameters
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    interpMethod = cfgSetup['interpMethod']
    resType = cfgSetup['resType']

    log.info('Prepare data for post-ptocessing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = aT.makeDomainTransfo(pathDict, cfgSetup)

    ###########################################################################
    # visualisation
    # TODO: needs to be moved somewhere else
    # read reference file
    nRef = pathDict['referenceFile']
    rasterSource = pathDict[resType][nRef]

    anaRaster = IOf.readRaster(rasterSource)
    slRaster = aT.transform(rasterSource, rasterTransfo, interpMethod)
    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = anaRaster['rasterData']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict, cfgFlags)
    #################################################################

    # transform resType_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning data to deskewed raster")
    newRasters['newRaster' + resType.upper()] = aT.assignData(pathDict[resType],
                                                 rasterTransfo, interpMethod)
    # assign dem data
    log.debug("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([pathDict['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.debug('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMECIndi(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')

    outAimec.visuRunoutStat(rasterTransfo, resAnalysis, newRasters, cfgSetup, pathDict, cfgFlags)
    outAimec.resultVisu(cfgSetup, pathDict, cfgFlags, rasterTransfo, resAnalysis)

    return rasterTransfo, newRasters, resAnalysis


# -----------------------------------------------------------
# Aimec analysis tools
# -----------------------------------------------------------
def postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags):
    """ Analyse pressure and depth transformed data

    Analyse pressure depth and speed.
    Calculate runout, Max Peak Pressure, Average PP...
    Get mass and entrainement

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgSetup: configParser objec
        parameters for data analysis
    pathDict: dict
        path to data to analyse
    cfgFlags: congfigParser object
        configuration settings flagMass used

    Returns
    -------
    resAnalysis: dict
        resAnalysis dictionary containing all results:
            -MMPPR: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak pressure
            -PPRCrossMax: 2D numpy array
                    containing for each simulation analyzed the
                    max peak pressure in each cross section
            -PPRCrossMean: 2D numpy array
                    containing for each simulation analyzed the
                    mean peak pressure in each cross section
            -MMPFD: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak flow depth
            -PFDCrossMax: 2D numpy array
                    containing for each simulation analyzed the
                    max peak flow depth in each cross section
            -PFDCrossMean: 2D numpy array
                    containing for each simulation analyzed the
                    mean peak flow depth in each cross section
            -MMPFV: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak flow velocity
            -PFVCrossMax: 2D numpy array
                    containing for each simulation analyzed the
                    max peak flow velocity in each cross section
            -PFVCrossMean: 2D numpy array
                    containing for each simulation analyzed the
                    mean peak flow velocity in each cross section

            -thresholdValue: float
                    threshold value for runout analysis
            -startOfRunoutAreaAngle: float
                    angle of the slope at the beginning of the run-out
                    area (given in input)

            -runout: 2D numpy array
                    containing for each simulation analyzed the x and
                    y coord of the runout point as well as the runout distance
                    measured from the begining of the path. run-out
                    calculated with the MAX pressure in each cross section
            -runoutMean: 2D numpy array
                    containing for each simulation analyzed the x
                    and y coord of the runout point as well as the runout
                    distance measured from the begining of the path.
                    run-out calculated with the MEAN pressure in each cross
                    section
            -elevRel: 1D numpy array
                    containing for each simulation analyzed the
                    elevation of the release area (based on first point with
                    peak field > thresholdValue)
            -deltaH: 1D numpy array
                    containing for each simulation analyzed the
                    elevation fall difference between elevRel and altitude of
                    run-out point

            -TP: float
                    ref = True sim2 = True
            -FN: float
                    ref = False sim2 = True
            -FP: float
                    ref = True sim2 = False
            -TN: float
                    ref = False sim2 = False

            if mass analysis is performed
            -relMass: 1D numpy array
                    containing for each simulation analyzed the
                    release mass
            -entMass: 1D numpy array
                    containing for each simulation analyzed the
                    entrained mass
            -finalMass: 1D numpy array
                    containing for each simulation analyzed the
                    final mass
            -relativMassDiff: 1D numpy array
                    containing for each simulation analyzed
                    the final mass diff with ref (in %)
            -growthIndex: 1D numpy array
                    containing for each simulation analyzed the
                    growth index
            -growthGrad: 1D numpy array
                    containing for each simulation analyzed the
                    growth gradient
    """
    # read inputs
    flagMass = cfgFlags.getboolean('flagMass')
    dataPressure = newRasters['newRasterPPR']
    dataDepth = newRasters['newRasterPFD']
    dataSpeed = newRasters['newRasterPFV']
    transformedDEMRasters = newRasters['newRasterDEM']

    resType = cfgSetup['resType']
    thresholdValue = cfgSetup.getfloat('thresholdValue')

    resultsAreaAnalysis = {}
    resultsAreaAnalysis['resType'] = resType
    # analyize all fields
    resultsAreaAnalysis = aT.analyzeField(rasterTransfo, dataPressure, 'ppr', resultsAreaAnalysis)
    resultsAreaAnalysis = aT.analyzeField(rasterTransfo, dataDepth, 'pfd', resultsAreaAnalysis)
    resultsAreaAnalysis = aT.analyzeField(rasterTransfo, dataSpeed, 'pfv', resultsAreaAnalysis)

    # compute runout based on resType
    runout, runoutMean, elevRel, deltaH = aT.computeRunOut(rasterTransfo, thresholdValue, resultsAreaAnalysis, transformedDEMRasters)

    runoutLength = runout[0]
    TP, FN, FP, TN, compPlotPath = aT.analyzeArea(rasterTransfo, runoutLength, resultsAreaAnalysis[resType]['transformedRasters'], cfgSetup, pathDict, cfgFlags)

    # affect values to output dictionary
    resAnalysis = {}
    resAnalysis['resType'] = resType
    resAnalysis['runout'] = runout
    resAnalysis['runoutMean'] = runoutMean
    resAnalysis['MMPPR'] = resultsAreaAnalysis['ppr']['maxaCrossMax']
    resAnalysis['MMPFD'] = resultsAreaAnalysis['pfd']['maxaCrossMax']
    resAnalysis['MMPFV'] = resultsAreaAnalysis['pfv']['maxaCrossMax']
    resAnalysis['elevRel'] = elevRel
    resAnalysis['deltaH'] = deltaH
    resAnalysis['PPRCrossMax'] = resultsAreaAnalysis['ppr']['aCrossMax']
    resAnalysis['PPRCrossMean'] = resultsAreaAnalysis['ppr']['aCrossMean']
    resAnalysis['PFDCrossMax'] = resultsAreaAnalysis['pfd']['aCrossMax']
    resAnalysis['PFDCrossMean'] = resultsAreaAnalysis['pfd']['aCrossMean']
    resAnalysis['PFVCrossMax'] = resultsAreaAnalysis['pfv']['aCrossMax']
    resAnalysis['PFVCrossMean'] = resultsAreaAnalysis['pfv']['aCrossMean']
    resAnalysis['thresholdValue'] = thresholdValue
    resAnalysis['startOfRunoutAreaAngle'] = rasterTransfo['startOfRunoutAreaAngle']
    resAnalysis['TP'] = TP
    resAnalysis['FN'] = FN
    resAnalysis['FP'] = FP
    resAnalysis['TN'] = TN
    resAnalysis['areasPlot'] = {'Aimec area analysis': compPlotPath}

    if flagMass:
        # perform mass analysis
        fnameMass = pathDict['massBal']
        resAnalysis = aT.analyzeMass(fnameMass, resAnalysis)

    return resAnalysis


def postProcessAIMECIndi(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags):
    """ Analyse one peak field transformed data

    Analyse peak field set in resType e.g. ppr
    Calculate runout, Max Peak Value, Average PP...
    Mass is currently not included

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgSetup: dict
        parameters for data analysis
    pathDict: dict
        path to data to analyse
    cfgFlags: congfigParser object
        configuration settings flagMass used

    Returns
    -------
    resAnalysis: dict
        resAnalysis dictionary containing all results:

    """
    # read inputs
    resType = cfgSetup['resType']
    thresholdValue = cfgSetup.getfloat('thresholdValue')
    dataResType = newRasters['newRaster' + resType.upper()]
    transformedDEMRasters = newRasters['newRasterDEM']

    resultsAreaAnalysis = {}
    resultsAreaAnalysis['resType'] = resType
    # get max and mean values along path for cross profiles
    resultsAreaAnalysis = aT.analyzeField(rasterTransfo, dataResType, resType, resultsAreaAnalysis)

    # compute runout based on resType
    runout, runoutMean, elevRel, deltaH = aT.computeRunOut(rasterTransfo, thresholdValue, resultsAreaAnalysis, transformedDEMRasters)

    runoutLength = runout[0]
    TP, FN, FP, TN, compPlotPath = aT.analyzeArea(rasterTransfo, runoutLength, dataResType, cfgSetup, pathDict, cfgFlags)

    # affect values to output dictionary
    resAnalysis = {}
    resAnalysis['resType'] = resType
    resAnalysis['runout'] = runout
    resAnalysis['runoutMean'] = runoutMean
    resAnalysis['MM' + resType.upper()] = resultsAreaAnalysis[resType]['maxaCrossMax']
    resAnalysis['elevRel'] = elevRel
    resAnalysis['deltaH'] = deltaH
    resAnalysis[resType.upper() + 'CrossMax'] = resultsAreaAnalysis[resType]['aCrossMax']
    resAnalysis[resType.upper() + 'CrossMean'] = resultsAreaAnalysis[resType]['aCrossMean']
    resAnalysis['thresholdValue'] = thresholdValue
    resAnalysis['startOfRunoutAreaAngle'] = rasterTransfo['startOfRunoutAreaAngle']
    resAnalysis['TP'] = TP
    resAnalysis['FN'] = FN
    resAnalysis['FP'] = FP
    resAnalysis['TN'] = TN
    resAnalysis['areasPlot'] = {'Aimec area analysis': compPlotPath}

    return resAnalysis


def aimecRes2ReportDict(resAnalysis, reportD, benchD, refSim):
    """ gather aimec results and append them to report dictionary """

    if refSim == 0:
        compSim = 1
    elif refSim == 1:
        compSim = 0
    else:
        log.warning('Reference File out of range - needs to be 0 or 1')

    reportD['Aimec analysis'] ={'type': 'list', 'runout [m]': resAnalysis['runout'][0][compSim],
    'max peak pressure [kPa]': resAnalysis['MMPPR'][compSim], 'max peak flow depth [m]': resAnalysis['MMPFD'][compSim],
    'max peak flow velocity [ms-1]': resAnalysis['MMPFV'][compSim]}

    benchD['Aimec analysis'] ={'type': 'list', 'runout [m]': resAnalysis['runout'][0][refSim],
    'max peak pressure [kPa]': resAnalysis['MMPPR'][refSim], 'max peak flow depth [m]': resAnalysis['MMPFD'][refSim],
    'max peak flow velocity [ms-1]': resAnalysis['MMPFV'][refSim]}
    # add mass info
    if "relMass" in resAnalysis:
        reportD['Aimec analysis'].update({'Initial mass [kg]': resAnalysis['relMass'][compSim]})
        reportD['Aimec analysis'].update({'Final mass [kg]': resAnalysis['finalMass'][compSim]})
        reportD['Aimec analysis'].update({'Entrained mass [kg]': resAnalysis['entMass'][compSim]})
        benchD['Aimec analysis'].update({'Initial mass [kg]': resAnalysis['relMass'][refSim]})
        benchD['Aimec analysis'].update({'Final mass [kg]': resAnalysis['finalMass'][refSim]})
        benchD['Aimec analysis'].update({'Entrained mass [kg]': resAnalysis['entMass'][refSim]})

    return reportD, benchD
