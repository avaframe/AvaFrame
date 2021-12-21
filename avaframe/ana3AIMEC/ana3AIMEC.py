"""
    AIMEC post processing workflow
"""

import logging
import numpy as np

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

    log.info('Prepare data for post-processing')
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
    inputData['xyHeader'] = raster['header']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict)
    #################################################################

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning pressure data to deskewed raster")
    newRasters['newRasterPPR'] = aT.assignData(pathDict['ppr'],  rasterTransfo, interpMethod)
    # assign depth data
    log.debug("Assigning depth data to deskewed raster")
    newRasters['newRasterPFD'] = aT.assignData(pathDict['pfd'], rasterTransfo, interpMethod)
    # assign speed data
    log.debug("Assigning speed data to deskewed raster")
    newRasters['newRasterPFV'] = aT.assignData(pathDict['pfv'], rasterTransfo, interpMethod)

    # assign dem data
    log.debug("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([pathDict['demSource']], rasterTransfo, interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.debug('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMEC(rasterTransfo, newRasters, cfgSetup, pathDict, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')

    plotName = outAimec.visuRunoutComp(rasterTransfo, resAnalysis, cfgSetup, pathDict)
    resAnalysis['slCompPlot'] = {'Aimec comparison of mean and max values along path': plotName}
    if cfgFlags.getboolean('flagMass') and (pathDict['numSim'] > 1):
        plotName = outAimec.visuMass(resAnalysis, pathDict)
        resAnalysis['massAnalysisPlot'] = {'Aimec mass analysis': plotName}

    return rasterTransfo, newRasters, resAnalysis


def mainAIMEC(pathDict, inputsDF, cfg):
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
    refSimulationName = pathDict['refSimulation']

    log.info('Prepare data for post-processing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = aT.makeDomainTransfo(pathDict, inputsDF, cfgSetup)

    # read reference file
    refSimulationName = pathDict['refSimulation']
    refSimulation = inputsDF[inputsDF['simName'] == refSimulationName]

    # ####################################################
    # visualisation
    # TODO: needs to be moved somewhere else

    rasterSource = refSimulation.iloc[0][cfgSetup['resType']]
    raster = IOf.readRaster(rasterSource)
    slRaster = aT.transform(rasterSource, rasterTransfo, interpMethod)
    newRasters = {}
    log.info("Assigning dem data to deskewed raster")
    newRasters['newRasterDEM'] = aT.transform(pathDict['demSource'], rasterTransfo, interpMethod)

    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = raster['rasterData']
    inputData['xyHeader'] = raster['header']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict)
    # ###########################################################

    # postprocess reference
    inputsDFrow = inputsDF.loc[inputsDF['simName'] == refSimulationName].squeeze()
    timeMass = 0
    resAnalysisDF = inputsDF[['simName']].copy()
    resAnalysisDF, newRasters, timeMass = postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDFrow, newRasters, timeMass, refSimulationName, resAnalysisDF)

    # postprocess other simulations
    for index, inputsDFrow in inputsDF.iterrows():
        simName = inputsDFrow['simName']
        if simName != refSimulationName:
            resAnalysisDF, newRasters, timeMass = postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDFrow, newRasters, timeMass, simName, resAnalysisDF)

    # -----------------------------------------------------------
    # result visualisation + report
    # ToDo: should we move this somewere else, this is now just plotting, it should be accessible from outside
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')
    outAimec.visuSimple(cfgSetup, rasterTransfo, resAnalysisDF, newRasters, pathDict)
    if len(resAnalysisDF.index) == 2:
        outAimec.visuRunoutComp(rasterTransfo, resAnalysisDF, cfgSetup, pathDict)
    else:
        outAimec.visuRunoutStat(rasterTransfo, resAnalysisDF, newRasters, cfgSetup, pathDict)

    outAimec.resultVisu(cfgSetup, inputsDF, pathDict, cfgFlags, rasterTransfo, resAnalysisDF)

    # -----------------------------------------------------------
    # write results to file
    # -----------------------------------------------------------
    log.info('Writing results to file')
    outAimec.resultWrite(pathDict, cfg, rasterTransfo, resAnalysisDF)

    return rasterTransfo, newRasters, resAnalysisDF


def postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDFrow, newRasters, timeMass, simName, resAnalysisDF):
    """ Apply domain transformation and analyse pressure, thickness and velocity data

    Apply the domain tranformation to peak results
    Analyse pressure thickness and speed.
    Calculate runout, Max Peak Pressure, Average PP...
    Get mass and entrainement

    Parameters
    ----------
    cfg: configParser objec
        parameters for AIMEC analysis
    rasterTransfo: dict
        transformation information
    pathDict: dict
        path to dem, lines...
    inputsDF: dataFrame
        path to simulation data to analyze
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation for the reference and the current simulation
    timeMass: 1D numpy array
        time array for mass analysis (if flagMass=True, otherwise None)
    simName: str
        name of the curent simulation to analyze
    resAnalysisDF: dataFrame
        results from Aimec Analysis

    Returns
    -------
    resAnalysisDF: dataFrame
        results from Aimec Analysis updated with results from curent simulation:
            -maxpprCrossMax: float
                    max max peak pressure
            -pprCrossMax: 1D numpy array
                    max peak pressure in each cross section
            -pprCrossMean: 1D numpy array
                    mean peak pressure in each cross section
            -maxpfdCrossMax: float
                    max max peak flow depth
            -pfdCrossMax: 1D numpy array
                    max peak flow depth in each cross section
            -pfdCrossMean: 1D numpy array
                    mean peak flow depth in each cross section
            -maxpfvCrossMax: float
                    max max peak flow velocity
            -pfvCrossMax: 1D numpy array
                    max peak flow velocity in each cross section
            -pfvCrossMean: 1D numpy array
                    mean peak flow velocity in each cross section
            -xRunout: float
                    x coord of the runout point calculated from the
                    MAX peak result in each cross section (resType provided in the ini file)
            -yRunout: float
                    y coord of the runout point calculated from the
                    MAX peak result in each cross section (resType provided in the ini file)
            -sRunout: float
                    projected runout distance calculated from the
                    MAX peak result in each cross section (resType provided in the ini file)
            -xMeanRunout: float
                    x coord of the runout point calculated from the
                    MEAN peak result in each cross section (resType provided in the ini file)
            -yMeanRunout: float
                    y coord of the runout point calculated from the
                    MEAN peak result in each cross section (resType provided in the ini file)
            -sMeanRunout: float
                    projected runout distance calculated from the
                    MEAN peak result in each cross section (resType provided in the ini file)
            -elevRel: float
                    elevation of the release area (based on first point with
                    peak field > thresholdValue)
            -deltaH: float
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
            -relMass: float
                    release mass
            -entMass: float
                    entrained mass
            -finalMass: float
                    final mass
            -relativMassDiff: float
                    the final mass diff with ref (in %)
            -growthIndex: float
                    growth index
            -growthGrad: float
                    growth gradient
    """
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    interpMethod = cfgSetup['interpMethod']
    flagMass = cfgFlags.getboolean('flagMass')
    refSimName = pathDict['refSimulation']
    # apply domain transformation

    refSimulationName = pathDict['refSimulation']
    log.info('Analyzing data in path coordinate system')
    log.debug("Assigning pressure data to deskewed raster")
    inputFiles = inputsDFrow['ppr']
    print(inputFiles)
    newRasterPPR = aT.transform(inputFiles, rasterTransfo, interpMethod)
    newRasters['newRasterPPR'] = newRasterPPR
    log.debug("Assigning thickness data to deskewed raster")
    inputFiles = inputsDFrow['pfd']
    newRasterPFD = aT.transform(inputFiles, rasterTransfo, interpMethod)
    newRasters['newRefRasterPFD'] = newRasterPFD
    log.debug("Assigning velocity data to deskewed raster")
    inputFiles = inputsDFrow['pfv']
    newRasterPFV = aT.transform(inputFiles, rasterTransfo, interpMethod)
    newRasters['newRasterPFV'] = newRasterPFV

    if simName == refSimName:
        newRasters['newRefRasterPPR'] = newRasterPPR
        newRasters['newRasterPFD'] = newRasterPFD
        newRasters['newRefRasterPFV'] = newRasterPFV
        dataTypeList = ['ppr', 'pfd', 'pfv']
        for dataType in dataTypeList:
            resAnalysisDF[dataType + 'CrossMax'] = np.nan
            resAnalysisDF[dataType + 'CrossMax'] = resAnalysisDF[dataType + 'CrossMax'].astype(object)
            resAnalysisDF[dataType + 'CrossMean'] = np.nan
            resAnalysisDF[dataType + 'CrossMean'] = resAnalysisDF[dataType + 'CrossMax'].astype(object)

    # analyize all fields
    dataPressure = newRasters['newRasterPPR']
    dataThickness = newRasters['newRasterPFD']
    dataSpeed = newRasters['newRasterPFV']
    resAnalysisDF = aT.analyzeField(simName, rasterTransfo, dataPressure, 'ppr', resAnalysisDF)
    resAnalysisDF = aT.analyzeField(simName, rasterTransfo, dataThickness, 'pfd', resAnalysisDF)
    resAnalysisDF = aT.analyzeField(simName, rasterTransfo, dataSpeed, 'pfv', resAnalysisDF)

    # compute runout based on resType
    resAnalysisDF = aT.computeRunOut(cfgSetup, rasterTransfo, resAnalysisDF, newRasters['newRasterDEM'], simName)

    if flagMass:
        # perform mass analysis
        fnameMass = inputsDFrow['massBal']
        if simName == refSimName:
            timeMass = 0
        resAnalysisDF, timeMass = aT.analyzeMass(fnameMass, simName, refSimulationName, resAnalysisDF, time=timeMass)
        outAimec.visuMass(resAnalysisDF, pathDict, simName, refSimulationName, timeMass)
    else:
        timeMass = None

    resAnalysisDF, compPlotPath = aT.analyzeArea(rasterTransfo, resAnalysisDF, simName, newRasters, cfgSetup, pathDict)

    resAnalysisDF.loc[simName, 'areasPlot'] = compPlotPath
    return resAnalysisDF, newRasters, timeMass


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
    outAimec.visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict)
    #################################################################

    # transform resType_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning data to deskewed raster")
    newRasters['newRaster' + resType.upper()] = aT.assignData(pathDict[resType], rasterTransfo, interpMethod)
    # assign dem data
    log.debug("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([pathDict['demSource']], rasterTransfo, interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.debug('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMECIndi(rasterTransfo, newRasters, cfgSetup, pathDict)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')

    outAimec.visuRunoutStat(rasterTransfo, resAnalysis, newRasters, cfgSetup, pathDict)
    outAimec.resultVisu(cfgSetup, pathDict, cfgFlags, rasterTransfo, resAnalysis)

    return rasterTransfo, newRasters, resAnalysis


def postProcessAIMECIndi(rasterTransfo, newRasters, cfgSetup, pathDict):
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
        path to dem, lines...

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
    runout, runoutMean, elevRel, deltaH = aT.computeRunOut(rasterTransfo, thresholdValue, resultsAreaAnalysis,
                                                           transformedDEMRasters)

    runoutLength = runout[0]
    TP, FN, FP, TN, compPlotPath = aT.analyzeArea(rasterTransfo, runoutLength, dataResType, cfgSetup, pathDict)

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

    reportD['Aimec analysis'] = {'type': 'list', 'runout [m]': resAnalysis['runout'][0][compSim],
                                 'max peak pressure [kPa]': resAnalysis['MMPPR'][compSim],
                                 'max peak flow depth [m]': resAnalysis['MMPFD'][compSim],
                                 'max peak flow velocity [ms-1]': resAnalysis['MMPFV'][compSim]}

    benchD['Aimec analysis'] = {'type': 'list', 'runout [m]': resAnalysis['runout'][0][refSim],
                                'max peak pressure [kPa]': resAnalysis['MMPPR'][refSim],
                                'max peak flow depth [m]': resAnalysis['MMPFD'][refSim],
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
