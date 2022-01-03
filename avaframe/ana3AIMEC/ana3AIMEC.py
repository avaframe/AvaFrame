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


def AIMEC2Report(pathDict, inputsDF, cfg):
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
    timeMass = None
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
    plotDict = {}
    log.info('Visualisation of AIMEC results')
    outAimec.visuSimple(cfgSetup, rasterTransfo, resAnalysisDF, newRasters, pathDict)
    plotName = outAimec.visuRunoutComp(rasterTransfo, resAnalysisDF, cfgSetup, pathDict)
    plotDict['slCompPlot'] = {'Aimec comparison of mean and max values along path': plotName}
    plotDict['areasPlot'] = {'Aimec area analysis': resAnalysisDF['areasPlot'][1]}
    if cfgFlags.getboolean('flagMass'):
        plotDict['massAnalysisPlot'] = {'Aimec mass analysis': resAnalysisDF['massPlotName'][1]}

    return rasterTransfo, resAnalysisDF, plotDict


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
        dictionary with paths to dem and lines for Aimec analysis
    inputsDF : dataFrame
        dataframe with simulations to analyze and associated path to raster data
    cfg : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini

    Returns
    -------
    rasterTransfo: dict
        domain transformation information
    resAnalysisDF: dataFrame
        results of ana3AIMEC analysis
    """

    # Extract input config parameters
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    interpMethod = cfgSetup['interpMethod']

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
    timeMass = None
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

    return rasterTransfo, resAnalysisDF


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
    pathDict : dict
        dictionary with paths to dem and lines for Aimec analysis
    inputsDF : dataFrame
        dataframe with simulations to analyze and associated path to raster data
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
    refSimulationName = pathDict['refSimulation']
    # apply domain transformation

    log.info('Analyzing data in path coordinate system')
    log.debug("Assigning pressure data to deskewed raster")
    inputFiles = inputsDFrow['ppr']
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

    if simName == refSimulationName:
        newRasters['newRefRasterPPR'] = newRasterPPR
        newRasters['newRasterPFD'] = newRasterPFD
        newRasters['newRefRasterPFV'] = newRasterPFV
        dataTypeList = ['ppr', 'pfd', 'pfv']
        for dataType in dataTypeList:
            resAnalysisDF[dataType + 'CrossMax'] = np.nan
            resAnalysisDF[dataType + 'CrossMax'] = resAnalysisDF[dataType + 'CrossMax'].astype(object)
            resAnalysisDF[dataType + 'CrossMean'] = np.nan
            resAnalysisDF[dataType + 'CrossMean'] = resAnalysisDF[dataType + 'CrossMax'].astype(object)

    # analyze all fields
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
        resAnalysisDF, timeMass = aT.analyzeMass(fnameMass, simName, refSimulationName, resAnalysisDF, time=timeMass)
        if simName != refSimulationName:
            massPlotName = outAimec.visuMass(resAnalysisDF, pathDict, simName, refSimulationName, timeMass)
            resAnalysisDF.loc[simName, 'massPlotName'] = massPlotName
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


def aimecRes2ReportDict(resAnalysisDF, reportD, benchD, refSimName):
    """ gather aimec results and append them to report dictionary """
    resAnalysisDFRef = resAnalysisDF[resAnalysisDF['simName']==refSimName]
    resAnalysisDFComp = resAnalysisDF[resAnalysisDF['simName']!=refSimName]
    reportD['Aimec analysis'] = {'type': 'list', 'runout [m]': resAnalysisDFComp['sRunout'][0],
                                 'max peak pressure [kPa]': resAnalysisDFComp['maxpprCrossMax'][0],
                                 'max peak flow depth [m]': resAnalysisDFComp['maxpfdCrossMax'][0],
                                 'max peak flow velocity [ms-1]': resAnalysisDFComp['maxpfvCrossMax'][0]}

    benchD['Aimec analysis'] = {'type': 'list', 'runout [m]': resAnalysisDFRef['sRunout'][0],
                                'max peak pressure [kPa]': resAnalysisDFRef['maxpprCrossMax'][0],
                                'max peak flow depth [m]': resAnalysisDFRef['maxpfdCrossMax'][0],
                                'max peak flow velocity [ms-1]': resAnalysisDFRef['maxpfvCrossMax'][0]}
    # add mass info
    if "relMass" in resAnalysisDF.columns:
        reportD['Aimec analysis'].update({'Initial mass [kg]': resAnalysisDFComp['relMass'][0]})
        reportD['Aimec analysis'].update({'Final mass [kg]': resAnalysisDFComp['finalMass'][0]})
        reportD['Aimec analysis'].update({'Entrained mass [kg]': resAnalysisDFComp['entMass'][0]})
        benchD['Aimec analysis'].update({'Initial mass [kg]': resAnalysisDFRef['relMass'][0]})
        benchD['Aimec analysis'].update({'Final mass [kg]': resAnalysisDFRef['finalMass'][0]})
        benchD['Aimec analysis'].update({'Entrained mass [kg]': resAnalysisDFRef['entMass'][0]})

    return reportD, benchD
