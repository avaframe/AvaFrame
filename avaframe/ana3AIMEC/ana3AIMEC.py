"""
    AIMEC post processing workflow
"""

import logging
import numpy as np
import pathlib

# Local imports
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.outAIMEC as outAimec

# create local logger
log = logging.getLogger(__name__)


def mainAIMEC(pathDict, inputsDF, cfg):
    """ Main logic for AIMEC postprocessing

    Reads the required files location for ana3AIMEC postprocessing
    given an avalanche directory
    Make the domain transformation corresponding to the input avalanche path
    Transform 2D fields (pressure, flow thickness ...)
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
    plotDict: dict
        dictionary for report
    """

    # Extract input config parameters
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    interpMethod = cfgSetup['interpMethod']

    log.info('Prepare data for post-processing')
    # Read input dem
    demSource = pathDict['demSource']
    dem = IOf.readRaster(demSource)
    # get hash of the reference
    refSimRowHash = pathDict['refSimRowHash']
    # read reference file and raster and config
    refResultSource = inputsDF.loc[refSimRowHash, cfgSetup['runoutResType']]
    refRaster = IOf.readRaster(refResultSource)
    refRasterData = refRaster['rasterData']
    refHeader = refRaster['header']

    log.debug('Data-file %s analysed and domain transformation done' % demSource)

    # Make domain transformation
    log.debug("Creating new deskewed raster and preparing new raster assignment function")
    log.debug('Analyze and make domain transformation on Data-file %s' % demSource)
    rasterTransfo = aT.makeDomainTransfo(pathDict, dem, refHeader['cellsize'], cfgSetup)

    # ####################################################
    # visualisation
    # TODO: needs to be moved somewhere else
    slRaster = aT.transform(refRaster, refResultSource, rasterTransfo, interpMethod)
    newRasters = {}
    log.debug("Assigning dem data to deskewed raster")
    newRasters['newRasterDEM'] = aT.transform(dem, demSource, rasterTransfo, interpMethod)

    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = refRasterData
    inputData['xyHeader'] = refHeader
    outAimec.visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict)

    # postprocess reference...
    timeMass = None
    resAnalysisDF, newRasters, timeMass = postProcessAIMEC(cfg, rasterTransfo, pathDict, inputsDF, newRasters,
                                                           timeMass, refSimRowHash)
    # postprocess other simulations
    for simRowHash, inputsDFrow in inputsDF.iterrows():
        if simRowHash != refSimRowHash:
            resAnalysisDF, newRasters, timeMass = postProcessAIMEC(cfg, rasterTransfo, pathDict, resAnalysisDF,
                                                                   newRasters, timeMass, simRowHash)
            pathDict['simRowHash'] = simRowHash

    # -----------------------------------------------------------
    # result visualisation + report
    # ToDo: should we move this somewere else, this is now just plotting, it should be accessible from outside
    # -----------------------------------------------------------
    plotDict = {}
    log.info('Visualisation of AIMEC results')
    if sorted(pathDict['resTypeList']) == sorted(['ppr', 'pft', 'pfv']):
        outAimec.visuSimple(cfgSetup, rasterTransfo, resAnalysisDF, newRasters, pathDict)
    if len(resAnalysisDF.index) == 2:
        if sorted(pathDict['resTypeList']) == sorted(['ppr', 'pft', 'pfv']):
            plotName = outAimec.visuRunoutComp(rasterTransfo, resAnalysisDF, cfgSetup, pathDict)
    else:
        plotName = outAimec.visuRunoutStat(rasterTransfo, inputsDF, resAnalysisDF, newRasters, cfgSetup, pathDict)
    areaDict = {('area comparison plot ' + resAnalysisDF.loc[k, 'simName']):
                v for k, v in resAnalysisDF['areasPlot'].to_dict().items()}
    plotDict['areasPlot'] = {'Aimec area analysis': areaDict}
    if cfgFlags.getboolean('flagMass'):
        massDict = {('mass comparison plot ' + resAnalysisDF.loc[k, 'simName']):
                    v for k, v in resAnalysisDF['massPlotName'].to_dict().items()}
        plotDict['massAnalysisPlot'] = {'Aimec mass analysis': massDict}

    outAimec.resultVisu(cfgSetup, inputsDF, pathDict, cfgFlags, rasterTransfo, resAnalysisDF)

    plotDict['slCompPlot'] = {'Aimec comparison of mean and max values along path': {'AIMEC '
                                                                                     + cfgSetup['runoutResType']
                                                                                     + ' analysis': plotName}}

    log.info('Writing results to file')
    outAimec.resultWrite(pathDict, cfg, rasterTransfo, resAnalysisDF)

    return rasterTransfo, resAnalysisDF, plotDict


def postProcessAIMEC(cfg, rasterTransfo, pathDict, resAnalysisDF, newRasters, timeMass, simRowHash):
    """ Apply domain transformation and analyse result data (for example pressure, thickness, velocity...)

    Apply the domain tranformation to peak results
    Analyse them.
    Calculate runout, Max Peak Values, Average Peak Values...
    Get mass and entrainement

    Parameters
    ----------
    cfg: configParser objec
        parameters for AIMEC analysis
    rasterTransfo: dict
        transformation information
    pathDict : dict
        dictionary with paths to dem and lines for Aimec analysis
    resAnalysisDF : dataFrame
        dataframe with simulations to analyze and associated path to raster data
    newRasters: dict
        dictionary containing pressure, velocity and flow thickness rasters after
        transformation for the reference and the current simulation
    timeMass: 1D numpy array
        time array for mass analysis (if flagMass=True, otherwise None)
    simRowHash: str
        dataframe hash of the current simulation to analyze

    Returns
    -------
    resAnalysisDF: dataFrame
        input DF with results from Aimec Analysis updated with results from curent simulation:
            -maxACrossMax: float
                    max max A (A depends on what is in resTypes)
            -ACrossMax: 1D numpy array
                    max A in each cross section (A depends on what is in resTypes)
            -ACrossMean: 1D numpy array
                    mean A in each cross section (A depends on what is in resTypes)
            -xRunout: float
                    x coord of the runout point calculated from the
                    MAX peak result in each cross section (runoutResType provided in the ini file)
            -yRunout: float
                    y coord of the runout point calculated from the
                    MAX peak result in each cross section (runoutResType provided in the ini file)
            -sRunout: float
                    projected runout distance calculated from the
                    MAX peak result in each cross section (runoutResType provided in the ini file)
            -xMeanRunout: float
                    x coord of the runout point calculated from the
                    MEAN peak result in each cross section (runoutResType provided in the ini file)
            -yMeanRunout: float
                    y coord of the runout point calculated from the
                    MEAN peak result in each cross section (runoutResType provided in the ini file)
            -sMeanRunout: float
                    projected runout distance calculated from the
                    MEAN peak result in each cross section (runoutResType provided in the ini file)
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
    refSimRowHash = pathDict['refSimRowHash']
    resTypeList = pathDict['resTypeList']

    # apply domain transformation
    log.info('Analyzing data in path coordinate system')

    for resType in resTypeList:
        log.debug("Assigning %s data to deskewed raster" % resType)
        inputFiles = resAnalysisDF.loc[simRowHash, resType]
        if isinstance(inputFiles, pathlib.PurePath):
            rasterData = IOf.readRaster(inputFiles)
            newRaster = aT.transform(rasterData, inputFiles, rasterTransfo, interpMethod)
            newRasters['newRaster' + resType.upper()] = newRaster
            if simRowHash == refSimRowHash:
                newRasters['newRefRaster' + resType.upper()] = newRaster
                resAnalysisDF[resType + 'CrossMax'] = np.nan
                resAnalysisDF[resType + 'CrossMax'] = resAnalysisDF[resType + 'CrossMax'].astype(object)
                resAnalysisDF[resType + 'CrossMean'] = np.nan
                resAnalysisDF[resType + 'CrossMean'] = resAnalysisDF[resType + 'CrossMax'].astype(object)

            # add max, min and std values of result fields
            resAnalysisDF.at[simRowHash, resType + 'FieldMax'] = np.nanmax(rasterData['rasterData'])
            resAnalysisDF.at[simRowHash, resType + 'FieldMean'] = np.nanmean(rasterData['rasterData'])
            resAnalysisDF.at[simRowHash, resType + 'FieldStd'] = np.nanstd(rasterData['rasterData'])

            # analyze all fields
            resAnalysisDF = aT.analyzeField(simRowHash, rasterTransfo, newRaster, resType, resAnalysisDF)

    # compute runout based on runoutResType
    resAnalysisDF = aT.computeRunOut(cfgSetup, rasterTransfo, resAnalysisDF, newRasters, simRowHash)
    if flagMass:
        # perform mass analysis
        fnameMass = resAnalysisDF.loc[simRowHash, 'massBal']
        resAnalysisDF, timeMass = aT.analyzeMass(fnameMass, simRowHash, refSimRowHash, resAnalysisDF, time=timeMass)

        if simRowHash != refSimRowHash:
            massPlotName = outAimec.visuMass(resAnalysisDF, pathDict, simRowHash, refSimRowHash, timeMass)
            resAnalysisDF.loc[simRowHash, 'massPlotName'] = massPlotName
    else:
        timeMass = None

    resAnalysisDF, compPlotPath = aT.analyzeArea(rasterTransfo, resAnalysisDF, simRowHash, newRasters, cfgSetup,
                                                 pathDict)

    resAnalysisDF.loc[simRowHash, 'areasPlot'] = compPlotPath
    return resAnalysisDF, newRasters, timeMass


def aimecRes2ReportDict(resAnalysisDF, reportD, benchD, pathDict):
    """ Gather aimec results and append them to report the dictionary

    Parameters
    ----------
    resAnalysisDF : dataFrame
        dataframe with aimec analysis results
    reportD: dict
        report dictionary to be updated
    benchD: dict
        benchmark dictionary to be updated
    pathDict : dict
        dictionary with refSimRowHash

    Returns
    --------
    reportD: dict
        report dictionary updated with the 'Aimec analysis' section
    benchD: dict
        benchmark dictionary updated with the 'Aimec analysis' section
    """
    refSimRowHash = pathDict['refSimRowHash']
    resAnalysisDFRef = resAnalysisDF.loc[refSimRowHash].squeeze().to_dict()
    resAnalysisDFComp = resAnalysisDF.loc[resAnalysisDF.index != refSimRowHash].squeeze().to_dict()
    reportD['Aimec analysis'] = {'type': 'list', 'runout [m]': resAnalysisDFComp['sRunout'],
                                 'max peak pressure [kPa]': resAnalysisDFComp['maxpprCrossMax'],
                                 'max peak flow thickness [m]': resAnalysisDFComp['maxpftCrossMax'],
                                 'max peak flow velocity [ms-1]': resAnalysisDFComp['maxpfvCrossMax']}

    benchD['Aimec analysis'] = {'type': 'list', 'runout [m]': resAnalysisDFRef['sRunout'],
                                'max peak pressure [kPa]': resAnalysisDFRef['maxpprCrossMax'],
                                'max peak flow thickness [m]': resAnalysisDFRef['maxpftCrossMax'],
                                'max peak flow velocity [ms-1]': resAnalysisDFRef['maxpfvCrossMax']}
    # add mass info
    if "relMass" in resAnalysisDF.columns:
        reportD['Aimec analysis'].update({'Initial mass [kg]': resAnalysisDFComp['relMass']})
        reportD['Aimec analysis'].update({'Final mass [kg]': resAnalysisDFComp['finalMass']})
        reportD['Aimec analysis'].update({'Entrained mass [kg]': resAnalysisDFComp['entMass']})
        benchD['Aimec analysis'].update({'Initial mass [kg]': resAnalysisDFRef['relMass']})
        benchD['Aimec analysis'].update({'Final mass [kg]': resAnalysisDFRef['finalMass']})
        benchD['Aimec analysis'].update({'Entrained mass [kg]': resAnalysisDFRef['entMass']})

    return reportD, benchD
