"""
    Main logic for AIMEC post processing
"""

import os
import logging
import glob
import numpy as np
import copy


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.outAIMEC as outAimec

# create local logger
log = logging.getLogger(__name__)



def AIMEC2Report(cfgPath, cfg):
    """ perform AIMEC analysis and generate plots for reports

    Reads the required files location for ana3AIMEC postpocessing
    given a path dictionary to the input files

    Parameters
    ----------
    cfgPath : dict
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
    pressureLimit = cfgSetup.getfloat('pressureLimit')
    interpMethod = cfgSetup['interpMethod']

    log.info('Prepare data for post-ptocessing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = aT.makeDomainTransfo(cfgPath, cfgSetup)

    ###########################################################################
    # visualisation
    # TODO: needs to be moved somewhere else
    # read reference file
    nRef = cfgPath['referenceFile']
    rasterSource = cfgPath['ppr'][nRef]

    pressureRaster = IOf.readRaster(rasterSource)
    slRaster = aT.transform(rasterSource, rasterTransfo, interpMethod)
    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = pressureRaster['rasterData']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgPath, cfgFlags)
    #################################################################

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning pressure data to deskewed raster")
    newRasters['newRasterPressure'] = aT.assignData(cfgPath['ppr'],
                                                 rasterTransfo, interpMethod)
    # assign depth data
    log.debug("Assigning depth data to deskewed raster")
    newRasters['newRasterDepth'] = aT.assignData(cfgPath['pfd'],
                                              rasterTransfo, interpMethod)
    # assign speed data
    log.debug("Assigning speed data to deskewed raster")
    newRasters['newRasterSpeed'] = aT.assignData(cfgPath['pfv'],
                                              rasterTransfo, interpMethod)

    # assign dem data
    log.debug("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([cfgPath['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.debug('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMECReport(rasterTransfo, pressureLimit, newRasters, cfgPath, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')
    plotName = outAimec.visuRunoutComp(rasterTransfo, resAnalysis, pressureLimit, newRasters, cfgPath, cfgFlags)
    resAnalysis['slCompPlot'] = {'Aimec comparison of mean and max values along path': plotName}

    return rasterTransfo, newRasters, resAnalysis


def mainAIMEC(cfgPath, cfg):
    """ Main logic for AIMEC postprocessing

    Reads the required files location for ana3AIMEC postpocessing
    given an avalanche directory

    Parameters
    ----------
    cfgPath : dict
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
    pressureLimit = cfgSetup.getfloat('pressureLimit')
    interpMethod = cfgSetup['interpMethod']

    log.info('Prepare data for post-ptocessing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = aT.makeDomainTransfo(cfgPath, cfgSetup)

    ###########################################################################
    # visualisation
    # TODO: needs to be moved somewhere else
    # read reference file
    nRef = cfgPath['referenceFile']
    rasterSource = cfgPath['ppr'][nRef]

    pressureRaster = IOf.readRaster(rasterSource)
    slRaster = aT.transform(rasterSource, rasterTransfo, interpMethod)
    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = pressureRaster['rasterData']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgPath, cfgFlags)
    #################################################################

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning pressure data to deskewed raster")
    newRasters['newRasterPressure'] = aT.assignData(cfgPath['ppr'],
                                                 rasterTransfo, interpMethod)
    # assign depth data
    log.debug("Assigning depth data to deskewed raster")
    newRasters['newRasterDepth'] = aT.assignData(cfgPath['pfd'],
                                              rasterTransfo, interpMethod)
    # assign speed data
    if cfgPath['pfv']:
        log.debug("Assigning speed data to deskewed raster")
        newRasters['newRasterSpeed'] = aT.assignData(cfgPath['pfv'],
                                                  rasterTransfo, interpMethod)

    # assign dem data
    log.info("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([cfgPath['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.info('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMEC(rasterTransfo, pressureLimit, newRasters, cfgPath, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')
    outAimec.visuSimple(rasterTransfo, resAnalysis, newRasters, cfgPath, cfgFlags)
    if cfgPath['numSim']==2:
        outAimec.visuRunoutComp(rasterTransfo, resAnalysis, pressureLimit, newRasters, cfgPath, cfgFlags)
        outAimec.visuMass(resAnalysis, cfgPath, cfgFlags)
    else:
        outAimec.visuRunoutStat(rasterTransfo, resAnalysis, pressureLimit, newRasters, cfgPath, cfgFlags)
    outAimec.resultVisu(cfgPath, cfgFlags, rasterTransfo, resAnalysis,
                        pressureLimit)

    # -----------------------------------------------------------
    # write results to file
    # -----------------------------------------------------------
    log.info('Writing results to file')

    outAimec.resultWrite(cfgPath, cfgSetup, rasterTransfo, resAnalysis)

    return rasterTransfo, newRasters, resAnalysis


def AIMECIndi(cfgPath, cfg):
    """ perform AIMEC analysis and generate plots for reports

    Reads the required files location for ana3AIMEC postpocessing
    given a path dictionary to the input files

    Parameters
    ----------
    cfgPath : dict
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
    anaLimit = cfgSetup.getfloat('pressureLimit')
    interpMethod = cfgSetup['interpMethod']
    resType = cfgPath['resType']

    log.info('Prepare data for post-ptocessing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = aT.makeDomainTransfo(cfgPath, cfgSetup)

    ###########################################################################
    # visualisation
    # TODO: needs to be moved somewhere else
    # read reference file
    nRef = cfgPath['referenceFile']
    rasterSource = cfgPath[resType][nRef]

    anaRaster = IOf.readRaster(rasterSource)
    slRaster = aT.transform(rasterSource, rasterTransfo, interpMethod)
    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = anaRaster['rasterData']
    outAimec.visuTransfo(rasterTransfo, inputData, cfgPath, cfgFlags)
    #################################################################

    # transform resType_data in new raster
    newRasters = {}
    # assign pressure data
    log.debug("Assigning data to deskewed raster")
    newRasters['newRasterResType'] = aT.assignData(cfgPath[resType],
                                                 rasterTransfo, interpMethod)
    # assign dem data
    log.debug("Assigning dem data to deskewed raster")
    newRasterDEM = aT.assignData([cfgPath['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.debug('Analyzing data in path coordinate system')
    resAnalysis = postProcessAIMECIndi(rasterTransfo, anaLimit, newRasters, cfgPath, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of AIMEC results')

    #outAimec.visuRunoutStat(rasterTransfo, resAnalysis, anaLimit, newRasters, cfgPath, cfgFlags)
    outAimec.resultVisu(cfgPath, cfgFlags, rasterTransfo, resAnalysis,
                        anaLimit, flagIndi=True)

    return rasterTransfo, newRasters, resAnalysis




# -----------------------------------------------------------
# Aimec analysis tools
# -----------------------------------------------------------
def postProcessAIMEC(rasterTransfo, pLim, newRasters, cfgPath, cfgFlags):
    """ Analyse pressure and depth transformed data

    Analyse pressure depth and speed.
    Calculate runout, Max Peak Pressure, Average PP...
    Get mass and entrainement

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    pLim: float
        numerical value of the pressure limit to use
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgPath: dict
        path to data to analyse

    Returns
    -------
    resAnalysis: dict
        resAnalysis dictionnary containing all results:
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
            -MMPPR: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak pressure
            -MMPFD: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak flow depth
            -MMPFV: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak flow velocity
            -elevRel: 1D numpy array
                    containing for each simulation analyzed the
                    elevation of the release area (based on first point with
                    peak pressure > pLim)
            -deltaH: 1D numpy array
                    containing for each simulation analyzed the
                    elevation fall difference between elevRel and altitude of
                    run-out point
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
            -PPRCrossMax: 2D numpy array
                    containing for each simulation analyzed the
                    max peak pressure in each cross section
            -PPRCrossMean: 2D numpy array
                    containing for each simulation analyzed the
                    mean peak pressure in each cross section
            -PFDCrossMax: 2D numpy array
                    containing for each simulation analyzed the
                    max peak flow depth in each cross section
            -PFDCrossMean: 2D numpy array
                    containing for each simulation analyzed the
                    mean peak flow depth in each cross section
            -PFVCrossMax: 2D numpy array
                    containing for each simulation analyzed the
                    max peak flow velocity in each cross section
            -PFVCrossMean: 2D numpy array
                    containing for each simulation analyzed the
                    mean peak flow velocity in each cross section
            -pressureLimit: float
                    pressure threshold
            -startOfRunoutAreaAngle: float
                    angle of the slope at the beginning of the run-out
                    area (given in input)
            -TP: float
                    ref = True sim2 = True
            -FN: float
                    ref = False sim2 = True
            -FP: float
                    ref = True sim2 = False
            -TN: float
                    ref = False sim2 = False
    """
    # read inputs
    fnameMass = cfgPath['mb']
    dataPressure = newRasters['newRasterPressure']
    dataDepth = newRasters['newRasterDepth']
    dataSpeed = newRasters['newRasterSpeed']
    transformedDEMRasters = newRasters['newRasterDEM']

    maxPPRCrossMax, PPRCrossMax, PPRCrossMean = aT.analyzeField(rasterTransfo, dataPressure, 'ppr')
    maxPFDCrossMax, PFDCrossMax, PFDCrossMean = aT.analyzeField(rasterTransfo, dataDepth, 'pfd')
    maxPFVCrossMax, PFVCrossMax, PFVCrossMean = aT.analyzeField(rasterTransfo, dataSpeed, 'pfv')

    runout, runoutMean, elevRel, deltaH = aT.computeRunOut(rasterTransfo, pLim, PPRCrossMax, PPRCrossMean, transformedDEMRasters)

    releaseMass, entrainedMass, entMassArray, totalMassArray, finalMass, relativMassDiff, grIndex, grGrad, time = aT.analyzeMass(fnameMass)

    runoutLength = runout[0]
    TP, FN, FP, TN, _ = aT.analyzeArea(rasterTransfo, runoutLength, pLim, dataPressure, cfgPath, cfgFlags)

    # affect values to output dictionary
    resAnalysis = {}
    resAnalysis['runout'] = runout
    resAnalysis['runoutMean'] = runoutMean
    resAnalysis['MMPPR'] = maxPPRCrossMax
    resAnalysis['MMPFD'] = maxPFDCrossMax
    resAnalysis['MMPFV'] = maxPFVCrossMax
    resAnalysis['elevRel'] = elevRel
    resAnalysis['deltaH'] = deltaH
    resAnalysis['relMass'] = releaseMass
    resAnalysis['entMass'] = entrainedMass
    resAnalysis['entMassArray'] = entMassArray
    resAnalysis['totalMassArray'] = totalMassArray
    resAnalysis['time'] = time
    resAnalysis['finalMass'] = finalMass
    resAnalysis['relativMassDiff'] = relativMassDiff
    resAnalysis['growthIndex'] = grIndex
    resAnalysis['growthGrad'] = grGrad
    resAnalysis['PPRCrossMax'] = PPRCrossMax
    resAnalysis['PPRCrossMean'] = PPRCrossMean
    resAnalysis['PFDCrossMax'] = PFDCrossMax
    resAnalysis['PFDCrossMean'] = PFDCrossMean
    resAnalysis['PFVCrossMax'] = PFVCrossMax
    resAnalysis['PFVCrossMean'] = PFVCrossMean
    resAnalysis['pressureLimit'] = pLim
    resAnalysis['startOfRunoutAreaAngle'] = rasterTransfo['startOfRunoutAreaAngle']
    resAnalysis['TP'] = TP
    resAnalysis['FN'] = FN
    resAnalysis['FP'] = FP
    resAnalysis['TN'] = TN

    return resAnalysis

def postProcessAIMECReport(rasterTransfo, pLim, newRasters, cfgPath, cfgFlags):
    """ Analyse pressure, depth and speed transformed data

    Analyse pressure depth and speed.
    Calculate runout, Max Peak Pressure, Average PP...
    Mass is currently not included

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    pLim: float
        numerical value of the pressure limit to use
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgPath: dict
        path to data to analyse

    Returns
    -------
    resAnalysis: dict
        resAnalysis dictionary containing all results:

    """
    # read inputs
    dataPressure = newRasters['newRasterPressure']
    dataDepth = newRasters['newRasterDepth']
    dataSpeed = newRasters['newRasterSpeed']
    transformedDEMRasters = newRasters['newRasterDEM']

    # get max and mean values along path for cross profiles
    maxPPRCrossMax, PPRCrossMax, PPRCrossMean = aT.analyzeField(rasterTransfo, dataPressure, 'ppr')
    maxPFDCrossMax, PFDCrossMax, PFDCrossMean = aT.analyzeField(rasterTransfo, dataDepth, 'pfd')
    maxPFVCrossMax, PFVCrossMax, PFVCrossMean = aT.analyzeField(rasterTransfo, dataSpeed, 'pfv')

    # compute runout based on peak pressure
    runout, runoutMean, elevRel, deltaH = aT.computeRunOut(rasterTransfo, pLim, PPRCrossMax, PPRCrossMean, transformedDEMRasters)

    runoutLength = runout[0]
    TP, FN, FP, TN, compPlotPath = aT.analyzeArea(rasterTransfo, runoutLength, pLim, dataPressure, cfgPath, cfgFlags)

    # affect values to output dictionary
    resAnalysis = {}
    resAnalysis['runout'] = runout
    resAnalysis['runoutMean'] = runoutMean
    resAnalysis['MMPPR'] = maxPPRCrossMax
    resAnalysis['MMPFD'] = maxPFDCrossMax
    resAnalysis['MMPFV'] = maxPFVCrossMax
    resAnalysis['elevRel'] = elevRel
    resAnalysis['deltaH'] = deltaH
    resAnalysis['PPRCrossMax'] = PPRCrossMax
    resAnalysis['PPRCrossMean'] = PPRCrossMean
    resAnalysis['PFDCrossMax'] = PFDCrossMax
    resAnalysis['PFDCrossMean'] = PFDCrossMean
    resAnalysis['PFVCrossMax'] = PFVCrossMax
    resAnalysis['PFVCrossMean'] = PFVCrossMean
    resAnalysis['pressureLimit'] = pLim
    resAnalysis['startOfRunoutAreaAngle'] = rasterTransfo['startOfRunoutAreaAngle']
    resAnalysis['TP'] = TP
    resAnalysis['FN'] = FN
    resAnalysis['FP'] = FP
    resAnalysis['TN'] = TN
    resAnalysis['areasPlot'] = {'Aimec area analysis': compPlotPath}

    return resAnalysis


def postProcessAIMECIndi(rasterTransfo, anaLim, newRasters, cfgPath, cfgFlags):
    """ Analyse pressure, depth and speed transformed data

    Analyse pressure depth and speed.
    Calculate runout, Max Peak Pressure, Average PP...
    Mass is currently not included

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    pLim: float
        numerical value of the pressure limit to use
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgPath: dict
        path to data to analyse

    Returns
    -------
    resAnalysis: dict
        resAnalysis dictionary containing all results:

    """
    # read inputs
    resType = cfgPath['resType']
    dataResType = newRasters['newRasterResType']
    transformedDEMRasters = newRasters['newRasterDEM']

    # get max and mean values along path for cross profiles
    maxPCrossMax, PCrossMax, PCrossMean = aT.analyzeField(rasterTransfo, dataResType, resType)

    # compute runout based on peak pressure
    runout, runoutMean, elevRel, deltaH = aT.computeRunOut(rasterTransfo, anaLim, PCrossMax, PCrossMean, transformedDEMRasters)

    runoutLength = runout[0]
    TP, FN, FP, TN, compPlotPath = aT.analyzeArea(rasterTransfo, runoutLength, anaLim, dataResType, cfgPath, cfgFlags)

    # affect values to output dictionary
    resAnalysis = {}
    resAnalysis['runout'] = runout
    resAnalysis['runoutMean'] = runoutMean
    resAnalysis['MMP'] = maxPCrossMax
    resAnalysis['elevRel'] = elevRel
    resAnalysis['deltaH'] = deltaH
    resAnalysis['PCrossMax'] = PCrossMax
    resAnalysis['PCrossMean'] = PCrossMean
    resAnalysis['anaLimit'] = anaLim
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

    return reportD, benchD
