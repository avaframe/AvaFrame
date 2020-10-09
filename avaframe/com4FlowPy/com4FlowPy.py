#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:23:00 2018

@author: Michael Neuhauser
"""

# Load modules
import pathlib
import glob
import numpy as np
from datetime import datetime
from multiprocessing import cpu_count
import multiprocessing as mp
import logging
import pickle


# Local imports
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.initialiseDirs as inDirs
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as gT
# Flow-Py Libraries
from avaframe.com4FlowPy import com4FlowPy
import avaframe.com4FlowPy.raster_io as io
import avaframe.com4FlowPy.flow_core as fc
import avaframe.com4FlowPy.split_and_merge as SPAM

# create local logger
log = logging.getLogger(__name__)


def readFlowPyinputs(avalancheDir, cfgFlowPy):
    cfgPath = {}
    avalancheDir = pathlib.Path(avalancheDir)
    # read release area
    releaseDir = avalancheDir / 'Inputs' / 'REL'
    # from shapefile
    relFiles = sorted(list(releaseDir.glob('*.shp')))
    tryTif = False
    if len(relFiles) == 0:
        log.info('Found not *.shp file containing the release area in %s, trying with *.tif' % releaseDir)
        tryTif = True
    elif len(relFiles) > 1:
        message = 'There should be exactly one *.shp file containing the release area in %s' % releaseDir
        log.error(message)
        raise AssertionError(message)
    else:
        log.info('Release area file is: %s' % relFiles[0])
        cfgPath['releasePath'] = relFiles[0]

    # from tif
    if tryTif:
        relFiles = sorted(list(releaseDir.glob('*.tif')))
        if len(relFiles) == 0:
            message = 'You need to provide one *.shp file or one *.tif file containing the release area in %s' % releaseDir
            log.error(message)
            raise FileNotFoundError(message)
        elif len(relFiles) > 1:
            message = 'There should be exactly one *.tif file containing the release area in %s' % releaseDir
            log.error(message)
            raise AssertionError(message)
        else:
            log.info('Release area file is: %s' % relFiles[0])
            cfgPath['releasePath'] = relFiles[0]

    # read infra area
    infraDir = avalancheDir / 'Inputs' / 'INFRA'
    infraPath = sorted(list(infraDir.glob('*.tif')))
    if len(infraPath) == 0 or cfgFlowPy.getboolean('SETUP', 'infra') is False:
        infraPath = ''
    elif len(infraPath) > 1:
        message = 'More than one Infrastructure file .%s file in %s not allowed' % (infraDir)
        log.error(message)
        raise AssertionError(message)
    else:
        infraPath = infraPath[0]
        log.info('Infrastructure area file is: %s' % infraPath)
    cfgPath['infraPath'] = infraPath

    # read DEM
    demPath = gI.getDEMPath(avalancheDir)
    log.info('DEM file is: %s' % demPath)
    cfgPath['demPath'] = demPath

    # make output path
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, 'com4FlowPy', False)
    cfgPath['outDir'] = outDir
    cfgPath['workDir'] = workDir

    return cfgPath


def com4FlowPyMain(cfgPath, cfgSetup):
    # Input Parameters
    alpha = float(cfgSetup['alpha'])
    exp = float(cfgSetup['exp'])
    flux_threshold = float(cfgSetup['flux_threshold'])
    max_z = float(cfgSetup['max_z'])
    # Input Paths
    outDir = cfgPath['outDir']
    workDir = cfgPath['workDir']
    demPath = cfgPath['demPath']
    releasePath = cfgPath['releasePath']
    infraPath = cfgPath['infraPath']

    print("Starting...")
    print("...")

    start = datetime.now().replace(microsecond=0)
    infraBool = False
    # Create result directory
    timeString = datetime.now().strftime("%Y%m%d_%H%M%S")
    resDir = outDir / 'res_{}'.format(timeString)
    fU.makeADir(resDir)
    tempDir = resDir / 'temp'
    fU.makeADir(tempDir)

    # # Setup logger
    # for handler in logging.root.handlers[:]:
    #     logging.root.removeHandler(handler)
    #
    # logging.basicConfig(level=logging.INFO,
    #                     format='%(asctime)s %(levelname)-8s %(message)s',
    #                     datefmt='%Y-%m-%d %H:%M:%S',
    #                     filename=(resDir / 'log_{}.txt'.format(timeString)),
    #                     filemode='w')

    # Start of Calculation
    log.info('Start Calculation')
    log.info('Alpha Angle: {}'.format(alpha))
    log.info('Exponent: {}'.format(exp))
    log.info('Flux Threshold: {}'.format(flux_threshold))
    log.info('Max Z_delta: {}'.format(max_z))

    # ToDo: this is a kind of inputs check, we should put it somewere else in a sub function
    # Read in raster files
    # ToDo: we actualy only need to read the header
    # demDict = IOf.readRaster(demPath)
    # demHeader = demDict['header']
    demHeader = IOf.readASCheader(demPath)
    # dem, header = io.read_raster(demPath)

    # read the release area
    if releasePath.suffix == '.shp':
        # the release is a shp polygon, we need to convert it to a raster
        # releaseLine = shpConv.readLine(releasePath, 'releasePolygon', demDict)
        releaseLine = shpConv.SHP2Array(releasePath, 'releasePolygon')
        thresholdPointInPoly = 0.01
        releaseLine = gT.prepareArea(releaseLine, demHeader, thresholdPointInPoly, combine=True, checkOverlap=False)
        # give the same header as the dem
        releaseAreaHeader = demHeader
        releaseArea = np.flipud(releaseLine['rasterData'])
        releasePathWork = workDir / "release.tif"
        io.output_raster(demPath, workDir / "release.asc", releaseArea)
        io.output_raster(demPath, workDir / "release.tif", releaseArea)
        del releaseLine
    else:
        releasePathWork = releasePath
        releaseArea, releaseAreaHeader = io.read_raster(releasePath)

    # Check if Layers have same size!!!
    if demHeader['ncols'] == releaseAreaHeader['ncols'] and demHeader['nrows'] == releaseAreaHeader['nrows']:
        print("DEM and Release Layer ok!")
    else:
        print("Error: Release Layer doesn't match DEM!")
        return

    if infraPath != '':
        infraArea, infraAreaHeader = io.read_raster(infraPath)
        if demHeader['ncols'] == infraAreaHeader['ncols'] and demHeader['nrows'] == infraAreaHeader['nrows']:
            print("Infra Layer ok!")
            infraBool = True
            log.info('Infrastructure File: %s' % (infraPath))
        else:
            print("Error: Infra Layer doesn't match DEM!")
            return
    else:
        infraArea = np.zeros((demHeader['nrows'], demHeader['ncols']))

    # ToDo: why? is it in case it was too big?
    del releaseArea, infraArea
    log.info('Files read in')

    cellsize = demHeader["cellsize"]
    nodata = demHeader["noDataValue"]

    if nodata is None:
        print("Your DEM Layer has no defined No Data Value!!!")

    # Here we split the computation domain in sub parts
    tileCOLS = int(15000 / cellsize)
    tileROWS = int(15000 / cellsize)
    U = int(5000 / cellsize)  # 5km overlap

    log.info("Start Tiling.")
    print("Start Tiling...")

    SPAM.tileRaster(demPath, "dem", tempDir, tileCOLS, tileROWS, U)
    SPAM.tileRaster(releasePathWork, "init", tempDir, tileCOLS, tileROWS, U, isInit=True)
    if infraBool:
        SPAM.tileRaster(infraPath, "infra", tempDir, tileCOLS, tileROWS, U)

    print("Finished Tiling...")
    nTiles = pickle.load(open(tempDir / "nTiles", "rb"))

    optList = []
    # das hier ist die batch-liste, die von mulitprocessing
    # abgearbeitet werden muss - sieht so aus:
    # [(0,0,alpha,exp,cellsize,-9999.),
    # (0,1,alpha,exp,cellsize,-9999.),
    # etc.]

    for i in range(nTiles[0]+1):
        for j in range(nTiles[1]+1):
            optList.append((i, j, alpha, exp, cellsize, nodata, flux_threshold, max_z, tempDir))

    # Calculation
    log.info('Multiprocessing starts, used cores: %i' % (cpu_count() - 1))
    print("%i Processes started and %i calculations to perform." % (mp.cpu_count() - 1, len(optList)))
    pool = mp.Pool(mp.cpu_count() - 1)

    if infraBool:
        pool.map(fc.calculation, optList)
    else:
        pool.map(fc.calculation_effect, optList)

    pool.close()
    pool.join()

    log.info('Calculation finished, merging results.')

    # Merge calculated tiles
    z_delta = SPAM.MergeRaster(tempDir, "res_z_delta")
    flux = SPAM.MergeRaster(tempDir, "res_flux")
    cell_counts = SPAM.MergeRaster(tempDir, "res_count")
    z_delta_sum = SPAM.MergeRaster(tempDir, "res_z_delta_sum")
    fp_ta = SPAM.MergeRaster(tempDir, "res_fp")
    sl_ta = SPAM.MergeRaster(tempDir, "res_sl")
    if infraBool:
        backcalc = SPAM.MergeRaster(tempDir, "res_backcalc")

    # time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
    log.info('Writing Output Files')
    output_format = '.tif'
    io.output_raster(demPath,
                     resDir / ("flux%s" % (output_format)),
                     flux)
    io.output_raster(demPath,
                     resDir / ("z_delta%s" % (output_format)),
                     z_delta)
    io.output_raster(demPath,
                     resDir / ("FP_travel_angle%s" % (output_format)),
                     fp_ta)
    io.output_raster(demPath,
                     resDir / ("SL_travel_angle%s" % (output_format)),
                     sl_ta)
    if not infraBool:  # if no infra
        io.output_raster(demPath,
                         resDir / ("cell_counts%s" % (output_format)),
                         cell_counts)
        io.output_raster(demPath,
                         resDir / ("z_delta_sum%s" % (output_format)),
                         z_delta_sum)
    if infraBool:  # if infra
        io.output_raster(demPath,
                         resDir / ("backcalculation%s" % (output_format)),
                         backcalc)

    print("Calculation finished")
    print("...")
    end = datetime.now().replace(microsecond=0)
    log.info('Calculation needed: ' + str(end - start) + ' seconds')
