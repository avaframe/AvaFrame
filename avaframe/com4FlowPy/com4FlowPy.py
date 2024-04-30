#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
com4FlowPy main function
mainly handling input of data, model params
and output of model results
"""
# Load modules
import pathlib
import numpy as np
from   datetime import datetime
import logging
import pickle
import shutil
import os, sys

# Local imports (avaFrame API)
from   avaframe.in1Data import getInput as gI
import avaframe.in3Utils.initialiseDirs as inDirs
from   avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as gT

# com4FlowPy Libraries
import avaframe.com4FlowPy.rasterIo as io
import avaframe.com4FlowPy.flowCore as fc
import avaframe.com4FlowPy.splitAndMerge as SPAM

# create local logger
log = logging.getLogger(__name__)

def com4FlowPyMain(cfgPath, cfgSetup):
    """ com4FlowPy main function - handles:
            * reading of input data and model Parameters
            * calculation
            * writing of result files

        the function assumes that all necessary directories inside cfgPath have
        already been created (workDir, tempDir)
        
        Parameters
        ----------
        cfgPath:  dictionary with paths from .ini file
        cfgSetup: dictionary with other model configs from .ini file

        Returns
        ----------
        nothing - performs model run and writes results to disk
    """
    _startTime = datetime.now().replace(microsecond=0) #used for timing model runtime

    modelParameters = {}
    # Flow-Py parameters
    modelParameters["alpha"]          = float(cfgSetup["alpha"])
    modelParameters["exp"]            = float(cfgSetup["exp"])
    modelParameters["flux_threshold"] = float(cfgSetup["flux_threshold"])
    modelParameters["max_z"]          = float(cfgSetup["max_z"])

    # Flags for use of Forest and/or Infrastructure
    modelParameters["infra"]  = cfgSetup["infra"]
    modelParameters["forest"] = cfgSetup["forest"]
    
    # Tiling Parameters used for calculation of large model-domains
    tilingParameters = {}
    tilingParameters["tileSize"]    = float(cfgSetup["tileSize"])
    tilingParameters["tileOverlap"] = float(cfgSetup["tileOverlap"])

    # Paths
    modelPaths = {}
    modelPaths["outDir"]        = cfgPath["outDir"]
    modelPaths["workDir"]       = cfgPath["workDir"]
    modelPaths["demPath"]       = cfgPath["demPath"]
    modelPaths["releasePath"]   = cfgPath["releasePath"]

    modelPaths["resDir"]        = cfgPath["resDir"]
    modelPaths["tempDir"]       = cfgPath["tempDir"]

    # check if 'customDirs' are used - alternative is 'default' AvaFrame Folder Structure
    modelPaths["useCustomDirs"]    = True if cfgPath["customDirs"]=='True' else False
    # check if the temp folder, where intermediate results are stored, should be deleted after writing output files
    modelPaths["deleteTempFolder"] = True if cfgPath["deleteTemp"]=="True" else False
    
    # Multiprocessing Options
    MPOptions = {}
    MPOptions["nCPU"]       = int(cfgSetup["cpuCount"]) #number of CPUs to use
    MPOptions["procPerCPU"] = int(cfgSetup["procPerCPUCore"]) #processes per core
    MPOptions["chunkSize"]  = int(cfgSetup["chunkSize"]) # default task size for MP
    MPOptions["maxChunks"]  = int(cfgSetup["maxChunks"]) # max number of tasks for MP
    
    # check if calculation with infrastructure
    if modelParameters["infra"]=='True':
        modelParameters["infraBool"] = True
        modelPaths["infraPath"]      = cfgPath["infraPath"]
    else:
        modelParameters["infraBool"] = False
        modelPaths["infraPath"]      = ""

    forestParams = {}
    # check if calculation with forest
    if modelParameters["forest"]=='True':
        modelParameters["forestBool"] = True

        forestParams["forestModule"] = cfgSetup["forestModule"]
        modelPaths["forestPath"]        = cfgPath["forestPath"]
        # 'forestFriction' and 'forestDetrainment' parameters
        forestParams["maxAddedFriction"] = float(cfgSetup["maxAddedFrictionFor"])
        forestParams["minAddedFriction"] = float(cfgSetup["minAddedFrictionFor"])
        forestParams["velThForFriction"] = float(cfgSetup["velThForFriction"])
        forestParams["maxDetrainment"]   = float(cfgSetup["maxDetrainmentFor"])
        forestParams["minDetrainment"]   = float(cfgSetup["minDetrainmentFor"])
        forestParams["velThForDetrain"]  = float(cfgSetup["velThForDetrain"])
        # 'forestFrictionLayer' parameter
        forestParams["fFrLayerType"]     = cfgSetup["forestFrictionLayerType"]

    else:
        modelParameters["forestBool"] = False
        modelPaths["forestPath"]      = ""

    # TODO: provide some kind of check for the model Parameters
    #       i.e. * sensible value ranges
    #            * contradicting options ...

    # write model parameters paths, etc. to logfile

    startLogging(modelParameters, forestParams, modelPaths, MPOptions)

    # check if release file is given als .shp and convert to .tif/.asc in that case
    # NOTE-TODO: maybe also handle this in ../runCom4FlowPy.py
    modelPaths = checkConvertReleaseShp2Tif(modelPaths)

    # check if input layers have same x,y dimensions
    checkInputLayerDimensions(modelParameters, modelPaths)

    # get information on cellsize and nodata value from demHeader
    demHeader = IOf.readASCheader(modelPaths["demPath"])
    rasterAttributes = {}
    rasterAttributes["cellsize"] = demHeader["cellsize"]
    rasterAttributes["nodata"]   = demHeader["nodata_value"]

    # tile input layers and write tiles (pickled np.arrays) to temp Folder
    nTiles = tileInputLayers(modelParameters, modelPaths, rasterAttributes, tilingParameters)

    # now run the model for all tiles and save the results for each tile to the temp Folder
    performModelCalculation(nTiles, modelParameters, modelPaths, rasterAttributes, forestParams, MPOptions)

    # merge results for the tiles stored in Temp Folder and write Output files
    mergeAndWriteResults(modelPaths, modelParameters)

    _endTime = datetime.now().replace(microsecond=0)
    log.info("==================================")
    log.info(":::> Total time needed: " + str(_endTime - _startTime) + " <:::")
    log.info("==================================")
    
    if (modelPaths["useCustomDirs"]==True) and (modelPaths["deleteTempFolder"] == True):
        deleteTempFolder(modelPaths["tempDir"])

def startLogging(modelParameters, forestParams, modelPaths, MPOptions):
    """ just used to move this chunk of code out of the main function
        only performs logging at the start of the simulation
    """
    # Start of Calculation (logging...)
    log.info("Starting...")
    log.info("========================")
    log.info(f"{'Alpha Angle:' : <20}{modelParameters['alpha'] : <5}")
    log.info(f"{'Exponent:' : <20}{modelParameters['exp'] : <5}")
    log.info(f"{'Flux Threshold:' : <20}{modelParameters['flux_threshold'] : <5}")
    log.info(f"{'Max Z_delta:' : <20}{modelParameters['max_z'] : <5}")
    log.info("------------------------")
    # Also log the used input-files
    log.info(f"{'DEM:' : <5}{'%s'%modelPaths['demPath'] : <5}")
    log.info(f"{'REL:' : <5}{'%s'%modelPaths['releasePath'] : <5}")
    #log.info("DEM: {}".format(modelPaths["demPath"]))
    #log.info("REL: {}".format(modelPaths["releasePath"]))
    log.info("------------------------")
    if modelParameters["forestBool"]:
        log.info("Calculation using forestModule: {}".format(forestParams["forestModule"]))
        log.info(f"{'FOREST LAYER:' : <14}{'%s'%modelPaths['forestPath'] : <5}")
        log.info("-----")
        for param, value in forestParams.items():
            log.info(f"{'%s:'%param : <20}{value : <5}")
        log.info("------------------------")
    if modelParameters["infraBool"]:
        log.info("calculation with Infrastructure")
        log.info(f"{'INFRA LAYER:' : <14}{'%s'%modelPaths['infraPath'] : <5}")
        log.info("------------------------")
    for param, value in MPOptions.items():
        log.info(f"{'%s:'%param : <20}{value : <5}")
        #log.info("{}:\t{}".format(param,value))
    log.info("------------------------")    
    log.info(f"{'WorkDir:' : <12}{'%s'%modelPaths['workDir'] : <5}")
    log.info(f"{'ResultsDir:' : <12}{'%s'%modelPaths['resDir'] : <5}")
    #log.info("WorkDir: {}".format(modelPaths["workDir"]))
    #log.info("ResultsDir: {}".format(modelPaths["resDir"]))
    log.info("========================")

def checkInputLayerDimensions(modelParameters, modelPaths):
    """
    check if all layers have the same size 
    and can be read from the provided paths
    """
    # Check if Layers have same size!!!
    try:
        log.info("checking input layer alignment ...")

        _demHeader = IOf.readASCheader(modelPaths["demPath"])
        _relHeader = io.read_header(modelPaths["releasePath"])

        if _demHeader["ncols"] == _relHeader["ncols"] and _demHeader["nrows"] == _relHeader["nrows"]:
            log.info("DEM and Release Layer ok!")
        else:
            log.error("Error: Release Layer doesn't match DEM!")
            sys.exit(1)

        if modelParameters["infraBool"]:
            _infraHeader = io.read_header(modelPaths["infraPath"])
            if _demHeader["ncols"] == _infraHeader["ncols"] and _demHeader["nrows"] == _infraHeader["nrows"]:
                log.info("Infra Layer ok!")
            else:
                log.error("Error: Infra Layer doesn't match DEM!")
                sys.exit(1)

        if modelParameters["forestBool"]:
            _forestHeader = io.read_header(modelPaths["forestPath"])
            if _demHeader["ncols"] == _forestHeader["ncols"] and _demHeader["nrows"] == _forestHeader["nrows"]:
                log.info("Forest Layer ok!")
            else:
                log.error("Error: Infra Layer doesn't match DEM!")
                sys.exit(1)

        log.info("========================")
    
    except:
        log.error("could not read all required Input Layers, please re-check files and paths provided in .ini files")
        #return
        sys.exit(1)

def tileInputLayers(modelParameters,modelPaths,rasterAttributes,tilingParameters):
    
    _tileCOLS = int(tilingParameters["tileSize"] / rasterAttributes["cellsize"])
    _tileROWS = int(tilingParameters["tileSize"] / rasterAttributes["cellsize"])
    _U = int(tilingParameters["tileOverlap"] / rasterAttributes["cellsize"])  

    log.info("Start Tiling...")
    log.info("---------------------")

    SPAM.tileRaster(modelPaths["demPath"], "dem", modelPaths["tempDir"], _tileCOLS, _tileROWS, _U)
    SPAM.tileRaster(modelPaths["releasePathWork"], "init", modelPaths["tempDir"], _tileCOLS, _tileROWS, _U, isInit=True)
    
    if modelParameters["infraBool"]:
        SPAM.tileRaster(modelPaths["infraPath"], "infra", modelPaths["tempDir"], _tileCOLS, _tileROWS, _U)
    if modelParameters["forestBool"]:
        SPAM.tileRaster(modelPaths["forestPath"], "forest", modelPaths["tempDir"], _tileCOLS, _tileROWS, _U)
    log.info("Finished Tiling All Input Rasters.")
    log.info("==================================")

    nTiles = pickle.load(open(modelPaths["tempDir"] / "nTiles", "rb"))

    return nTiles

def performModelCalculation(nTiles, modelParameters, modelPaths, rasterAttributes, forestParams, MPOptions):
    """ wrapper around fc.run()
        handles passing of model paths, configurations to fc.run()
        also responsible for processing input-data tiles in sequence
    """
    
    optList = []

    for i in range(nTiles[0] + 1):
        for j in range(nTiles[1] + 1):
            optList.append(
                (i,j,
                 modelParameters,
                 modelPaths,
                 rasterAttributes,
                 forestParams,
                 MPOptions
                )
            )

    log.info(' >> Start Calculation << ')
    log.info("-------------------------")
    # Calculation, i.e. iterating over the list of Tiles which have to be processed with fc.run()
    for i,optTuple in enumerate(optList):
        log.info("processing tile %i of %i"%(i+1,len(optList)))
        fc.run(optTuple)
        log.info("finished   tile %i of %i"%(i+1,len(optList)))
        log.info("-------------------------")

    log.info(' >> Calculation Finished << ')
    log.info("==================================")

def mergeAndWriteResults(modelPaths, modelOptions):
    """ function handles merging of results for all tiles inside the temp Folder
        and also writing result files to the resultDir
    """

    log.info(" merging results ...")
    log.info("-------------------------")

    # Merge calculated tiles
    z_delta     = SPAM.MergeRaster(modelPaths["tempDir"], "res_z_delta")
    flux        = SPAM.MergeRaster(modelPaths["tempDir"], "res_flux")
    cell_counts = SPAM.MergeRaster(modelPaths["tempDir"], "res_count")
    z_delta_sum = SPAM.MergeRaster(modelPaths["tempDir"], "res_z_delta_sum")
    fp_ta       = SPAM.MergeRaster(modelPaths["tempDir"], "res_fp")
    sl_ta       = SPAM.MergeRaster(modelPaths["tempDir"], "res_sl")

    if modelOptions["infraBool"]:
        backcalc = SPAM.MergeRaster(modelPaths["tempDir"], "res_backcalc")
    
    # Write Output Files to Disk
    log.info("-------------------------")
    log.info(" writing output files ...")
    log.info("-------------------------")
    output_format = ".tif"
    io.output_raster(modelPaths["demPath"], modelPaths["resDir"] / ("flux%s" % (output_format)), flux)
    io.output_raster(modelPaths["demPath"], modelPaths["resDir"] / ("z_delta%s" % (output_format)), z_delta)
    io.output_raster(modelPaths["demPath"], modelPaths["resDir"] / ("FP_travel_angle%s" % (output_format)), fp_ta)
    io.output_raster(modelPaths["demPath"], modelPaths["resDir"] / ("SL_travel_angle%s" % (output_format)), sl_ta)
    
    #TODO: List of result files, which are produced should be specified also in the .ini file!!!!
    #NOTE: Probably good to have "default" output files (z_delta,FP_travel_angle,cell_counts) 
    #      and only write other output files if set accordingly
    if not modelOptions["infraBool"]:  # if no infra
        io.output_raster(modelPaths["demPath"], modelPaths["resDir"] / ("cell_counts%s" % (output_format)), cell_counts)
        io.output_raster(modelPaths["demPath"], modelPaths["resDir"] / ("z_delta_sum%s" % (output_format)), z_delta_sum)
    if modelOptions["infraBool"]:  # if infra
        io.output_raster(modelPaths["demPath"], modelPaths["resDir"] / ("backcalculation%s" % (output_format)), backcalc)
    
def checkConvertReleaseShp2Tif(modelPaths):
    """ function checks if release area is a .shp file and tries to convert to tif in that case
        
        Parameters:
        ---------------
        modelPaths: {} - dict with modelPaths
        
        Returns:
        ---------------
        modelPaths: {} - dict with added ["releasePathWork"]

    """
    # the release is a shp polygon, we need to convert it to a raster
    # releaseLine = shpConv.readLine(releasePath, 'releasePolygon', demDict)
    if modelPaths["releasePath"].suffix == ".shp":
        
        dem = IOf.readRaster(demPath)
        demHeader = IOf.readASCheader(demPath)

        releaseLine = shpConv.SHP2Array(modelPaths["releasePath"], "releasePolygon")
        thresholdPointInPoly = 0.01
        releaseLine = gT.prepareArea(
            releaseLine, dem, thresholdPointInPoly, combine=True, checkOverlap=False
        )
        # give the same header as the dem
        releaseAreaHeader = demHeader
        releaseArea = np.flipud(releaseLine["rasterData"])
        modelPaths["releasePathWork"] = workDir / "release.tif"
        io.output_raster(modelPaths["demPath"], modelPaths["workDir"] / "release.asc", releaseArea)
        io.output_raster(modelPaths["demPath"], modelPaths["workDir"] / "release.tif", releaseArea)
        del releaseLine
    else:
        modelPaths["releasePathWork"] = modelPaths["releasePath"]
    
    return modelPaths

def deleteTempFolder(tempFolderPath):
    """ delete tempFolder containing the pickled np.arrays of the input data and output data tiles.
            - should be called after all merged model results are written to disk.
        performs a few checks to make sure the folder is indeed a com4FlowPy tempFolder, i.e.
            - does not contain subfolders
            - no other file-extensions than '.npy' and ''
    """
    
    log.info('+++++++++++++++++++++++')
    log.info("deleteTempFolder = True in (local_)com4FlowPyCfg.ini")

    log.info("... checking if folder is a com4FlowPy temp Folder")
    # check if path exists and is directory
    isDir = os.path.isdir(tempFolderPath)
    validTemp = True
    
    for f in os.scandir(tempFolderPath):
        # check if there's a nested folder inside tempDir
        if f.is_dir():
            validTemp=False
            break
        # check if all files are either .npy or start with "ext, "nTi"
        elif f.is_file():
            if not f.path.endswith(".npy"):
                if not f.name[:3] in ["ext","nTi"]:
                    validTemp=False
                    break
    
    if isDir and validTemp:
        log.info("Tempfolder checked: isDir:{} isTemp:{}".format(isDir,validTemp))
        try:
            shutil.rmtree(tempFolderPath)
            log.info("Deleted temp folder {}".format(tempFolderPath))
        except OSError as e:
            log.info("deletion of temp folder {} failed".format(tempFolderPath))
            print ("Error: %s : %s" %(tempFolderPath, e.strerror))
    else:
        log.info("deletion of temp folder {} failed".format(tempFolderPath))
        log.info(" isDir:{} isTemp:{}}".format(isDir,validTemp))
    
    log.info('+++++++++++++++++++++++')

