"""
Fetch input data for avalanche simulations
"""

import logging

# Load modules
import os
import pathlib

import numpy as np
import pandas as pd
import shapely as shp

import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.com1DFA.deriveParameterSet as dP
import avaframe.in2Trans.rasterUtils as IOf
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in3Utils.geoTrans as geoTrans

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.out3Plot import in1DataPlots

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def readDEM(avaDir):
    """read the DEM (ascii or tif) file from a provided avalanche directory

    Parameters
    ----------
    avaDir : str
        path to avalanche directory

    Returns
    -------
    dem : dict
        dict with header and raster data
    """

    # get dem file name
    demSource = getDEMPath(avaDir)

    log.debug("Read DEM: %s" % demSource)

    dem = IOf.readRaster(demSource)

    return dem


def getDEMPath(avaDir):
    """get the DEM file path from a provided avalanche directory

    Parameters
    ----------
    avaDir : str
        path to avalanche directory

    Returns
    -------
    demFile : str (first element of list)
        full path to DEM .asc/.tif file
    """

    # if more than one .asc / .tif file found throw error
    inputDir = pathlib.Path(avaDir, "Inputs")
    demFile = list(inputDir.glob("*.tif")) + list(inputDir.glob("*.asc"))

    if len(demFile) > 1:
        message = "There should be exactly one topography .asc/.tif file in %s/Inputs/" % avaDir
        log.error(message)
        raise AssertionError(message)

    elif len(demFile) == 0:
        message = "No topography .asc / .tif file in %s/Inputs/" % avaDir
        log.error(message)
        raise FileNotFoundError(message)

    return demFile[0]


def getDEMFromConfig(avaDir, fileName=""):
    """get dem file path in avaDir/Inputs

    Parameters
    -----------
    avaDir: str or pathlib path
        to avalancheDir
    fileName: pathlib path
        to dem file with filename in avaDir/Inputs

    Returns
    --------
    demFile: pathlib path
        to dem file
    """

    inputDir = pathlib.Path(avaDir, "Inputs")
    demFile = inputDir / fileName
    if demFile.is_file() is False:
        message = "Dem file: %s does not exist" % (str(demFile))
        log.error(message)
        raise FileNotFoundError(message)

    return demFile


def getInputData(avaDir, cfg):
    """Fetch input datasets required for simulation, duplicated function because
        simulation type set differently in com1DFAOrig compared to com1DFA:
        TODO: remove duplicate once it is not required anymore

    Parameters
    ----------
    avaDir : str
        path to avalanche directory
    cfg : dict
        configuration read from com1DFA simulation ini file

    Returns
    -------
    demFile[0] : str (first element of list)
        list of full path to DEM .asc file
    relFiles : list
        list of full path to release area scenario .shp files
    entFile : str
        full path to entrainment area .shp file
    resFile : str
        full path to resistance area .shp file
    wallFile: str
        full path to wall line .shp file
    entResInfo : flag dict
        flag if Yes entrainment and/or resistance areas found and used for simulation
    """

    # Set directories for inputs, outputs and current work
    avaDir = pathlib.Path(avaDir)
    inputDir = avaDir / "Inputs"

    # Set flag if there is an entrainment or resistance area
    entResInfo = {"flagEnt": "No", "flagRes": "No"}

    # Initialise release areas, default is to look for shapefiles
    releaseDir = "REL"
    if cfg["releaseScenario"] != "":
        relFiles = []
        releaseFiles = cfg["releaseScenario"].split("|")
        for rel in releaseFiles:
            if ".shp" in rel:
                relf = inputDir / releaseDir / rel
            else:
                relf = inputDir / releaseDir / ("%s.shp" % rel)
            if not os.path.isfile(relf):
                message = "No release scenario called: %s" % relf
                log.error(message)
                raise FileNotFoundError(message)
            relFiles.append(relf)
        log.debug("Release area file is specified to be: %s" % relFiles)
    else:
        relDir = inputDir / releaseDir
        relFiles = sorted(relDir.glob("*.shp"))

    log.info("Release area files are: %s" % relFiles)

    # Initialise resistance areas
    resFile, entResInfo["flagRes"] = getAndCheckInputFiles(inputDir, "RES", "Resistance", fileExt="shp")
    if resFile is None:
        resFile = ""
    # Initialise entrainment areas
    entFile, entResInfo["flagEnt"] = getAndCheckInputFiles(inputDir, "ENT", "Entrainment", fileExt="shp")
    if entFile is None:
        entFile = ""

    # Initialise dam line
    wallFile, entResInfo["flagWall"] = getAndCheckInputFiles(inputDir, "DAM", "Dam", fileExt="shp")
    # Initialise DEM
    demFile = getDEMPath(avaDir)

    return demFile, relFiles, entFile, resFile, wallFile, entResInfo


def getInputDataCom1DFA(avaDir):
    """Fetch input datasets required for simulation, duplicated function because
    now fetch all available files simulation type set differently in com1DFA compared
    to com1DFAOrig: TODO: remove duplicate once it is not required anymore

    Parameters
    ----------
    avaDir : str or pathlib object
        path to avalanche directory

    Returns
    -------
    inputSimFiles: dict
        dictionary with all the input files

        - demFile : str (first element of list), list of full path to DEM .asc file
        - relFiles : list, list of full path to release area scenario .shp files
        - secondaryReleaseFile : str, full path to secondary release area .shp file
        - entFile : str, full path to entrainment area .shp file
        - resFile : str, full path to resistance area .shp file
        - entResInfo : flag dict
        flag if Yes entrainment and/or resistance areas found and used for simulation
        flag True if a Secondary Release file found and activated

    """

    # Set directories for inputs, outputs and current work
    inputDir = pathlib.Path(avaDir, "Inputs")

    # Set flag if there is an entrainment or resistance area
    entResInfo = {}

    releaseDir = inputDir / "REL"
    relFiles = sorted(list(releaseDir.glob("*.shp")))
    log.info("Release area files are: %s" % [str(relFilestr) for relFilestr in relFiles])

    # check if relThFile is available
    relThFile, entResInfo["releaseThicknessFile"] = getAndCheckInputFiles(
        inputDir, "RELTH", "release thickness data", fileExt="raster"
    )

    # Initialise secondary release areas
    secondaryReleaseFile, entResInfo["flagSecondaryRelease"] = getAndCheckInputFiles(
        inputDir, "SECREL", "Secondary release", fileExt="shp"
    )
    if secondaryReleaseFile:
        log.info("Secondary release file is: %s" % secondaryReleaseFile)

    # Initialise resistance areas
    resFile, entResInfo["flagRes"] = getAndCheckInputFiles(inputDir, "RES", "Resistance", fileExt="shp")
    if resFile:
        log.info("Resistance file is: %s" % resFile)

    # Initialise entrainment areas
    entFile, entResInfo["flagEnt"] = getAndCheckInputFiles(inputDir, "ENT", "Entrainment", fileExt="shp")
    if entFile:
        log.info("Entrainment file is: %s" % entFile)

    # Initialise dam line
    damFile, entResInfo["dam"] = getAndCheckInputFiles(inputDir, "DAM", "Dam", fileExt="shp")
    if damFile:
        log.info("Dam file is: %s" % damFile)

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    # check if frictionParameter file  is available
    muFile, entResInfo["mu"] = getAndCheckInputFiles(
        inputDir, "RASTERS", "mu parameter data", fileExt="raster", fileSuffix="_mu"
    )

    # check if frictionParameter file  is available
    xiFile, entResInfo["xi"] = getAndCheckInputFiles(
        inputDir, "RASTERS", "xi parameter data", fileExt="raster", fileSuffix="_xi"
    )

    # return DEM, first item of release, entrainment and resistance areas
    inputSimFiles = {
        "demFile": demFile,
        "relFiles": relFiles,
        "secondaryReleaseFile": secondaryReleaseFile,
        "entFile": entFile,
        "resFile": resFile,
        "damFile": damFile,
        "entResInfo": entResInfo,
        "relThFile": relThFile,
        "muFile": muFile,
        "xiFile": xiFile,
    }

    return inputSimFiles


def getAndCheckInputFiles(inputDir, folder, inputType, fileExt="shp", fileSuffix=""):
    """Fetch fileExt files and check if they exist and if it is not more than one. If fileExt is set to "raster" both
    .asc and .tif files will be searched.

    Raises error if there is more than one fileExt file.

    If fileExt is empty, only suffix will be checked.

    If a file is found, it is checked to be of an allowed extension; currently .shp .asc and .tif are supported.

    Parameters
    ----------
    inputDir : pathlib object or str
        path to avalanche input directory
    folder : str
        subfolder name where the shape file should be located (SECREL, ENT or RES)
    inputType : str
        type of input (used for the logging messages).
    fileExt: str
        file extension e.g. shp, asc, tif - optional; default is shp
    fileSuffix: str
        file name part before extension

    Returns
    -------
    OutputFile: str
        path to file checked
    available: str
        Yes or No depending on if there is a  file available (if No, OutputFile is None)
    """
    available = "No"

    supportedFileFormats = [".shp", ".asc", ".tif"]

    # Define the directory to search and the extensions
    if fileExt == "":
        extensions = [""]
    elif fileExt.lower() == "raster":
        extensions = ["asc", "tif"]
    else:
        extensions = [fileExt]

    inDir = pathlib.Path(inputDir, folder)

    if fileSuffix == "":
        OutputFile = list([file for ext in extensions for file in inDir.glob(f"*.{ext}")])
    else:
        OutputFile = list([file for ext in extensions for file in inDir.glob(f"*{fileSuffix}.{ext}")])

    # check for number of files
    if len(OutputFile) < 1:
        OutputFile = None
    elif len(OutputFile) > 1:
        message = "More than one %s .%s file in %s/%s/ not allowed" % (
            inputType,
            fileExt,
            inputDir,
            folder,
        )
        log.error(message)
        raise AssertionError(message)
    else:
        available = "Yes"
        OutputFile = OutputFile[0]

        if OutputFile.suffix not in supportedFileFormats:
            message = (
                "Unsupported file format found for OutputFile %s; shp, asc, tif are allowed" % OutputFile
            )
            log.error(message)
            raise AssertionError(message)

    return OutputFile, available


def getThicknessInputSimFiles(inputSimFiles):
    """add thickness of shapefiles to dictionary

    Parameters
    -----------
    inputSimFiles: dict
        dictionary with info on release and entrainment file paths

    Returns
    --------
    inputSimFiles: dict
        updated dictionary with thickness info read from shapefile attributes
        now includes one separate dictionary for each release, entrainment or secondary release
        scenario with a thickness and id value for each feature (given as list)
    """

    # fetch thickness attribute of entrainment area and secondary release
    for thType in ["entFile", "secondaryReleaseFile"]:
        if inputSimFiles[thType] is not None:
            thicknessList, idList, ci95List = shpConv.readThickness(inputSimFiles[thType])
            inputSimFiles[inputSimFiles[thType].stem] = {
                "thickness": thicknessList,
                "id": idList,
                "ci95": ci95List,
            }

    # initialize release scenario list
    releaseScenarioList = []

    # fetch thickness attribute of release areas and add info to input dict
    for releaseA in inputSimFiles["relFiles"]:
        # fetch thickness and id info from input data
        thicknessList, idList, ci95List = shpConv.readThickness(releaseA)
        inputSimFiles[releaseA.stem] = {
            "thickness": thicknessList,
            "id": idList,
            "ci95": ci95List,
        }
        # append release scenario name to list
        releaseScenarioList.append(releaseA.stem)

    # append release scenario names
    inputSimFiles["releaseScenarioList"] = releaseScenarioList

    return inputSimFiles


def updateThicknessCfg(inputSimFiles, cfgInitial):
    """add available release scenarios to ini file and
    set thickness values in ini files

    Parameters
    -----------
    inputSimFiles: dict
        dictionary with info on release and entrainment file paths
    cfgInitial: configParser object
        with the current (and possibly overridden) configuration

    Returns
    --------
    inputSimFiles: dict
        updated dictionary with thickness info read from shapefile attributes
        now includes one separate dictionary for each release, entrainment or secondary release
        scenario with a thickness and id value for each feature (given as list)
    cfgInitial: configparser object
        updated config object with release scenario, thickness info, etc.

    """

    # check if thickness info is required from entrainment and secondary release according to simType
    simTypeList = cfgInitial["GENERAL"]["simTypeList"].split("|")
    thTypeList = []
    if any(simType in ["ent", "entres", "available"] for simType in simTypeList):
        thTypeList.append("entFile")
    if cfgInitial["GENERAL"].getboolean("secRelArea"):
        thTypeList.append("secondaryReleaseFile")

    # initialize release scenario list
    releaseScenarioIni = cfgInitial["INPUT"]["releaseScenario"]
    if releaseScenarioIni == "":
        releaseScenarioList = inputSimFiles["releaseScenarioList"]
    else:
        releaseScenarioList = cfgInitial["INPUT"]["releaseScenario"].split("|")

    # add input data info to cfg object
    # fetch thickness attribute of release areas and add info to input dict
    for releaseA in releaseScenarioList:
        # update configuration with thickness value to be used for simulations
        cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, releaseA, "relTh")
        if cfgInitial["GENERAL"].getboolean("relThFromFile"):
            if inputSimFiles["relThFile"] is None:
                message = "relThFromFile set to True but no relTh file found"
                log.error(message)
                raise FileNotFoundError(message)
            else:
                cfgInitial["INPUT"]["relThFile"] = str(
                    pathlib.Path("RELTH", inputSimFiles["relThFile"].name)
                )

    # add entrainment and secondary release thickness in input data info and in cfg object
    if inputSimFiles["entFile"] != None and "entFile" in thTypeList:
        cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, inputSimFiles["entFile"].stem, "entTh")
        cfgInitial["INPUT"]["entrainmentScenario"] = inputSimFiles["entFile"].stem
    if inputSimFiles["secondaryReleaseFile"] != None and "secondaryReleaseFile" in thTypeList:
        cfgInitial = dP.getThicknessValue(
            cfgInitial,
            inputSimFiles,
            inputSimFiles["secondaryReleaseFile"].stem,
            "secondaryRelTh",
        )
        cfgInitial["INPUT"]["secondaryReleaseScenario"] = inputSimFiles["secondaryReleaseFile"].stem

    # create cfg string from release scenario list and add to cfg object
    releaseScenarioName = cfgUtils.convertToCfgList(releaseScenarioList)
    if cfgInitial["INPUT"]["releaseScenario"] == "":
        cfgInitial["INPUT"]["releaseScenario"] = releaseScenarioName
    else:
        for relIniFileName in cfgInitial["INPUT"]["releaseScenario"].split("|"):
            if relIniFileName not in releaseScenarioList:
                message = "Chosen release scenario: %s not available" % relIniFileName
                log.error(message)
                raise FileNotFoundError(message)
        else:
            log.info("Chosen release scenarios: %s" % cfgInitial["INPUT"]["releaseScenario"])

    return cfgInitial


def initializeRelTh(cfg, dOHeader):
    """check for relThFile and load to dict (if relTh is read from file)

    Parameters
    -----------
    cfg: configparser
        used are avalanchDir and relThFile
    demHeader:
        header of dem to check for matching rows and numbers

    Returns
    --------
    relThFieldData: ndarray or str
        with release thickness field data
        if not relThFromFile an empty string is returned
    relThFile: path
        updated path to relThFile (needed in case of remeshing)
    """

    avaDir = cfg["GENERAL"]["avalancheDir"]
    relThFile = pathlib.Path(avaDir, "Inputs", cfg["INPUT"]["relThFile"])

    if relThFile != None and cfg["GENERAL"].getboolean("relThFromFile"):
        relThField = IOf.readRaster(relThFile)
        relThFieldData = relThField["rasterData"]
        if (
            dOHeader["ncols"] != relThField["header"]["ncols"]
            or dOHeader["nrows"] != relThField["header"]["nrows"]
        ):
            message = (
                "Release thickness field read from %s does not match the number of rows and columns of the dem"
                % relThFile
            )
            log.error(message)
            raise AssertionError(message)
        elif np.isnan(relThFieldData).any() == True:
            message = (
                "Release thickness field contains nans - not allowed no release thickness must be set to 0"
            )
            log.error(message)
            raise AssertionError(message)
    else:
        relThFieldData = ""

    return relThFieldData, relThFile


def initializeDEM(avaDir, demPath=""):
    """check for dem and load to dict

    Parameters
    -----------
    avaDir: str or pathlib path
        path to avalanche directory
    demPath: str or pathlib Path
        path to dem relative to Inputs - optional if not provided read DEM from Inputs

    Returns
    --------
    demOri: dict
        dem dictionary with header and data
    """

    if demPath == "":
        dem = readDEM(avaDir)
    else:
        # build full path and load data to dict
        demFile = pathlib.Path(avaDir, "Inputs", demPath)
        dem = IOf.readRaster(demFile, noDataToNan=True)

    return dem


def selectReleaseFile(inputSimFiles, releaseScenario):
    """select release scenario

    Parameters
    -----------
    inputSimFiles: dict
        dictionary with info on input data
    releaseScenario: str
        name of release scenario


    Returns
    -------
    inputSimFiles: dict
        dictionary with info on input data updated with releaseScenario
    """

    # fetch release file path for scenario
    relFiles = inputSimFiles["relFiles"]
    for relF in relFiles:
        if relF.stem == releaseScenario:
            releaseScenarioPath = relF

    inputSimFiles["releaseScenario"] = releaseScenarioPath

    return inputSimFiles


def fetchReleaseFile(inputSimFiles, releaseScenario, cfgSim, releaseList):
    """select release scenario, update configuration to only include thickness info
    of current scenario and return file path

    Parameters
    -----------
    inputSimFiles: dict
        dictionary with info on input data
    releaseScenario: str
        name of release scenario
    cfgSim: conigparser object
        configuration of simulation
    releaseList: list
        list of available release scenarios

    Returns
    -------
    releaseScenarioPath: pathlib path
        file path to release scenario shp file
    cfgSim: configparser object
        updated cfg object, removed thickness info from not other release scenarios than used
        one and rename thickness values of chosen scenario to relThThickness, relThId, ...
    """

    # fetch release files paths
    relFiles = inputSimFiles["relFiles"]

    foundScenario = False
    for relF in relFiles:
        if relF.stem == releaseScenario:
            releaseScenarioPath = relF
            foundScenario = True

    if foundScenario is False:
        message = "Release area scenario %s not found - check input data" % releaseScenario
        log.error(message)
        raise FileNotFoundError(message)

    # update config entry for release scenario, thickness and id
    cfgSim["INPUT"]["releaseScenario"] = str(releaseScenario)
    if cfgSim["GENERAL"]["relThFromShp"] == "True":
        for scenario in releaseList:
            if scenario == releaseScenario:
                cfgSim["INPUT"]["relThId"] = cfgSim["INPUT"][scenario + "_" + "relThId"]
                cfgSim["INPUT"]["relThThickness"] = cfgSim["INPUT"][scenario + "_" + "relThThickness"]
                cfgSim["INPUT"]["relThCi95"] = cfgSim["INPUT"][scenario + "_" + "relThCi95"]
            # remove thickness, id and ci95 values specified by releaseScenario
            cfgSim["INPUT"].pop(scenario + "_" + "relThId")
            cfgSim["INPUT"].pop(scenario + "_" + "relThThickness")
            cfgSim["INPUT"].pop(scenario + "_" + "relThCi95")

    return releaseScenarioPath, cfgSim


def createReleaseStats(avaDir, cfg):
    """create a csv file with info on release shp file on:
    max, mean and min elevation, slope and projected and real area

    Parameters
    -------------

    Returns
    --------
    fPath: pathlib Path
        file path
    """

    # fetch input data (dem and release area shp file)
    inputData = getInputDataCom1DFA(avaDir)

    # get line from release area polygon
    relFiles = inputData["relFiles"]

    # initialize dem and get real area of dem
    dem = IOf.readRaster(inputData["demFile"], noDataToNan=True)
    dem["originalHeader"] = dem["header"].copy()
    methodMeshNormal = cfg.getfloat("GENERAL", "methodMeshNormal")
    # get normal vector of the grid mesh
    dem = geoTrans.getNormalMesh(dem, methodMeshNormal)
    dem = DFAtls.getAreaMesh(dem, methodMeshNormal)

    # loop over all relFiles and compute information saved to dictionary of dataframes
    relPath = pathlib.Path(avaDir, "Outputs", "com1DFA", "releaseInfoFiles")
    fU.makeADir(relPath)
    relDFDict = {}
    # loop over all relFiles
    for relFile in relFiles:
        # for relF in relFiles:
        releaseLine = shpConv.readLine(relFile, "release1", dem)
        releaseLine["file"] = relFile
        releaseLine["type"] = "Release"
        releaseLine["thicknessSource"] = ["artificial"] * len(releaseLine["id"])
        releaseLine["thickness"] = ["1."] * len(releaseLine["id"])
        # convert release line to a raster with 1 set inside of the release line
        # and compute projected and actual areas
        areaActualList, areaProjectedList, releaseLine = computeAreasFromRasterAndLine(releaseLine, dem)

        # compute max, min elevation, mean slope and save all to dict
        relInfo = computeRelStats(releaseLine, dem)

        # create dict and dataFrame
        relD = {
            "release feature": relInfo["featureNames"],
            "MaxZ [m]": relInfo["zMax"],
            "MinZ [m]": relInfo["zMin"],
            "slope [deg]": relInfo["meanSlope"],
            "projected area [ha]": [areaP / 10000.0 for areaP in areaProjectedList],
            "actual area [ha]": [area / 10000.0 for area in areaActualList],
        }
        relDF = pd.DataFrame(data=relD, index=np.arange(len(relInfo["featureNames"])))
        relInfo = relPath / ("%s.csv" % relFile.stem)
        relDF.to_csv(relInfo, index=False, sep=";", float_format="%.4f")

        log.info("Written relInfo file for %s to %s" % (relFile.stem, str(relPath)))
        relDFDict[relFile.stem] = relDF

    return relDFDict


def computeAreasFromRasterAndLine(line, dem):
    """compute the area covered by a polygon by creating a raster from polygon
    projected area and actual area using a dem info

    Parameters
    -----------
    line: dict
        dictionary with info on line
        x, y coordinates start, end of each line feature
    dem: dict
        dictionary with dem data, header and areaRaster

    Returns
    --------
    areaActual: float
        actual area taking slope into account by using dem area
    areaProjected: float
        projected area in xy plane
    """

    line = geoTrans.prepareArea(line, dem, 0.01, combine=False, checkOverlap=False)

    csz = dem["header"]["cellsize"]
    # create dict for raster data for each feature
    areaProjectedList = []
    areaActualList = []
    for index, lineRaster in enumerate(line["rasterData"]):
        lineRasterOnes = np.where(lineRaster > 0, 1.0, 0.0)
        areaActualList.append(np.nansum(lineRasterOnes * dem["areaRaster"]))
        areaProjectedList.append(np.sum(csz * csz * lineRasterOnes))

    return areaActualList, areaProjectedList, line


def computeRelStats(line, dem):
    """compute stats of a polygon and a dem
    actual area (taking slope into account), projected area,
    max, mean, min elevation, mean slope

    Parameters
    -----------
    line: dict
        dictionary with info on line (x, y, Start, Length, rasterData, ...)
    dem: dict
        dictionary with info on dem (header, rasterData, normals, areaRaster)

    Returns
    --------
    lineDict: dict
        dictionary with stats info:
        featureNames, zMax, zMin, Slope, Area, AreaP
    """

    # compute projected areas using shapely
    projectedAreas = computeAreasFromLines(line)

    # create dict for raster data for each feature
    lineDict = {}
    for index, relRaster in enumerate(line["rasterData"]):
        zArray = np.where(relRaster > 0, dem["rasterData"], np.nan)

        # compute slope of release area
        _, _, NzNormed = DFAtls.normalize(dem["Nx"], dem["Ny"], dem["Nz"])
        nzArray = np.where(relRaster > 0, NzNormed, np.nan)
        slopeArray = np.rad2deg(np.arccos(nzArray))
        lineDict.setdefault("meanSlope", []).append(np.nanmean(slopeArray))

        # compute elevation stats
        lineDict.setdefault("zMax", []).append(np.nanmax(zArray))
        lineDict.setdefault("zMean", []).append(np.nanmean(zArray))
        lineDict.setdefault("zMin", []).append(np.nanmin(zArray))

        # create a dataframe
        lineDict.setdefault("featureNames", []).append(line["Name"][index])

        # print info
        log.debug("++++ Release feature %s ++++++++++" % line["Name"][index])
        log.debug("maximum elevation: %.2f" % np.nanmax(zArray))
        log.debug("mean elevation: %.2f" % np.nanmean(zArray))
        log.debug("minimum elevation: %.2f" % np.nanmin(zArray))
        log.debug("mean slope: %.2f" % np.nanmean(slopeArray))

    return lineDict


def computeAreasFromLines(line):
    """compute the area of a polygon in xy using shapely

    Parameters
    -----------
    line: dict
        dictionary with info on line
        x, y coordinates start, end of each line feature

    Returns
    --------
    projectedAreas: list
        list of projected area for each polygon in line dict
    """

    projectedAreas = []
    name = line["Name"]
    start = line["Start"]
    Length = line["Length"]
    # fetch individual polygons
    for i in range(len(name)):
        end = start[i] + Length[i]
        x = line["x"][int(start[i]) : int(end)]
        y = line["y"][int(start[i]) : int(end)]

        # create shapely polygon
        for m in range(len(x)):
            avaPoly = shp.Polygon(list(zip(x, y)))
        projectedAreas.append(shp.area(avaPoly))

    return projectedAreas


def getInputPaths(avaDir):
    """Fetch paths to dem and first release area shp file found

    Parameters
    ----------
    avaDir : str or pathlib object
        path to avalanche directory


    Returns
    -------
    demFile : pathlib path
        full path to DEM .asc file
    relFiles : list
        list of full paths to release area scenario .shp files found in avaDir/Inputs/REL
    relFieldFiles : list
        list of full paths to release area thickness .asc files found in avaDir/Inputs/RELTH

    """

    # Set directories for inputs, outputs and current work
    inputDir = pathlib.Path(avaDir, "Inputs")

    # fetch release area shp files
    releaseShpDir = inputDir / "REL"
    relFiles = sorted(list(releaseShpDir.glob("*.shp")))
    log.info("Release area files are: %s" % [str(relFilestr) for relFilestr in relFiles])

    # fetch release thickness fields
    releaseFieldDir = inputDir / "RELTH"
    relFieldFiles = sorted(list(releaseFieldDir.glob("*.asc")) + list(releaseFieldDir.glob("*.tif")))
    if len(relFieldFiles) > 0:
        log.info("Release area files are: %s" % [str(relFFilestr) for relFFilestr in relFieldFiles])
    else:
        relFieldFiles = None

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    return demFile, relFiles, relFieldFiles


def checkForMultiplePartsShpArea(avaDir, lineDict, modName, type=""):
    """check if in polygon read from shape file holes are present, if so error and save a plot to Outputs/com1DFA
    procedure: check if polygon has several parts

    Parameters
    -----------
    avaDir: str
        path to avalanche directory
    lineDict: dict
        dictionary with info read from shape file
        used: x, y, Start, Length, nParts, nFeatures
    modName: str
        name of computational module where to save error plots to Outputs subfolder
    type: str
        type of shp file area (release, secondary release, entrainment, resistance)

    Returns
    --------
    error if number of parts is bigger 1
    save plot showing all parts of each polygon feature in Outputs/modName

    """

    foundMultipleParts = False
    # loop over all polygons in scenario
    for lineFeature in range(lineDict["nFeatures"]):
        if len(lineDict["nParts"][lineFeature]) > 2:
            # fetch number of parts for each lineFeature
            nParts = lineDict["nParts"][lineFeature]
            # create x, y coordinates of each lineFeature
            xFeat = lineDict["x"][
                int(lineDict["Start"][lineFeature]) : int(
                    lineDict["Start"][lineFeature] + lineDict["Length"][lineFeature]
                )
            ]
            yFeat = lineDict["y"][
                int(lineDict["Start"][lineFeature]) : int(
                    lineDict["Start"][lineFeature] + lineDict["Length"][lineFeature]
                )
            ]

            # check for type of area
            if type.lower() in ["entrainment", "resistance", "secondary release"]:
                lineFileName = lineDict["fileName"]
            else:
                lineFileName = lineDict["file"]

            # create plot of parts of feature for analysis
            title = "Parts of lineFeature %d of %s" % (lineFeature, lineFileName.stem)
            outDir = pathlib.Path(avaDir, "Outputs", modName)
            pathDict = {
                "pathResult": outDir,
                "title": title,
                "outFileName": ("%s feature%d_errorPlot" % (type, lineFeature)),
            }
            in1DataPlots.plotAreaShpError(xFeat, yFeat, nParts, pathDict)

            # set flag
            foundMultipleParts = True

    # if polygon has multiple parts - error
    if foundMultipleParts:
        message = "One or more %s features in %s have holes - check error plots in %s" % (
            type,
            lineFileName.name,
            outDir,
        )
        log.error(message)
        raise AssertionError(message)
