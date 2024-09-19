"""
Run script for module com4FlowPy
"""

# Load modules
import pathlib
import os
import sys
from datetime import datetime
import logging
import json

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.com4FlowPy import com4FlowPy
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initialiseDirs as inDirs

# Local imports (avaFrame API)
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.initialiseDirs as inDirs
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as gT


def main():
    """this is a wrapper around com4FlowPy.py that handles the following tasks:
    * reading inputs from (local_)avaframeCfg.ini and (local_)com4FlowPyCfg.ini
    * constructing cfgPath and cfgSetup dictionaries for passing to com4FlowPy.com4FlowPyMain()
    * handling creation of working directories ...

    NOTE-TODO:
        * This function needs clean-up!
    """
    # log file name; leave empty to use default runLog.log
    logName = "runcom4FlowPy"
    # Read main Config and com4FlowPy config from files
    cfgMain = cfgUtils.getGeneralConfig()
    cfg = cfgUtils.getModuleConfig(com4FlowPy)

    # check and handle outputFiles list provided in (local_)com4FlowPyCfg.ini
    cfg["PATHS"]["outputFiles"] = checkOutputFilesFormat(cfg["PATHS"]["outputFiles"])

    cfgSetup = cfg["GENERAL"]
    cfgFlags = cfg["FLAGS"]
    cfgCustomPaths = cfg["PATHS"]

    # if customPaths == False --> use AvaFrame Folder structure
    if cfgCustomPaths["useCustomPaths"] == "False":
        # if "useCustomPaths" == False, we also use the AvaDir Info for the
        # creation of the simulaiton uid
        cfg['GENERAL']['avaDir'] = cfgMain["MAIN"]["avalancheDir"]
        uid = cfgUtils.cfgHash(cfg)
        # Load avalanche directory from general configuration file
        avalancheDir = cfgMain["MAIN"]["avalancheDir"]
        # Clean input directory of old work and output files from module
        initProj.cleanModuleFiles(avalancheDir, com4FlowPy, deleteOutput=False)

        # Start logging
        log = logUtils.initiateLogger(avalancheDir, logName)
        log.info("==================================")
        log.info("MAIN SCRIPT")
        log.info("Current avalanche: %s", avalancheDir)
        log.info("==================================")

        # Extract input file locations
        cfgPath = readFlowPyinputs(avalancheDir, cfg, log)

        # IMPORTANT!! - this is a quick'n'dirty hack to set nSims to 9999, so that
        # min(nCPU,nSims) in cfgUtils.getNumberOfProcesses returns nCPU
        cfgSetup["cpuCount"] = str(cfgUtils.getNumberOfProcesses(cfgMain, 9999))
        cfgPath["customDirs"] = cfgCustomPaths["useCustomPaths"]

        # Create result directory
        # NOTE-TODO: maybe move into separate function as well ...
        timeString = datetime.now().strftime("%Y%m%d_%H%M%S")

        cfgPath["resDir"] = cfgPath["outDir"] / "peakFiles" / "res_{}".format(uid)  # (timeString)
        # check if simulation with same uid already has results folder
        if os.path.isdir(cfgPath["resDir"]):
            log.info("folder with same name already exists - aborting")
            log.info("simulation results folder with same .ini parameters already exists: simulation {}".format(uid))
            sys.exit(1)
        else:
            fU.makeADir(cfgPath["resDir"])
            cfgPath["tempDir"] = cfgPath["workDir"] / "temp"
            fU.makeADir(cfgPath["tempDir"])

        # writing config to .json file
        successToJSON = writeCfgJSON(cfg, uid, cfgPath['outDir'])

        if successToJSON is True:
            log.info('wrote config to {}/{}.json'.format(cfgPath['outDir'], uid))
        else:
            log.info('could not write  config to {}/{}.json'.format(cfgPath['outDir'], uid))
            log.error("Exception occurred: %s", str(successToJSON), exc_info=True)

        cfgPath["deleteTemp"] = "False"

        cfgPath["uid"] = uid
        cfgPath["timeString"] = timeString
        cfgPath["outputFiles"] = cfgCustomPaths["outputFiles"]

        com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)

    # if customPaths == True --> check
    elif cfgCustomPaths["useCustomPaths"] == "True":
        # if "useCustomPaths" == True, we don't need the AvaDir Info for the
        # creation of the simulaiton uid
        uid = cfgUtils.cfgHash(cfg)
        cfgPath = {}

        # Handling Custom directory creation
        workDir = pathlib.Path(cfgCustomPaths["workDir"])

        timeString = datetime.now().strftime("%Y%m%d_%H%M%S")
        try:
            os.makedirs(workDir / "res_{}".format(uid))  # (time_string))
            res_dir = workDir / "res_{}".format(uid)   # (time_string)
        except FileExistsError:
            log.info("folder with same name already exists - aborting")
            log.info("simulation results folder with same .ini parameters already exists: simulation {}".format(uid))
            sys.exit(1)
        try:
            os.makedirs(workDir / res_dir / "temp")
            temp_dir = workDir / res_dir / "temp"
        except FileExistsError:
            log.info("temp folder for simualtion {} already exists - aborting".format(uid))
            sys.exit(1)
        log = logUtils.initiateLogger(res_dir, logName)

        # writing config to .json file
        successToJSON = writeCfgJSON(cfg, uid, workDir)

        if successToJSON is True:
            log.info('wrote config to {}/{}.json'.format(workDir, uid))
        else:
            log.info('could not write  config to {}/{}.json'.format(workDir, uid))
            log.error("Exception occurred: %s", str(successToJSON), exc_info=True)

        cfgPath["workDir"] = pathlib.Path(workDir)
        cfgPath["outDir"] = pathlib.Path(res_dir)
        cfgPath["resDir"] = cfgPath["outDir"]
        cfgPath["tempDir"] = pathlib.Path(temp_dir)
        cfgPath["demPath"] = pathlib.Path(cfgCustomPaths["demPath"])
        cfgPath["releasePath"] = pathlib.Path(cfgCustomPaths["releasePath"])
        cfgPath["infraPath"] = pathlib.Path(cfgCustomPaths["infraPath"])
        cfgPath["forestPath"] = pathlib.Path(cfgCustomPaths["forestPath"])
        cfgPath["varUmaxPath"] = pathlib.Path(cfgCustomPaths["varUmaxPath"])
        cfgPath["varAlphaPath"] = pathlib.Path(cfgCustomPaths["varAlphaPath"])
        cfgPath["varExponentPath"] = pathlib.Path(cfgCustomPaths["varExponentPath"])
        cfgPath["deleteTemp"] = cfgCustomPaths["deleteTempFolder"]
        cfgPath["outputFileFormat"] = cfgCustomPaths["outputFileFormat"]
        cfgPath["outputFiles"] = cfgCustomPaths["outputFiles"]

        cfgSetup["cpuCount"] = str(cfgUtils.getNumberOfProcesses(cfgMain, 9999))
        cfgPath["customDirs"] = cfgCustomPaths["useCustomPaths"]

        cfgPath["uid"] = uid
        cfgPath["timeString"] = timeString

        com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)

    else:
        print(
            "INPUT SETTINGS incorrect - please check (local_)avaframeCfg.ini and (local_)com4FlowPyCfg.ini"
        )
        sys.exit(1)


def readFlowPyinputs(avalancheDir, cfgFlowPy, log):
    """function is used to read necessary flowPy Inputs
    function is only called from '../runCom4FlowPy.py',
    which in turn creates the necessary
    'cfgPath' dictionary, that is passed to
    the com4FlowPyMain() function in this file...
    TODO-20240418: check this function - this is still from MTo Flow-Py --> avaFrame port
    """

    cfgPath = {}
    avalancheDir = pathlib.Path(avalancheDir)
    # read release area
    releaseDir = avalancheDir / "Inputs" / "REL"

    # from shapefile
    relFiles = sorted(list(releaseDir.glob("*.shp")))
    tryTif = False
    if len(relFiles) == 0:
        log.info("Found no *.shp file containing the release area in %s, trying with *.tif" % releaseDir)
        tryTif = True
    elif len(relFiles) > 1:
        message = "There should be exactly one *.shp file containing the release area in %s" % releaseDir
        log.error(message)
        raise AssertionError(message)
    else:
        log.info("Release area file is: %s" % relFiles[0])
        cfgPath["releasePath"] = relFiles[0]

    # from tif
    if tryTif:
        relFiles = sorted(list(releaseDir.glob("*.tif")))
        if len(relFiles) == 0:
            message = (
                "You need to provide one *.shp file or one *.tif file containing the release area in %s"
                % releaseDir
            )
            log.error(message)
            raise FileNotFoundError(message)
        elif len(relFiles) > 1:
            message = "There should be exactly one *.tif file containing the release area in %s" % releaseDir
            log.error(message)
            raise AssertionError(message)
        else:
            log.info("Release area file is: %s" % relFiles[0])
            cfgPath["releasePath"] = relFiles[0]

    # read infra area
    infraDir = avalancheDir / "Inputs" / "INFRA"
    infraPath = sorted(list(infraDir.glob("*.tif")))
    if len(infraPath) == 0 or cfgFlowPy.getboolean("GENERAL", "infra") is False:
        infraPath = ""
    elif len(infraPath) > 1:
        message = "More than one Infrastructure file .%s file in %s not allowed" % (infraDir)
        log.error(message)
        raise AssertionError(message)
    else:
        infraPath = infraPath[0]
        log.info("Infrastructure area file is: %s" % infraPath)
    cfgPath["infraPath"] = infraPath

    # read uMax Limit Raster
    varUmaxDir = avalancheDir / "Inputs" / "UMAX"
    varUmaxPath = sorted(list(varUmaxDir.glob("*.tif")))
    if len(varUmaxPath) == 0 or cfgFlowPy.getboolean("GENERAL", "variableUmaxLim") is False:
        varUmaxPath = ""
    elif len(varUmaxPath) > 1:
        message = "More than one uMax Limit file .%s file in %s not allowed" % (varUmaxDir)
        log.error(message)
        raise AssertionError(message)
    else:
        varUmaxPath = varUmaxPath[0]
        log.info("uMax Limit file is: %s" % varUmaxPath)
    cfgPath["varUmaxPath"] = varUmaxPath

    # read variable Alpha Angle Raster
    varAlphaDir = avalancheDir / "Inputs" / "ALPHA"
    varAlphaPath = sorted(list(varAlphaDir.glob("*.tif")))
    if len(varAlphaPath) == 0 or cfgFlowPy.getboolean("GENERAL", "variableAlpha") is False:
        varAlphaPath = ""
    elif len(varAlphaPath) > 1:
        message = "More than one variable alpha file .%s file in %s not allowed" % (varAlphaDir)
        log.error(message)
        raise AssertionError(message)
    else:
        varAlphaPath = varAlphaPath[0]
        log.info("variable Alpha file is: %s" % varAlphaPath)
    cfgPath["varAlphaPath"] = varAlphaPath

    # read variable Exponent Raster
    varExponentDir = avalancheDir / "Inputs" / "EXP"
    varExponentPath = sorted(list(varExponentDir.glob("*.tif")))
    if len(varExponentPath) == 0 or cfgFlowPy.getboolean("GENERAL", "variableExponent") is False:
        varExponentPath = ""
    elif len(varExponentPath) > 1:
        message = "More than one variable exponent file .%s file in %s not allowed" % (varExponentDir)
        log.error(message)
        raise AssertionError(message)
    else:
        varExponentPath = varExponentPath[0]
        log.info("variable Exponent file is: %s" % varExponentPath)
    cfgPath["varExponentPath"] = varExponentPath

    # check if forest should be used (assumed to be in the RES - 'RESISTANCE' directory)

    forestDir = avalancheDir / "Inputs" / "RES"  # use directory for Resistance Layer for the forest layer

    if cfgFlowPy.getboolean("GENERAL", "forest") is False:
        forestPath = ""
    else:
        patterns = ("*.tif", "*.asc", "*.TIF", "*.tiff", "*.TIFF", "*.ASC")
        forestPath = [f for f in forestDir.iterdir() if any(f.match(p) for p in patterns)]

        if len(forestPath) == 0:
            message = (
                "Please provide a Forest file in %s or set 'forest = False' in the .ini file" % forestDir
            )
            log.error(message)
            raise AssertionError(message)
        elif len(forestPath) > 1:
            message = (
                "Please provide exactly one Forest file in %s or set 'forest = False' in the .ini file"
                % forestDir
            )
            log.error(message)
            raise AssertionError(message)
        else:
            forestPath = forestPath[0]
            log.info("Forest File file is: %s" % forestPath)

    cfgPath["forestPath"] = forestPath

    # read DEM
    demDir = avalancheDir / "Inputs"
    patterns = ("*.tif", "*.asc", "*.TIF", "*.tiff", "*.TIFF", "*.ASC")

    _demPath = [f for f in demDir.iterdir() if any(f.match(p) for p in patterns)]

    if len(_demPath) == 0:
        message = (
                "Please provide a DEM file in %s" % demDir
                  )
        log.error(message)
        raise AssertionError(message)
    elif len(_demPath) > 1:
        message = (
                "Please provide exactly 1 (One!) DEM file in %s" % demDir
                  )
        log.error(message)
        raise AssertionError(message)
    else:
        if os.path.splitext(_demPath[0])[1] in [".asc", ".ASC"]:
            demPath = gI.getDEMPath(avalancheDir)
        else:
            demPath = _demPath[0]

    log.info("DEM file is: %s" % demPath)
    cfgPath["demPath"] = demPath

    # make output path
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, "com4FlowPy", False)
    cfgPath["outDir"] = outDir
    cfgPath["workDir"] = workDir

    cfgPath["outputFileFormat"] = cfgFlowPy["PATHS"]["outputFileFormat"]

    return cfgPath


def checkOutputFilesFormat(strOutputFiles):
    """check if outputFiles option is provided in proper format, else return
    default string in right format

    Parameters:
    ---------------
    strOutputFiles: string - outputFiles string from (local_)com4FlowPy.ini

    Returns:
    ---------------
    strOutputFiles: string - returns the input if the format is ok, else default
                    value string is returned
    """

    try:
        setA = set(strOutputFiles.split('|'))
        setB = set(['zDelta', 'cellCounts', 'fpTravelAngle', 'travelLength',
                    'slTravelAngle', 'flux', 'zDeltaSum'])
        # if there is at least 1 correct ouputfile defined, we use the string provided in the .ini file
        if (setA & setB):
            return strOutputFiles
        else:
            raise ValueError('outputFiles defined in .ini have wrong format - using default settings')
    except ValueError:
        # else we return the default options
        return 'zDelta|cellCounts|travelLength|fpTravelAngle'


def writeCfgJSON(cfg, uid, workDir):
    """
    writes a JSON file containing all the input parameters from the configFile
    using the same uid as the simulation results

    Parameters:
    ---------------
    cfg: configParser object - all the model configs are in here
    uid: string - UID created based on the cfg object
    workDir: string - workDirectory (place to write the .json file to)

    Returns:
    ---------------
    success: boolean/Exception - True if file is written successfully, else Exception

    """

    cfgDict = cfgUtils.convertConfigParserToDict(cfg)

    try:
        with open(workDir / "{}.json".format(uid), 'w') as outfile:
            jsonDict = json.dumps(cfgDict, sort_keys=True, ensure_ascii=True)
            outfile.write(jsonDict)

        return True

    except Exception as e:
        return e


if __name__ == "__main__":
    main()
