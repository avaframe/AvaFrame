"""
Run script for module com4FlowPy
"""
# Load modules
import pathlib
import os
from datetime import datetime
import logging

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
    """ this is a wrapper around com4FlowPy.py that handles the following tasks:
            * reading inputs from (local_)avaframeCfg.ini and (local_)com4FlowPyCfg.ini
            * constructing cfgPath and cfgSetup dictionaries for passing to com4FlowPy.com4FlowPyMain()
            * handling creation of working directories ...

            NOTE-TODO:
                * Forest Interaction currently only works with customDirs set to True
                  Currently this is not implemented in the with the AvaFrame File structure ...
                * This function needs clean-up!
    """
    # log file name; leave empty to use default runLog.log
    logName = 'runcom4FlowPy'
    # Read main Config and com4FlowPy config from files
    cfgMain = cfgUtils.getGeneralConfig()
    cfg     = cfgUtils.getModuleConfig(com4FlowPy)

    cfgSetup = cfg['SETUP']
    cfgFlags = cfg['FLAGS']
    cfgCustomPaths = cfg['PATHS']

    # if customPaths == False --> use AvaFrame Folder structure
    if cfgCustomPaths["useCustomPaths"]=='False':
        # Load avalanche directory from general configuration file
        avalancheDir = cfgMain['MAIN']['avalancheDir']
        # Clean input directory of old work and output files from module
        initProj.cleanModuleFiles(avalancheDir, com4FlowPy, deleteOutput=False)
        # Extract input file locations
        cfgPath = readFlowPyinputs(avalancheDir, cfg)
        # Start logging
        
        log = logUtils.initiateLogger(avalancheDir, logName)
        log.info('MAIN SCRIPT')
        log.info('Current avalanche: %s', avalancheDir)
        
        # IMPORTANT!! - this is a quick'n'dirty hack to set nSims to 9999, so that
        # min(nCPU,nSims) in cfgUtils.getNumberOfProcesses returns nCPU
        cfgSetup['cpuCount']  = str(cfgUtils.getNumberOfProcesses(cfgMain,9999))
        cfgPath["customDirs"] = cfgCustomPaths["useCustomPaths"]

        # Create result directory
        # NOTE-TODO: maybe move into separate function as well ...
        timeString = datetime.now().strftime("%Y%m%d_%H%M%S")
        cfgPath["resDir"] = cfgPath["outDir"] / "res_{}".format(timeString)
        fU.makeADir(cfgPath["resDir"])
        cfgPath["tempDir"] = cfgPath["workDir"] / "temp"
        fU.makeADir(cfgPath["tempDir"])

        cfgPath["deleteTemp"]='False'

        com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)

    # if customPaths == True --> check
    elif cfgCustomPaths["useCustomPaths"]=='True':
        
        cfgPath = {}

        # Handling Custom directory creation
        workDir = pathlib.Path(cfgCustomPaths['workDir'])
        
        time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        try:
            os.makedirs(workDir / 'res_{}'.format(time_string))
            res_dir = (workDir / 'res_{}'.format(time_string))
        except FileExistsError:
            print('folder with same name already exists - aborting')
            sys.exit(1)
        try:
            os.makedirs(workDir / res_dir / 'temp')
            temp_dir = (workDir / res_dir / 'temp')
        except FileExistsError:
            print('temp folder already exists - aborting')
            sys.exit(1)

        cfgPath["workDir"]=pathlib.Path(workDir)
        cfgPath["outDir"]=pathlib.Path(res_dir)
        cfgPath["resDir"]=cfgPath["outDir"]
        cfgPath["tempDir"]=pathlib.Path(temp_dir)
        cfgPath["demPath"]=pathlib.Path(cfgCustomPaths["demPath"])
        cfgPath["releasePath"]=pathlib.Path(cfgCustomPaths["releasePath"])
        cfgPath["infraPath"]=pathlib.Path(cfgCustomPaths["infraPath"])
        cfgPath["forestPath"]=pathlib.Path(cfgCustomPaths["forestPath"])
        cfgPath["deleteTemp"]=cfgCustomPaths["deleteTempFolder"]

        log = logUtils.initiateLogger(cfgPath["outDir"], logName)

        cfgSetup["cpuCount"]  = str(cfgUtils.getNumberOfProcesses(cfgMain,9999))
        cfgPath["customDirs"] = cfgCustomPaths["useCustomPaths"]

        com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)

    else:
        print('INPUT SETTINGS incorrect - please check (local_)avaframeCfg.ini and (local_)com4FlowPyCfg.ini')
        sys.exit(1)

def readFlowPyinputs(avalancheDir, cfgFlowPy):
    """ function is used to read necessary flowPy Inputs 
        function is only called from '../runCom4FlowPy.py', 
        which in turn creates the necessary
        'cfgPath' dictionary, that is passed to
        the com4FlowPyMain() function in this file...
        TODO-20240418: check this function - this is still from MTo Flow-Py --> avaFrame port
    """
    log = logging.getLogger(__name__)
    cfgPath = {}
    avalancheDir = pathlib.Path(avalancheDir)
    # read release area
    releaseDir = avalancheDir / "Inputs" / "REL"
    # from shapefile
    relFiles = sorted(list(releaseDir.glob("*.shp")))
    tryTif = False
    if len(relFiles) == 0:
        log.info("Found not *.shp file containing the release area in %s, trying with *.tif" % releaseDir)
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
    if len(infraPath) == 0 or cfgFlowPy.getboolean("SETUP", "infra") is False:
        infraPath = ""
    elif len(infraPath) > 1:
        message = "More than one Infrastructure file .%s file in %s not allowed" % (infraDir)
        log.error(message)
        raise AssertionError(message)
    else:
        infraPath = infraPath[0]
        log.info("Infrastructure area file is: %s" % infraPath)
    cfgPath["infraPath"] = infraPath

    # read DEM
    demPath = gI.getDEMPath(avalancheDir)
    log.info("DEM file is: %s" % demPath)
    cfgPath["demPath"] = demPath

    # make output path
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, "com4FlowPy", False)
    cfgPath["outDir"] = outDir
    cfgPath["workDir"] = workDir

    return cfgPath

if __name__ == '__main__':
    main()

