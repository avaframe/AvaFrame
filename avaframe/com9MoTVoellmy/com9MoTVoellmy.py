import os
import platform
import logging
import numpy as np
import pathlib
import time
import shutil
import sys


if os.name == "nt":
    from multiprocessing.pool import ThreadPool as Pool
elif platform.system() == "Darwin":
    from multiprocessing.pool import ThreadPool as Pool
else:
    from multiprocessing import Pool

import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in2Trans import rasterUtils as rU
from avaframe.com1DFA import particleInitialisation as pI
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils.cfgUtils import cfgToRcf
from avaframe.in3Utils.MoTUtils import rewriteDEMtoZeroValues, runAndCheckMoT, MoTGenerateConfigs, copyMoTFiles


# create a local logger
log = logging.getLogger(__name__)


def com9MoTVoellmyMain(cfgMain, cfgInfo=None):
    """Run MoT-Voellmy simulations using specified configurations.

    This function executes MoT-Voellmy simulations by handling
    preprocessing, parallel execution, and postprocessing steps.

    Parameters
    ----------
    cfgMain : configparser.ConfigParser
        Main configuration settings for the simulation
    cfgInfo : str or pathlib.Path or configparser.ConfigParser, optional
        Additional configuration information, by default None
        Can be:
        - Path to configuration file for overrides
        - ConfigParser object with initial configuration
        - None to use local/default configuration

    Returns
    -------
    None
        Results are written to the output directory structure

    Notes
    -----
    The function performs several key steps:
    - Generates configuration settings for all simulations
    - Preprocesses input data and creates RCF configuration files
    - Executes simulations in parallel using multiple processes
    - Handles DEM preparation by converting NaN values to zeros
    """

    # Get all necessary information from the configuration files
    currentModule = sys.modules[__name__] # As if you would to import com9MoTVoellmy
    simDict, inputSimFiles = MoTGenerateConfigs(cfgMain, cfgInfo, currentModule)

    # convert DEM from nan to 0 values
    # TODO: suggest MoT-PSA to handle nan values
    rewriteDEMtoZeroValues(inputSimFiles["demFile"])

    log.info("The following simulations will be performed")
    for key in simDict:
        log.info("Simulation: %s" % key)

    # Preprocess the simulations, mainly creating the rcf files
    rcfFiles = com9MoTVoellmyPreprocess(simDict, inputSimFiles, cfgMain)

    # And now we run the simulations
    startTime = time.time()

    log.info("--- STARTING (potential) PARALLEL PART ----")

    # Get number of CPU Cores wanted
    nCPU = cfgUtils.getNumberOfProcesses(cfgMain, len(rcfFiles))

    # Create parallel pool and run
    # with multiprocessing.Pool(processes=nCPU) as pool:
    with Pool(processes=nCPU) as pool:
        results = pool.map(com9MoTVoellmyTask, rcfFiles)
        pool.close()
        pool.join()

    timeNeeded = "%.2f" % (time.time() - startTime)
    log.info("Overall (parallel) com9MoTVoellmy computation took: %s s " % timeNeeded)
    log.info("--- ENDING (potential) PARALLEL PART ----")

    # Postprocess the simulations
    com9MoTVoellmyPostprocess(simDict, cfgMain, inputSimFiles)


def com9MoTVoellmyPostprocess(simDict, cfgMain, inputSimFiles):
    """Post-process MoT-Voellmy simulation results.

    This function handles post-processing tasks after MoT-Voellmy simulations complete,
    including copying result files to output directories and generating visualization plots.

    Parameters
    ----------
    simDict : dict
        Dictionary containing simulation configurations, with one entry per simulation
    cfgMain : configparser.ConfigParser
        Main configuration settings for the simulation
    inputSimFiles : dict
        Dictionary containing paths to input files (DEM, release areas, etc.)

    Returns
    -------
    None
        Results are written to the output directory structure

    Notes
    -----
    The function performs several key tasks:
    - Creates output directory structure
    - Copies simulation result files (DataTime.txt, ppr, pfd, pfv files)
    - Renames files according to conventions
    - Generates visualization plots of peak fields
    """
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    # Copy max files to output directory
    outputDirPeakFile = pathlib.Path(avalancheDir) / "Outputs" / "com9MoTVoellmy" / "peakFiles"
    fU.makeADir(outputDirPeakFile)

    for key in simDict:
        workDir = pathlib.Path(avalancheDir) / "Work" / "com9MoTVoellmy" / str(key)

        # Copy ppr files
        copyMoTFiles(workDir, outputDirPeakFile, "p_max", "ppr")

        # Copy pfd files
        copyMoTFiles(workDir, outputDirPeakFile, "h_max", "pfd")

        # Copy pfv files
        copyMoTFiles(workDir, outputDirPeakFile, "s_max", "pfv")

    # create plots and report
    modName = __name__.split(".")[-1]
    reportDir = pathlib.Path(avalancheDir, "Outputs", modName, "reports")
    fU.makeADir(reportDir)

    dem = rU.readRaster(inputSimFiles["demFile"])

    # Generate plots for all peakFiles
    oP.plotAllPeakFields(avalancheDir, cfgMain["FLAGS"], modName, demData=dem)


def com9MoTVoellmyTask(rcfFile):
    """Execute a single MoT-PSA simulation using the provided configuration file.

    Parameters
    ----------
    rcfFile : pathlib.Path
        Path to the RCF configuration file for the simulation

    Returns
    -------
    list
        The command that was executed as a list containing the executable path
        and configuration file path

    Notes
    -----
    Changes to the directory containing this module before executing the simulation.
    Uses runAndCheckMoT to execute and monitor the simulation process.
    """
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    if os.name == "nt":
        exeName = "MoT-Voellmy_win.exe"
    elif platform.system() == "Darwin":
        message = "MoT-Voellmy does not support MacOS at the moment"
        log.error(message)
        raise OSError(message)
    else:
        exeName = "./MoT-Voellmy_linux.exe"

    command = [exeName, rcfFile]
    log.info("Run simulation: %s" % rcfFile)
    runAndCheckMoT(command)
    return command


def com9MoTVoellmyPreprocess(simDict, inputSimFiles, cfgMain):
    """Preprocess data for MoT-PSA simulations.

    This function prepares the input data and configuration files needed to run
    MoT-PSA simulations. It processes release areas, creates required directories,
    and generates configuration files for each simulation.

    Parameters
    ----------
    simDict : dict
        Dictionary containing simulation configurations with one entry per simulation
    inputSimFiles : dict
        Dictionary containing paths to input files (DEM, release areas, etc.)
    cfgMain : configparser.ConfigParser
        Main configuration settings for the simulation

    Returns
    -------
    list
        List of pathlib.Path objects pointing to generated RCF configuration files
        that will be used to run the simulations

    Notes
    -----
    The function performs several key steps:
    - Creates working directories for each simulation
    - Processes release areas and converts them to raster format 
    - Generates configuration files in RCF format
    - Handles both single and multiple release scenarios
    """
    # Load avalanche directory from general configuration file
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    workDir = pathlib.Path(avalancheDir) / "Work" / "com9MoTVoellmy"
    cfgFileDir = pathlib.Path(avalancheDir) / "Outputs" / "com9MoTVoellmy" / "configurationFiles"
    fU.makeADir(cfgFileDir)
    rcfFiles = list()

    for key in simDict:
        # Generate command and run via subprocess.run
        # Configuration that needs adjustment

        # load configuration object for current sim
        cfg = simDict[key]["cfgSim"]

        # convert release shape to raster with values for current sim
        # select release area input data according to a chosen release scenario
        inputSimFiles = gI.selectReleaseFile(inputSimFiles, cfg["INPUT"]["releaseScenario"])
        # create the required input from input files
        demOri, inputSimLines = com1DFA.prepareInputData(inputSimFiles, cfg)

        if cfg["GENERAL"].getboolean("iniStep"):
            # append buffered release Area
            inputSimLines = pI.createReleaseBuffer(cfg, inputSimLines)

        # set thickness values for the release area, entrainment and secondary release areas
        relName, inputSimLines, badName = com1DFA.prepareReleaseEntrainment(
            cfg, inputSimFiles["releaseScenario"], inputSimLines
        )

        releaseLine = inputSimLines["releaseLine"]
        # check if release features overlap between features
        # TODO: split releaseheight -> question NGI
        dem = rU.readRaster(inputSimFiles["demFile"])
        dem["originalHeader"] = dem["header"].copy()
        # releaseLine = geoTrans.prepareArea(releaseLine, dem, np.sqrt(2), combine=True, checkOverlap=False)
        if len(inputSimLines["relThField"]) == 0:
            # if no release thickness field or function - set release according to shapefile or ini file
            # this is a list of release rasters that we want to combine
            releaseLine = geoTrans.prepareArea(
                releaseLine,
                dem,
                np.sqrt(2),
                thList=releaseLine["thickness"],
                combine=True,
                checkOverlap=False,
            )
            releaseField = releaseLine["rasterData"]
        else:
            # if relTh provided - set release thickness with field or function
            releaseLine = geoTrans.prepareArea(
                releaseLine, dem, np.sqrt(2), combine=True, checkOverlap=False
            )
            relRasterPoly = releaseLine["rasterData"].copy()
            releaseRelThCombined = np.where(relRasterPoly > 0, inputSimLines["relThField"], 0)
            releaseField = releaseRelThCombined

        # Generate the work and data dirs for the current simHash
        cuWorkDir = workDir / key
        workInputDir = cuWorkDir / "Input"
        workOutputDir = cuWorkDir / key

        fU.makeADir(cuWorkDir)
        fU.makeADir(workInputDir)

        zeroRaster = np.full_like(releaseLine["rasterData"], 0)

        release = workInputDir / "release"
        bedDepth = workInputDir / "dummyBedDepth"
        bedShear = workInputDir / "dummyBedShear"
        rU.writeResultToRaster(dem["header"], releaseField, release, flip=True)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedDepth)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedShear)

        # set configuration for MoT-Voellmy
        cfg["Run information"]["Area of Interest"] = cfgMain["MAIN"]["avalancheDir"]
        cfg["Run information"]["UTM zone"] = "32N"
        cfg["Run information"]["EPSG geodetic datum code"] = "31287"
        cfg["Run information"]["Run name"] = cfgMain["MAIN"]["avalancheDir"]
        cfg["File names"]["Grid filename"] = str(inputSimFiles["demFile"])
        cfg["File names"]["Release depth filename"] = str(release) + ".asc"
        cfg["File names"]["Bed depth filename"] = str(bedDepth) + ".asc"
        cfg["File names"]["Bed shear strength filename"] = str(bedShear) + ".asc"
        cfg["File names"]["Output filename root"] = str(workOutputDir)

        rcfFileName = cfgFileDir / (str(key) + ".rcf")

        currentModule = sys.modules[__name__]
        cfgUtils.writeCfgFile(avalancheDir, currentModule, cfg, str(key))
        cfgToRcf(cfg, rcfFileName)
        rcfFiles.append(rcfFileName)
    return rcfFiles


