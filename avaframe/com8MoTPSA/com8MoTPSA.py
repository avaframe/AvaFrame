import os
import subprocess
import platform
import logging
import pandas as pd
import numpy as np

import avaframe.com1DFA.com1DFATools as com1DFATools
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in2Trans import rasterUtils as rU
from avaframe.com1DFA import particleInitialisation as pI
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.geoTrans as geoTrans

# create local logger
log = logging.getLogger(__name__)


# TODO move to whereever?
def _runAndCheck(command):
    if os.name == "nt":
        useShell = True
    elif platform.system() == "Darwin":
        useShell = False
    else:
        useShell = False

    # This starts the subprocess
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=useShell,
        encoding="utf-8",
        errors="replace",
        universal_newlines=True,
    )

    printCounter = 0
    counter = 1

    while True:
        realtimeOutput = process.stdout.readline()

        if realtimeOutput == "" and process.poll() is not None:
            break

        if realtimeOutput:
            line = realtimeOutput.strip()

            # do not pollute output window with time step prints
            # TODO: hacky for now
            if "Step" in line:
                counter = counter + 1
                printCounter = printCounter + 1
                if printCounter > 100:
                    # print('\r' + line, flush=True, end='')
                    msg = (
                            "Process is running. Reported time steps (all sims): "
                            + str(counter)
                    )
                    print(msg)
                    printCounter = 0

            elif "find_dt" in line:
                continue
            elif "h1" in line:
                continue
            elif "h2" in line:
                continue
            else:
                print(line)


# TODO move to utils
def cfgToRcf(cfg, fileName):
    with open(fileName, "w") as f:
        for section in cfg.sections():
            if section in ("FOREST_EFFECTS", "ENTRAINMENT"):
                pass
            elif section in ("GENERAL", "INPUT"):
                continue
            else:
                f.write(f"# {section.replace('_', ' ')}\n")
                f.write("#\n")
            for key, value in cfg.items(section):
                # key = key.replace('_', ' ')
                key = key.strip()
                f.write(f"{key:<40}{value}\n")
            f.write("#\n")


def rewriteDEMtoZeroValues(demFile):
    demData = rU.readRaster(demFile)
    demData["rasterData"][np.isnan(demData["rasterData"])] = 0.0
    demData["header"]["nodata_value"] = 0.0
    newFileName = demFile.parent / demFile.stem
    rU.writeResultToRaster(demData["header"], demData["rasterData"], newFileName)




def com8MoTPSAMain(cfgMain, cfgInfo=None):
    # Get all necessary information from the configuration files

    # Load avalanche directory from general configuration file
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    # fetch type of cfgInfo
    typeCfgInfo = com1DFATools.checkCfgInfoType(cfgInfo)

    # preprocessing to create configuration objects for all simulations to run
    simDict, outDir, inputSimFiles, simDFExisting = com1DFA.com1DFAPreprocess(cfgMain, typeCfgInfo, cfgInfo)

    # convert DEM from nan to 0 values
    # TODO: suggest MoT-PSA to handle nan values
    rewriteDEMtoZeroValues(inputSimFiles["demFile"])


    log.info("The following simulations will be performed")
    for key in simDict:
        log.info("Simulation: %s" % key)
        exportFlag = simDict[key]["cfgSim"]["EXPORTS"].getboolean("exportData")

    for key in simDict:
        # Generate command and run via subprocess.run
        # Configuration that needs adjustment

        simDF = pd.DataFrame()

        # load configuration object for current sim
        cfg = simDict[key]["cfgSim"]

        # fetch simHash for current sim
        simHash = simDict[key]["simHash"]

        # append configuration to dataframe
        simDF = cfgUtils.appendCgf2DF(simHash, key, cfg, simDF)

        # convert release shape to raster with values for current sim
        relFile = simDict[key]["relFile"]

        # select release area input data according to chosen release scenario
        inputSimFiles = gI.selectReleaseFile(inputSimFiles, cfg["INPUT"]["releaseScenario"])
        # create required input from input files
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
        releaseLine = geoTrans.prepareArea(releaseLine, dem, np.sqrt(2), combine=True, checkOverlap=False)
        releaseL1 = inputSimFiles["demFile"].parent / "releaseLine1"
        releaseL2 = inputSimFiles["demFile"].parent / "releaseLine2"
        bedDepth = inputSimFiles["demFile"].parent / "dummyBedDepth"
        bedDepo = inputSimFiles["demFile"].parent / "dummyBedDepo"
        bedShear = inputSimFiles["demFile"].parent / "dummyBedShear"
        rU.writeResultToRaster(dem["header"], releaseLine["rasterData"], releaseL1)
        rU.writeResultToRaster(dem["header"], releaseLine["rasterData"], releaseL2)
        emptyRaster = np.full_like(releaseLine["rasterData"], 0)
        rU.writeResultToRaster(dem["header"], emptyRaster, bedDepth)
        rU.writeResultToRaster(dem["header"], emptyRaster, bedDepo)
        rU.writeResultToRaster(dem["header"], emptyRaster, bedShear)

        # set configuration for MoT-PSA
        cfgInfo["Run information"]["Area of Interest"] = cfgMain["MAIN"]["avalancheDir"]
        cfgInfo["Run information"]["UTM zone"] = "32N"
        cfgInfo["Run information"]["EPSG geodetic datum code"] = "31287"
        cfgInfo["Run information"]["Run name"] = cfgMain["MAIN"]["avalancheDir"]
        cfgInfo["File names"]["Grid filename"] = "./" + str(inputSimFiles["demFile"])
        cfgInfo["File names"]["Release depth 1 filename"] = "./" + str(releaseL1) + ".asc"
        cfgInfo["File names"]["Release depth 2 filename"] = "./" + str(releaseL2) + ".asc"
        cfgInfo["File names"]["Bed depth filename"] = "./" + str(bedDepth) + ".asc"
        cfgInfo["File names"]["Bed deposition filename"] = "./" + str(bedDepo) + ".asc"
        cfgInfo["File names"]["Bed shear strength filename"] = "./" + str(bedShear) + ".asc"

        # select release area input data according to chosen release scenario
        # inputSimFiles = gI.selectReleaseFile(inputSimFiles, cfg["INPUT"]["releaseScenario"])

        # create required input from input files
        # demOri, inputSimLines = prepareInputData(inputSimFiles, cfg)

        rcfFile = key + ".rcf"

        cfgToRcf(cfgInfo, rcfFile)

        command = ['./MoT-PSA', rcfFile]
        _runAndCheck(command)

        sys.exit()
