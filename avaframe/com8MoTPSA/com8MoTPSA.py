import os
import subprocess
import platform
import logging
import pandas as pd
import numpy as np
import pathlib
import time
import shutil

if os.name == "nt":
    from multiprocessing.pool import ThreadPool as Pool
elif platform.system() == "Darwin":
    from multiprocessing.pool import ThreadPool as Pool
else:
    from multiprocessing import Pool

import avaframe.com1DFA.com1DFATools as com1DFATools
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in2Trans import rasterUtils as rU
from avaframe.com1DFA import particleInitialisation as pI
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.out1Peak import outPlotAllPeak as oP

import avaframe.com8MoTPSA.com8MoTPSA as com8MoTPSA

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
                            "Process is running. Reported time steps: "
                            + str(counter)
                    )
                    log.info(msg)
                    printCounter = 0

            elif "find_dt" in line:
                continue
            elif "h1" in line:
                continue
            elif "h2" in line:
                continue
            elif "write_data" in line:
                continue
            else:
                log.info(line)


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
    rU.writeResultToRaster(demData["header"], demData["rasterData"], newFileName, flip=True)




def com8MoTPSAMain(cfgMain, cfgInfo=None):
    # Get all necessary information from the configuration files


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

    # Preprocess the simulations, mainly creating the rcf files
    rcfFiles = com8MoTPSAPreprocess(simDict, inputSimFiles, cfgMain, cfgInfo)

    # And now we run the simulations
    startTime = time.time()

    log.info("--- STARTING (potential) PARALLEL PART ----")

    # Get number of CPU Cores wanted
    nCPU = cfgUtils.getNumberOfProcesses(cfgMain, len(rcfFiles))

    # Create parallel pool and run
    # with multiprocessing.Pool(processes=nCPU) as pool:
    with Pool(processes=nCPU) as pool:
        results = pool.map(com8MoTPSATask, rcfFiles)
        pool.close()
        pool.join()

    timeNeeded = "%.2f" % (time.time() - startTime)
    log.info("Overall (parallel) com8MoTPSA computation took: %s s " % timeNeeded)
    log.info("--- ENDING (potential) PARALLEL PART ----")

    # Postprocess the simulations
    com8MoTPSAPostprocess(simDict, cfgMain, inputSimFiles)


def com8MoTPSAPostprocess(simDict, cfgMain, inputSimFiles):
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]
    # Copy max files to output directory

    outputDir = pathlib.Path(avalancheDir) / "Outputs" / "com8MoTPSA"
    outputDirPeakFile = pathlib.Path(avalancheDir) / "Outputs" / "com8MoTPSA" / "peakFiles"
    fU.makeADir(outputDirPeakFile)

    for key in simDict:
        workDir = pathlib.Path(avalancheDir) / "Work" / "com8MoTPSA" / str(key)

        # Copy DataTime.txt
        dataTimeFile = workDir / "DataTime.txt"
        shutil.copy2(dataTimeFile, outputDir / (str(key) + "_DataTime.txt"))

        # TODO: functionize it
        # Copy ppr files
        pprFiles = list(workDir.glob("*p?_max*"))
        targetFiles = [pathlib.Path(str(f.name).replace('null_psa_p1_max', 'lay1_psa_ppr')) for f in pprFiles]
        targetFiles = [pathlib.Path(str(f).replace('null_psa_p2_max', 'lay2_psa_ppr')) for f in targetFiles]
        targetFiles = [outputDirPeakFile / f for f in targetFiles]
        for source, target in zip(pprFiles, targetFiles):
            shutil.copy2(source, target)

        # Copy pfd files
        pfdFiles = list(workDir.glob("*h?_max*"))
        targetFiles = [pathlib.Path(str(f.name).replace('null_psa_h1_max', 'lay1_psa_pfd')) for f in pfdFiles]
        targetFiles = [pathlib.Path(str(f).replace('null_psa_h2_max', 'lay2_psa_pfd')) for f in targetFiles]
        targetFiles = [outputDirPeakFile / f for f in targetFiles]
        for source, target in zip(pfdFiles, targetFiles):
            shutil.copy2(source, target)

        # Copy pfv files
        pfvFiles = list(workDir.glob("*s?_max*"))
        targetFiles = [pathlib.Path(str(f.name).replace('null_psa_s1_max', 'lay1_psa_pfv')) for f in pfvFiles]
        targetFiles = [pathlib.Path(str(f).replace('null_psa_s2_max', 'lay2_psa_pfv')) for f in targetFiles]
        targetFiles = [outputDirPeakFile / f for f in targetFiles]
        for source, target in zip(pfvFiles, targetFiles):
            shutil.copy2(source, target)

    # create plots and report
    modName = __name__.split('.')[-1]
    reportDir = pathlib.Path(avalancheDir, "Outputs", modName, "reports")
    fU.makeADir(reportDir)
    print(inputSimFiles["demFile"])

    dem = rU.readRaster(inputSimFiles["demFile"])
    # Generate plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avalancheDir, cfgMain["FLAGS"], modName, demData=dem)


def com8MoTPSATask(rcfFile):
    command = ['./MoT-PSA', rcfFile]
    log.info("Run simulation: %s" % rcfFile)
    _runAndCheck(command)
    return command


def com8MoTPSAPreprocess(simDict, inputSimFiles, cfgMain, cfgInfo):
    # Load avalanche directory from general configuration file
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    workDir = pathlib.Path(avalancheDir) / "Work" / "com8MoTPSA"
    cfgFileDir = pathlib.Path(avalancheDir) / "Outputs" / "com8MoTPSA" / "configurationFiles"
    fU.makeADir(cfgFileDir)
    rcfFiles = list()

    for key in simDict:
        # Generate command and run via subprocess.run
        # Configuration that needs adjustment

        simDF = pd.DataFrame()

        # load configuration object for current sim
        cfg = simDict[key]["cfgSim"]

        # fetch simHash for current sim
        simHash = simDict[key]["simHash"]

        # # append configuration to dataframe
        # simDF = cfgUtils.appendCgf2DF(simHash, key, cfg, simDF)

        # convert release shape to raster with values for current sim
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

        # Generate the work and data dirs for the current simHash

        cuWorkDir = workDir / key
        workInputDir = cuWorkDir / "Input"
        workOutputDir = cuWorkDir / key

        fU.makeADir(cuWorkDir)
        fU.makeADir(workInputDir)

        zeroRaster = np.full_like(releaseLine["rasterData"], 0)

        releaseL1 = workInputDir / "releaseLayer1"
        releaseL2 = workInputDir / "releaseLayer2"
        bedDepth = workInputDir / "dummyBedDepth"
        bedDepo = workInputDir / "dummyBedDepo"
        bedShear = workInputDir / "dummyBedShear"
        rU.writeResultToRaster(dem["header"], releaseLine["rasterData"], releaseL1, flip=True)
        rU.writeResultToRaster(dem["header"], releaseLine["rasterData"], releaseL2, flip=True)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedDepth)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedDepo)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedShear)

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
        cfgInfo["File names"]["Output filename root"] = "./" + str(workOutputDir)

        rcfFileName = cfgFileDir / (str(key) + ".rcf")

        cfgUtils.writeCfgFile(avalancheDir, com8MoTPSA, cfgInfo, str(key))

        cfgToRcf(cfgInfo, rcfFileName)
        rcfFiles.append(rcfFileName)

    return rcfFiles
