import os
import sys
import platform
import logging
import numpy as np
import pathlib
import time
import shutil

from avaframe.in3Utils.cfgUtils import cfgToRcf

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
from avaframe.in3Utils.MoTUtils import rewriteDEMtoZeroValues, runAndCheckMoT, MoTGenerateConfigs
from avaframe.in3Utils.initializeProject import _checkForFolderAndDelete

# create local logger
log = logging.getLogger(__name__)


def com8MoTPSAMain(cfgMain, cfgInfo=None):
    # Get all necessary information from the configuration files
    currentModule = sys.modules[__name__]
    simDict, inputSimFiles = MoTGenerateConfigs(cfgMain, cfgInfo, currentModule)

    # convert DEM from nan to 0 values
    # TODO: suggest MoT-PSA to handle nan values
    rewriteDEMtoZeroValues(inputSimFiles["demFile"])

    log.info("The following simulations will be performed")
    for key in simDict:
        log.info("Simulation: %s" % key)

    # Preprocess the simulations, mainly creating the rcf files
    rcfFiles = com8MoTPSAPreprocess(simDict, inputSimFiles, cfgMain)

    # And now we run the simulations
    startTime = time.time()

    log.info("--- STARTING (potential) PARALLEL PART ----")

    # if length of rcfFiles is too long, split it into chunks, this is easier to handle
    # else: run all simulations at once
    chunkSize = 8
    if len(rcfFiles) > chunkSize:
        for i in range(0, len(rcfFiles), chunkSize):
            # splits rcfFiles into segments
            rcfFilesChunk = rcfFiles[i:i + chunkSize]
            simNamesChunk = [p.stem for p in rcfFilesChunk]

            # Get number of CPU Cores wanted
            nCPU = cfgUtils.getNumberOfProcesses(cfgMain, len(rcfFilesChunk))

            # check if there is any simulation to run
            if bool(simNamesChunk):
                # Create parallel pool and run
                # with multiprocessing.Pool(processes=nCPU) as pool:
                with Pool(processes=nCPU) as pool:
                    results = pool.map(com8MoTPSATask, rcfFilesChunk)
                    pool.close()
                    pool.join()

                timeNeeded = "%.2f" % (time.time() - startTime)
                log.info("Overall (parallel) com8MoTPSA computation took: %s s " % timeNeeded)
                log.info("--- ENDING (potential) PARALLEL PART ----")

                # Postprocess the simulations
                com8MoTPSAPostprocess(simNamesChunk, cfgMain, inputSimFiles)

                # Delete folder in Work directory after postprocessing to reduce memory costs
                avaDir = cfgMain['MAIN']['avalancheDir']
                for sim in simNamesChunk:
                    folderName = 'Work/com8MoTPSA/' + sim
                    _checkForFolderAndDelete(avaDir, folderName)
            else:
                log.warning("There is no simulation to be performed for releaseScenario")
    else:
        # Get number of CPU Cores wanted
        nCPU = cfgUtils.getNumberOfProcesses(cfgMain, len(rcfFiles))

        simNames = [p.stem for p in rcfFiles]
        # check if there is any simulation to run
        if bool(simNames):
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
            com8MoTPSAPostprocess(simNames, cfgMain, inputSimFiles)
        else:
            log.warning("There is no simulation to be performed for releaseScenario")


def com8MoTPSAPostprocess(simNames, cfgMain, inputSimFiles):
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]
    # Copy max files to output directory

    outputDir = pathlib.Path(avalancheDir) / "Outputs" / "com8MoTPSA"
    outputDirPeakFile = pathlib.Path(avalancheDir) / "Outputs" / "com8MoTPSA" / "peakFiles"
    fU.makeADir(outputDirPeakFile)

    for key in simNames:
        workDir = pathlib.Path(avalancheDir) / "Work" / "com8MoTPSA" / str(key)

        # Copy DataTime.txt
        dataTimeFile = workDir / "DataTime.txt"
        shutil.copy2(dataTimeFile, outputDir / (str(key) + "_DataTime.txt"))

        # TODO: functionize it
        # Copy ppr files
        pprFiles = list(workDir.glob("*p?_max*"))
        targetFiles = [
            pathlib.Path(str(f.name).replace("null_psa_p1_max", "null_dfa_ppr")) for f in pprFiles
        ]
        targetFiles = [pathlib.Path(str(f).replace("null_psa_p2_max", "null_psa_ppr")) for f in targetFiles]
        targetFiles = [outputDirPeakFile / f for f in targetFiles]
        for source, target in zip(pprFiles, targetFiles):
            shutil.copy2(source, target)

        # Copy pfd files
        pfdFiles = list(workDir.glob("*h?_max*"))
        targetFiles = [
            pathlib.Path(str(f.name).replace("null_psa_h1_max", "null_dfa_pfd")) for f in pfdFiles
        ]
        targetFiles = [pathlib.Path(str(f).replace("null_psa_h2_max", "null_psa_pfd")) for f in targetFiles]
        targetFiles = [outputDirPeakFile / f for f in targetFiles]
        for source, target in zip(pfdFiles, targetFiles):
            shutil.copy2(source, target)

        # Copy pfv files
        pfvFiles = list(workDir.glob("*s?_max*"))
        targetFiles = [
            pathlib.Path(str(f.name).replace("null_psa_s1_max", "null_dfa_pfv")) for f in pfvFiles
        ]
        targetFiles = [pathlib.Path(str(f).replace("null_psa_s2_max", "null_psa_pfv")) for f in targetFiles]
        targetFiles = [outputDirPeakFile / f for f in targetFiles]
        for source, target in zip(pfvFiles, targetFiles):
            shutil.copy2(source, target)

        # write text file to Outputs/com8MoTPSA/configurationFilesDone to indicate that this simulation has been performed
        configFileName = "%s.ini" % key
        for saveDir in ["configurationFilesDone", "configurationFilesLatest"]:
            configDir = pathlib.Path(avalancheDir, "Outputs", "com8MoTPSA", "configurationFiles", saveDir)
            with open((configDir / configFileName), "w") as fi:
                fi.write("see directory configurationFiles for info on config")
            fi.close()

    # create plots and report
    modName = __name__.split(".")[-1]
    reportDir = pathlib.Path(avalancheDir, "Outputs", modName, "reports")
    fU.makeADir(reportDir)
    print(inputSimFiles["demFile"])

    dem = rU.readRaster(inputSimFiles["demFile"])
    # Generate plots for all peakFiles
    oP.plotAllPeakFields(avalancheDir, cfgMain["FLAGS"], modName, demData=dem)


def com8MoTPSATask(rcfFile):
    # TODO: Obvious...
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    command = ["./MoT-PSA", rcfFile]
    # command = ['/home/felix/Versioning/AvaFrame/avaframe/com8MoTPSA/MoT-PSA', rcfFile]
    log.info("Run simulation: %s" % rcfFile)
    runAndCheckMoT(command)
    return command


def com8MoTPSAPreprocess(simDict, inputSimFiles, cfgMain):
    # Load avalanche directory from general configuration file
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    workDir = pathlib.Path(avalancheDir) / "Work" / "com8MoTPSA"
    cfgFileDir = pathlib.Path(avalancheDir) / "Outputs" / "com8MoTPSA" / "configurationFiles"
    fU.makeADir(cfgFileDir)
    rcfFiles = list()

    for key in simDict:
        # Generate command and run via subprocess.run
        # Configuration that needs adjustment

        # load configuration object for current sim
        cfg = simDict[key]["cfgSim"]

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

        releaseL1 = workInputDir / "releaseLayer1"
        releaseL2 = workInputDir / "releaseLayer2"
        bedDepth = workInputDir / "dummyBedDepth"
        bedDepo = workInputDir / "dummyBedDepo"
        bedShear = workInputDir / "dummyBedShear"
        rU.writeResultToRaster(dem["header"], releaseField, releaseL1, flip=True)
        rU.writeResultToRaster(dem["header"], zeroRaster, releaseL2, flip=True)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedDepth)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedDepo)
        rU.writeResultToRaster(dem["header"], zeroRaster, bedShear)

        # set configuration for MoT-PSA
        cfg["Run information"]["Area of Interest"] = cfgMain["MAIN"]["avalancheDir"]
        cfg["Run information"]["UTM zone"] = "32N"
        cfg["Run information"]["EPSG geodetic datum code"] = "31287"
        cfg["Run information"]["Run name"] = cfgMain["MAIN"]["avalancheDir"]
        cfg["File names"]["Grid filename"] = str(inputSimFiles["demFile"])
        cfg["File names"]["Release depth 1 filename"] = str(releaseL1) + ".asc"
        cfg["File names"]["Release depth 2 filename"] = str(releaseL2) + ".asc"
        cfg["File names"]["Bed depth filename"] = str(bedDepth) + ".asc"
        cfg["File names"]["Bed deposition filename"] = str(bedDepo) + ".asc"
        cfg["File names"]["Bed shear strength filename"] = str(bedShear) + ".asc"
        cfg["File names"]["Output filename root"] = str(workOutputDir)

        rcfFileName = cfgFileDir / (str(key) + ".rcf")
        currentModule = sys.modules[__name__]
        cfgUtils.writeCfgFile(avalancheDir, currentModule, cfg, str(key))
        cfgToRcf(cfg, rcfFileName)
        rcfFiles.append(rcfFileName)
    return rcfFiles


