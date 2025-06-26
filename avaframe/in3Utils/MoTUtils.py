import os
import platform
import subprocess
import logging
import pathlib
import shutil

import numpy as np

from avaframe.com1DFA import com1DFATools as com1DFATools, com1DFA as com1DFA
from avaframe.in2Trans import rasterUtils as rU

log = logging.getLogger(__name__)

def rewriteDEMtoZeroValues(demFile):
    """Set all NaN values in a DEM raster to zero and update the nodata value.

    This function reads a DEM raster file, replaces all NaN values with 0.0,
    updates the nodata value in the header to 0.0, and writes the modified
    raster back to a new file.

    Parameters
    ----------
    demFile : pathlib.Path
        Path to the input DEM raster file

    Returns
    -------
    None
        Writes a new raster file with zero values instead of NaN values.
        The output file is saved in the same directory as the input file,
        using the same stem name.

    Notes
    -----
    The function uses the rasterUtils module for reading and writing raster data.
    The output raster is flipped during writing.
    """
    demData = rU.readRaster(demFile)
    demData["rasterData"][np.isnan(demData["rasterData"])] = 0.0
    demData["header"]["nodata_value"] = 0.0
    newFileName = demFile.parent / demFile.stem
    rU.writeResultToRaster(demData["header"], demData["rasterData"], newFileName, flip=True)


def runAndCheckMoT(command):
    """Execute MoT command and monitor its output.

    This function runs a MoT command as a subprocess and monitors its output,
    filtering and logging specific messages while tracking time steps.

    Parameters
    ----------
    command : str or list
        The command to execute. Can be a string or list of arguments.

    Returns
    -------
    None
        Function runs the command and logs output but does not return a value.

    Notes
    -----
    - Uses different shell settings based on operating system
    - Filters output to reduce noise from common status messages
    - Logs time step progress every 100 steps
    - Handles UTF-8 encoding with replacement of invalid characters
    """
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
                    msg = "Process is running. Reported time steps: " + str(counter)
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
            elif "update_boundaries" in line:
                continue
            else:
                log.info(line)


def MoTGenerateConfigs(cfgMain, cfgInfo, currentModule):
    """
    Creates configuration objects for com8MoTPSA.

    Parameters
    ------------
    cfgMain: configparser object
        main configuration of AvaFrame
    cfgInfo: str or pathlib Path or configparser object
        path to configuration file if overwrite is desired - optional
        if not local (if available) or default configuration will be loaded
        if cfgInfo is a configparser object take this as initial config
    currentModule: module object
        is being passed to cfgUtils.writeCfgFile to create the correct the cfg


    Returns
    --------
    simDict: dict
        dictionary with one key per simulation to perform including its config object
    inputSimFiles: dict
        dictionary with input files info
    """

    # fetch type of cfgInfo
    typeCfgInfo = com1DFATools.checkCfgInfoType(cfgInfo)

    if typeCfgInfo == "cfgFromDir":
        # preprocessing to create configuration objects for all simulations to run by reading multiple cfg files
        simDict, inputSimFiles, simDFExisting, outDir = com1DFATools.createSimDictFromCfgs(
            cfgMain, cfgInfo, module=currentModule
        )
    else:
        # preprocessing to create configuration objects for all simulations to run
        simDict, outDir, inputSimFiles, simDFExisting = com1DFA.com1DFAPreprocess(
            cfgMain, typeCfgInfo, cfgInfo, module=currentModule
        )
    return simDict, inputSimFiles

def copyMoTFiles(workDir, outputDir, searchString, replaceString):
    """
    Copy and rename MoT result files from work directory to output directory.

    Parameters
    ----------
    workDir : pathlib.Path
        Source directory containing the original p_max files
    outputDir : pathlib.Path
        Destination directory where renamed ppr files will be copied to
    searchString : str
        String pattern to search for in the original filenames
    replaceString : str
        String to replace the searchString with in the new filenames

    Returns
    -------
    None
        Files are copied to the destination directory with renamed extensions
    """
    varFiles = list(workDir.glob("*"+searchString+"*"))
    targetFiles = [
        pathlib.Path(str(f.name).replace(searchString, replaceString))
        for f in varFiles
    ]
    targetFiles = [outputDir / f for f in targetFiles]

    for source, target in zip(varFiles, targetFiles):
        shutil.copy2(source, target)