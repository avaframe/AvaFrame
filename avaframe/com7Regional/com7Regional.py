"""Module for handling regional avalanche simulations."""

import pathlib
import shutil
import logging

import avaframe.in3Utils.initializeProject as initProj
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils, cfgHandling
from avaframe.in3Utils import logUtils

# create local logger
log = logging.getLogger(__name__)

def processAvaDirCom1Regional(cfgMain, cfgCom7, avalancheDir):
    """Run com1DFA simulation in a specific avalanche directory with regional settings.

    Parameters
    ----------
    cfgMain : configparser.ConfigParser
        Main configuration settings
    cfgCom7 : configparser.ConfigParser
        Regional configuration settings with potential overrides
    avalancheDir : pathlib.Path or str
        Path to the avalanche directory to process

    Returns
    -------
    tuple
        (avalancheDir, status) where status is "Success" if simulation completed
    """

    # Initialize log for each process
    log = logUtils.initiateLogger(avalancheDir, logName='runCom1DFA')
    log.info('COM1DFA PROCESS CALLED BY REGIONAL RUN')
    log.info('Current avalanche: %s', avalancheDir)

    # Update cfgMain setting to reflect the current avalancheDir
    cfgMain['MAIN']['avalancheDir'] = str(avalancheDir)

    # Clean input directory of old work and output files from module
    initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)

    # Create com1DFA configuration for the current avalanche directory and override with regional settings
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False,
                                          onlyDefault=cfgCom7['com1DFA_com1DFA_override'].getboolean('defaultConfig'))
    cfgCom1DFA, cfgCom7 = cfgHandling.applyCfgOverride(cfgCom1DFA, cfgCom7, com1DFA, addModValues=False)

    # Run com1DFA in the current avalanche directory
    com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgCom1DFA)

    return avalancheDir, "Success"

def findAvaDirs(Dir):
    """Find all valid avalanche directories within a given directory.

    A directory is considered a valid avalanche directory if it contains an "Inputs" folder.

    Parameters
    ----------
    Dir : pathlib.Path or str
        Path to the directory to search in

    Returns
    -------
    list
        List of pathlib.Path objects pointing to valid avalanche directories

    Notes
    -----
    Logs the total number and names of found avalanche directories
    """
    avaDirs = [pathlib.Path(p).parent for p in pathlib.Path(Dir).glob("*/Inputs")]
    log.info(f"Found a total of '{len(avaDirs)}' avalanche directories in: {Dir}:")
    for avaDir in avaDirs:
        log.info(f"'{avaDir.name}'")

    return avaDirs

def moveOrCopyFile(src, dst, copy=False):
    """Move or copy a file from source to destination location.

    Parameters
    ----------
    src : pathlib.Path or str
        Source file location
    dst : pathlib.Path or str
        Destination file location
    copy : bool, optional
        If True, copy the file; if False, move it (default: False)
    """
    if copy:
        shutil.copy(str(src), str(dst))
        log.debug(f"Copied {src} to {dst}")
    else:
        shutil.move(str(src), str(dst))
        log.debug(f"Moved {src} to {dst}")

def moveOrCopyPeakFiles(cfg, avalancheDir, avaDirs):
    """Consolidate peak files from multiple avalanche directories.

    Creates two directories:
    1. allPeakFiles: Contains peak files from all avalanche directories
    2. allPeakFiles/allTimeSteps: Contains time step files from all avalanche directories

    Parameters
    ----------
    cfg : configparser.ConfigParser
        Configuration containing GENERAL.copyPeakFiles setting
    avalancheDir : pathlib.Path or str
        Base directory where allPeakFiles will be created
    avaDirs : list
        List of avalanche directories to process

    Returns
    -------
    tuple
        (allPeakFilesDir, allTimeStepsDir) paths to the created directories
    """

    # Get setting from cfg
    copyPeakFiles = cfg['GENERAL'].getboolean('copyPeakFiles')

    fileOperation = moveOrCopyFile
    fileOperationKwargs = {'copy': copyPeakFiles}

    # Set up allPeakFilesDir and allTimeStepsDir
    allPeakFilesDir = pathlib.Path(avalancheDir, 'allPeakFiles')
    if allPeakFilesDir.exists():
        shutil.rmtree(str(allPeakFilesDir)) #remove it if it already exists
    allPeakFilesDir.mkdir(parents=True, exist_ok=True)
    allTimeStepsDir = allPeakFilesDir / 'allTimeSteps'
    if allTimeStepsDir.exists():
        shutil.rmtree(str(allTimeStepsDir)) #remove it if it already exists
    allTimeStepsDir.mkdir(parents=True, exist_ok=True)

    # Move or copy the peak files and time steps
    for avadir in avaDirs:
        peakFilesDir = pathlib.Path(avadir, 'Outputs', 'com1DFA', 'peakFiles')
        if peakFilesDir.is_dir():
            for file in peakFilesDir.glob('*.asc'):
                fileOperation(file, allPeakFilesDir, **fileOperationKwargs)

            timeStepsDir = peakFilesDir / 'timeSteps'
            if timeStepsDir.is_dir():
                for file in timeStepsDir.glob('*.asc'):
                    fileOperation(file, allTimeStepsDir, **fileOperationKwargs)

    return allPeakFilesDir, allTimeStepsDir
