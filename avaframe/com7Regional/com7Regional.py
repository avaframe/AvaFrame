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
    """
    avaDirs = [pathlib.Path(p).parent for p in pathlib.Path(Dir).glob("*/Inputs")]
    log.info(f"Found a total of '{len(avaDirs)}' avalanche directories in: {Dir}:")
    for avaDir in avaDirs:
        log.info(f"'{avaDir.name}'")

    return avaDirs

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

def moveOrCopyPeakFiles(cfg, avalancheDir, avaDirs):
    """Consolidate peak files from multiple avalanche directories.

    Creates two directories:
    1. allPeakFiles: Contains peak files from all avalanche directories
    2. allPeakFiles/allTimeSteps: Contains time step files from all avalanche directories

    Parameters
    ----------
    cfg : configparser.ConfigParser
        Configuration containing GENERAL settings:
        - copyPeakFiles: If True, copy/move files; if False, do nothing
        - moveInsteadOfCopy: If True, move files instead of copying
    avalancheDir : pathlib.Path or str
        Base directory where allPeakFiles will be created
    avaDirs : list
        List of avalanche directories to process

    Returns
    -------
    tuple
        (allPeakFilesDir, allTimeStepsDir) paths to the created directories
    """
    # Get settings from cfg
    copyPeakFiles = cfg['GENERAL'].getboolean('copyPeakFiles')
    moveInsteadOfCopy = cfg['GENERAL'].getboolean('moveInsteadOfCopy', False)

    if not copyPeakFiles:
        log.info("copyPeakFiles is False - no files will be copied or moved")
        return None, None

    # Set up directories
    allPeakFilesDir = pathlib.Path(avalancheDir, 'allPeakFiles')
    if allPeakFilesDir.exists():
        shutil.rmtree(str(allPeakFilesDir))  #remove it if it already exists #todo: developer convenience - remove
    allPeakFilesDir.mkdir(parents=True, exist_ok=True)

    allTimeStepsDir = allPeakFilesDir / 'allTimeSteps'
    if allTimeStepsDir.exists():
        shutil.rmtree(str(allTimeStepsDir))  #remove it if it already exists #todo: developer convenience - remove
    allTimeStepsDir.mkdir(parents=True, exist_ok=True)

    # Gather and process files
    peakFiles = []
    timeStepFiles = []
    for avaDir in avaDirs:
        peakFilesDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'peakFiles')
        if peakFilesDir.is_dir():
            peakFiles.extend(list(peakFilesDir.glob('*.asc')))
            
            timeStepsDir = peakFilesDir / 'timeSteps'
            if timeStepsDir.is_dir():
                timeStepFiles.extend(list(timeStepsDir.glob('*.asc')))
    
    fileOp = shutil.move if moveInsteadOfCopy else shutil.copy
    
    for file in peakFiles:
        fileOp(str(file), str(allPeakFilesDir))
        log.debug(f"{'Moved' if moveInsteadOfCopy else 'Copied'} {file} to {allPeakFilesDir}")

    for file in timeStepFiles:
        fileOp(str(file), str(allTimeStepsDir))
        log.debug(f"{'Moved' if moveInsteadOfCopy else 'Copied'} {file} to {allTimeStepsDir}")

    action = "Moving" if moveInsteadOfCopy else "Copying"
    log.debug(f"{action} completed: {len(peakFiles)} peak files and {len(timeStepFiles)} timestep files processed")

    return allPeakFilesDir, allTimeStepsDir
