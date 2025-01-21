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

    Note: This function calls com1DFA within each avalanche directory within the input directory. 
    If wanted it may be used as a template to call another operation within each directory, such as com2AB, ana5Utils, etc.

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
    log.info('COM1DFA PROCESS CALLED BY COM7REGIONAL RUN')
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
    if not cfg['GENERAL'].getboolean('copyPeakFiles'):
        log.info("copyPeakFiles is False - no files will be copied or moved")
        return None, None

    # Set up dirs
    allPeakFilesDir = pathlib.Path(avalancheDir, 'allPeakFiles')
    allTimeStepsDir = allPeakFilesDir / 'allTimeSteps'
    
    # Create fresh dirs, remove old ones
    for dirPath in [allPeakFilesDir, allTimeStepsDir]:
        if dirPath.exists():
            shutil.rmtree(str(dirPath))
        dirPath.mkdir(parents=True, exist_ok=True)

    # Set file operation based on settings (move or copy)
    fileOp = shutil.move if cfg['GENERAL'].getboolean('moveInsteadOfCopy', False) else shutil.copy
    opName = 'Moving' if fileOp == shutil.move else 'Copying'

    # Process files
    nFiles = {'peak': 0, 'timestep': 0}
    for avaDir in avaDirs:
        peakFilesDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'peakFiles')
        if not peakFilesDir.is_dir():
            continue
            
        for file in peakFilesDir.glob('*.asc'):
            fileOp(str(file), str(allPeakFilesDir))
            nFiles['peak'] += 1
            log.debug(f"{opName} {file} to {allPeakFilesDir}")
        
        timeStepsDir = peakFilesDir / 'timeSteps'
        if timeStepsDir.is_dir():
            for file in timeStepsDir.glob('*.asc'):
                fileOp(str(file), str(allTimeStepsDir))
                nFiles['timestep'] += 1
                log.debug(f"{opName} {file} to {allTimeStepsDir}")

    log.debug(f"{opName} completed: {nFiles['peak']} peak files and {nFiles['timestep']} timestep files processed")
    return allPeakFilesDir, allTimeStepsDir
