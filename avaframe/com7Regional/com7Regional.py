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
    """Function to call com1DFA in each avalanche directory with regional override settings."""

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
    """Function to find all valid avalanche directories within a directory based on if there is an "Inputs" folder."""

    avaDirs = [avaDir for avaDir in pathlib.Path(Dir).iterdir() if avaDir.is_dir() and
               (avaDir / "Inputs").is_dir()]
    log.info(f"Found a total of '{len(avaDirs)}' avalanche directories in '{Dir}':")
    for avaDir in avaDirs:
        log.info(f"'{avaDir.name}'")

    return avaDirs

def moveOrCopyFile(src, dst, copy=False):
    """Function to move or copy a file from a source location to a destination location

    Parameters
    ----------
    src: path object
        the source location
    dst: path object
        destination location
    copy: bool
        whether to copy or move the file
    """
    if copy:
        shutil.copy(str(src), str(dst))
        log.debug(f"Copied {src} to {dst}")
    else:
        shutil.move(str(src), str(dst))
        log.debug(f"Moved {src} to {dst}")

def moveOrCopyPeakFiles(cfg, avalancheDir, avaDirs):
    """Function to move or copy peak files from each avalanche directory to a directory called allPeakFiles.
    Also copy all timeSteps to an allTimeSteps directory"""

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
