"""Run script for running simulations in parallel based on input from a regional folder containing
 multiple avaFolders"""

import time
import pathlib
import argparse
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

from avaframe.com7Regional import com7Regional
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

def runCom7Regional(avalancheDir=''):
    """
    Parameters
    ----------
    avalancheDir:
        input directory which should contain multiple valid avalanche directories, i.e. the regional directory!

    Returns
    -------
    outDir: str
        path to output directory
    """

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName='runCom7Regional')
    log.info('MAIN SCRIPT')

    # Time the whole routine
    startTime = time.time()

    # Load the (REGIONAL) avalanche directory from general configuration file
    # More information about the configuration can be found
    # on the Configuration page in the documentation
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    #Load regional configuration
    cfgCom7 = cfgUtils.getModuleConfig(com7Regional, fileOverride='', toPrint=False, onlyDefault=False)

    # List avalanche directories that contain an 'Inputs' folder - logic is that they are valid avaFolders
    avaDirs = [avadir for avadir in pathlib.Path(avalancheDir).iterdir() if avadir.is_dir() and
               (avadir / "Inputs").is_dir()]
    log.info(f"Found a total of '{len(avaDirs)}' avalanche directories in '{avalancheDir}'.")
    log.info(f"Processing the following directories in this order:")
    for avadir in avaDirs:
        log.info(f"'{avadir.name}'")

    # Run processes in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor() as executor:
        # Submit each avalanche directory to the executor
        futures = {executor.submit(com7Regional.processAvaDirCom1DFA, cfgMain, cfgCom7, avadir):
                       avadir for avadir in avaDirs}
        # Log results as each future completes
        for future in as_completed(futures):
            avadir = futures[future]
            try:
                resultDir, status = future.result()
                log.info(f"{status} in directory: {resultDir.relative_to(pathlib.Path(avalancheDir))}"
                         f" at {time.time() - startTime:.1f} s")
            except Exception as e:
                log.error(f"Error processing {avadir}: {e}")

    # Keep in mind, this is only com1DFA specific
    # Move all files from the 'Outputs/com1DFA/peakFiles' folder from each of the subfolders (avadir) to a folder
    # called (allPeakFiles) in the main regional folder
    if cfgCom7['GENERAL'].getboolean('movePeakFiles'):
        allPeakFilesDir = pathlib.Path(avalancheDir, 'allPeakFiles')
        if allPeakFilesDir.exists():
            shutil.rmtree(str(allPeakFilesDir))
        allPeakFilesDir.mkdir(parents=True, exist_ok=True)

        allTimeStepsDir = allPeakFilesDir / 'allTimeSteps'
        allTimeStepsDir.mkdir(parents=True, exist_ok=True)

        for avadir in avaDirs:
            peakFilesDir = pathlib.Path(avadir, 'Outputs', 'com1DFA', 'peakFiles')
            if peakFilesDir.is_dir():
                for file in peakFilesDir.iterdir():
                    if file.is_file():  # Only move files, not directories
                        shutil.move(str(file), str(allPeakFilesDir))
                        log.info(f"Moved {file.name} to {allPeakFilesDir}")
                timeStepsDir = peakFilesDir / 'timeSteps'
                if timeStepsDir.is_dir():
                    for file in timeStepsDir.iterdir():
                        shutil.move(str(file), str(allTimeStepsDir))
                        log.debug(f"Moved {file.name} to {allTimeStepsDir}")

        log.info(f"Moved all peak files to {allPeakFilesDir}")

    # Print time needed
    endTime = time.time()
    log.info('Regional process finished after %.1f seconds.' % (endTime - startTime))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run regional workflow')
    parser.add_argument('avadir', metavar='avalancheDir', type=str, nargs='?', default='',
                        help='the avalanche directory')

    args = parser.parse_args()
    runCom7Regional(str(args.avadir))
