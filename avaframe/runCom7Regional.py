"""Run script for running simulations in parallel based on input from a regional folder containing
 multiple avaFolders"""

import time
import pathlib
import argparse
import os
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

    # Time the whole routine
    startTime = time.time()

    # Load the avalanche directory from command input or the general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Set the regional directory
    regionalDir = os.path.join(avalancheDir, 'Regional')

    # Start logging
    log = logUtils.initiateLogger(regionalDir, logName='runCom7Regional')
    log.info('MAIN SCRIPT')

    #Load regional configuration
    cfgCom7 = cfgUtils.getModuleConfig(com7Regional, fileOverride='', toPrint=False, onlyDefault=False)

    # List valid avalanche directories within the regional directory
    avaDirs = com7Regional.findAvaDirs(regionalDir)

    # Get number of processes to perform #ToDo: get nVariations as well, adjust log output. need to adjust function?
    nProcesses = cfgUtils.getNumberOfProcesses(cfgMain, len(avaDirs))

    # Set nCPU for com1 to 1 to avoid dual parallelization, i.e. each subAvaDir variation is processed sequentially
    cfgMain['MAIN']['nCPU'] = '1'

    # Process each avalanche directory within the regional folder in parallel
    with ProcessPoolExecutor(max_workers=nProcesses) as executor:
        # Submit each avalanche directory to the executor
        futures = {executor.submit(com7Regional.processAvaDirCom1Regional, cfgMain, cfgCom7, avaDir):
                       avaDir for avaDir in avaDirs}
        # Log results as each future completes
        for future in as_completed(futures):
            avaDir = futures[future]
            try:
                resultDir, status = future.result()
                log.info(f"{status} in directory: {resultDir.relative_to(pathlib.Path(regionalDir))}"
                         f" at {time.time() - startTime:.1f} s")
            except Exception as e:
                log.error(f"Error processing {avaDir}: {e}")

    # Move or copy files from the 'Outputs/com1DFA/peakFiles' folder from each of the subfolders (avaDir) to a folder
    # called (allPeakFiles) in the main regional folder. Clear them if they already exist
    # Keep in mind, this is only com1DFA specific for now, also keep in mind that the information where it comes from
    # is lost when it's moved (i.e. things like layer rename is no longer possible when they are imported to QGIS)
    # preliminary feature until an "import output rasters" is added to QGIS
    if cfgCom7['GENERAL'].getboolean('movePeakFiles'):
        allPeakFilesDir, allTimeStepsDir = com7Regional.moveOrCopyPeakFiles(cfgCom7, regionalDir, avaDirs)

        copyPeakFiles = cfgCom7['GENERAL'].getboolean('copyPeakFiles')
        log.info(f"{'Copied' if copyPeakFiles else 'Moved'} peakFiles to "
                  f"{allPeakFilesDir}")
        log.info(f"{'Copied' if copyPeakFiles else 'Moved'} timeSteps to "
                  f"{allTimeStepsDir}")

    # Print time needed
    endTime = time.time()
    log.info('Regional process finished after %.1f seconds.' % (endTime - startTime))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run regional workflow')
    parser.add_argument('avadir', metavar='avalancheDir', type=str, nargs='?', default='',
                        help='the avalanche directory')

    args = parser.parse_args()
    runCom7Regional(str(args.avadir))
