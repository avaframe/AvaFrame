"""
Run script for running simulations in parallel based on input from a regional folder containing
multiple avaFolders
"""

import time
import pathlib
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

from avaframe.com7Regional import com7Regional as com7
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

def runCom7Regional(avalancheDir=''):
    """
    Parameters
    ----------
    avalancheDir: str
        input directory which should contain multiple valid avalanche directories within the regionalDir

    Returns
    -------
    outDir: str
        path to output directory

    Notes
    -----
    The function expects the following directory structure:
    avalancheDir/
    └── regionalDir/
        ├── avalanche1/
        ├── avalanche2/
        └── ...

    Where the 'regionalDir' contains multiple individual avalanche directories
    that will be processed in parallel.
    """

    # Time the whole routine
    startTime = time.time()

    # Load the avalanche directory from command input or the general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Start logging
    log = logUtils.initiateLogger(str(avalancheDir), logName='runCom7Regional')
    log.info('MAIN SCRIPT')

    # Load module configuration
    cfg = cfgUtils.getModuleConfig(com7, fileOverride='', toPrint=False, onlyDefault=False)

    # Define the regional directory in relation to the avalanche directory
    regionalDirFromCfg = str(cfg['GENERAL']['regionalDir'])
    regionalDir = pathlib.Path(avalancheDir) / regionalDirFromCfg

    # List valid avalanche directories within the regional directory
    avaDirs = com7.findAvaDirs(regionalDir)

    # Get number of processes to perform #ToDo: somehow get nVariations as well so we can display total amount of sims
    nProcesses = cfgUtils.getNumberOfProcesses(cfgMain, len(avaDirs))

    # Set nCPU for com1 to 1 to avoid dual parallelization, i.e. each subAvaDir variation is
    # processed sequentially. Preliminary solution for now, may be optimized.
    cfgMain['MAIN']['nCPU'] = '1'

    # Process each avalanche directory within the regional folder in parallel
    with ProcessPoolExecutor(max_workers=nProcesses) as executor:
        # Submit each avalanche directory to the executor
        futures = {executor.submit(com7.processAvaDirCom1Regional, cfgMain, cfg, avaDir):
                       avaDir for avaDir in avaDirs}
        # Log results as each future completes
        nSuccesses = 0
        for future in as_completed(futures):
            avaDir = futures[future]
            try:
                resultDir, status = future.result()
                log.info(f"{status} in directory: {resultDir.relative_to(pathlib.Path(regionalDir))}"
                         f" at {time.time() - startTime:.1f} s")
                if status == "Success":
                    nSuccesses += 1
            except Exception as e:
                log.error(f"Error processing {avaDir}: {e}")
    log.info(f"Processing complete. Success in '{nSuccesses}' out of '{len(avaDirs)}' directories.")

    # Copy (or move) files from the 'Outputs/com1DFA/peakFiles' directory from each of the subdirectories (avaDirs) to a directory
    # called 'allPeakFiles' in the main regional directory. Also copy (or move) 'timeSteps' to allPeakFiles/allTimeSteps
    if cfg['GENERAL'].getboolean('copyPeakFiles'):
        allPeakFilesDir, allTimeStepsDir = com7.moveOrCopyPeakFiles(cfg, regionalDir, avaDirs)

        moved = cfg['GENERAL'].getboolean('moveInsteadOfCopy')
        log.info(f"{'Moved' if moved else 'Copied'} peakFiles to "
                  f"{allPeakFilesDir}")
        log.info(f"{'Moved' if moved else 'Copied'} timeSteps to "
                  f"{allTimeStepsDir}")

    # Merge output rasters from all avalanche directories
    mergedRastersDir = com7.mergeOutputRasters(cfg, regionalDir)
    if mergedRastersDir:
        log.info(f"Merged output rasters saved to {mergedRastersDir}")

    # Print time needed
    endTime = time.time()
    log.info('Regional process finished after %.1f seconds.' % (endTime - startTime))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run regional workflow')
    parser.add_argument('avadir', metavar='avalancheDir', type=str, nargs='?', default='',
                        help='the avalanche directory')

    args = parser.parse_args()
    runCom7Regional(str(args.avadir))
