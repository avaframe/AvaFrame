"""Run script for running simulations in parallel based on input from multiple avaFolders
(and optionally copying all results to one output folder)"""

import time
import pathlib
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initializeProject as initProj
from avaframe.com1DFA import com1DFA

def process_avadir_com1DFA(cfgMain, avalancheDir):
    """Function to process com1DFA in each avalanche directory in parallel."""
    avalancheDir = pathlib.Path(avalancheDir)

    # Initialize log for each process
    log = logUtils.initiateLogger(avalancheDir, logName='runCom1DFA')
    log.info('PROCESS CALLED BY REGIONAL RUN')
    log.info('Current avalanche: %s', avalancheDir)

    # Update cfgMain to reflect the current avalancheDir
    cfgMain['MAIN']['avalancheDir'] = str(avalancheDir)

    # Clean input directory of old work and output files from module
    initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)

    # Create com1DFA configuration for the current avalanche directory
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False, onlyDefault=False)

    # Run com1DFA in the current avalanche directory
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgCom1DFA)
    return avalancheDir, "Success"

def runRegional(avalancheDir='', outDir=''):
    """
    Parameters
    ----------
    avalancheDir:
        input directory which should contain multiple valid avalanche directories, i.e. the regional directory!
    outDir:
        here the output directory may be set

    Returns
    -------
    outDir: str
        path to output directory
    """

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName='runRegional')
    log.info('MAIN SCRIPT')

    # Time the whole routine
    startTime = time.time()

    # Load the (REGIONAL) avalanche directory from general configuration file
    # More information about the configuration can be found here
    # on the Configuration page in the documentation
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # List avalanche directories that contain an 'Inputs' folder - logic is that they are valid avaFolders
    avaDirs = [avadir for avadir in pathlib.Path(avalancheDir).iterdir() if avadir.is_dir() and (avadir / "Inputs").is_dir()]
    log.info(f"Found a total of '{len(avaDirs)}' avalanche directories in '{avalancheDir}'.")
    log.info(f"Processing the following directories:")
    for avadir in avaDirs:
        log.info(f"'{avadir.name}'")

    # Prepare arguments for each process
    # Run processes in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor() as executor:
        # Submit each avalanche directory to the executor
        futures = {executor.submit(process_avadir_com1DFA, cfgMain, avadir): avadir for avadir in avaDirs}

        # Log results as each future completes
        for future in as_completed(futures):
            avadir = futures[future]
            try:
                result_dir, status = future.result()
                log.info(f"{status} in directory: {result_dir.relative_to(pathlib.Path(avalancheDir))}")
            except Exception as e:
                log.error(f"Error processing {avadir}: {e}")

    #todo:  gather the outputs and copy them to the outDir - then clean original output location? alternatively,
    #todo:  copy them directly to new OutDir/'Outputs' in the Regional folder, should introduce some options here

    #todo: maybe down the line we could write a small report, i.e. which avalanche took the longest to calculate, how long did it take on average, etc...

    # Print time needed
    endTime = time.time()
    log.info('Took %6.1f seconds to calculate.' % (endTime - startTime))

    return str(outDir)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run regional workflow')
    parser.add_argument('avadir', metavar='avalancheDir', type=str, nargs='?', default='',
                        help='the avalanche directory')

    args = parser.parse_args()
    runRegional(str(args.avadir))
