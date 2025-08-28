"""runScript for splitting avalanche input data into multiple folders based on individual release areas"""

import time
import pathlib
import argparse

from avaframe.in3Utils import cfgUtils, logUtils
from avaframe.com7Regional import splitInputs

def runSplitInputs(avalancheDir=''):
    """Main function to split input data for avalanche scenarios.

    Parameters
    ----------
    avalancheDir : str
        Path to the avalanche directory. If not provided, the function will
        retrieve it from the general configuration file.
    
    Returns
    -------
    None

    """
    # Time the whole routine
    startTime = time.time()

    # Load the avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName='runSplitInputs')
    log.info('MAIN SCRIPT')

    # Define the output directory
    outputDir = pathlib.Path(avalancheDir) / 'com7Regional'

    # Load module configuration
    cfg = cfgUtils.getModuleConfig(splitInputs)

    # Run splitting process
    splitInputs.splitInputsMain(pathlib.Path(avalancheDir), outputDir, cfg, cfgMain)

    # Print time needed
    endTime = time.time()
    log.info(f"Completed splitting input data after "f"{endTime - startTime:.1f} seconds.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run split inputs for avalanche directories")
    parser.add_argument('avalancheDir', type=str, nargs='?', default='',
                        help="Directory containing the main avalanche data")

    args = parser.parse_args()
    runSplitInputs(args.avalancheDir)
