"""Run script for running simulations in parallel based on input from multiple avaFolders
(and optionally copying all results to one output folder)"""

import time
import pathlib
import argparse

from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initializeProject as initProj
from avaframe.com1DFA import com1DFA

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

    # Time the whole routine
    startTime = time.time()

    # Load avalanche directory from general configuration file
    # More information about the configuration can be found here
    # on the Configuration page in the documentation
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName='runRegional')
    log.info('MAIN SCRIPT')

    # List avalanche directories that contain an 'Inputs' folder
    avaDirs = [dir for dir in pathlib.Path(avalancheDir).iterdir() if dir.is_dir() and (dir / "Inputs").is_dir()]
    log.info(f"Found a total of '{len(avaDirs)}' avalanche directories in '{avalancheDir}'.")
    log.info(f"Processing the following directories:")
    for dir in avaDirs:
        log.info(f"'{dir.name}'")

    ### iterate through each avadir and run the e.g. com1DFA process
    for dir in avaDirs:
        avalancheDir = dir
        log.info(f"PROCESSING DIRECTORY: {avalancheDir.relative_to(pathlib.Path(avalancheDir).parent)}")

        # Update the avalancheDir in the configuration for com1DFA to use
        cfgMain['MAIN']['avalancheDir'] = str(avalancheDir)

        # Clean input directory of old work and output files from module
        initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)
        # create com1DFA configuration
        cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False, onlyDefault=False)
        # Run com1DFA in the current avalanche directory
        dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgCom1DFA)

    #todo:gather the outputs and copy them to the outDir - then clean original output location? alternatively,
    #copy them directly to new OutDir/'Outputs' in the Regional folder

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