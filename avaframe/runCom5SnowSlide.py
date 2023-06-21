"""
Run the snow slide tool of com1DFA
"""

import pathlib
import time
import argparse

# Local imports
# import config and init tools
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.log2Report import generateReport as gR

# import computation modules
from avaframe.com5SnowSlide import com5SnowSlide


def runSnowSlide(avalancheDir=''):
    """ Run com1DFA with snow slide parameters with only an
    avalanche directory as input

    Parameters
    ----------
    avalancheDir: str
        path to avalanche directory (setup eg. with init scipts)

    Returns
    -------
    peakFilesDF: pandas dataframe
        dataframe with info about com1DFA peak file locations
    """
    # Time the whole routine
    startTime = time.time()

    # log file name; leave empty to use default runLog.log
    logName = 'runCom5SnowSlide'

    # Load avalanche directory from general configuration file
    # More information about the configuration can be found here
    # on the Configuration page in the documentation
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    # ----------------
    # Clean input directory(ies) of old work files
    initProj.cleanSingleAvaDir(avalancheDir, deleteOutput=False)

    # load snow slide tool config
    snowSlideCfg = cfgUtils.getModuleConfig(com5SnowSlide)

    # perform com1DFA simulation with snow slide settings
    _, plotDict, reportDictList, _ = com5SnowSlide.runSnowSlide(cfgMain, snowSlideCfg)

    # Get peakfiles to return to QGIS
    avaDir = pathlib.Path(avalancheDir)
    inputDir = avaDir / 'Outputs' / 'com1DFA' / 'peakFiles'
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    # ----------------
    # Collect results/plots/report  to a single directory
    # Set directory for report
    reportDir = avaDir / 'Outputs' / 'reports'
    fU.makeADir(reportDir)
    # write report and copy plots to report dir
    gR.writeReport(reportDir, reportDictList,
                   cfgMain['FLAGS'].getboolean('reportOneFile'),
                   plotDict=plotDict,
                   standaloneReport=True)

    # Print time needed
    endTime = time.time()
    log.info('Took %6.1f seconds to calculate.' % (endTime - startTime))

    return peakFilesDF


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run snow slide workflow')
    parser.add_argument('avadir', metavar='avadir', type=str, nargs='+',
                        help='the avalanche directory')

    args = parser.parse_args()
    runSnowSlide(str(args.avadir[0]))
