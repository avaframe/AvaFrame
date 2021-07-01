"""
    Run script for running the operational workflow
"""

# Load modules
# importing general python modules 
import time
import logging
import pathlib
import pandas

# Import the avaframe modules you want to use. This list will
# grow as you use more avaframe modules. You can refer to the different
# computational modules documentation to know which imports are required

# first import the always useful tools
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import fileHandlerUtils as fU

# then depending on which computational module you want to use
from avaframe import runCom1DFA
from avaframe.com2AB import com2AB
from avaframe.out3Plot import outAB


def runOperational(avalancheDir=''):
    """ Run com1DFA and com2AB in the default configuration with only an
    avalanche directory as input

    Parameters
    ----------
    avalancheDir: str
        path to avalanche directory (setup eg. with init scipts)

    Returns
    -------
    abShpFile: str
        path to com2AB results as shapefile
    peakFilesDF: pandas dataframe
        dataframe with info about com1DFA peak file locations
    """
    # Time the whole routine
    startTime = time.time()

    # log file name; leave empty to use default runLog.log
    logName = 'runOperational'

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
    # Clean input directory(ies) of old work and output files
    # If you just created the ``avalancheDir`` this one should be clean but if you
    # already did some calculations you might want to clean it::
    initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

    # ----------------
    # Run dense flow
    _, _, _, _, _, reportDictList = runCom1DFA.runCom1DFA(avalancheDir)

    # Get peakfiles to return to QGIS
    avaDir = pathlib.Path(avalancheDir)
    inputDir = avaDir / 'Outputs' / 'com1DFA' / 'peakFiles'
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    # Run Alpha Beta
    cfgAB = cfgUtils.getModuleConfig(com2AB)
    resAB = com2AB.com2ABMain(cfgAB, avalancheDir)
    abShpFile = outAB.writeABtoSHP(resAB)

    # ----------------
    # TODO: make report and feed info back (for QGis)
    # Collect results/plots/report  to a single directory
    # make simple plots (com1DFA, com2AB)
    # peak file plot

    # # Generata plots for all peakFiles
    # plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'])
    # reportDictList = []
    # reportDictList, _, _ = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

    # # # Set directory for report
    # reportDir = os.path.join(avalancheDir, 'Outputs')
    # # # write report
    # gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    # Print time needed
    endTime = time.time()
    log.info('Took %6.1f seconds to calculate.' % (endTime - startTime))

    return str(abShpFile), peakFilesDF


if __name__ == '__main__':
    runOperational()
