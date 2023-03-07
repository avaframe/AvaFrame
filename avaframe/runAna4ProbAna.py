"""
    Run script for performing an com1DFA avalanche simulation with parameter variation and performing a
    probability analysis with the simulation results
    Define settings in ana4Stats/probAnaCfg.ini or your local copy - local_probAnaCfg.ini
"""

# Load modules
import pathlib
import shutil
import argparse

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.out3Plot import statsPlots as sP
from avaframe.ana4Stats import probAna
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.out3Plot import outQuickPlot as oP


def runProbAna(avalancheDir=''):
    """ Run a com1DFA probability analysis with parameters and only an
    avalanche directory as input

    Parameters
    ----------
    avalancheDir: str
        path to avalanche directory (setup eg. with init scipts)

    Returns
    -------
    """

    # log file name; leave empty to use default runLog.log
    logName = 'runAna4ProbAna'

    # Load general configuration filee
    cfgMain = cfgUtils.getGeneralConfig()

    # Load avalanche directory from general configuration file
    # More information about the configuration can be found here
    # on the Configuration page in the documentation
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    avalancheDir = pathlib.Path(avalancheDir)

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    # Clean input directory(ies) of old work files
    initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

    # Load configuration file for probabilistic run and analysis
    cfgProb = cfgUtils.getModuleConfig(probAna)

    # create configuration files for com1DFA simulations including parameter variation - defined in the probabilistic config
    cfgFiles, cfgPath = probAna.createComModConfig(cfgProb, avalancheDir, com1DFA, cfgFileMod='')

    # perform com1DFA simulations
    outDir = pathlib.Path(avalancheDir, 'Outputs')
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgPath)

    # check if sampling strategy is from full sample - then only one configuration is possible
    probabilityConfigurations = probAna.fetchProbConfigs(cfgProb)

    # perform pobability analysis
    for probConf in probabilityConfigurations:

        # filter simulations according to probability configurations
        cfgProb['FILTER'] = probabilityConfigurations[probConf]
        log.info('Perform proba analysis for configuration: %s' % probConf)
        # provide optional filter criteria for simulations
        parametersDict = fU.getFilterDict(cfgProb, 'FILTER')

        # perform probability analysis
        analysisPerformed, contourDict = probAna.probAnalysis(avalancheDir, cfgProb, com1DFA, parametersDict=parametersDict)
        if analysisPerformed is False:
            log.warning('No files found for configuration: %s' % probConf)
        # make a plot of the map
        inputDir = pathlib.Path(avalancheDir, 'Outputs', 'ana4Stats')
        sP.plotProbMap(avalancheDir, inputDir, cfgProb, demPlot=True)
        # make a plot of the contours
        pathDict = {'pathResult': str(inputDir), 'avaDir': str(avalancheDir), 'plotScenario': probConf}
        oP.plotContours(contourDict, cfgProb['GENERAL']['peakVar'], cfgProb['GENERAL']['peakLim'], pathDict)

        # copy outputs to folder called like probability configurations
        outputFiles = avalancheDir / 'Outputs' / 'ana4Stats'
        saveFiles = avalancheDir / 'Outputs' / ('ana4Stats_' + probConf)
        shutil.move(outputFiles, saveFiles)

    return


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Run ana4ProbAna workflow')
    parser.add_argument('avadir', metavar='a', type=str, nargs='+',
                        help='the avalanche directory')

    args = parser.parse_args()
    runProbAna(str(args.avadir[0]))
