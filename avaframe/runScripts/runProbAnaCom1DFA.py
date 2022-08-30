"""
    Run script for performing an avalanche simulation with parameter variation and performing a
    probability analysis with the simulation results
    Define settings in ana4Stats/probAnaCfg.ini or your local copy - local_probAnaCfg.ini
"""

# Load modules
import time
import glob
import pathlib
import shutil
import numpy as np

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.out3Plot import statsPlots as sP
from avaframe.ana4Stats import probAna
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.fileHandlerUtils as fU


# -------USER INPUT ------------
# probability configurations are used to filter simulations for creating probability maps:
# set for mu and relTh your reference values (from your local_com1DFACfg.ini file)
# also other parameters can be used for filtering, here for example only null simulations are taken into account
# for testRelTh and testMu
probabilityConfigurations = {'testAll': {}, 'testRelTh': {'scenario': 'relTh',
    'simTypeActual': 'null'}, 'testMu': {'scenario': 'mu', 'simTypeActual': 'null'}}
# -----------------------------


# log file name; leave empty to use default runLog.log
logName = 'runCom1DFAandProbAna'

# Load general configuration filee
cfgMain = cfgUtils.getGeneralConfig()

# Read avalanche directory
avaDir = cfgMain['MAIN']['avalancheDir']
avaDir = pathlib.Path(avaDir)
avaName = avaDir.name

# Start logging
log = logUtils.initiateLogger(avaDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avaDir)

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avaDir, keep=logName)

# Load configuration file for probabilistic run and analysis
cfgProb = cfgUtils.getModuleConfig(probAna)

# create configuration files for com1DFA simulations including parameter variation - defined in the probabilistic config
cfgFiles = probAna.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod='')

# perform com1DFA simulations
for varPar in cfgFiles:
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avaDir, cfgMain,
    cfgFile=cfgFiles[varPar]['cfgFile'])

    # Clean input directory(ies) of old work files but keep outputs
    initProj.cleanSingleAvaDir(avaDir, keep=logName, deleteOutput=False)

# perform pobability analysis
for probConf in probabilityConfigurations:

    # filter simulations according to probability configurations
    cfgProb['FILTER'] = probabilityConfigurations[probConf]
    log.info('Perform proba analysis for configuration: %s' % probConf)
    # provide optional filter criteria for simulations
    parametersDict = fU.getFilterDict(cfgProb, 'FILTER')

    # perform probability analysis
    analysisPerformed = probAna.probAnalysis(avaDir, cfgProb, com1DFA, parametersDict=parametersDict)
    if analysisPerformed is False:
        log.warning('No files found for configuration: %s' % probConf)
    # make a plot of the map
    inputDir = pathlib.Path(avaDir, 'Outputs', 'ana4Stats')
    sP.plotProbMap(avaDir, inputDir, cfgProb, demPlot=True)

    # copy outputs to folder called like probability configurations
    outputFiles = avaDir / 'Outputs' / 'ana4Stats'
    saveFiles = avaDir / 'Outputs' / ('ana4Stats_' + probConf)
    shutil.move(outputFiles, saveFiles)
