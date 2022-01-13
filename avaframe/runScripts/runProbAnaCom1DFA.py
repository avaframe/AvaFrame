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


# log file name; leave empty to use default runLog.log
logName = 'runCom1DFAandProbAna'

# Load general configuration filee
cfgMain = cfgUtils.getGeneralConfig()

# Read avalanche directory
avaDir = cfgMain['MAIN']['avalancheDir']
avaDir = pathlib.Path(avaDir)
avaName = avaDir.name

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avaDir, keep=logName)

# Load configuration file for probabilistic run and analysis
cfgProb = cfgUtils.getModuleConfig(probAna)

# create configuration files for com1DFA simulations including parameter variation - defined in the probabilistic config
cfgFiles = probAna.createComModConfig(cfgProb, avaDir, com1DFA, cfgFileMod='')

# -------USER INPUT ------------
# probability configurations
probabilityConfigurations = {'testAll': {}, 'testRelTh': {'mu': '0.15500',
    'simTypeActual': 'null'}, 'testMu': {'relTh': '1.', 'simTypeActual': 'null'}}
# -----------------------------

# Start logging
log = logUtils.initiateLogger(avaDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avaDir)

# perform com1DFA simulations
for varPar in cfgFiles:
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avaDir, cfgMain,
    cfgFile=cfgFiles[varPar]['cfgFile'], relThField='')

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
    sP.plotProbMap(avaDir, inputDir, cfgProb)

    # copy outputs to folder called like probability configurations
    outputFiles = avaDir / 'Outputs' / 'ana4Stats'
    saveFiles = avaDir / 'Outputs' / ('ana4Stats_' + probConf)
    shutil.copytree(outputFiles, saveFiles)
