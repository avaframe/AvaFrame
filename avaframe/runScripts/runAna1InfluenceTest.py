"""
    Run script to run the script influenceTest.py from module ana1Test
"""

import pathlib
import configparser
config = configparser.RawConfigParser()
config.optionxform = str


# local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.com1DFA import com1DFA
from avaframe.com2AB import com2AB
import avaframe.runScripts.runAna3AIMEC as runAna3AIMEC
from avaframe.out3Plot import outAB
from avaframe.ana1Tests import influenceTest


# log file name; leave empty to use default runLog.log
logName = 'runAna1InfluenceTest'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Clean input directory of old work and output files from module

initProj.cleanSingleAvaDir(avalancheDir, deleteOutput=False)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for influence test
influenceTestCfg = pathlib.Path('ana1Tests', 'influenceTestCfg.ini')
# ABCfg = cfgUtils.getModuleConfig(com2AB)
ataCfg = pathlib.Path('ana1Tests', 'ataCfg.ini')

# Collect data from config file
config.read(influenceTestCfg)
visc2Study = config['INFLUENCETEST']['visc2Study']

# Generate an .ini file and run visc2Study simulation(s)
log.info("Generate " + visc2Study + " '.ini' file")
influenceTest.generateAutoIniFile(influenceTestCfg, visc2Study)
visc2StudyCfg = pathlib.Path('ana1Tests', visc2Study + 'Cfg.ini')
# Run simulations themselves
log.info('Run com1DFA with SAMOS viscosity')
_, _, _, _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=visc2StudyCfg, relThField='', variationDict='')

# Compute the runout reference
log.info('Compute the runout reference')
config.read(influenceTestCfg)
ref = config['ref']['ref']
if ref=='dfa':
    # Clean input directory of old work and output files from module
    initProj.cleanSingleAvaDir(avalancheDir, deleteOutput=False)
    log.info("Generate REF '.ini' file")
    influenceTest.generateAutoIniFile(influenceTestCfg, 'ref')
    refCfg = pathlib.Path('ana1Tests', 'refCfg.ini')
    log.info('Run Com1DFA with ref solution')
    _, _, _, _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=refCfg, relThField='', variationDict='')
# resAB = com2AB.com2ABMain(ABCfg, avalancheDir)
# # Analyse/ plot/ write results #
# reportDictList = []
# _, plotFile, writeFile = outAB.writeABpostOut(resAB, ABCfg, reportDictList)
# log.info('Plotted to: %s' % [str(plotFileName) for plotFileName in plotFile])
# log.info('Data written: %s' % [str(writeFileName) for writeFileName in writeFile])

# Get all simulation Data Frame
log.info('Get all simulation Data Frame')
simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False)

# Run AIMEC on all DFA simulations
log.info('Run AIMEC on all DFA simulations')
pathDict, _, _, resAIMEC = runAna3AIMEC.runAna3AIMEC(avalancheDir)

# Run influenceTest
simDF = influenceTest.mainInfluenceTest(avalancheDir, influenceTestCfg, simDF, resAIMEC, pathDict)
