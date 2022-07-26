"""
Run a the energy line test
Compare the DFA simulation result to the energy solution
"""
import pathlib

# Local imports
# import config and init tools
from avaframe.in3Utils import initializeProject as iP
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# import computation modules
from avaframe.com1DFA import com1DFA

# import analysis modules
from avaframe.ana1Tests import energyLineTest


# +++++++++REQUIRED+++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runEnergyLineTest'
# do you want to run the DFA module (all results in the Outputs/com1DFA forlder will be deleted)
runDFAModule = False
# ++++++++++++++++++++++++++++++

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# ----------------
# Clean input directory(ies) of old work and output files
iP.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

# get path to com1DFA configuration file used for the energy line test
energyLineTestCfgFile = pathlib.Path('ana1Tests', 'energyLineTest_com1DFACfg.ini')
energyLineTestCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=energyLineTestCfgFile, modInfo=True)
# run the com1DFA module or load the results from com1DFA
dem, simDF, _ = com1DFA.runOrLoadCom1DFA(avalancheDir, cfgMain, runDFAModule=runDFAModule,
                                         cfgFile=energyLineTestCfgFile)

# generate mass averaged path profile
for simName in simDF.index:
    # make analysis and generate plots
    resultEnergyTest, savePath = energyLineTest.mainEnergyLineTest(avalancheDir, energyLineTestCfg, simName, dem)
