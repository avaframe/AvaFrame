"""
Run the rotation test
Analyze the effect of the grid direction on DFA simulation results
"""
import pathlib

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.ana1Tests import rotationTest


# log file name; leave empty to use default runLog.log
logName = 'runenergyLineTest'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# ----------------
# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

# get path to com1DFA configuration file used for the energy line test
rotationTestCfgFile = pathlib.Path('ana1Tests', 'rotationTest_com1DFACfg.ini')
rotationTest.mainRotationTest(cfgMain, rotationTestCfgFile)
