"""
Run a the energy line test
Compare the DFA simulation result to the energy solution
"""
import time

# Local imports
# import config and init tools
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.ana1Tests import ana1EnergyLineTest


# Time the whole routine
startTime = time.time()

# log file name; leave empty to use default runLog.log
logName = 'runana1EnergyLineTest'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# ----------------
# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=True)

ana1EnergyLineTest.mainEnergyLineTest(cfgMain)
