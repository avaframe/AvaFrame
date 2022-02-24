"""
Run a the energy line test
Compare the DFA simulation result to the energy solution
"""
import time

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.ana1Tests import energyLineTest


# Time the whole routine
startTime = time.time()

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
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=True)

energyLineTest.mainEnergyLineTest(cfgMain)
