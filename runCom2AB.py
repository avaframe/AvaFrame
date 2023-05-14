"""
    Run script for module com2AB
"""

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3Plot import outAB
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runCom2AB'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
# write config to log file
cfg = cfgUtils.getModuleConfig(com2AB)

# Calculate ALPHABETA
pathDict, dem, splitPoint, eqParams, resAB = com2AB.com2ABMain(cfg, avalancheDir)
abShpFile = outAB.writeABtoSHP(pathDict, resAB)

# Analyse/ plot/ write results #
reportDictList = []
_, plotFile, writeFile = outAB.writeABpostOut(pathDict, dem, splitPoint, eqParams, resAB, cfgMain, reportDictList)

log.info('Plotted to: %s' % [str(plotFileName) for plotFileName in plotFile])
log.info('Data written: %s' % [str(writeFileName) for writeFileName in writeFile])
