"""
    Run script for running Standalone DFA
    This file is part of Avaframe.
"""

# Load modules
import os
import time

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'runCom1DFA'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load input parameters from configuration file
# write config to log file
cfg = cfgUtils.getModuleConfig(com1DFA)

startTime = time.time()

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

# Run Standalone DFA
reportDictList = com1DFA.com1DFAMain(cfg, avalancheDir)

# Print time needed
endTime = time.time()
log.info(('Took %s seconds to calculate.' % (endTime - startTime)))

# Generata plots for all peakFiles
plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'])

# Set directory for report
reportDir = os.path.join(avalancheDir, 'Outputs', 'com1DFA', 'reports')
# write report
gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)
