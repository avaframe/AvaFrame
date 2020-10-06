"""
    Run script for running Standalone DFA
    This file is part of Avaframe.
"""

# Load modules
import os

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.out3SimpPlot import outPlotAllPeak as oP
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import time

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
# Run Standalone DFA
reportDictList = com1DFA.runSamos(cfg, avalancheDir)

# Print time needed
endTime = time.time()
log.info(('Took %s seconds to calculate.' % (endTime - startTime)))

# Generata plots for all peakFiles
plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'], flagReport=True)

# write report
gR.writeReport(avalancheDir, reportDictList, plotDict, reportOneFile=True)
