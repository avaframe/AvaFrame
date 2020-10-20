"""
    Run script for running the operational workflow
    This file is part of Avaframe.
"""

# Load modules
import os
import time
import logging

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.com2AB import com2AB
from avaframe.log2Report import generateReport as gR
from avaframe.out3SimpPlot import outPlotAllPeak as oP
from avaframe.out3Plot import outAB
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import initializeProject as initProj

# Time the whole shebang
startTime = time.time()

# log file name; leave empty to use default runLog.log
logName = 'runOperational'

# Load avalanche directory from general configuration file
# v0.1 -> one avalanche is fine,
# v0.2 -> if given dir contains more than one avalancheDir -> run loop
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# ----------------
# Load input parameters from configuration files

# ----------------
# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

# ----------------
# Run dense flow
# cfg = cfgUtils.getModuleConfig(com1DFA)
# reportDictList = com1DFA.com1DFAMain(cfg, avalancheDir)

# ----------------
# Run Alpha Beta
cfgAB = cfgUtils.getModuleConfig(com2AB)
resAB = com2AB.com2ABMain(cfgAB, avalancheDir)

# ----------------
# Collect results/plots/report  to a single directory
# make simple plots (com1DFA, com2AB)
# peak file plot

# Generata plots for all peakFiles
reportDictList = []
plotDict = ''
# plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'])
reportDictList, _, _ = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

# Set directory for report
reportDir = os.path.join(avalancheDir, 'Outputs')
# write report
gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

# Print time needed
endTime = time.time()
log.info('Took %s seconds to calculate.' % (endTime - startTime))
