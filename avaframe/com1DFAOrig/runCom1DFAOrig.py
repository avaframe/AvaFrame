"""
    Run script for running Standalone DFA
"""

# Load modules
import os
import time

# Local imports
from avaframe.com1DFAOrig import com1DFAOrig
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'runCom1DFA'
modName = 'com1DFAOrig'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.debug('MAIN SCRIPT')
log.debug('Current avalanche: %s', avalancheDir)

# Load input parameters from configuration file
# write config to log file
cfg = cfgUtils.getModuleConfig(com1DFAOrig)

startTime = time.time()

# Clean input directory(ies) of module files
initProj.cleanModuleFiles(avalancheDir, com1DFAOrig, modName)

# Run Standalone DFA
reportDictList = com1DFAOrig.com1DFAOrigMain(cfg, avalancheDir)

# Print time needed
endTime = time.time()
log.info(('Took %s seconds to calculate.' % (endTime - startTime)))

# append parameters from logFile
for reportD in reportDictList:
    simName = reportD['simName']['name']
    parameterDict, reportD = fU.extractParameterInfo(avalancheDir, simName, reportD)


# Generate plots for all peakFiles
plotDict = oP.plotAllPeakFields(avalancheDir, cfgMain['FLAGS'], modName, demData='')

# Set directory for report
reportDir = os.path.join(avalancheDir, 'Outputs', 'com1DFAOrig', 'reports')
# write report
gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'].getboolean('reportOneFile'), plotDict)

# write configuration to file
cfgUtils.writeCfgFile(avalancheDir, com1DFAOrig, cfg)
