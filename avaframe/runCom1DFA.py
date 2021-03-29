"""
    Run script for running Standalone DFA
"""

# Load modules
import os
import time

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import fileHandlerUtils as fU
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
log.debug('MAIN SCRIPT')
log.debug('Current avalanche: %s', avalancheDir)

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


for reportD in reportDictList:
    simName = reportD['simName']['name']
    parameterDict = fU.extractParameterInfo(avalancheDir, simName)

    # get parameters from log
    reportD['Simulation Parameters'].update({'stop criterion': parameterDict['stop criterion']})
    reportD['Simulation Parameters'].update({'CPU time [s]': parameterDict['CPU time [s]']})

# Generata plots for all peakFiles
plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'])

# Set directory for report
reportDir = os.path.join(avalancheDir, 'Outputs', 'com1DFA', 'reports')
# write report
gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

# write configuration to file
cfgUtils.writeCfgFile(avalancheDir, com1DFA, cfg)
