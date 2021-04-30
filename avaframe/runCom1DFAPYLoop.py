"""
    Run script for running the tests
"""

# Load modules
import os
import time
import numpy as np

# Local imports
from avaframe.com1DFAPy import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'runCom1DFAPyLoop'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# name of avalanche
avaName = 'avaParabola'
avaDir = 'data' + os.sep + avaName

# release thickness values
releaseThickness = np.linspace(0.5, 1.5, 5)

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avaDir,  keep=logName)

# run Standard Tests sequentially
for relTh in releaseThickness:

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('Current avalanche: %s', avaDir)

    # Load input parameters from configuration file
    # write config to log file
    standardCfg = os.path.join(avaDir, 'Inputs', '%s_com1DFACfgPy.ini' % avaName)
    cfg = cfgUtils.getModuleConfig(com1DFA, standardCfg)

    # set release thickness
    cfg['GENERAL']['relTh'] = str(relTh)
    modName = 'com1DFAPy'

    # RUN com1DFAPy
    # Set timing
    startTime = time.time()
    # Run Standalone DFA
    Particles, Fields, Tsave, dem, reportDictList = com1DFA.com1DFAMain(cfg, avaDir, relThField='')

    # +++++++++EXPORT RESULTS AND PLOTS++++++++++++++++++++++++
    # Generate plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avaDir, cfg, cfgMain['FLAGS'], modName)

    # Print time needed
    endTime = time.time()
    timeNeeded = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeeded)))

    # Set directory for report
    reportDir = os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'reports')

    cfgMain['FLAGS']['reportOneFile'] = 'False'
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)
