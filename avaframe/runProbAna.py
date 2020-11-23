"""
    Run script for running Standalone DFA
    This file is part of Avaframe.
"""

# Load modules
import os
import time
import glob

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.ana4Prob import probAna
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import fileHandlerUtils as fU


# log file name; leave empty to use default runLog.log
logName = 'runCom1DFAandProbAna'

# Load general configuration filee
cfgMain = cfgUtils.getGeneralConfig()

# Define avalanche directories for tests
avalancheDirectories = ['data/avaHockey']

# run simulations sequentially
for avaDir in avalancheDirectories:

    # get path to executable
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
    com1Exe = cfgCom1DFA['GENERAL']['com1Exe']

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avaDir)
    #
    # Load input parameters from configuration file
    # write config to log file
    avaName = os.path.basename(avaDir)
    probSimCfg = os.path.join('..', 'benchmarks', avaName, '%sProbAna_com1DFACfg.ini' % avaName)
    cfg = cfgUtils.getModuleConfig(com1DFA, probSimCfg)
    cfg['GENERAL']['com1Exe'] = com1Exe

    startTime = time.time()

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # Run Standalone DFA
    reportDictList = com1DFA.com1DFAMain(cfg, avaDir)

    # Print time needed
    endTime = time.time()
    log.info(('Took %s seconds to calculate.' % (endTime - startTime)))

    # Generata plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avaDir, cfg, cfgMain['FLAGS'])

    # Set directory for report
    reportDir = os.path.join(avaDir, 'Outputs', 'com1DFA', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    # Load input parameters from configuration file
    cfgProb = cfgUtils.getModuleConfig(probAna)

    # perform probability analysis
    probAna.probAnalysis(avaDir, cfgProb, cfg)

    # make a plot of the map
    inputDir = os.path.join(avaDir, 'Outputs', 'ana4Prob')
    outputDir = os.path.join(avaDir, 'Outputs', 'ana4Prob', 'plots')
    oP.plotAllFields(avaDir, inputDir, outputDir, cfgProb)
