"""
    Run script for running python DFA kernel
"""

import time
import copy
import os

# Local imports
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in3Utils.initializeProject as initProj
import avaframe.com1DFAPy.com1DFA as com1DFA
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def runCom1DFAPy(avaDir='', cfgFile='', relThField=''):
    """ run com1DFAPy module """

    # +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
    # log file name; leave empty to use default runLog.log
    logName = 'runCom1DFAPy'

    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avaDir != '':
        avalancheDir = avaDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']
    # set module name, reqiured as long we are in dev phase
    # - because need to create e.g. Output folder for com1DFAPy to distinguish from
    # current com1DFA
    modName = 'com1DFAPy'

    # Clean input directory(ies) of old work and output files
    # initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)
    initProj.cleanModuleFiles(avalancheDir, com1DFA, modName)

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    # Load configuration
    if cfgFile != '':
        cfg = cfgUtils.getModuleConfig(com1DFA, cfgFile)
    else:
        cfg = cfgUtils.getModuleConfig(com1DFA)
    cfgGen = cfg['GENERAL']
    cfgGen['avalancheDir'] = avalancheDir

    # +++++++++++++++++++++++++++++++++
    # ------------------------
    Particles, Fields, Tsave, dem, reportDictList = com1DFA.com1DFAMain(cfg, avalancheDir, relThField)

    # +++++++++EXPORT RESULTS AND PLOTS++++++++++++++++++++++++
    # Generate plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'], modName)

    # Set directory for report
    reportDir = os.path.join(avalancheDir, 'Outputs', 'com1DFAPy', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    return Particles, Fields, Tsave, dem, plotDict, reportDictList


if __name__ == "__main__":
    runCom1DFAPy()
