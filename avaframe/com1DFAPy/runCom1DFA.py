"""
    Run script for running python DFA kernel
"""

import os
import json

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.com1DFAPy.com1DFA as com1DFA
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFAPy.com1DFA as com1DFAPy
import avaframe.com1DFAPy.deriveParameterSet as dP


def runCom1DFAPy(avaDir='', cfgFile='', relThField='', variationDict=''):

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

    # Load standard/ default configuration
    standardCfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    # add avalanche directory info
    standardCfg['GENERAL']['avalancheDir'] = avalancheDir

    # Create output and work directories
    # set module name, reqiured as long we are in dev phase
    # - because need to create e.g. Output folder for com1DFAPy to distinguish from
    modName = 'com1DFAPy'
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName)

    # fetch input data - dem, release-, entrainment- and resistance areas
    inputSimFiles = gI.getInputDataCom1DFAPy(avalancheDir, standardCfg['FLAGS'])

    # generate list of simulations from desired configuration
    if variationDict == '':
        variationDict = dP.getVariationDict(avalancheDir, com1DFA, standardCfg, cfgFile=cfgFile)
    else:
        log.info('Variations are performed for:')
        for key in variationDict:
            log.info('%s: %s' % (key, variationDict[key]))

    # create a list of simulations
    # if need to reproduce exactely the hash - need to be strings with exately the same number of digits!!
    simDict = com1DFA.prepareVarSimDict(standardCfg, inputSimFiles, variationDict, varyAll=True)

    log.info('The following simulations will be performed')
    for key in simDict:
        log.info('Simulation: %s' % key)

    reportDictList = []
    # loop over all simulations
    for cuSim in simDict:

        # load configuration dictionary for cuSim
        cfg = simDict[cuSim]['cfgSim']

        # save configuration settings for each simulation
        simHash = simDict[cuSim]['simHash']
        cfgUtils.writeCfgFile(avalancheDir, com1DFAPy, cfg, fileName=simHash)

        # log simulation name
        log.info('Run simulation: %s' % cuSim)

        # set release area scenario
        inputSimFiles['releaseScenario'] = simDict[cuSim]['relFile']

        # +++++++++++++++++++++++++++++++++
        # ------------------------
        particlesList, fieldsList, Tsave, dem, reportDict, cfgFinal = com1DFA.com1DFAMain(cfg, avalancheDir, cuSim, inputSimFiles, outDir, relThField)

        # +++++++++EXPORT RESULTS AND PLOTS++++++++++++++++++++++++
        # Generate plots for all peakFiles
        plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'], modName)

        reportDictList.append(reportDict)

        # create hash to check if config didnt change
        simHashFinal = cfgUtils.cfgHash(cfgFinal)
        if simHashFinal != simHash:
            log.warning('simulation configuration has been changed since start')
            cfgUtils.writeCfgFile(avalancheDir, com1DFAPy, cfg, fileName='%s_butModified' % simHash)


    # Set directory for report
    reportDir = os.path.join(avalancheDir, 'Outputs', 'com1DFAPy', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)
    # export for visulation
    if cfg['VISUALISATION'].getboolean('writePartToCSV'):
        outDir = os.path.join(avalancheDir, 'Outputs', modName)
        com1DFA.savePartToCsv(cfg['VISUALISATION']['particleProperties'], particlesList, outDir)

    # read all simulation configuration files and return dataFrame and write to csv
    simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg, writeCSV=True)

    return particlesList, fieldsList, Tsave, dem, plotDict, reportDictList


if __name__ == "__main__":
    runCom1DFAPy()
