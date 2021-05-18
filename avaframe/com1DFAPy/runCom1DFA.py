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

    # Create output and work directories
    # - because need to create e.g. Output folder for com1DFAPy to distinguish from
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName)

    # generate list of simulations from desired configuration
    if variationDict == '':
        # Load full configuration
        modCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=cfgFile, modInfo=True)
        variationDict = dP.getVariationDict(avalancheDir, modCfg, modInfo)
    else:
        # check if variationDict items exist and are provided in correct format
        # Load standard/ default configuration
        modCfg = cfgUtils.getDefaultModuleConfig(com1DFA)
        variationDict = dP.validateVarDict(variationDict, modCfg)
        log.info('Variations are performed for:')
        for key in variationDict:
            log.info('%s: %s' % (key, variationDict[key]))

    # add avalanche directory info to cfg
    modCfg['GENERAL']['avalancheDir'] = avalancheDir

    # fetch input data - dem, release-, entrainment- and resistance areas
    inputSimFiles = gI.getInputDataCom1DFAPy(avalancheDir, modCfg['FLAGS'])

    # write full configuration file to file
    cfgUtils.writeCfgFile(avalancheDir, com1DFAPy, modCfg, fileName='overallConfiguration')

    # create a list of simulations
    # if need to reproduce exactely the hash - need to be strings with exactely the same number of digits!!
    simDict = com1DFA.prepareVarSimDict(modCfg, inputSimFiles, variationDict['GENERAL'])

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

        # export for visulation
        if cfg['VISUALISATION'].getboolean('writePartToCSV'):
            outDir = os.path.join(avalancheDir, 'Outputs', modName)
            com1DFA.savePartToCsv(cfg['VISUALISATION']['particleProperties'], particlesList, outDir)

        # create hash to check if config didnt change
        simHashFinal = cfgUtils.cfgHash(cfgFinal)
        if simHashFinal != simHash:
            log.warning('simulation configuration has been changed since start')
            cfgUtils.writeCfgFile(avalancheDir, com1DFAPy, cfg, fileName='%s_butModified' % simHash)

    # Set directory for report
    reportDir = os.path.join(avalancheDir, 'Outputs', 'com1DFAPy', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    # read all simulation configuration files and return dataFrame and write to csv
    standardCfg = cfgUtils.getDefaultModuleConfig(com1DFA)
    simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg, writeCSV=True)

    return particlesList, fieldsList, Tsave, dem, plotDict, reportDictList


if __name__ == "__main__":
    runCom1DFAPy()
