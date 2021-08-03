"""
    Run script for running python DFA kernel
"""

import pathlib

# Local imports
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in3Utils.initializeProject as initProj
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in1Data import getInput as gI
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFA.deriveParameterSet as dP


def runCom1DFA(avaDir='', cfgFile='', relThField='', variationDict=''):

    """ run com1DFA module """

    # +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
    # log file name; leave empty to use default runLog.log
    logName = 'runCom1DFA'

    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avaDir != '':
        avalancheDir = avaDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # set module name, reqiured as long we are in dev phase
    # - because need to create e.g. Output folder for com1DFA to distinguish from
    # current com1DFA
    modName = 'com1DFA'

    # Clean input directory(ies) of old work and output files
    # initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)
    initProj.cleanModuleFiles(avalancheDir, com1DFA, modName)

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    # Create output and work directories
    # - because need to create e.g. Output folder for com1DFA to distinguish from
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
    inputSimFiles = gI.getInputDataCom1DFA(avalancheDir, modCfg['FLAGS'])

    # write full configuration file to file
    cfgUtils.writeCfgFile(avalancheDir, com1DFA, modCfg, fileName='sourceConfiguration')

    # create a list of simulations
    # if need to reproduce exactely the hash - need to be strings with exactely the same number of digits!!
    simDict = com1DFA.prepareVarSimDict(modCfg, inputSimFiles, variationDict)

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
        cfgUtils.writeCfgFile(avalancheDir, com1DFA, cfg, fileName=cuSim)

        # log simulation name
        log.info('Run simulation: %s' % cuSim)

        # set release area scenario
        inputSimFiles['releaseScenario'] = simDict[cuSim]['relFile']

        # +++++++++++++++++++++++++++++++++
        # ------------------------
        particlesList, fieldsList, tSave, dem, reportDict, cfgFinal = com1DFA.com1DFAMain(cfg, avalancheDir, cuSim, inputSimFiles, outDir, relThField)

        # +++++++++EXPORT RESULTS AND PLOTS++++++++++++++++++++++++

        reportDictList.append(reportDict)

        # export for visulation
        if cfg['VISUALISATION'].getboolean('writePartToCSV'):
            outDir = pathlib.Path(avalancheDir, 'Outputs', modName)
            com1DFA.savePartToCsv(cfg['VISUALISATION']['particleProperties'], particlesList, outDir)

        # create hash to check if config didnt change
        simHashFinal = cfgUtils.cfgHash(cfgFinal)
        if simHashFinal != simHash:
            log.warning('simulation configuration has been changed since start')
            cfgUtils.writeCfgFile(avalancheDir, com1DFA, cfg, fileName='%s_butModified' % simHash)

    # Set directory for report
    reportDir = pathlib.Path(avalancheDir, 'Outputs', modName, 'reports')
    # Generate plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avalancheDir, cfgMain['FLAGS'], modName, demData=dem)
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    # read all simulation configuration files and return dataFrame and write to csv
    standardCfg = cfgUtils.getDefaultModuleConfig(com1DFA, toPrint=False)
    simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg=standardCfg, writeCSV=True)

    return particlesList, fieldsList, tSave, dem, plotDict, reportDictList


if __name__ == "__main__":
    runCom1DFA()
