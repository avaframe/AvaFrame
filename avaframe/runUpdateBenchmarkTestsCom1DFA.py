"""
    Run script for running the standard tests with com1DFA
    in this test all the available tests tagged standardTest are performed
"""

# Load modules
import time
import pathlib
from configupdater import ConfigUpdater

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.ana1Tests import testUtilities as tU
from avaframe.log2Report import generateReport as gR
from avaframe.log2Report import generateCompareReport
from avaframe.ana3AIMEC import ana3AIMEC, dfa2Aimec, aimecTools
from avaframe.out3Plot import outQuickPlot
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from benchmarks import simParametersDict


# +++++++++REQUIRED+++++++++++++
# Which result types for comparison plots
outputVariable = ['ppr', 'pfd', 'pfv']
# ++++++++++++++++++++++++++++++

# log file name; leave empty to use default runLog.log
logName = 'runUpdateBenchmarkTestsCom1DFA'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False)

# filter benchmarks to extract only desiered ones
type = 'TAGS'
valuesList = ['standardTest', 'test']
testList = tU.filterBenchmarks(testDictList, type, valuesList, condition='and')

# Set directory for full standard test report
outDir = pathlib.Path.cwd() / 'tests' / 'reportsCom1DFA'
fU.makeADir(outDir)

# Start writing markdown style report for standard tests
reportFile = outDir / 'updateBenchmarkTestsCom1DFA.md'
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Updated benchmarks Report \n')

log = logUtils.initiateLogger(outDir, logName)
log.info('The following benchmark tests will be updated ')
for test in testList:
    log.info('%s' % test['NAME'])

# run Standard Tests sequentially
for test in testList:

    avaDir = test['AVADIR']

    # Fetch benchmark test info
    # toDo: this needs to be changed. We want to read this from the json file
    benchDict = simParametersDict.fetchBenchParameters(test['NAME'])
    simNameRef = test['simNameRef']
    refDir = pathlib.Path('..', 'benchmarks', test['NAME'])
    simType = benchDict['simType']
    rel = benchDict['Simulation Parameters']['Release Area Scenario']

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # Load input parameters from configuration file for standard tests
    avaName = pathlib.Path(avaDir).name
    standardCfg = refDir / ('%s_com1DFACfg.ini' % test['AVANAME'])
    modName = 'com1DFA'

    # change the ini file to force particle initialization
    updater = ConfigUpdater()
    updater.read(standardCfg)
    updater['GENERAL']['initialiseParticlesFromFile'].value = False

    updater.update_file()
    # Set timing
    startTime = time.time()
    # call com1DFA run
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avaDir, cfgMain, cfgFile=standardCfg)
    endTime = time.time()
    timeNeeded = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeeded)))

    # change the ini file to read particles from file
    updater = ConfigUpdater()
    updater.read(standardCfg)
    updater['GENERAL']['initialiseParticlesFromFile'].value = True
    updater['GENERAL']['particleFile'].value = '../benchmarks/' + test['NAME']

    print(reportDictList[0])
    print(test)
    rep = reportDictList[0]
    test['simName'] = rep['simName']
    test['Simulation Parameters'] = rep['Simulation Parameters']
    test['SimNameRef'] = rep['Simulation Parameters']
    tU.writeDesDicttoJson(test, 'testing', refDir)
    # # Fetch correct reportDict according to simType and release area scenario
    # # read all simulation configuration files and return dataFrame and write to csv
    # parametersDict = {'simTypeActual': simType, 'releaseScenario': rel}
    # simNameComp = cfgUtils.filterSims(avaDir, parametersDict)
    # if len(simNameComp) > 1:
    #     log.error('more than one matching simulation found for criteria! ')
    # else:
    #     simNameComp = simNameComp[0]

    # find report dictionary corresponding to comparison simulation
    # reportD = {}
    # for dict in reportDictList:
    #     if simNameComp in dict['simName']['name']:
    #         reportD = dict
    # if reportD == {}:
    #     message = 'No matching simulation found for reference simulation: %s' % simNameRef
    #     log.error(message)
    #     raise ValueError(message)
    # log.info('Reference simulation %s and comparison simulation %s ' % (simNameRef, simNameComp))
    #
    # # set result files directory
    # compDir = pathlib.Path(avaDir, 'Outputs', modName, 'peakFiles')
    #
    # # Add info on run time
    # reportD['runTime'] = timeNeeded
    #
    # # add plot info to general report Dict
    # reportD['Simulation Results'] = plotPaths
    #
    # # write report
    # generateCompareReport.writeCompareReport(reportFile, reportD, benchDict, avaName, cfgRep)
