"""
    Run script for running the standard tests
"""

# Load modules
import os
import time
import pathlib
import glob

# Local imports
from avaframe import runCom1DFA
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
import avaframe.com1DFA.com1DFA as com1DFA


# log file name; leave empty to use default runLog.log
logName = 'runStandardTests'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# set desired filter criteria to filter individual simulations
parametersDict = {'simTypeActual': [], 'releaseScenario': []}

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False)

# filter benchmarks for a tag
type = 'TAGS'
valuesList = ['standardTest']
testList = tU.filterBenchmarks(testDictList, type, valuesList, condition='and')

# Set directory for full standard test report
outDir = os.path.join(os.getcwd(), 'tests', 'reportsCom1DFA')
fU.makeADir(outDir)

# Start writing markdown style report for standard tests
reportFile = os.path.join(outDir, 'standardTestsReportPy.md')
with open(reportFile, 'w') as pfile:
    # Write header
    pfile.write('# Standard Tests Report \n')
    pfile.write('## Compare com1DFA simulations to benchmark results \n')

log = logUtils.initiateLogger(outDir, logName)
log.info('The following benchmark tests will be fetched ')
for test in testList:
    log.info('%s' % test['NAME'])

for test in testList:

    # Start logging
    avaDir = test['AVADIR']

    # Fetch benchmark test info
    benchDict = simParametersDict.fetchBenchParameters(test['NAME'])
    simNameRef = test['simNameRef']
    refDir = pathlib.Path('..', 'benchmarks', test['NAME'])

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # Load input parameters from configuration file for standard tests
    # write config to log file
    avaName = os.path.basename(avaDir)
    standardCfg = os.path.join('..', 'benchmarks', test['NAME'], '%s_com1DFACfg.ini' % test['AVANAME'])

    particlesList, fieldsList, Tsave, dem, plotDict, reportDictList = runCom1DFA.runCom1DFA(avaDir=avaDir, cfgFile=standardCfg, relThField='', variationDict='')

    modName = 'com1DFA'

    # Fetch correct reportDict according to parametersDict
    simNameComp = cfgUtils.filterSims(avaDir, parametersDict)
    if len(simNameComp) > 1:
        log.error('more than one matching simulation found for criteria! ')
    else:
        simNameComp = simNameComp[0]

    for dict in reportDictList:
        if simNameComp == dict['simName']['name']:
            reportD = dict

    # set result files directory
    compDir = pathlib.Path(avaDir, 'Outputs', modName, 'peakFiles')

    # +++++++Aimec analysis
    # load configuration
    aimecCfg = os.path.join('..', 'benchmarks', test['NAME'], '%s_AIMECCfg.ini' % test['AVANAME'])
    cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC, aimecCfg)
    cfgAimec['AIMECSETUP']['resType'] = 'ppr'
    cfgAimec['AIMECSETUP']['thresholdValue'] = '1'
    cfgAimec['AIMECSETUP']['diffLim'] = '5'
    cfgAimec['AIMECSETUP']['contourLevels'] = '1|3|5|10'
    cfgAimec['FLAGS']['flagMass'] = 'False'
    cfgAimec['AIMECSETUP']['comModules'] = 'benchmarkReference|com1DFA'
    cfgAimec['AIMECSETUP']['testName'] = test['NAME']

    # Setup input from com1DFA and reference
    pathDict = []
    pathDict = dfa2Aimec.dfaBench2Aimec(avaDir, cfgAimec, simNameRef, simNameComp)
    pathDict['numSim'] = len(pathDict['ppr'])
    log.info('reference file comes from: %s' % pathDict['compType'][1])

    # Extract input file locations
    pathDict = aimecTools.readAIMECinputs(avaDir, pathDict, dirName=simNameComp)

    # perform analysis
    rasterTransfo, newRasters, resAnalysis = ana3AIMEC.AIMEC2Report(pathDict, cfgAimec)

    # add aimec results to report dictionary
    reportD, benchDict = ana3AIMEC.aimecRes2ReportDict(resAnalysis, reportD, benchDict, pathDict['referenceFile'])
    # +++++++++++Aimec analysis

    # Create plots for report
    # Load input parameters from configuration file
    cfgRep = cfgUtils.getModuleConfig(generateCompareReport)

    # REQUIRED+++++++++++++++++++
    # Which parameter to filter data, e.g. varPar = 'simType', values = ['null'] or
    # varPar = 'Mu', values = ['0.055', '0.155']; values need to be given as list, also if only one value
    outputVariable = ['ppr', 'pfd', 'pfv']
    plotListRep = {}
    reportD['Simulation Difference'] = {}
    reportD['Simulation Stats'] = {}
    # ++++++++++++++++++++++++++++

    # Plot data comparison for all output variables defined in suffix
    for var in outputVariable:
        plotDict = outQuickPlot.quickPlotBench(avaDir, simNameRef, simNameComp, refDir, compDir, cfgMain, var)
        for plot in plotDict['plots']:
            plotListRep.update({var: plot})
            reportD['Simulation Difference'].update({var: plotDict['difference']})
            reportD['Simulation Stats'].update({var: plotDict['stats']})

    # copy files to report directory
    plotPaths = generateCompareReport.copyQuickPlots(avaName, test['NAME'], outDir, plotListRep)
    aimecPlots = [resAnalysis['slCompPlot'], resAnalysis['areasPlot']]
    plotPaths = generateCompareReport.copyAimecPlots(aimecPlots, test['NAME'], outDir, plotPaths)

    # add plot info to general report Dict
    reportD['Simulation Results'] = plotPaths

    # write report
    generateCompareReport.writeCompareReport(reportFile, reportD, benchDict, avaName, cfgRep)
