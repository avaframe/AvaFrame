"""
    Run script for running the variations tests
"""

# Load modules
import time
import pathlib

# Local imports
from avaframe.com1DFAOrig import com1DFAOrig
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
from benchmarks import simParametersVar

# +++++++++REQUIRED+++++++++++++
# Which result types for comparison plots
outputVariable = ['ppr', 'pfd', 'pfv']
# aimec settings
aimecResType = 'ppr'
aimecThresholdValue = '1'
aimecDiffLim = '5'
aimecContourLevels = '1|3|5|10'
aimecFlagMass = 'False'
aimecComModules = 'benchmarkReference|com1DFAOrig'
# ++++++++++++++++++++++++++++++

# log file name; leave empty to use default runLog.log
logName = 'runVariationsTestsCom1DFAOrig'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False)

# filter benchmarks for a tag
type = 'TAGS'
valuesList = ['varParTest']
testList = tU.filterBenchmarks(testDictList, type, valuesList, condition='and')

# Set directory for full standard test report
outDir = pathlib.Path.cwd() / 'tests' / 'reportsVariationsCom1DFAOrig'
fU.makeADir(outDir)

# Start writing markdown style report for standard tests
reportFile = outDir / 'variationTestsReportCom1DFAOrig.md'
with open(reportFile, 'w') as pfile:
    # Write header
    pfile.write('# Variation Tests Report \n')
    pfile.write('## Compare com1DFAOrig simulations to benchmark results \n')

log = logUtils.initiateLogger(outDir, logName)
log.info('The following benchmark tests will be fetched ')
for test in testList:
    log.info('%s' % test['NAME'])

for test in testList:

    # Start logging
    avaDir = test['AVADIR']

    # Fetch benchmark test info
    benchDict = simParametersVar.fetchBenchParameters(test['NAME'])
    simNameRef = test['simNameRef']
    simNameRefTest = simNameRef.replace('ref', 'dfa')
    refDir = pathlib.Path('..', 'benchmarks', test['NAME'])
    simType = benchDict['simType']
    rel = benchDict['Simulation Parameters']['Release Area Scenario']

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # get path to executable
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFAOrig)
    com1Exe = cfgCom1DFA['GENERAL']['com1Exe']

    # Load input parameters from configuration file for standard tests
    # write config to log file
    avaName = pathlib.Path(avaDir).name
    standardCfg = refDir / ('%sVarPar_com1DFAOrigCfg.ini' % test['AVANAME'])
    cfg = cfgUtils.getModuleConfig(com1DFAOrig, standardCfg)
    cfg['GENERAL']['com1Exe'] = com1Exe

    # Set timing
    startTime = time.time()
    # Run Standalone DFA
    reportDictList = com1DFAOrig.com1DFAOrigMain(cfg, avaDir)

    modName = 'com1DFAOrig'

    # Print time needed
    endTime = time.time()
    timeNeeded = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeeded)))

    # Generata plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avaDir, cfgMain['FLAGS'], modName)

    # Set directory for report
    reportDir = pathlib.Path(avaDir, 'Outputs', modName, 'reports')
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    # Fetch correct reportDict according to flagEntRes
    reportD = {}
    for dict in reportDictList:
        if simNameRefTest == dict['simName']['name']:
            reportD = dict

    simNameComp = reportD['simName']['name']

    # set result files directory
    compDir = pathlib.Path(avaDir, 'Outputs', modName, 'peakFiles')

    # Add info on run time
    reportD['runTime'] = timeNeeded

    # +++++++Aimec analysis
    # load configuration
    aimecCfg = refDir / ('%s_AIMECCfg.ini' % test['AVANAME'])
    cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC, aimecCfg)
    cfgAimec['AIMECSETUP']['resType'] = aimecResType
    cfgAimec['AIMECSETUP']['thresholdValue'] = aimecThresholdValue
    cfgAimec['AIMECSETUP']['diffLim'] = aimecDiffLim
    cfgAimec['AIMECSETUP']['contourLevels'] = aimecContourLevels
    cfgAimec['FLAGS']['flagMass'] = aimecFlagMass
    cfgAimec['AIMECSETUP']['comModules'] = aimecComModules
    cfgAimec['AIMECSETUP']['testName'] = test['NAME']

    # Setup input from com1DFA and reference
    pathDict = {}
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

    plotListRep = {}
    reportD['Simulation Difference'] = {}
    reportD['Simulation Stats'] = {}

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