"""
    Run script for running the standard tests
"""

# Load modules
import os
import time

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.ana1Tests import testUtilities as tU
from avaframe.log2Report import generateReport as gR
from avaframe.log2Report import generateCompareReport
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.ana3AIMEC import ana3AIMEC, dfa2Aimec
from avaframe.out3Plot import outQuickPlot
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from benchmarks import simParameters

# log file name; leave empty to use default runLog.log
logName = 'runStandardTests'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False)
# filter benchmarks for tag standardTest
type = 'TAGS'
valuesList = ['standardTest']
testList = tU.filterBenchmarks(testDictList, type, valuesList, condition='and')

# Set directory for full standard test report
outDir = os.path.join(os.getcwd(), 'tests', 'reports')
fU.makeADir(outDir)

# Start writing markdown style report for standard tests
reportFile = os.path.join(outDir, 'standardTestsReport.md')
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Standard Tests Report \n')
    pfile.write('## Compare com1DFA simulations to benchmark results \n')

# run Standard Tests sequentially
for test in testList:

    # avalanche directory
    avaDir = 'data' + os.sep + test['AVANAME']

    # get path to executable
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
    com1Exe = cfgCom1DFA['GENERAL']['com1Exe']

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('Current avalanche: %s', avaDir)

    # Load input parameters from configuration file for standard tests
    # write config to log file
    avaName = os.path.basename(avaDir)
    standardCfg = os.path.join('..', 'benchmarks', test['NAME'], '%s_com1DFACfg.ini' % avaName)
    cfg = cfgUtils.getModuleConfig(com1DFA, standardCfg)
    cfg['GENERAL']['com1Exe'] = com1Exe

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir,  keep=logName)

    # Set timing
    startTime = time.time()
    # Run Standalone DFA
    reportDictList = com1DFA.com1DFAMain(cfg, avaDir)

    # Print time needed
    endTime = time.time()
    timeNeeded = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeeded)))

    # Generata plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avaDir, cfg, cfgMain['FLAGS'])

    # Set directory for report
    reportDir = os.path.join(avaDir, 'Outputs', 'com1DFA', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    # get release area scenarios
    relArea = []
    for dict in reportDictList:
        relArea.append(dict['Simulation Parameters']['Release Area Scenario'])
    relAreaSet = sorted(set(relArea))

    for rel in relAreaSet:

        # -----------Compare to benchmark results
        # Fetch simulation info from benchmark results
        benchDictList = simParameters.fetchBenchParameters(avaDir)
        benchDict = ''
        for bDict in benchDictList:
            if bDict['testName'] ==  test['NAME'] and rel == bDict['Simulation Parameters']['Release Area Scenario']:
                    benchDict = bDict
        benchSimName = benchDict['simName']
        # Check if simulation with entrainment and/or resistance or standard simulation
        simType = 'null'
        if 'entres' in benchSimName:
            simType = 'entres'

        # Fetch correct reportDict according to flagEntRes
        for dict in reportDictList:
            if simType in dict['simName']['name'] and dict['Simulation Parameters']['Release Area Scenario'] == rel:
                reportD = dict

        # Add info on run time
        reportD['runTime'] = timeNeeded

        # +++++++Aimec analysis
        # load configuration
        aimecCfg = os.path.join('..', 'benchmarks', test['NAME'], '%s_AIMECCfg.ini' % avaName)
        cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC, aimecCfg)
        cfgAimecSetup = cfgAimec['AIMECSETUP']
        cfgAimecSetup['testName'] = test['NAME']

        # Setup input from com1DFA and reference
        pathDictList = dfa2Aimec.dfaComp2Aimec(avaDir, cfgAimecSetup)

        for pathD in pathDictList:
            if pathD == reportD['simName']['name']:
                pathDict = pathDictList[pathD]

        # TODO: define referenceFile
        pathDict['numSim'] = len(pathDict['ppr'])
        pathDict['referenceFile'] = 0

        # Extract input file locations
        pathDict = ana3AIMEC.readAIMECinputs(avaDir, pathDict, dirName=reportD['simName']['name'])

        # perform analysis
        rasterTransfo, newRasters, resAnalysis = ana3AIMEC.AIMEC2Report(pathDict, cfgAimec)

        # +++++++++++Aimec analysis


        # Create plots for report
        # Load input parameters from configuration file
        cfgRep = cfgUtils.getModuleConfig(generateCompareReport)

        # REQUIRED+++++++++++++++++++
        # Which parameter to filter data, e.g. varPar = 'simType', values = ['null'] or
        # varPar = 'Mu', values = ['0.055', '0.155']; values need to be given as list, also if only one value
        outputVariable = ['ppr', 'pfd', 'pfv']
        values = simType
        parameter = 'simType'
        plotListRep = {}
        reportD['Simulation Difference'] = {}
        reportD['Simulation Stats'] = {}
        # ++++++++++++++++++++++++++++

        # Plot data comparison for all output variables defined in suffix
        for var in outputVariable:
            plotList = outQuickPlot.quickPlot(avaDir, test['NAME'], var, values, parameter, cfgMain, cfgRep, rel, simType=simType)
            for pDict in plotList:
                if rel in pDict['relArea']:
                    plotDict = pDict
            for plot in plotDict['plots']:
                plotListRep.update({var: plot})
                reportD['Simulation Difference'].update({var: plotDict['difference']})
                reportD['Simulation Stats'].update({var: plotDict['stats']})

        # copy files to report directory
        plotPaths = generateCompareReport.copyPlots(avaName, test['NAME'], outDir, plotListRep, rel)

        # add plot info to general report Dict
        reportD['Simulation Results'] = plotPaths

        # write report
        generateCompareReport.writeCompareReport(reportFile, reportD, benchDict, avaName, cfgRep)
