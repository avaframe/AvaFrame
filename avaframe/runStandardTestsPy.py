"""
    Run script for running the standard tests
"""

# Load modules
import os
import time

# Local imports
from avaframe.com1DFAPy import com1DFA
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
from benchmarks import simParameters

# log file name; leave empty to use default runLog.log
logName = 'runStandardTestsPy'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False)
# filter benchmarks for tag standardTest
type = 'TAGS'
valuesList = ['standardTest', 'null']
testList = tU.filterBenchmarks(testDictList, type, valuesList, condition='and')

# Set directory for full standard test report
outDir = os.path.join(os.getcwd(), 'tests', 'reportsPy')
fU.makeADir(outDir)

# Start writing markdown style report for standard tests
reportFile = os.path.join(outDir, 'standardTestsReportPy.md')
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Standard Tests Report \n')
    pfile.write('## Compare com1DFAPy simulations to benchmark results \n')

# run Standard Tests sequentially
for test in testList:

    avaDir = 'data' + os.sep + test['AVANAME']

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('Current avalanche: %s', avaDir)

    # Load input parameters from configuration file for standard tests
    # write config to log file
    avaName = os.path.basename(avaDir)
    standardCfg = os.path.join('..', 'benchmarks', test['NAME'], '%s_com1DFACfgPy.ini' % test['AVANAME'])
    cfg = cfgUtils.getModuleConfig(com1DFA, standardCfg)
    cfg['GENERAL']['frictModel'] = 'samosAT'
    cfg['GENERAL']['sphKernelRadius'] = '5.'
    cfg['GENERAL']['massPerParticleDeterminationMethod'] = 'MPPDH'
    cfg['GENERAL']['curvature'] = '1'
    cfg['GENERAL']['meshCellSize'] = '5'
    cfg['GENERAL']['meshCellSizeThreshold'] = '0.001'
    cfg['GENERAL']['sphOption'] = '1'
    cfg['GENERAL']['interpOption'] = '2'
    cfg['GENERAL']['methodMeshNormal'] = '1'
    cfg['GENERAL']['resType'] = 'ppr_pfd_pfv'
    cfg['GENERAL']['tSteps'] = '0'
    cfg['GENERAL']['seed'] = '12345'
    modName = 'com1DFAPy'

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir,  keep=logName)

    # Set timing
    startTime = time.time()

    # Run Standalone DFA
    # call com1DFAPy to perform simulation - provide configuration
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
        benchDictList = simParameters.fetchBenchParameters(test['AVANAME'])
        benchDict = ''
        for bDict in benchDictList:
            if test['NAME'] == bDict['testName'] and rel == bDict['Simulation Parameters']['Release Area Scenario']:
                benchDict = bDict
        benchSimName = benchDict['simName']['name']
        # Check if simulation with entrainment and/or resistance or standard simulation
        simType = benchDict['simType']

        # Fetch correct reportDict according to flagEntRes
        for dict in reportDictList:
            if simType in dict['simName']['name'] and dict['Simulation Parameters']['Release Area Scenario'] == rel:
                reportD = dict

        # Add info on run time
        reportD['runTime'] = timeNeeded

        # +++++++Aimec analysis
        # load configuration
        aimecCfg = os.path.join('..', 'benchmarks', test['NAME'], '%s_AIMECPyCfg.ini' % test['AVANAME'])
        cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC, aimecCfg)
        cfgAimecSetup = cfgAimec['AIMECSETUP']
        cfgAimecSetup['testName'] = test['NAME']
        cfgAimec['AIMECSETUP']['resType'] = 'ppr'
        cfgAimec['AIMECSETUP']['thresholdValue'] = '1'
        cfgAimec['AIMECSETUP']['diffLim'] = '5'
        cfgAimec['AIMECSETUP']['contourLevels'] = '1|3|5|10'
        cfgAimec['FLAGS']['flagMass'] = 'False'
        cfgAimec['AIMECSETUP']['comModules'] = 'benchmarkReference|com1DFAPy'

        # Setup input from com1DFA and reference
        pathDictList = dfa2Aimec.dfaComp2Aimec(avaDir, cfgAimecSetup)
        for pathD in pathDictList:
            if pathD == reportD['simName']['name']:
                pathDict = pathDictList[pathD]

        #
        pathDict['numSim'] = len(pathDict['ppr'])
        log.info('reference file comes from: %s' % pathDict['compType'][1])

        # Extract input file locations
        pathDict = aimecTools.readAIMECinputs(avaDir, pathDict, dirName=reportD['simName']['name'])

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
        values = simType
        parameter = 'simType'
        plotListRep = {}
        reportD['Simulation Difference'] = {}
        reportD['Simulation Stats'] = {}
        # ++++++++++++++++++++++++++++

        # Plot data comparison for all output variables defined in suffix
        for var in outputVariable:
            plotList = outQuickPlot.quickPlot(avaDir, test['NAME'], var, values, parameter, cfgMain, cfgRep, rel, simType=simType, comModule='com1DFAPy')
            for pDict in plotList:
                if rel in pDict['relArea']:
                    plotDict = pDict
            for plot in plotDict['plots']:
                plotListRep.update({var: plot})
                reportD['Simulation Difference'].update({var: plotDict['difference']})
                reportD['Simulation Stats'].update({var: plotDict['stats']})

        # copy files to report directory
        plotPaths = generateCompareReport.copyQuickPlots(avaName, test['NAME'], outDir, plotListRep, rel)
        aimecPlots = [resAnalysis['slCompPlot'], resAnalysis['areasPlot']]
        plotPaths = generateCompareReport.copyAimecPlots(aimecPlots, test['NAME'], outDir, rel, plotPaths)

        # add plot info to general report Dict
        reportD['Simulation Results'] = plotPaths

        # write report
        generateCompareReport.writeCompareReport(reportFile, reportD, benchDict, avaName, cfgRep)
