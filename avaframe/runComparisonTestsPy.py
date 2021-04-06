"""
    Run script for running the standard tests
"""

# Load modules
import os
import time

# Local imports
from avaframe.com1DFAPy import runCom1DFA
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
from benchmarks import simParameters

# log file name; leave empty to use default runLog.log
logName = 'runComparisonTestsPy'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
outNew = 'OutputsFloatSphArtifdt001'
testList = ['avaInclinedPlane', 'avaInclinedPlaneDiag', 'avaParabola', 'avaHelix', 'avaHelixChannel', 'avaWog', 'avaKot']
simType = 'null'
# Set directory for full standard test report
outDir = os.path.join(os.getcwd(), 'tests', 'reportsComparisonSamosPy')
fU.makeADir(outDir)

# Start writing markdown style report for standard tests
reportFile = os.path.join(outDir, 'ComparisonSamosPy.md')
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Standard Tests Report \n')
    pfile.write('## Compare com1DFAPy simulations to benchmark results \n')

# run Standard Tests sequentially
for test in testList:

    avaDir = 'data' + os.sep + test


    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('Current avalanche: %s', avaDir)

    # Load input parameters from configuration file for standard tests
    # write config to log file
    avaName = os.path.basename(avaDir)

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir,  keep=logName)

    #####################################################################
    ########################## Run SamosAT ##############################

    # get path to executable
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
    # Set timing
    startTime = time.time()
    # Run Standalone DFA
    reportDictListSamos = com1DFA.com1DFAMain(cfgCom1DFA, avaDir)

    # Print time needed
    endTime = time.time()
    timeNeededSamos = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeededSamos)))

    # Generata plots for all peakFiles
    plotDictSamos = oP.plotAllPeakFields(avaDir, cfgCom1DFA, cfgMain['FLAGS'])

    # Set directory for report
    outDirOld = os.path.join(avaDir, 'Outputs')
    outDirNew = os.path.join(avaDir, outNew)
    reportDir = os.path.join(outDirOld, 'com1DFA', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictListSamos, cfgMain['FLAGS'], plotDictSamos)


    #####################################################################
    ########################## Run Com1DFAPy #############################
    # Set timing
    startTime = time.time()
    # Run Standalone DFA
    # call com1DFAPy to perform simulation - provide configuration file and release thickness function
    _, _, _, _, plotDictPy, reportDictListPy = runCom1DFA.runCom1DFAPy(avaDir=avaDir, cfgFile='', relThField='')

    # Print time needed
    endTime = time.time()
    timeNeededPy = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeededPy)))

    # Set directory for report
    reportDir = os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictListPy, cfgMain['FLAGS'], plotDictPy)

    #######################################################
    ############ Analyze results ##########################
    reportSamos = reportDictListSamos[0]
    reportPy = reportDictListPy[0]
    # Add info on run time
    reportSamos['runTime'] = timeNeededSamos
    reportPy['runTime'] = timeNeededPy
    rel = reportSamos['Simulation Parameters']['Release Area Scenario']

    # +++++++Aimec analysis
    # load configuration
    cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC)
    initProj.cleanModuleFiles(avaDir, ana3AIMEC)

    # write configuration to file
    cfgUtils.writeCfgFile(avaDir, ana3AIMEC, cfgAimec)

    cfgAimecSetup = cfgAimec['AIMECSETUP']

    # Setup input from com1DFA and reference
    pathDictList = dfa2Aimec.dfaComp2Aimec(avaDir, cfgAimecSetup)
    pathDict = pathDictList[reportSamos['simName']['name']]
    pathDict['numSim'] = len(pathDict['ppr'])
    log.info('reference file comes from: %s' % pathDict['compType'][1])

    # Extract input file locations
    pathDict = aimecTools.readAIMECinputs(avaDir, pathDict, dirName=reportSamos['simName']['name'])
    # ana3AIMEC.mainAIMEC(pathDict, cfgAimec)

    # perform analysis
    rasterTransfo, newRasters, resAnalysis = ana3AIMEC.AIMEC2Report(pathDict, cfgAimec)

    # add aimec results to report dictionary
    reportSamos, reportPy = ana3AIMEC.aimecRes2ReportDict(resAnalysis, reportSamos, reportPy, pathDict['referenceFile'])
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
    reportPy['Simulation Difference'] = {}
    reportPy['Simulation Stats'] = {}
    # ++++++++++++++++++++++++++++

    # Plot data comparison for all output variables defined in suffix
    for var in outputVariable:
        plotList = outQuickPlot.quickPlot(avaDir, test, var, values, parameter, cfgMain, cfgRep, rel, simType=simType, comModule='com1DFA', comModule2='com1DFAPy')
        for pDict in plotList:
            if rel in pDict['relArea']:
                plotDict = pDict
        for plot in plotDict['plots']:
            plotListRep.update({var: plot})
            reportPy['Simulation Difference'].update({var: plotDict['difference']})
            reportPy['Simulation Stats'].update({var: plotDict['stats']})

    # copy files to report directory
    plotPaths = generateCompareReport.copyQuickPlots(avaName, test, outDir, plotListRep, rel)
    aimecPlots = [resAnalysis['slCompPlot'], resAnalysis['areasPlot']]
    plotPaths = generateCompareReport.copyAimecPlots(aimecPlots, test, outDir, rel, plotPaths)

    # add plot info to general report Dict
    reportSamos['Simulation Results'] = plotPaths

    reportPy['Test Info'] = {'type': 'text', 'Test Info': 'This test uses a bowl-shaped geometry.'}
    reportSamos['Test Info'] = {'type': 'text', 'Test Info': 'This test uses a bowl-shaped geometry.'}

    # write report
    generateCompareReport.writeCompareReport(reportFile, reportPy, reportSamos, avaName, cfgRep)

    # rename output folder
    # initProj.cleanModuleFiles(avaDir, com1DFA, alternativeName=outNew)
    os.rename(outDirOld, outDirNew)
