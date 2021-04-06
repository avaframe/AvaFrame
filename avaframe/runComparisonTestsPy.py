"""
    Run script for running the standard tests
"""

# Load modules
import os
import time

# Local imports
from avaframe.com1DFAPy import runCom1DFA
from avaframe.com1DFA import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.log2Report import generateCompareReport
from avaframe.ana3AIMEC import ana3AIMEC, dfa2Aimec, aimecTools
from avaframe.out3Plot import outQuickPlot
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'runComparisonTestsPy'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
outNew = 'OutputsFloatAllSamos'
testList = ['avaInclinedPlane', 'avaParabola', 'avaHelix', 'avaHelixChannel', 'avaWog', 'avaKot'] #, 'avaInclinedPlaneDiag'
simType = 'null'
# Set directory for full standard test report
outDir = os.path.join(os.getcwd(), 'tests', 'reportscom1DFAvsPy')
fU.makeADir(outDir)

# Start writing markdown style report for standard tests
reportFile = os.path.join(outDir, 'com1DFAvsPy.md')
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Standard Tests Report \n')
    pfile.write('## Compare com1DFA simulation to com1DFAPy results \n')

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
    ########################## Run com1DFA ##############################
    # use the local config file

    # get path to executable
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
    # Set timing
    startTime = time.time()
    # Run Standalone DFA
    reportDictListcom1DFA = com1DFA.com1DFAMain(cfgCom1DFA, avaDir)

    # Print time needed
    endTime = time.time()
    timeNeededcom1DFA = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeededcom1DFA)))

    # Generata plots for all peakFiles
    plotDictcom1DFA = oP.plotAllPeakFields(avaDir, cfgCom1DFA, cfgMain['FLAGS'])

    # Set directory for report
    outDirOld = os.path.join(avaDir, 'Outputs')
    # name needed to rename the outputs at the end
    outDirNew = os.path.join(avaDir, outNew)
    reportDir = os.path.join(outDirOld, 'com1DFA', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictListcom1DFA, cfgMain['FLAGS'], plotDictcom1DFA)


    #####################################################################
    ########################## Run Com1DFAPy #############################
    # Set timing
    startTime = time.time()
    # Run Standalone DFA
    # call com1DFAPy to perform simulation - provide configuration file and release thickness function
    _, _, _, _, plotDictcom1DFAPy, reportDictListcom1DFAPy = runCom1DFA.runCom1DFAPy(avaDir=avaDir, cfgFile='', relThField='')

    # Print time needed
    endTime = time.time()
    timeNeededcom1DFAPy = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeededcom1DFAPy)))

    # Set directory for report
    reportDir = os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictListcom1DFAPy, cfgMain['FLAGS'], plotDictcom1DFAPy)

    #######################################################
    ############ Analyze results ##########################
    reportcom1DFA = reportDictListcom1DFA[0]
    reportcom1DFAPy = reportDictListcom1DFAPy[0]
    # Add info on run time
    reportcom1DFA['runTime'] = timeNeededcom1DFA
    reportcom1DFAPy['runTime'] = timeNeededcom1DFAPy
    rel = reportcom1DFA['Simulation Parameters']['Release Area Scenario']

    # +++++++Aimec analysis
    # load configuration
    cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC)
    initProj.cleanModuleFiles(avaDir, ana3AIMEC)

    # write configuration to file
    cfgUtils.writeCfgFile(avaDir, ana3AIMEC, cfgAimec)

    cfgAimecSetup = cfgAimec['AIMECSETUP']

    # Setup input from com1DFA and reference
    pathDictList = dfa2Aimec.dfaComp2Aimec(avaDir, cfgAimecSetup)
    pathDict = pathDictList[reportcom1DFA['simName']['name']]
    pathDict['numSim'] = len(pathDict['ppr'])
    log.info('reference file comes from: %s' % pathDict['compType'][1])

    # Extract input file locations
    pathDict = aimecTools.readAIMECinputs(avaDir, pathDict, dirName=reportcom1DFA['simName']['name'])

    # perform analysis
    rasterTransfo, newRasters, resAnalysis = ana3AIMEC.AIMEC2Report(pathDict, cfgAimec)

    # add aimec results to report dictionary
    reportSamos, reportPy = ana3AIMEC.aimecRes2ReportDict(resAnalysis, reportcom1DFA, reportcom1DFAPy, pathDict['referenceFile'])
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

    reportcom1DFAPy['Test Info'] = {'type': 'text', 'Test Info': 'Compare com1DFA to com1DFAPy results.'}
    reportcom1DFA['Test Info'] = {'type': 'text', 'Test Info': 'Compare com1DFA to com1DFAPy results.'}

    # write report
    generateCompareReport.writeCompareReport(reportFile, reportcom1DFAPy, reportcom1DFA, avaName, cfgRep)

    # rename output folder
    # os.rename(outDirOld, outDirNew)
