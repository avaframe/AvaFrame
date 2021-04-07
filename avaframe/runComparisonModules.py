"""
    Run script for running the standard tests
"""

# Load modules
import os

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
logName = 'runComparisonModules'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testList = ['avaInclinedPlane', 'avaParabola', 'avaHelix', 'avaHelixChannel', 'avaWog', 'avaKot']
# Set directory for full standard test report
outDirReport = os.path.join(os.getcwd(), 'tests', 'reportscom1DFAvsPy')
fU.makeADir(outDirReport)

# Start writing markdown style report for standard tests
reportFile = os.path.join(outDirReport, 'com1DFAvsPy.md')
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Standard Tests Report \n')
    pfile.write('## Compare com1DFA simulation to com1DFAPy simulation results \n')

# run Standard Tests sequentially
for avaName in testList:

    avaDir = 'data' + os.sep + avaName

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('Current avalanche: %s', avaDir)
    outDir = os.path.join(avaDir, 'Outputs')

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir,  keep=logName)

    #####################################################################
    # ######################### Run com1DFA ##############################
    # get module configuration (path to executable...)
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
    # Run Standalone DFA
    reportDictListcom1DFA = com1DFA.com1DFAMain(cfgCom1DFA, avaDir)

    for reportD1 in reportDictListcom1DFA:
        simName1 = reportD1['simName']['name']
        parameterDict = fU.extractParameterInfo(avaDir, simName1, reportD1)

    # Generata plots for all peakFiles
    plotDictcom1DFA = oP.plotAllPeakFields(avaDir, cfgCom1DFA, cfgMain['FLAGS'])

    # Set directory for com1DFA report
    reportDir = os.path.join(outDir, 'com1DFA', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictListcom1DFA, cfgMain['FLAGS'], plotDictcom1DFA)

    #####################################################################
    # ######################### Run Com1DFAPy #############################
    # Run python DFA
    # call com1DFAPy to perform simulation - provide configuration file and release thickness function
    _, _, _, _, plotDictcom1DFAPy, reportDictListcom1DFAPy = runCom1DFA.runCom1DFAPy(avaDir=avaDir, cfgFile='', relThField='')

    # Set directory for com1DFAPy report
    reportDir = os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictListcom1DFAPy, cfgMain['FLAGS'], plotDictcom1DFAPy)

    #######################################################
    # ########### Analyze results ###########################
    # get release area scenarios
    relArea = []
    for dict in reportDictListcom1DFAPy:
        relArea.append(dict['Simulation Parameters']['Release Area Scenario'])
    relAreaSet = sorted(set(relArea))

    for rel in relAreaSet:
        reportDcom1DFAentres = ''
        for dict in reportDictListcom1DFA:
            com1DFASimName = dict['simName']['name']
            if (rel == dict['Simulation Parameters']['Release Area Scenario']):
                # is it an entres or null simulation?
                if ('entres' in com1DFASimName):
                    reportDcom1DFAentres = dict
                else:
                    reportDcom1DFA = dict
                    simType = 'null'
        # take the entres sim if it exists
        if reportDcom1DFAentres:
            reportDcom1DFA = reportDcom1DFAentres
            simType = 'entres'
        com1DFASimName = reportDcom1DFA['simName']['name']

        # Fetch corresponding com1DFAPy
        for dict in reportDictListcom1DFAPy:
            if simType in dict['simName']['name'] and dict['Simulation Parameters']['Release Area Scenario'] == rel:
                reportDcom1DFAPy = dict

        # Aimec analysis
        # load configuration
        cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC)
        initProj.cleanModuleFiles(avaDir, ana3AIMEC)

        # write configuration to file
        cfgUtils.writeCfgFile(avaDir, ana3AIMEC, cfgAimec)
        if simType == 'entres':
            cfgAimec['FLAGS']['analyzeMass'] = 'True'
        cfgAimecSetup = cfgAimec['AIMECSETUP']

        # Setup input from com1DFA and com1DFAPy
        pathDictList = dfa2Aimec.dfaComp2Aimec(avaDir, cfgAimecSetup)
        for pathD in pathDictList:
            if pathD == reportDcom1DFA['simName']['name']:
                pathDict = pathDictList[pathD]

        pathDict['numSim'] = len(pathDict['ppr'])
        log.info('reference file comes from: %s' % pathDict['compType'][1])

        # Extract input file locations
        pathDict = aimecTools.readAIMECinputs(avaDir, pathDict, dirName=reportDcom1DFA['simName']['name'])

        # perform analysis
        rasterTransfo, newRasters, resAnalysis = ana3AIMEC.AIMEC2Report(pathDict, cfgAimec)

        # add aimec results to report dictionary
        reportDcom1DFA, reportDcom1DFAPy = ana3AIMEC.aimecRes2ReportDict(resAnalysis, reportDcom1DFA, reportDcom1DFAPy, pathDict['referenceFile'])

        # Create plots for report
        # Load input parameters from configuration file
        cfgRep = cfgUtils.getModuleConfig(generateCompareReport)
        cfgRep['PLOT']['refModel'] = 'dfa'

        # REQUIRED+++++++++++++++++++
        # Which parameter to filter data, e.g. varPar = 'simType', values = ['null'] or
        # varPar = 'Mu', values = ['0.055', '0.155']; values need to be given as list, also if only one value
        outputVariable = ['ppr', 'pfd', 'pfv']
        values = simType
        parameter = 'simType'
        plotListRep = {}
        reportDcom1DFAPy['Simulation Difference'] = {}
        reportDcom1DFAPy['Simulation Stats'] = {}
        # ++++++++++++++++++++++++++++

        # Plot data comparison for all output variables defined in suffix
        for var in outputVariable:
            plotList = outQuickPlot.quickPlot(avaDir, avaName, var, values, parameter, cfgMain, cfgRep, rel, simType=simType, comModule='com1DFA', comModule2='com1DFAPy')
            for pDict in plotList:
                if rel in pDict['relArea']:
                    plotDict = pDict
            for plot in plotDict['plots']:
                plotListRep.update({var: plot})
                reportDcom1DFAPy['Simulation Difference'].update({var: plotDict['difference']})
                reportDcom1DFAPy['Simulation Stats'].update({var: plotDict['stats']})

        # copy files to report directory
        plotPaths = generateCompareReport.copyQuickPlots(avaName, avaName, outDir, plotListRep, rel)
        aimecPlots = [resAnalysis['slCompPlot'], resAnalysis['areasPlot'], resAnalysis['massAnalysisPlot']]
        plotPaths = generateCompareReport.copyAimecPlots(aimecPlots, avaName, outDir, rel, plotPaths)

        # add plot info to general report Dict
        reportDcom1DFAPy['Simulation Results'] = plotPaths

        reportDcom1DFA['Test Info'] = {'type': 'text', 'Test Info': 'Compare com1DFA (Reference) to com1DFAPy (Simulation) results.'}

        # write report
        generateCompareReport.writeCompareReport(reportFile, reportDcom1DFAPy, reportDcom1DFA, avaName, cfgRep)
