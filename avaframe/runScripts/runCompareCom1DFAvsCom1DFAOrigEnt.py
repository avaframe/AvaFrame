"""
    Run script for comparing com1DFAOrig to com1DFA results
    configuration for com1DFAOrig and com1DFA is read from com1DFAOrigCfg.ini and com1DFACfg.ini or if available
    your local copy of them; for aimec the default configuration is read and some parameters are directly set
    within this runScript
    avalancheDir has to be set as well as the simType that should be compared in the top section: 'required settings'
"""

# Load modules
import os
import pathlib

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.com1DFAOrig import com1DFAOrig
from avaframe.log2Report import generateReport as gR
from avaframe.log2Report import generateCompareReport
from avaframe.ana3AIMEC import ana3AIMEC, dfa2Aimec, aimecTools
from avaframe.out3Plot import outQuickPlot
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# +++++++++REQUIRED+++++++++++++
# name of avalanche directory as list, multiple possible
avaList = ['avaKot1', 'avaAlr1', 'avaHit1', 'avaMal1', 'avaWog1', 'avaHelixChannel']
# simType that should be compared (options: null, ent, entres, res) -
# must also be set in the ini files for the computational modules
simType = 'ent'
values = simType
parameter = 'simType'
# Which result types for comparison plots
outputVariable = ['ppr', 'pft', 'pfv']
# aimec settings
aimecResType = 'ppr'
aimecThresholdValue = '1'
aimecDiffLim = '5'
aimecContourLevels = '1|3|5|10'
aimecFlagMass = 'True'
aimecComModules = 'com1DFAOrig|com1DFA'
startOfRunoutAreaAngle = '20'
# ++++++++++++++++++++++++++++++

# define simTypeString for finding simulations
simTypeString = '_' + simType + '_'

# log file name; leave empty to use default runLog.log
logName = 'runCompareCom1DFAOrigvsCom1DFAEnt'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# Set directory for full standard test report
outDirReport = pathlib.Path.cwd() / 'tests' / 'reportscom1DFAOrigvsCom1DFAEnt'
fU.makeADir(outDirReport)

# Start writing markdown style report for standard tests
reportFile = outDirReport / 'com1DFAOrigvsCom1DFAEnt.md'
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Comparison Report \n')
    pfile.write('## Compare com1DFAOrig simulation to com1DFA simulation results \n')

# run Standard Tests sequentially
for avaName in avaList:

    # set avaDir
    avaDir = 'data' + os.sep + avaName

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('Current avalanche: %s', avaDir)
    outDir = pathlib.Path(avaDir, 'Outputs')

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    #####################################################################
    # ######################### Run com1DFAOrig ##############################
    # get module configuration (path to executable...)
    cfgCom1DFAOrig = cfgUtils.getModuleConfig(com1DFAOrig)
    # Run Standalone DFA
    reportDictListcom1DFAOrig = com1DFAOrig.com1DFAOrigMain(cfgCom1DFAOrig, avaDir)

    for reportD1 in reportDictListcom1DFAOrig:
        simName1 = reportD1['simName']['name']
        parameterDict, reportD1 = fU.extractParameterInfo(avaDir, simName1, reportD1)

    # Generata plots for all peakFiles
    modNameOrig = 'com1DFAOrig'
    plotDictcom1DFAOrig = oP.plotAllPeakFields(avaDir, cfgMain['FLAGS'], modNameOrig)

    # Set directory for com1DFA report
    reportDirOrig = outDir / 'com1DFAOrig' / 'reports'
    # write report
    gR.writeReport(reportDirOrig, reportDictListcom1DFAOrig, cfgMain['FLAGS'].getboolean('reportOneFile'),
                   plotDictcom1DFAOrig)

    #####################################################################
    # ######################### Run Com1DFA #############################
    # Run python DFA
    # call com1DFA to perform simulation - provide configuration file and release thickness function
    _, plotDictcom1DFA, reportDictListcom1DFA, simDF = com1DFA.com1DFAMain(avaDir, cfgMain, cfgFile='')

    # Set directory for com1DFA report
    reportDir = outDir / 'com1DFA' / 'reports'
    # write report
    gR.writeReport(reportDir, reportDictListcom1DFA, cfgMain['FLAGS'].getboolean('reportOneFile'), plotDictcom1DFA)

    #######################################################
    # ########### Analyze results ##########################
    # Aimec analysis
    initProj.cleanModuleFiles(avaDir, ana3AIMEC)

    # get release area scenarios
    relArea = []
    for dict in reportDictListcom1DFA:
        relArea.append(dict['Simulation Parameters']['Release Area Scenario'])
    relAreaSet = sorted(set(relArea))

    for rel in relAreaSet:
        reportDcom1DFAOrig = ''
        for dict in reportDictListcom1DFAOrig:
            com1DFASimName = dict['simName']['name']
            if (rel == dict['Simulation Parameters']['Release Area Scenario']):
                if simTypeString in com1DFASimName:
                    reportDcom1DFAOrig = dict
                    log.info('Comparison based on releaseScenario: %s and simType: %s' % (rel, simType))
                    log.info('Reference simulation: %s' % com1DFASimName)
                    break
                else:
                    log.error('No reference simulation found based on releaseScenario: %s and simType: %s' %
                        (rel, simType))

        simNameRef = reportDcom1DFAOrig['simName']['name']
        refDir = pathlib.Path(avaDir, 'Outputs', 'com1DFAOrig', 'peakFiles')
        if reportDcom1DFAOrig:
            com1DFASimName = reportDcom1DFAOrig['simName']['name']
            # Fetch corresponding com1DFA
            for dict in reportDictListcom1DFA:
                if (simTypeString in dict['simName']['name'] and
                        dict['Simulation Parameters']['Release Area Scenario'] == rel):
                    reportDcom1DFA = dict
                    log.info('Comparison simulation: %s' % dict['simName']['name'])
                    break
                else:
                    log.error('No matching simulation found based on releaseScenario: %s and simType: %s' %
                        (rel, simType))

            simNameComp = reportDcom1DFA['simName']['name']
            compDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'peakFiles')

            # +++++++Aimec analysis
            # load configuration
            cfgAimec = cfgUtils.getDefaultModuleConfig(ana3AIMEC)
            cfgAimec['AIMECSETUP']['resType'] = aimecResType
            cfgAimec['AIMECSETUP']['thresholdValue'] = aimecThresholdValue
            cfgAimec['AIMECSETUP']['diffLim'] = aimecDiffLim
            cfgAimec['AIMECSETUP']['contourLevels'] = aimecContourLevels
            cfgAimec['FLAGS']['flagMass'] = aimecFlagMass
            cfgAimec['AIMECSETUP']['comModules'] = aimecComModules
            cfgAimec['AIMECSETUP']['startOfRunoutAreaAngle'] = startOfRunoutAreaAngle

            # write configuration to file
            cfgUtils.writeCfgFile(avaDir, ana3AIMEC, cfgAimec)
            # Setup input from com1DFA and com1DFAPy
            inputsDF, pathDict = dfa2Aimec.dfaBench2Aimec(avaDir, cfgAimec, simNameRef, simNameComp)
            pathDict['refSimulation'] = inputsDF.index[0]
            log.info('reference file comes from: %s' % pathDict['refSimulation'])

            # Extract input file locations
            pathDict = aimecTools.readAIMECinputs(avaDir, pathDict, dirName=reportDcom1DFAOrig['simName']['name'])

            # perform analysis
            rasterTransfo, resAnalysisDF, aimecPlotDict = ana3AIMEC.mainAIMEC(pathDict, inputsDF, cfgAimec)

            # add aimec results to report dictionary
            reportDcom1DFA, reportDcom1DFAOrig = ana3AIMEC.aimecRes2ReportDict(resAnalysisDF, reportDcom1DFA,
                reportDcom1DFAOrig, pathDict['refSimulation'])

            # Create plots for report
            # Load input parameters from configuration file
            cfgRep = cfgUtils.getModuleConfig(generateCompareReport)
            cfgRep['PLOT']['refModel'] = 'dfa'

            # initialize plot dict
            plotListRep = {}
            reportDcom1DFA['Simulation Difference'] = {}
            reportDcom1DFA['Simulation Stats'] = {}

            # Plot data comparison for all output variables defined in suffix
            for var in outputVariable:
                plotDict = outQuickPlot.quickPlotBench(avaDir, simNameRef, simNameComp, refDir, compDir, cfgMain, var)
                for plot in plotDict['plots']:
                    plotListRep.update({var: plot})
                    reportDcom1DFA['Simulation Difference'].update({var: plotDict['difference']})
                    reportDcom1DFA['Simulation Stats'].update({var: plotDict['stats']})

            # copy files to report directory
            plotPaths = generateCompareReport.copyQuickPlots(avaName, avaName, outDirReport, plotListRep, rel=rel)
            aimecPlots = [aimecPlotDict['slCompPlot'], aimecPlotDict['areasPlot']]
            if cfgAimec.getboolean('FLAGS', 'flagMass'):
                aimecPlots.append(aimecPlotDict['massAnalysisPlot'])
            plotPaths = generateCompareReport.copyAimecPlots(aimecPlots, avaName, outDirReport, plotPaths, rel=rel)

            # add plot info to general report Dict
            reportDcom1DFA['Simulation Results'] = plotPaths

            reportDcom1DFAOrig['Test Info'] = {'type': 'text',
                           'Test Info': 'Compare com1DFAOrig (Reference) to com1DFA (Simulation) results.'}

            # write report
            generateCompareReport.writeCompareReport(reportFile, reportDcom1DFA, reportDcom1DFAOrig, avaName, cfgRep)
