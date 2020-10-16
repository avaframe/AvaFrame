"""
    Run script for running Standalone DFA for the standard tests
    This file is part of Avaframe.
"""

# Load modules
import os

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.log2Report import generateCompareReport
from avaframe.out3SimpPlot import outPlotAllPeak as oP
from avaframe.out3Plot import outQuickPlot
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import time

# log file name; leave empty to use default runLog.log
logName = 'runStandardTests'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
#
standardNames = ['data/avaBowl',
                                'data/avaFlatPlane',
                                'data/avaHelix',
                                'data/avaHelixChannel',
                                'data/avaHockey',
                                'data/avaHockeySmoothChannel',
                                'data/avaHockeySmoothSmall',
                                'data/avaInclinedPlane']

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
for avaDir in standardNames:

    # Start logging
    log = logUtils.initiateLogger(avaDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avaDir)

    # Load input parameters from configuration file
    # write config to log file
    cfg = cfgUtils.getModuleConfig(com1DFA)

    # set timing
    startTime = time.time()
    # Run Standalone DFA
    reportDictList = com1DFA.runSamos(cfg, avaDir)
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

    # -----------Compare to benchmark results
    # Fetch simulation info from benchmark results
    benchDict = generateCompareReport.fetchBenchParameters(avaDir)
    benchSimName = benchDict['simName']
    # Check if simulation with entrainment and/or resistance or standard simulation
    simType = 'null'
    if 'entres' in benchSimName:
        simType = 'entres'

    # Fetch correct reportDict according to flagEntRes
    for dict in reportDictList:
        if simType in dict['simName']['name']:
            reportD = dict

    # Add info on run time
    reportD['runTime'] = timeNeeded

    # Create plots for report
    # Load input parameters from configuration file
    cfgRep = cfgUtils.getModuleConfig(generateCompareReport)

    # REQUIRED+++++++++++++++++++
    # Which parameter to filter data, e.g. varPar = 'simType', values = ['null'] or
    # varPar = 'Mu', values = ['0.055', '0.155']; values need to be given as list, also if only one value
    outputVariable = ['ppr', 'pfd', 'pv']
    values = simType
    parameter = 'simType'
    plotListRep = {}
    # ++++++++++++++++++++++++++++

    # Plot data comparison for all output variables defined in suffix
    for var in outputVariable:
        plotList = outQuickPlot.quickPlot(avaDir, var, values, parameter, cfgMain, cfgRep)
        for plot in plotList:
            plotListRep.update({var: plot})

    # copy files to report directory
    avaName = os.path.basename(avaDir)
    plotPaths = generateCompareReport.copyPlots(avaName, outDir, plotListRep)

    # add plot info to general report Dict
    reportD['Simulation Results'] = plotPaths

    # write report
    generateCompareReport.writeCompareReport(reportFile, reportD, benchDict, avaName)
