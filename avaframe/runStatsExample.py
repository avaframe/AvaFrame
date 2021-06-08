"""
    Run script for computing statistics of results including the runout derived from aimec
"""

# Load modules
import os
import time
import glob
import configparser
import shutil
import numpy as np

# Local imports
from avaframe.out3Plot import statsPlots as sPlot
from avaframe.ana4Stats import getStats
from avaframe.com1DFAOrig import com1DFAOrig
from avaframe.ana3AIMEC import ana3AIMEC, dfa2Aimec
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.log2Report import generateReport as gR
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.out3Plot import plotUtils
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runGetStats'
modName = 'com1DFAOrig'

# Load general configuration filee
cfgMain = cfgUtils.getGeneralConfig()
flagShow = cfgMain['FLAGS'].getboolean('showPlot')
# get path to executable
cfgCom1DFA = cfgUtils.getModuleConfig(com1DFAOrig)
com1Exe = cfgCom1DFA['GENERAL']['com1Exe']

peakDictList = []
avalancheDirs = ['data/avaHockeyChannel','data/avaHockeyChannel']
cfgFull = cfgUtils.getModuleConfig(getStats)
cfg = cfgFull['GENERAL']
# set output directory, first ava in list
outDir = os.path.join(avalancheDirs[0], 'Outputs', 'ana4Stats')
cfg['outDir'] = outDir
# Specify where you want the results to be stored
fU.makeADir(outDir)

# Start logging
log = logUtils.initiateLogger(outDir, logName)

count = 0
for avaDir in avalancheDirs:

    # --- run Com1DFA -----
    avaName = os.path.basename(avaDir)
    avaNameTest = avaName + 'StatsTest'
    statsSimCfg = os.path.join('..', 'benchmarks', avaNameTest, '%sStats%d_com1DFACfg.ini' % (avaName, count))
    cfgDFA = cfgUtils.getModuleConfig(com1DFAOrig, statsSimCfg)
    cfgDFA['GENERAL']['com1Exe'] = com1Exe

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # Run Standalone DFA
    reportDictList = com1DFAOrig.com1DFAOrigMain(cfgDFA, avaDir)
    fU.makeADir(outDir)

    # Generata plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avaDir, cfgDFA, cfgMain['FLAGS'], modName)

    # Set directory for report
    reportDir = os.path.join(avaDir, 'Outputs', 'com1DFAOrig', 'reports')
    # write report
    gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)

    if cfg.getboolean('aimec') == True:
        # run aimec
        statsAimecCfg = os.path.join('..', 'benchmarks', avaNameTest, '%sStats_ana3AIMECCfg.ini' % (avaName))
        cfgAIMEC = cfgUtils.getModuleConfig(ana3AIMEC, statsAimecCfg)
        cfgAimecSetup = cfgAIMEC['AIMECSETUP']

        # Setup input from com1DFA
        dfa2Aimec.mainDfa2Aimec(avaDir, cfgAimecSetup, modName)
        # Extract input file locations
        cfgPath = ana3AIMEC.readAIMECinputs(avaDir, dirName='com1DFA')
        # Run AIMEC postprocessing
        ana3AIMEC.mainAIMEC(cfgPath, cfgAIMEC)

    # copy Outputs to scenario directory
    dirOutput = os.path.join(avaDir, 'Outputs')
    outputNew = os.path.join(avaDir, 'Outputs%s' % cfgDFA['GENERAL']['releaseScenario'])
    shutil.copytree(dirOutput, outputNew, dirs_exist_ok=True)

    # ----- determine max values of peak fields
    # set directory of peak files
    inputDir = os.path.join(avaDir, 'Outputs%s' % cfgDFA['GENERAL']['releaseScenario'], 'com1DFAOrig', 'peakFiles')

    # get max values of peak files
    peakValues = getStats.extractMaxValues(inputDir, cfgDFA, avaDir, nameScenario=cfgDFA['GENERAL']['releaseScenario'])

    # log to screen
    for key in peakValues:
        print('peakValues:', key, peakValues[key])

    peakDictList.append(peakValues)
    count = count + 1

 #++++++++++++++ Plot max values +++++++++++++++++
sPlot.plotValuesScatter(peakDictList, 'pfd', 'pfv', cfgDFA['PARAMETERVAR']['varPar'], cfg, avalancheDirs, flagShow)
sPlot.plotValuesScatterHist(peakDictList, 'pfd', 'pfv', cfgDFA['PARAMETERVAR']['varPar'], cfg, avalancheDirs, flagShow, flagHue=True)

log.info('Plots have been saved to: %s' % outDir)
