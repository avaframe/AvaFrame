"""
    Run com1DFA kernel and compare it to the similarity solution
    This script computes the similarity solution for a gliding avalanche on
    a inclined plane according to similarity solution from :
    Hutter, K., Siegel, M., Savage, S.B. et al.
    Two-dimensional spreading of a granular avalanche down an inclined plane
    Part I. theory. Acta Mechanica 100, 37–68 (1993).
    https://doi.org/10.1007/BF01176861
    and compares it to the DFA kernel com1DFA
"""

import pathlib
import pandas as pd
from configupdater import ConfigUpdater

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import logUtils
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in1Data import getInput as gI
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.outAna1Plots as outAna1Plots


# +++++++++ REQUIRED ++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# if left empty, use the simiSolTestCfg.ini and local_simiSolTestCfg.ini configuration files
# use 'ana1Tests/figConvergence_simiSolTestCfg.ini' to produce the similarity solution plots from the Theory Paper
fileOverride = ''
# ++++++++++++++++++++++++++++++

# Load general configuration
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaSimilaritySol'

# Clean input directory(ies) of old work files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)
# setup work folder
workPath = pathlib.Path(avalancheDir, 'Work', 'ana1Tests', 'simiSolTest')
fU.makeADir(workPath)
# create output directory for test result plots
outDirTest = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests')
fU.makeADir(outDirTest)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for similarity solution test
simiSolCfg = cfgUtils.getModuleConfig(simiSolTest, fileOverride=fileOverride)
# ++++++++++ set configurations for all the used modules and override ++++++++++++
# get comDFA configuration and save to file
com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
                                      onlyDefault=simiSolCfg['com1DFA_override']['defaultConfig'])
com1DFACfg, simiSolCfg = cfgHandling.applyCfgOverride(com1DFACfg, simiSolCfg, com1DFA, addModValues=False)
com1DFACfgFile = cfgUtils.writeCfgFile(avalancheDir, com1DFA, com1DFACfg, fileName='com1DFA_settings',
                                       filePath=workPath)

# run DFA simulations
sphKernelRadiusList = simiSolCfg['SIMISOL']['sphKernelRadius'].split('|')
sphKernelRadiusList = [float(i) for i in sphKernelRadiusList]
for sphKernelRadius in sphKernelRadiusList:
    updater = ConfigUpdater()
    updater.read(com1DFACfgFile)
    updater['GENERAL']['sphKernelRadius'].value = sphKernelRadius
    updater['GENERAL']['meshCellSize'].value = sphKernelRadius
    updater.update_file()

    # Define release thickness distribution
    demFile = gI.getDEMPath(avalancheDir)
    relDict = simiSolTest.getReleaseThickness(avalancheDir, simiSolCfg['SIMISOL'], demFile, sphKernelRadius)
    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=False)
    _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=com1DFACfgFile)


# analyze simulations
# Load configuration info of all com1DFA simulations
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# compute the similartiy solution (this corresponds to our reference)
log.info('Computing similarity solution')
solSimi = simiSolTest.mainSimilaritySol(simiSolCfg['SIMISOL'], com1DFACfg)

# if the analysis already exists and you only want to replot uncomment this (and put the good result name)
# pathToResults = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests', 'simiSolTest', 'results20.p')
# if pathToResults.is_file():
#     simDF = pd.read_pickle(pathToResults)
simDF = simiSolTest.postProcessSimiSol(avalancheDir, cfgMain, simiSolCfg, simDF, solSimi, outDirTest)

# do some filtering before plotting
# simDF = simDF[simDF['nPPK0'] == 30]

# make convergence plot (if you add the fitting lines, make sure only the coloredBy and sizedBy parameters are varied)
fig1, ax1 = outAna1Plots.plotErrorConvergence(simDF, outDirTest, simiSolCfg['SIMISOL'], 'nPart', 'vhErrorL2',
                                              'aPPK', 'sphKernelRadius', logScale=True, fit=True)

# # make cpu time plot
# outAna1Plots.plotTimeCPULog(simDF, outDirTest, simiSolCfg['SIMISOL'], 'nPart', 'sphKernelRadius', 'sphKernelRadius')
#
# # make convergence plot (if you add the fiting lines, make sure only the coloredBy and sizedBy parameters are varied)
# # same as plotErrorConvergence but adds the points corresponding to different coloredBy values one after the others
# # and saves itermediate plots
# fig1, ax1 = outAna1Plots.plotPresentation(simDF, outDirTest, simiSolCfg['SIMISOL'], 'nPart', 'hErrorL2',
#                                           'aPPK', 'sphKernelRadius', logScale=True, fit=True)
