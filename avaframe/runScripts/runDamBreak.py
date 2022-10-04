"""
    Run script for running the dam break test DFA simulations on an inclined plane and compare to simulation results
"""

# Load modules
import pathlib
from configupdater import ConfigUpdater

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import logUtils
import avaframe.ana1Tests.damBreakTest as damBreakTest
import avaframe.out3Plot.outAna1Plots as outAna1Plots


# +++++++++ REQUIRED ++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runDamBreakTest'
# if left empty, use the damBreakTestCfg.ini and local_damBreakTestCfg.ini configuration files
# use 'ana1Tests/figConvergence_damBreakTestCfg.ini' to produce the similarity solution plots from the Theory Paper
fileOverride = ''
# ++++++++++++++++++++++++++++++

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaDamBreak'
cfgMain['MAIN']['avalancheDir'] = avalancheDir

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)
# setup work folder
workPath = pathlib.Path(avalancheDir, 'Work', 'ana1Tests', 'damBreakTest')
fU.makeADir(workPath)
# create output directory for test result plots
outDirTest = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests', 'damBreakTest')
fU.makeADir(outDirTest)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for dam break test
damBreakCfg = cfgUtils.getModuleConfig(damBreakTest, fileOverride=fileOverride)
# ++++++++++ set configurations for all the used modules and override ++++++++++++
# get comDFA configuration and save to file
com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
                                      onlyDefault=damBreakCfg['com1DFA_override']['defaultConfig'])
com1DFACfg, damBreakCfg = cfgHandling.applyCfgOverride(com1DFACfg, damBreakCfg, com1DFA, addModValues=False)
com1DFACfgFile = cfgUtils.writeCfgFile(avalancheDir, com1DFA, com1DFACfg, fileName='com1DFA_settings',
                                       filePath=workPath)

# run DFA simulations
sphKernelRadiusList = damBreakCfg['DAMBREAK']['sphKernelRadius'].split('|')
sphKernelRadiusList = [float(i) for i in sphKernelRadiusList]
for sphKernelRadius in sphKernelRadiusList:
    updater = ConfigUpdater()
    updater.read(com1DFACfgFile)
    updater['GENERAL']['sphKernelRadius'].value = sphKernelRadius
    updater['GENERAL']['meshCellSize'].value = sphKernelRadius
    updater.update_file()

    # call com1DFA to perform simulation - provide configuration file and release thickness function
    initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=False)
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=com1DFACfgFile)

# analyze simulations
# Load configuration info of all com1DFA simulations
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# Load flow thickness from analytical solution
solDam = damBreakTest.damBreakSol(avalancheDir, cfgMain, damBreakCfg['DAMBREAK'], com1DFACfg, outDirTest)


# if the analysis already exists and you only want to replot uncomment this (and put the good result name)
# pathToResults = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests', 'damBreakTest', 'results5.p')
# if pathToResults.is_file():
#     simDF = pd.read_pickle(pathToResults)
simDF = damBreakTest.postProcessDamBreak(avalancheDir, cfgMain, damBreakCfg, simDF, solDam, outDirTest)

# do some filtering for the presentation plot
# simDF = simDF[simDF['sphKernelRadius']==3]

# make convergence plot (if you add the fiting lines, make sure only the coloredBy and sizedBy parameters are varied)
fig1, ax1 = outAna1Plots.plotErrorConvergence(simDF, outDirTest, damBreakCfg['DAMBREAK'], 'nPart', 'vhErrorL2',
                                              'aPPK', 'sphKernelRadius', logScale=True, fit=True)

# # make convergence plot
# outAna1Plots.plotTimeCPULog(simDF, outDirTest, damBreakCfg['DAMBREAK'], 'nPart', 'aPPK', 'sphKernelRadius')
#
# # make convergence plot (if you add the fiting lines, make sure only the coloredBy and sizedBy parameters are varied)
# # same as plotErrorConvergence but adds the points corresponding to different coloredBy values one after the others
# # and saves itermediate plots
# fig1, ax1 = outAna1Plots.plotPresentation(simDF, outDirTest, damBreakCfg['DAMBREAK'], 'nPart', 'vhErrorL2',
#                                           'aPPK', 'sphKernelRadius', logScale=True, fit=True)
