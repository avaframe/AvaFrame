"""
    Run com1DFA kernel and compare it tothe similarity solution
    This script computes the similarity solution for a gliding avalanche on
    a inclined plane according to similarity solution from :
    Hutter, K., Siegel, M., Savage, S.B. et al.
    Two-dimensional spreading of a granular avalanche down an inclined plane
    Part I. theory. Acta Mechanica 100, 37–68 (1993).
    https://doi.org/10.1007/BF01176861
    and compares it to the DFA kernel com1DFA
"""

import pathlib
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


# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# Load general configuration
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaSimilaritySol'

# Clean input directory(ies) of old work and output files
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
simiSolCfg = cfgUtils.getModuleConfig(simiSolTest)
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
