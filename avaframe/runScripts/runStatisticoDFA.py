"""
Run a combination of the DFA kernel to get an alphaBeta path
to run alphaBeta model to get the alpha angle
and run the DFA kernel again
"""

import pathlib
import time
from configupdater import ConfigUpdater
import matplotlib.pyplot as plt
import numpy as np

# Local imports
# import config and init tools
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in1Data import getInput
from avaframe.ana5Hybrid import hybridTools
import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import geoTrans

# import computation modules
from avaframe.com1DFA import com1DFA
from avaframe.com2AB import com2AB

# import analysis tools
from avaframe.out3Plot import outAB
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools

# import plotting tools
from avaframe.out3Plot import outCom1DFA


# Time the whole routine
startTime = time.time()

# log file name; leave empty to use default runLog.log
logName = 'runSuperSmartANDFastModel'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# ----------------
# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

# Load configuration for hybrid model
hybridModelDFACfg = pathlib.Path('ana5Hybrid', 'hybridModel_com1DFACfg.ini')
updater = ConfigUpdater()
updater.read(hybridModelDFACfg)
updater['GENERAL']['mu'].value = str(0.6)
# updater['GENERAL']['mu'].value = str(np.tan(np.radians(20)))
updater.update_file()
# ----------------
# Run dense flow with coulomb friction
particlesList, fieldsList, tSaveList, dem, plotDict, reportDictList, simDF0 = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=hybridModelDFACfg)
# postprocess to extract path and energy line# read DEM
demOri = getInput.readDEM(avalancheDir)
avaProfilePart, avaProfileMass, avaProfileKE, _, _ = hybridTools.getCom1DFAPath(particlesList, tSaveList, demOri)
avaProfilePart, avaProfileMass, avaProfileKE = hybridTools.elongateCom1DFAPath(demOri, particlesList[0], avaProfilePart, avaProfileMass, avaProfileKE)
# save profile as AB profile in Inputs

pathAB = pathlib.Path(avalancheDir, 'Inputs', 'LINES', 'pathAB_aimec')
name = 'massAvaPath'
hybridTools.writeLine2SHPfile(avaProfileMass, name, pathAB)

# Run Alpha Beta
hybridModelABCfg = pathlib.Path('ana5Hybrid', 'hybridModel_com2ABCfg.ini')
cfgAB = cfgUtils.getModuleConfig(com2AB, fileOverride=hybridModelABCfg)
# take the path extracted from the DFA model as input
resAB = com2AB.com2ABMain(cfgAB, avalancheDir)


reportDictList = []
_, plotFile, writeFile = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

log.info('Plotted to: %s', plotFile)

# ----------------
# Run dense flow with new coulomb friction parameter
alpha = resAB[name]['alpha']
updater = ConfigUpdater()
updater.read(hybridModelDFACfg)
updater['GENERAL']['mu'].value = str(np.tan(np.radians(alpha)))
updater.update_file()

particlesList, fieldsList, tSaveList, _, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=hybridModelDFACfg)
simID = simDF.index[0]
avaProfilePartNew, avaProfileMassNew, avaProfileKENew, V2Path, EkinPath = hybridTools.getCom1DFAPath(particlesList, tSaveList, demOri)
_, avaProfileMassNewExtended, _ = hybridTools.elongateCom1DFAPath(demOri, particlesList[0], avaProfilePartNew.copy(),
                                                                  avaProfileMassNew.copy(), avaProfileKENew.copy())
# save profile as AB profile in Inputs

pathAB = pathlib.Path(avalancheDir, 'Inputs', 'LINES', 'pathAB_aimec')
name = 'massAvaPath'
hybridTools.writeLine2SHPfile(avaProfileMassNewExtended, name, pathAB)

# Run Alpha Beta
hybridModelABCfg = pathlib.Path('ana5Hybrid', 'hybridModel_com2ABCfg.ini')
cfgAB = cfgUtils.getModuleConfig(com2AB, fileOverride=hybridModelABCfg)
# take the path extracted from the DFA model as input
resABNew = com2AB.com2ABMain(cfgAB, avalancheDir)


reportDictList = []
_, plotFile, writeFile = outAB.writeABpostOut(resABNew, cfgAB, reportDictList)

log.info('Plotted to: %s', plotFile)
# Load all input Parameters from config file
# get the configuration of an already imported module
cfgAIMEC = cfgUtils.getModuleConfig(ana3AIMEC)

initProj.cleanModuleFiles(avalancheDir, ana3AIMEC)

# write configuration to file
cfgUtils.writeCfgFile(avalancheDir, ana3AIMEC, cfgAIMEC)

cfgSetup = cfgAIMEC['AIMECSETUP']
anaMod = cfgSetup['anaMod']

# Setup input from com1DFA
pathDict = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfgAIMEC)

# TODO: define referenceFile
pathDict['numSim'] = len(pathDict['pfd'])

# define reference simulation
pathDict = aimecTools.fetchReferenceSimNo(pathDict, cfgSetup)

pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)

log.info("Running ana3AIMEC model on test case DEM: \n %s \n with profile: \n %s ",
         pathDict['demSource'], pathDict['profileLayer'])
# Run AIMEC postprocessing
rasterTransfo, newRasters, resAnalysis = ana3AIMEC.mainAIMEC(pathDict, cfgAIMEC)
hybridTools.plotHybridRes(avalancheDir, resAB, resABNew, name, pathDict, simID, rasterTransfo, resAnalysis, dem, demOri,
                          particlesList, fieldsList, avaProfileMass, avaProfileMassNew, V2Path, EkinPath)
