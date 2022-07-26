"""
Run the rotation test
Analyze the effect of the grid orientation/alignment on DFA simulation results
"""
import pathlib

# Local imports
# import config and init tools
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import initializeProject as iP
# import computation modules
from avaframe.com1DFA import com1DFA

# import analysis modules
from avaframe.runScripts import runAna3AIMEC
from avaframe.ana3AIMEC import ana3AIMEC, aimecTools
from avaframe.ana1Tests import rotationTest

# +++++++++REQUIRED+++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runRotationTest'
# for the report
comModule = 'com1DFA'
# get path to com1DFA configuration file used for the rotation test
rotationTestCfgFile = ''
# if you want to specify a different file and disregard the local one
# rotationTestCfgFile = pathlib.Path('ana1Tests', 'rotationTest_com1DFACfg.ini')
# do you want to run the DFA module (all results in the Outputs/com1DFA folder will be deleted)
runDFAModule = True
# for aimec analysis (overwrites the parameters in the config files, default and local)
anaMod = 'com1DFARotated'
referenceSimName = 'rel0'
# start with this, will be changed to True if the simulations have entrainment
flagMass = False
# ++++++++++++++++++++++++++++++

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']
avalancheDir = pathlib.Path(avalancheDir)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# ----------------
# Clean input directory(ies) of old work and output files
iP.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

# run the com1DFA module or load the results from com1DFA
dem, simDF, resTypeList = com1DFA.runOrLoadCom1DFA(avalancheDir, cfgMain, runDFAModule=runDFAModule,
                                                   cfgFile=rotationTestCfgFile)
# fetch reference simulation based in referenceSimName
refSimRowHash, refSimName = aimecTools.defineRefOnSimName(referenceSimName, simDF)
# initialize report
reportRotationTest = rotationTest.initializeRotationTestReport(avalancheDir, resTypeList, comModule, refSimName)
# Rotate raster results for aimec analysis and run the energy line test
simDF, flagMass = rotationTest.mainRotationTest(avalancheDir, dem, simDF, resTypeList, flagMass, refSimRowHash)

# proceed to aimec analysis on the rotated raster results
iP.cleanModuleFiles(avalancheDir, ana3AIMEC)
# prepare the configuration
cfgAimec = cfgUtils.getModuleConfig(ana3AIMEC)
cfgAimec['AIMECSETUP']['resTypeList'] = 'ppr|pfv|pft'
cfgAimec['AIMECSETUP']['varParList'] = ''
cfgAimec['AIMECSETUP']['referenceSimName'] = referenceSimName
cfgAimec['AIMECSETUP']['anaMod'] = anaMod
cfgAimec['FLAGS']['flagMass'] = str(flagMass)
# write configuration to file
cfgUtils.writeCfgFile(avalancheDir, ana3AIMEC, cfgAimec)

pathDict, rasterTransfo, resAnalysisDF, aimecPlotDict = runAna3AIMEC.runAna3AIMEC(avalancheDir, cfgAimec)
resAnalysisDF = resAnalysisDF.reset_index().merge(simDF[['simName', 'relAngle']], on=['simName']).set_index('index')

# finalize rotation test report
reportRotationTest = rotationTest.buildRotationTestReport(avalancheDir, reportRotationTest, simDF, resAnalysisDF,
                                                          aimecPlotDict, flagMass)
