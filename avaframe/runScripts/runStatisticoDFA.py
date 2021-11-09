"""
Run a combination of the DFA kernel to get an alphaBeta path
to run alphaBeta model to get the alpha angle
and run the DFA kernel again
"""

import pathlib
import time

# Local imports
# import config and init tools
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in1Data import getInput
from avaframe.ana5hybrid import hybridTools

# import computation modules
from avaframe.com1DFA import com1DFA
from avaframe.com2AB import com2AB

# import analysis tools
from avaframe.out3Plot import outAB


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

# ----------------
# Run dense flow with coulomb friction
particlesList, fieldssList, tSaveList, _, plotDict, reportDictList, _ = com1DFA.com1DFAMain(avalancheDir, cfgMain)

# postprocess to extract path and energy line

# read DEM
dem = getInput.readDEM(avalancheDir)
avaProfilePart, avaProfileMass, avaProfileKE = hybridTools.getCom1DFAPath(particlesList, tSaveList, dem)
avaProfilePart, avaProfileMass, avaProfileKE = hybridTools.elongateCom1DFAPath(dem, avaProfilePart, avaProfileMass, avaProfileKE)
# save profile as AB profile in Inputs

pathAB = pathlib.Path(avalancheDir, 'Inputs', 'LINES', 'pathAB_aimec')
hybridTools.writeLine2SHPfile(avaProfileMass, 'massAvePath', pathAB)

# Run Alpha Beta
cfgAB = cfgUtils.getModuleConfig(com2AB)
# take the path extracted from the DFA model as input
resAB = com2AB.com2ABMain(cfgAB, avalancheDir)


reportDictList = []
_, plotFile, writeFile = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

log.info('Plotted to: %s', plotFile)

# ----------------
# Run dense flow with new coulomb friction parameter
