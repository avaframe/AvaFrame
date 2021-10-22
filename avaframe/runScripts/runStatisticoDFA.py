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
_, _, _, _, plotDict, reportDictList, _ = com1DFA.com1DFAMain(avalancheDir, cfgMain)

# Get peakfiles to return to QGIS
avaDir = pathlib.Path(avalancheDir)
inputDir = avaDir / 'Outputs' / 'com1DFA' / 'peakFiles'
peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

# Run Alpha Beta
cfgAB = cfgUtils.getModuleConfig(com2AB)
resAB = com2AB.com2ABMain(cfgAB, avalancheDir)
