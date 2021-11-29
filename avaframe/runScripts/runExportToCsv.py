"""
    prepare data for visualisation
"""

# Load modules
import os

# Local imports
from avaframe.com1DFAPy import com1DFA
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in1Data.getInput as gI
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFA.particleTools as particleTools

# -----------Required settings-----------------
# log file name; leave empty to use default runLog.log
logName = 'runPrepareforVis'
# which particle properties shall be exported
particleProperties = 'velocityMagnitude|m'
# ---------------------------------------------

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)

# load particles and DEM Info
particlesList, TimeStepInfo = particleTools.readPartFromPickle(avalancheDir, flagAvaDir=True)
demInfo = gI.readDEM(avalancheDir)

# save particle properties to csv
outDir = os.path.join(avalancheDir, 'Outputs', 'com1DFAPy')
fU.makeADir(outDir)
com1DFA.savePartToCsv(particleProperties, particlesList, outDir)
