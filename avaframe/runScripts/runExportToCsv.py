"""
    prepare data for visualisation
"""

# Load modules
import pathlib

# Local imports
from avaframe.com1DFA import com1DFA
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
outDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA')
fU.makeADir(outDir)
particleTools.savePartToCsv(particleProperties, particlesList, outDir)
