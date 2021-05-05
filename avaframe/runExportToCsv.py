"""
    prepare data for visualisation
"""

# Load modules
import time
import os

# Local imports
from avaframe.com1DFAPy import com1DFA
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in1Data.getInput as gI
from avaframe.log2Report import generateReport as gR
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# -----------Required settings-----------------
# log file name; leave empty to use default runLog.log
logName = 'runPrepareforVis'

# ---------------------------------------------
# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)

# load particles and DEM Info
particlesList, TimeStepInfo = com1DFA.readPartFromPickle(avalancheDir, flagAvaDir=True)
demInfo = gI.readDEM(avalancheDir)

# save particle properties to csv
particleProperties = 'velocityMagnitude|m'
outDir = os.path.join(avalancheDir, 'Outputs', 'com1DFAPy')
fU.makeADir(outDir)
com1DFA.savePartToCsv(particleProperties, particlesList, outDir)
