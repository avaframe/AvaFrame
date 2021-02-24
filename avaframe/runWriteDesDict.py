"""
    Run script for generating a description json file for a test
"""

# Load modules
import os
import time
import glob

# Local imports
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.ana1Tests import standardTests as sT

# log file name; leave empty to use default runLog.log
logName = 'runDescriptionDict'

# set directory of benchmarks
avaTest = 'avaMyTest'
testName = 'nullTest'
inDir = os.path.join('..', 'benchmarks', avaTest)
fU.makeADir(inDir)

# Start logging
log = logUtils.initiateLogger(inDir, logName)

# create empty description dictionary template
desDict = sT.createDesDictTemplate()

# fill this empty dictionary with test info
desDict['TAGS'] = ['null']
desDict['DESCRIPTION'] = " this is my test"
desDict['TYPE'] = "rasterFile"
desDict['FILES'] = ["release1HS_entres_ref_0.15500_pfv.asc", "release1HS_entres_ref_0.15500_ppr.asc"]

# write dictionary to json file
fileName = sT.writeDesDicttoJson(desDict, testName, inDir)

# read dictionary from json file
testDict = sT.readDesDictFromJson(fileName)

for key in testDict:
    log.info('%s: %s' % (key, testDict[key]))
