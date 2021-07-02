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
from avaframe.ana1Tests import testUtilities as tU

# log file name; leave empty to use default runLog.log
logName = 'runDescriptionDict'

# set directory of benchmark test and if not yet done, this directory is created
testName = 'avaMyTest'
inDir = os.path.join('..', 'benchmarks', testName)
fU.makeADir(inDir)

# Start logging
log = logUtils.initiateLogger(inDir, logName)

# create empty description dictionary template
desDict = tU.createDesDictTemplate()

# fill this empty dictionary with test info
desDict['TAGS'] = ['null']                              # in which category does this test fall
desDict['DESCRIPTION'] = " this is my null test"        # what is this test about
desDict['TYPE'] = ["2DPeak"]                            # what type of benchmark data does the test provide
desDict['FILES'] = ["mytest1.asc", "mytest2.asc"]       # which files does the test provide
desDict['AVANAME'] = 'avaInclinedPlane'                 # which avalanache does the test refer to
desDict['AVADIR'] = 'data/InclinedPlane'                # path to the avalanche directory
desDict['simNameRef'] = 'release1IP_null_ref_0.15500' # simulation name 

# write dictionary to json file
fileName = tU.writeDesDicttoJson(desDict, testName, inDir)


# read all the benchmark info
testDictsList = tU.readAllBenchmarkDesDicts(info=False)

# list all the tests and the respective tags
for key in testDictsList:
    log.info('Test: %s Tag: %s ' % (key['NAME'], key['TAGS']))
