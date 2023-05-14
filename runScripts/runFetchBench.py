"""
    Run script for running the standard tests
"""

# Load modules
import pathlib
import glob

# Local imports
from avaframe.ana1Tests import testUtilities as tU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from benchmarks import simParametersDict


# log file name; leave empty to use default runLog.log
logName = 'fetchBenchmarkTest'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False)

# filter benchmarks for a tag
type = 'TAGS'
valuesList = ['standardTest']
testList = tU.filterBenchmarks(testDictList, type, valuesList, condition='and')

for test in testList:

    # Start logging
    avaDir = test['AVADIR']
    log = logUtils.initiateLogger(avaDir, logName)

    # Fetch benchmark test results
    refFiles = tU.fetchBenchmarkResults(test['NAME'], resTypes=['ppr', 'pft', 'pfv'])
    refDir = pathlib.Path('..', 'benchmarks', test['NAME'])
    benchDict = simParametersDict.fetchBenchParameters(test['NAME'])
    simNameRef = test['simNameRef']
