"""
    test
"""

# Load modules
import os
import glob
import logging
import json

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def createDesDictTemplate():
    """ create an empty dictionary with all required test info keys """

    desDict = {'TAGS': [],
                'DESCRIPTION': '',
                'TYPE': [],
                'FILES': [],
                'AVANAME': ''}

    return desDict


def writeDesDicttoJson(desDict, testName, outDir):
    """ create a json file with all required test info from desdict """

    jsonDict = json.dumps(desDict)
    fileName = os.path.join(outDir,"%s_desDict.json" % testName)
    f = open(fileName, "w")
    f.write(jsonDict)
    f.close()

    return fileName


def readDesDictFromJson(fileName):
    """ read dict from json file """

    with open(fileName) as f:
        desDict = json.load(f)

    return desDict


def readAllBenchmarkDesDicts(info=False):
    """ get descritption dicts for all benchmark tests and add test name as key"""

    inDir = os.path.join('..', 'benchmarks')
    testDirs = glob.glob(inDir + os.sep + 'ava*')
    testDictList = []

    for testDir in testDirs:
        desDictFile = glob.glob(testDir + os.sep + '*desDict.json')
        if desDictFile != []:
            testName = os.path.basename(desDictFile[0]).split('_')[0]
            desDict = readDesDictFromJson(desDictFile[0])
            desDict.update({'NAME': testName})
            if info:
                log.info('%s: Tags: %s' % (desDict['NAME'], desDict['TAGS']))
                log.debug('%s: Full Info: %s' % (testDir, desDict))

            testDictList.append(desDict)
        else:
            testName = os.path.basename(testDir)
            log.debug('%s: No test description file provided! Test skipped' % testName)

    return testDictList


def filterBenchmarks(testDictList, type, valuesList, condition='or'):
    """ filter benchmarks according to characteristic

        Parameters
        -----------
        testDictList: list
            list of benchmark dictionaries
        type: str
            key which is used for filtering (options: TAGS, TYPE)
        valuesList: list
            list of values used for filtering
        condition: str
            if multiple values given in valuesList - if 'or' then all test dictionaries
            are returned that have either of the values, if 'and' then only those are
            returned that feature both values

        Returns
        ---------
        testList: list
            list of all the benchmark dicionaries that meet filter criterion

    """

    testList = []
    flag = False
    for testDict in testDictList:
        if type == 'AVANAME':
            if testDict['AVANAME'] == valuesList:
                testList.append(testDict)
        else:
            if condition == 'or':
                if any(values in testDict[type] for values in valuesList):
                    testList.append(testDict)
            elif condition == 'and':
                 if all(values in testDict[type] for values in valuesList):
                     testList.append(testDict)
            else:
                flag = True

    if flag:
        log.warning('Not a valid condition - please chose and or or as condition')

    return testList


def getTestAvaDirs(testList):
    """ get a list of avalanche directories of tests """

    avaDirs = []
    for tests in testList:
        avaDir = 'data' + os.sep + tests['AVANAME']
        avaDirs.append(avaDir)

    return avaDirs
