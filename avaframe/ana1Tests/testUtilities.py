"""
    test
"""

# Load modules
import os
import glob
import pathlib
import logging
import json

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFAPy.com1DFA as com1DFA


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
            testName = os.path.basename(testDir)
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
            if testDict['AVANAME'] in valuesList:
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

def fetchBenchmarkResults(testName, resTypes=[]):
    """ Fetch a list of all paths to result files in a benchmark test,
        default all result file paths, if resType list provided only for resTypes in list

        Parameters
        -----------
        refDir: str
            path to benchmark test directory
        resTypes: list
            list of all resTypes that shall be returned

        Returns
        -------
        refFiles: list
            list of paths to all result files found for benchmark test and resType
    """


    refDir = pathlib.Path('..', 'benchmarks', testName)
    if resTypes != []:
        refFiles = []
        for resType in resTypes:
            refFiles.extend(list(refDir.glob('*' + resType + '.asc')))
    else:
        refFiles = list(refDir.glob('*.asc'))

    for refFile in refFiles:
        log.info('Benchmarktest fetched result file: %s' % refFile)

    return refFiles
