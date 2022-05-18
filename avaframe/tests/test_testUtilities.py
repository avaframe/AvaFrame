"""Tests for module testUtilities"""

from avaframe.ana1Tests import testUtilities as tU
import pathlib
import json


def test_createDesDictTemplate():
    """ test create dictionary with specified keys """

    # call function to be tested
    testDict = tU.createDesDictTemplate()

    assert testDict['TAGS'] == []
    assert {'TAGS', 'DESCRIPTION', 'TYPE', 'FILES', 'AVANAME'} <= set(testDict)


def test_writeDesDicttoJson(tmp_path):
    """ test writing to json of dict """

    # setup required inputs
    desDict = {'testkey': 'testName', 'testList': [1, 0]}
    testName = 'testDict'
    outDir = pathlib.Path(tmp_path, 'testDir')
    outDir.mkdir(exist_ok=True)

    # call function to be tested
    fileName = tU.writeDesDicttoJson(desDict, testName, outDir)

    with open(fileName) as f:
        testDict = json.load(f)

    assert fileName.is_file()
    assert testDict['testkey'] == 'testName'
    assert testDict['testList'] == [1, 0]


def test_readDesDictFromJson():
    """ test reading des dict from json file """

    # get required input data
    dirName = pathlib.Path(__file__).parents[0]
    testFile = dirName / 'data' / 'testTestUtilities_desDict.json'

    # call function to be tested
    desDict = tU.readDesDictFromJson(testFile)

    assert {'TAGS', 'DESCRIPTION', 'TYPE', 'FILES', 'AVANAME'} <= set(desDict)
    assert desDict['AVADIR'] == "data/avaKot"
    assert ('TEST' in desDict) is False


def test_readAllBenchmarkDesDicts():
    """ test reading all benchmark desDicts """

    # call function to be tested
    dirName = pathlib.Path(__file__).parents[0]
    inDir = dirName / '..' / '..' / 'benchmarks'
    testDictList = tU.readAllBenchmarkDesDicts(info=True, inDir=inDir)

    print('testDictList', testDictList)

    assert isinstance(testDictList[0], dict)


def test_filterBenchmarks():
    """ test filtering benchmark tests """

    # set up required input data
    dict1 = {'TAGS': ['null', 'standard'], 'AVANAME': 'avaTest'}
    dict2 = {'TAGS': 'entres', 'AVANAME': 'avaTest2'}
    dict3 = {'TAGS': 'res', 'AVANAME': 'avaTest4'}
    testDictList = [dict1, dict2, dict3]

    # define filtering criteria
    type = 'TAGS'
    valuesList = ['null']
    type2 = 'AVANAME'
    valuesList2 = ['avaTest2', 'avaTest3']
    type3 = 'TAGS'
    valuesList3 = ['null', 'res']

    # call function to be tested
    testList1 = tU.filterBenchmarks(testDictList, type, valuesList, condition='or')
    testList2 = tU.filterBenchmarks(testDictList, type2, valuesList2, condition='or')
    testList3 = tU.filterBenchmarks(testDictList, type3, valuesList3, condition='or')
    testList4 = tU.filterBenchmarks(testDictList, type3, valuesList3, condition='and')

    assert testList1 == [dict1]
    assert testList2 == [dict2]
    assert testList3 == [dict1, dict3]
    assert testList4 == []


def test_getTestAvaDirs():
    """ test fetching avaDirs from tests """

    # set up required input data
    dict1 = {'TAGS': ['null', 'standard'], 'AVANAME': 'avaTest'}
    dict2 = {'TAGS': 'entres', 'AVANAME': 'avaTest2'}
    dict3 = {'TAGS': 'res', 'AVANAME': 'avaTest4'}
    testList = [dict1, dict2, dict3]

    # call function to be tested
    avaDirs = tU.getTestAvaDirs(testList)

    assert 'data/avaTest' in avaDirs
    assert 'data/avaTest2' in avaDirs
    assert 'data/avaTest4' in avaDirs
    assert ('data/avaTest3' in avaDirs) is False


def test_fetchBenchmarkResults():
    """ test fetching benchmark results path """

    # setup required inputs
    testName = 'avaAlrNullTest'

    dirName = pathlib.Path(__file__).parents[0]
    refDir = dirName / '..' / '..' / 'benchmarks' / testName
    refFiles = tU.fetchBenchmarkResults(testName, resTypes=[], refDir=refDir)

    assert 'relAlr_125e697996_null' in str(refFiles[0])
    print(refFiles)
    assert len(refFiles) == 3
