"""Tests for module cfgUtils"""

from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import pathlib
import pytest
import configparser
import sys
from pathlib import Path
import os
import multiprocessing
from types import ModuleType
import pandas as pd
import numpy as np



def test_getModuleConfig():
    '''Test for module getModuleConfig'''

    avalancheDir = pathlib.Path(__file__).parents[0]
    # test with both a default and local .ini file
    cfg = cfgUtils.getModuleConfig(test_logUtils)
    sections = cfg.sections()
    cfg.write(sys.stdout)

    assert sections == ['GENERAL', 'FLAGS', 'GOODSECTION1', 'GOODSECTION2', 'TESTSECTION1']
    assert cfg['GENERAL']['inputDir'] == 'path/to/avalanche'
    assert cfg['GENERAL']['fullLog'] == 'True'
    assert cfg['GOODSECTION1']['goodKey1'] == '1'
    assert cfg['GOODSECTION1']['goodKey2'] == 'False'
    hasKey = cfg.has_option('GOODSECTION1', 'testKey1')
    assert hasKey is True

    # test reading a different file
    filename = avalancheDir / 'local_test_logUtilsCfg.ini'
    cfg = cfgUtils.getModuleConfig(test_logUtils, fileOverride=filename)
    sections = cfg.sections()
    assert sections == ['GENERAL', 'FLAGS', 'GOODSECTION1', 'GOODSECTION2', 'TESTSECTION1']
    assert sections != ['GENERAL', 'FLAGS', 'GOODSECTION1', 'TESTSECTION1']
    assert cfg['GENERAL']['inputDir'] == 'path/to/avalanche'
    assert cfg['GENERAL']['fullLog'] == 'True'
    assert cfg['GOODSECTION1']['goodKey1'] == '1'
    assert cfg['GOODSECTION1']['goodKey2'] == 'False'
    assert cfg['GOODSECTION2']['goodKey3'] == '0'
    assert cfg['GOODSECTION2']['goodKey4'] == 'False'


def test_getGeneralConfig():
    '''Test for module getGeneralConfig'''


    # test with both a default and local .ini file
    cfg = cfgUtils.getGeneralConfig()
    sections = cfg.sections()

    assert sections == ['MAIN', 'FLAGS']
    assert cfg['MAIN']['avalancheDir'] == 'data/avaParabola'
    assert cfg['FLAGS']['savePlot'] == 'True'
    assert cfg.has_option('FLAGS', 'showPlot')
    assert cfg.has_option('FLAGS', 'reportOneFile')
    assert cfg.has_option('MAIN', 'CPUPercent')

    dirPath = pathlib.Path(__file__).parents[0]
    nameFile = dirPath / 'data' / 'testCfg.ini'
#    print(nameFile)
    # test with both a default and local .ini file
    cfg2 = cfgUtils.getGeneralConfig(nameFile=nameFile)
#    print(cfg2["MAIN"]["avalancheDir"])
    sections = cfg2.sections()

    assert sections == ['MAIN', 'FLAGS']
    assert cfg2['MAIN']['avalancheDir'] == 'data/avaHockeyChannel1'
    assert cfg2['FLAGS']['showPlot'] == 'True'
    assert cfg2.has_option('FLAGS', 'showPlot')
    assert cfg2.has_option('FLAGS', 'reportOneFile')
    assert cfg2.has_option('MAIN', 'CPUPercent')

    dirPath = pathlib.Path(__file__).parents[0]
    nameFile = dirPath / 'data' / 'test2Cfg.ini'
    # test with both a default and local .ini file
    cfg0 = cfgUtils.getGeneralConfig(nameFile=nameFile)
    sections = cfg0.sections()

    assert sections == ['MAIN', 'FLAGS']
    assert cfg0['MAIN']['avalancheDir'] == 'data/avaHockey'
    assert cfg0['FLAGS']['savePlot'] == 'False'
    assert cfg0.has_option('FLAGS', 'showPlot')
    assert cfg0.has_option('MAIN', 'CPUPercent')


def test_cfgHash():
    '''Test for the uid hash generation '''

    avalancheDir = pathlib.Path(__file__).parents[0]

    filename = avalancheDir / 'local_test_logUtilsCfg.ini'
    cfg = cfgUtils.getModuleConfig(test_logUtils, fileOverride=filename)

    uid = cfgUtils.cfgHash(cfg)

    # test for the correct uid
    assert uid == 'f9169a4908'

    # change and test again
    cfg['GOODSECTION1']['goodKey1'] = '1.5'
    uid = cfgUtils.cfgHash(cfg)
    # make sure it is not the same hash
    assert uid != '9d5f7aee27'


def test_convertConfigParserToDict():
    """ test conversion of cfg to dict """

    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'name': 'nameTest', 'relTh': '1.0'}
    cfg['TEST1'] = {'mu': '0.15500', 'id': 'name'}

    cfgDict = cfgUtils.convertConfigParserToDict(cfg)

    assert cfg['GENERAL'] == cfgDict['GENERAL']
    assert cfg['TEST1'] == cfgDict['TEST1']


def test_convertDictToConfigParser():
    """ test conversion to cfg """

    cfgDict = {}
    cfgDict['GENERAL'] = {'name': 'nameTest', 'relTh': 1.0}
    cfgDict['TEST1'] = {'mu': 0.15500, 'id': 'name'}

    cfg = cfgUtils.convertDictToConfigParser(cfgDict)

    assert cfgDict['GENERAL']['relTh'] == cfg['GENERAL'].getfloat('relTh')
    assert cfgDict['TEST1']['mu'] != cfg['TEST1']['mu']
    assert cfgDict['TEST1']['mu'] == cfg['TEST1'].getfloat('mu')


def test_createConfigurationInfo(tmp_path):
    """ test configuration info generation as DF """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir

    simDF = cfgUtils.createConfigurationInfo(avaDir, standardCfg='', writeCSV=False, specDir='')

    assert simDF.loc['0dcd58fc86']['releaseScenario'] == 'release1HS'
    assert simDF.loc['0dcd58fc86']['mu'] == 0.155
    assert simDF.loc['3d519adab0']['releaseScenario'] == 'release2HS'
    assert simDF.loc['3d519adab0']['relTh0'] == 1.0


def test_readAllConfigurationInfo():
    """ test readAllConfigurationInfo as DF """

    dirPath = pathlib.Path(__file__).parents[0]
    configDir = dirPath / 'data' / 'testCfgUtils'

    simDF, simDFName = cfgUtils.readAllConfigurationInfo(configDir, specDir=configDir)
    for simHash, simDFrow in simDF.iterrows():
        assert simHash in ['88cfde04b5', '1adb9373f1']


def test_appendCgf2DF(tmp_path):
    """ test appendCgf2DF """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir

    simDF = ''
    cFile = avaDir / 'Outputs' / 'com1DFA' / 'configurationFiles' / 'release1HS_0dcd58fc86_ent_dfa.ini'
    simName = pathlib.Path(cFile).stem
    simHash = 'd10bdc1e81'
    cfgObject = cfgUtils.readCfgFile(avaDir, fileName=cFile)
    simDF = cfgUtils.appendCgf2DF(simHash, simName, cfgObject, simDF)

    assert simDF.loc['d10bdc1e81']['releaseScenario'] == 'release1HS'
    assert simDF.loc['d10bdc1e81']['mu'] == '0.155'

    cFile = avaDir / 'Outputs' / 'com1DFA' / 'configurationFiles' / 'release2HS_3d519adab0_ent_dfa.ini'
    simName = pathlib.Path(cFile).stem
    simHash = 'e2145362b7'
    cfgObject = cfgUtils.readCfgFile(avaDir, fileName=cFile)
    simDF = cfgUtils.appendCgf2DF(simHash, simName, cfgObject, simDF)

    assert simDF.loc['d10bdc1e81']['releaseScenario'] == 'release1HS'
    assert simDF.loc['d10bdc1e81']['mu'] == '0.155'
    assert simDF.loc['e2145362b7']['releaseScenario'] != 'release1HS'
    assert simDF.loc['e2145362b7']['relTh0'] == '1.0'
    simDF = cfgUtils.convertDF2numerics(simDF)

    assert simDF.loc['d10bdc1e81']['releaseScenario'] == 'release1HS'
    assert simDF.loc['d10bdc1e81']['mu'] == 0.155
    assert simDF.loc['e2145362b7']['releaseScenario'] != 'release1HS'
    assert simDF.loc['e2145362b7']['relTh0'] == 1.0


def test_readConfigurationInfoFromDone(tmp_path):
    """ test checking for config files and which ones were perfromed and read simDF"""

    dirPath = pathlib.Path(__file__).parents[0]
    configDir = dirPath / 'data' / 'testReadConfig' / 'Outputs'/ 'com1DFA'
    configDir2 = dirPath / 'data' / 'testReadConfig'

    simDF, simNameExisting = cfgUtils.readConfigurationInfoFromDone(configDir, specDir=configDir, latest=False)

    assert len(simDF) == 9

    simDF, simNameExisting = cfgUtils.readConfigurationInfoFromDone(configDir2, specDir="", latest=True)

    assert len(simDF) == 3

def test_cfgToRcf(tmp_path):
    """Test the conversion of ConfigParser object to RCF format"""

    # Create a test ConfigParser object with various sections
    testCfg = configparser.ConfigParser()
    testCfg['GENERAL'] = {
        'parameter1': 'value1',
        'parameter2': 'value2'
    }
    testCfg['INPUT'] = {
        'input1': 'value1',
        'input2': 'value2'
    }
    testCfg['FOREST_EFFECTS'] = {
        'forest1': 'value1',
        'forest2': 'value2'
    }
    testCfg['ENTRAINMENT'] = {
        'entrainment1': 'value1',
        'entrainment2': 'value2'
    }
    testCfg['TEST_SECTION'] = {
        'test1': 'value1',
        'test2': 'value2'
    }

    # Create temporary output file
    outputFile = tmp_path / "test.rcf"

    print(outputFile)
    # Convert configuration to RCF format
    from avaframe.in3Utils.cfgUtils import cfgToRcf
    cfgToRcf(testCfg, outputFile)

    # Read the generated file
    with open(outputFile, 'r') as f:
        fileContent = f.read().splitlines()

    # Verify the content
    expectedContent = [
        'forest1                                 value1',
        'forest2                                 value2',
        '#',
        'entrainment1                            value1',
        'entrainment2                            value2',
        '#',
        '# TEST SECTION',
        '#',
        'test1                                   value1',
        'test2                                   value2',
        '#'
    ]

    # Check that GENERAL and INPUT sections are not included
    assert '# GENERAL' not in fileContent
    assert '# INPUT' not in fileContent

    # Check that FOREST_EFFECTS and ENTRAINMENT sections are handled correctly
    assert '# FOREST EFFECTS' not in fileContent
    assert '# ENTRAINMENT' not in fileContent

    # Check the formatting of the included section
    assert fileContent == expectedContent

def test_cfgToRcf_empty_config(tmp_path):
    """Test conversion of empty ConfigParser object"""
    emptyCfg = configparser.ConfigParser()
    outputFile = tmp_path / "empty.rcf"

    from avaframe.in3Utils.cfgUtils import cfgToRcf
    cfgToRcf(emptyCfg, outputFile)

    with open(outputFile, 'r') as f:
        fileContent = f.read()

    assert fileContent == ''


def test_cfgToRcf_special_characters(tmp_path):
    """Test handling of special characters in keys and values"""
    testCfg = configparser.ConfigParser()
    testCfg['TEST_SECTION'] = {
        'key with spaces': 'value with spaces',
        'key_with_underscore': 'value_with_underscore',
        'key.with.dots': 'value.with.dots'
    }

    outputFile = tmp_path / "special.rcf"

    from avaframe.in3Utils.cfgUtils import cfgToRcf
    cfgToRcf(testCfg, outputFile)

    with open(outputFile, 'r') as f:
        fileContent = f.read().splitlines()

    # Verify handling of special characters
    assert any('key with spaces' in line for line in fileContent)
    assert any('key_with_underscore' in line for line in fileContent)
    assert any('key.with.dots' in line for line in fileContent)


def test_cfgToRcf_file_handling(tmp_path):
    """Test file handling edge cases"""
    testCfg = configparser.ConfigParser()
    testCfg['TEST'] = {'key': 'value'}

    # Test with Path object
    pathObj = tmp_path / "test_path.rcf"
    from avaframe.in3Utils.cfgUtils import cfgToRcf
    cfgToRcf(testCfg, pathObj)
    assert pathObj.exists()

    # Test with string path
    strPath = str(tmp_path / "test_str.rcf")
    cfgToRcf(testCfg, strPath)
    assert os.path.exists(strPath)

def test_getNumberOfProcesses():
    """Test getting the number of CPU cores for parallel processing"""

    # Create test configuration
    testCfg = configparser.ConfigParser()
    testCfg['MAIN'] = {}
    maxCoresAvailable = multiprocessing.cpu_count()

    # Test case 1: Auto CPU with percentage
    testCfg['MAIN']['nCPU'] = 'auto'
    testCfg['MAIN']['CPUPercent'] = '50'
    numSims = 100

    from avaframe.in3Utils.cfgUtils import getNumberOfProcesses
    resultCores = getNumberOfProcesses(testCfg, numSims)
    expectedCores = int(maxCoresAvailable * 0.5)
    assert resultCores == expectedCores

    # Test case 2: Explicit CPU number
    testCfg['MAIN']['nCPU'] = '4'
    numSims = 100
    resultCores = getNumberOfProcesses(testCfg, numSims)
    assert resultCores == 4

    # Test case 3: Number of sims less than CPU cores
    testCfg['MAIN']['nCPU'] = '8'
    numSims = 3
    resultCores = getNumberOfProcesses(testCfg, numSims)
    assert resultCores == 3

    # Test case 4: Auto CPU with 100%
    testCfg['MAIN']['nCPU'] = 'auto'
    testCfg['MAIN']['CPUPercent'] = '100'
    numSims = 100
    resultCores = getNumberOfProcesses(testCfg, numSims)
    assert resultCores == maxCoresAvailable

    # Test case 5: Auto CPU with low percentage
    testCfg['MAIN']['nCPU'] = 'auto'
    testCfg['MAIN']['CPUPercent'] = '25'
    numSims = 100
    resultCores = getNumberOfProcesses(testCfg, numSims)
    expectedCores = int(maxCoresAvailable * 0.25)
    assert resultCores == expectedCores


def test_getNumberOfProcesses_edge_cases():
    """Test edge cases for getNumberOfProcesses"""

    testCfg = configparser.ConfigParser()
    testCfg['MAIN'] = {}

    # Test case 1: Zero simulations
    testCfg['MAIN']['nCPU'] = '4'
    numSims = 0

    from avaframe.in3Utils.cfgUtils import getNumberOfProcesses
    resultCores = getNumberOfProcesses(testCfg, numSims)
    assert resultCores == 0

    # Test case 2: Very high CPU number
    testCfg['MAIN']['nCPU'] = '999'
    numSims = 5
    resultCores = getNumberOfProcesses(testCfg, numSims)
    assert resultCores == 5

    # Test case 3: Auto with very high percentage
    testCfg['MAIN']['nCPU'] = 'auto'
    testCfg['MAIN']['CPUPercent'] = '200'
    numSims = 100
    resultCores = getNumberOfProcesses(testCfg, numSims)
    assert resultCores == multiprocessing.cpu_count() * 2


def test_getNumberOfProcesses_invalid_input():
    """Test invalid inputs for getNumberOfProcesses"""

    testCfg = configparser.ConfigParser()
    testCfg['MAIN'] = {}

    # Test case: Invalid nCPU value
    testCfg['MAIN']['nCPU'] = 'invalid'
    numSims = 10

    from avaframe.in3Utils.cfgUtils import getNumberOfProcesses
    with pytest.raises(ValueError):
        getNumberOfProcesses(testCfg, numSims)

import pytest


def test_convertToCfgList():
    """Test converting lists to configuration string format"""

    # Test case 1: Basic list with strings
    inputList = ['value1', 'value2', 'value3']

    from avaframe.in3Utils.cfgUtils import convertToCfgList
    resultString = convertToCfgList(inputList)
    assert resultString == 'value1|value2|value3'

    # Test case 2: Empty list
    emptyList = []
    resultString = convertToCfgList(emptyList)
    assert resultString == ''

    # Test case 3: Single item list
    singleList = ['value1']
    resultString = convertToCfgList(singleList)
    assert resultString == 'value1'

    # Test case 4: List with numeric values as strings
    numericList = ['1.0', '2.5', '3.7']
    resultString = convertToCfgList(numericList)
    assert resultString == '1.0|2.5|3.7'

    # Test case 5: List with mixed content
    mixedList = ['text', '123', 'special_chars!@#']
    resultString = convertToCfgList(mixedList)
    assert resultString == 'text|123|special_chars!@#'


def test_convertToCfgList_edge_cases():
    """Test edge cases for convertToCfgList"""

    from avaframe.in3Utils.cfgUtils import convertToCfgList

    # Test case 1: List with empty strings
    emptyStrList = ['', '', '']
    resultString = convertToCfgList(emptyStrList)
    assert resultString == '||'

    # Test case 2: List with spaces
    spacesList = ['item 1', 'item 2', 'item 3']
    resultString = convertToCfgList(spacesList)
    assert resultString == 'item 1|item 2|item 3'

    # Test case 3: List with pipe characters in items
    pipeList = ['value|1', 'value|2']
    resultString = convertToCfgList(pipeList)
    assert resultString == 'value|1|value|2'


def test_getDefaultModuleConfig(tmp_path):
    """Test getting default configuration for a module"""

    # Create a mock module
    testModule = ModuleType('testModule')
    testModule.__file__ = str(tmp_path / 'testModule.py')

    # Create default config file
    defaultConfigPath = tmp_path / 'testModuleCfg.ini'
    testConfig = configparser.ConfigParser()
    testConfig['GENERAL'] = {
        'parameter1': 'value1',
        'parameter2': 'value2'
    }
    testConfig['TEST'] = {
        'setting1': 'test1',
        'setting2': 'test2'
    }

    with open(defaultConfigPath, 'w') as configFile:
        testConfig.write(configFile)

    # Test getting default config
    from avaframe.in3Utils.cfgUtils import getDefaultModuleConfig
    resultConfig = getDefaultModuleConfig(testModule, toPrint=False)

    # Verify the configuration
    assert 'GENERAL' in resultConfig.sections()
    assert 'TEST' in resultConfig.sections()
    assert resultConfig['GENERAL']['parameter1'] == 'value1'
    assert resultConfig['GENERAL']['parameter2'] == 'value2'
    assert resultConfig['TEST']['setting1'] == 'test1'
    assert resultConfig['TEST']['setting2'] == 'test2'


def test_getDefaultModuleConfig_empty_file(tmp_path):
    """Test behavior with empty config file"""

    # Create a mock module
    testModule = ModuleType('testModule')
    testModule.__file__ = str(tmp_path / 'testModule.py')

    # Create empty config file
    defaultConfigPath = tmp_path / 'testModuleCfg.ini'
    with open(defaultConfigPath, 'w') as configFile:
        configFile.write('')

    from avaframe.in3Utils.cfgUtils import getDefaultModuleConfig
    resultConfig = getDefaultModuleConfig(testModule, toPrint=False)

    # Verify empty configuration
    assert len(resultConfig.sections()) == 0


def test_getDefaultModuleConfig_invalid_input():
    """Test handling of invalid inputs"""

    from avaframe.in3Utils.cfgUtils import getDefaultModuleConfig

    # Test with None
    with pytest.raises(AttributeError):
        getDefaultModuleConfig(None)

    # Test with non-module object
    with pytest.raises(AttributeError):
        getDefaultModuleConfig("not a module")

    # Test with module without __file__ attribute
    invalidModule = ModuleType('invalidModule')
    with pytest.raises(AttributeError):
        getDefaultModuleConfig(invalidModule)



def test_writeCfgFile(tmp_path):
    """Test writing configuration to a file with default settings"""

    # Create test avalanche directory
    avaDir = tmp_path / "avaTest"
    avaDir.mkdir()

    # Create mock module
    testModule = ModuleType('testModule')
    testModule.__file__ = str(tmp_path / 'testModule.py')

    # Create test configuration
    testConfig = configparser.ConfigParser()
    testConfig['GENERAL'] = {
        'parameter1': 'value1',
        'parameter2': 'value2'
    }

    # Write configuration file
    from avaframe.in3Utils.cfgUtils import writeCfgFile
    resultPath = writeCfgFile(avaDir, testModule, testConfig)

    # Check if file exists in correct location
    expectedPath = avaDir / "Outputs" / "testModule" / "configurationFiles" / "testModule.ini"
    assert resultPath == expectedPath
    assert resultPath.exists()

    # Verify file contents
    writtenConfig = configparser.ConfigParser()
    writtenConfig.read(resultPath)
    assert writtenConfig['GENERAL']['parameter1'] == 'value1'
    assert writtenConfig['GENERAL']['parameter2'] == 'value2'


def test_writeCfgFile_custom_path(tmp_path):
    """Test writing configuration to a custom path"""

    # Create custom directory
    customPath = tmp_path / "custom" / "path"
    customPath.mkdir(parents=True)

    # Create mock module
    testModule = ModuleType('testModule')
    testModule.__file__ = str(tmp_path / 'testModule.py')

    # Create test configuration
    testConfig = configparser.ConfigParser()
    testConfig['TEST'] = {'key': 'value'}

    # Write configuration file with custom path
    from avaframe.in3Utils.cfgUtils import writeCfgFile
    resultPath = writeCfgFile(tmp_path, testModule, testConfig, filePath=customPath)

    # Verify file location and contents
    assert resultPath == customPath / "testModule.ini"
    assert resultPath.exists()


def test_writeCfgFile_custom_filename(tmp_path):
    """Test writing configuration with custom filename"""

    # Create mock module
    testModule = ModuleType('testModule')
    testModule.__file__ = str(tmp_path / 'testModule.py')

    # Create test configuration
    testConfig = configparser.ConfigParser()
    testConfig['TEST'] = {'key': 'value'}

    # Write configuration file with custom filename
    customName = "customConfig"
    from avaframe.in3Utils.cfgUtils import writeCfgFile
    resultPath = writeCfgFile(tmp_path, testModule, testConfig, fileName=customName)

    # Verify custom filename
    expectedPath = tmp_path / "Outputs" / "testModule" / "configurationFiles" / "customConfig.ini"
    assert resultPath == expectedPath
    assert resultPath.exists()


def test_writeCfgFile_invalid_path(tmp_path):
    """Test writing configuration with invalid path"""

    # Create mock module
    testModule = ModuleType('testModule')
    testModule.__file__ = str(tmp_path / 'testModule.py')

    # Create test configuration
    testConfig = configparser.ConfigParser()
    testConfig['TEST'] = {'key': 'value'}

    # Create a file instead of a directory
    invalidPath = tmp_path / "invalid"
    invalidPath.touch()

    # Test with invalid path
    from avaframe.in3Utils.cfgUtils import writeCfgFile
    with pytest.raises(NotADirectoryError):
        writeCfgFile(tmp_path, testModule, testConfig, filePath=invalidPath)


def test_writeCfgFile_empty_config(tmp_path):
    """Test writing empty configuration"""

    # Create mock module
    testModule = ModuleType('testModule')
    testModule.__file__ = str(tmp_path / 'testModule.py')

    # Create empty configuration
    emptyConfig = configparser.ConfigParser()

    # Write empty configuration
    from avaframe.in3Utils.cfgUtils import writeCfgFile
    resultPath = writeCfgFile(tmp_path, testModule, emptyConfig)

    # Verify empty file was created
    assert resultPath.exists()
    assert resultPath.stat().st_size == 0



def test_setStrnanToNan_basic():
    """Test basic functionality of converting string 'nan' to np.nan"""

    # Create test dataframe
    testData = {
        'column1': ['nan', '1.0', 'NaN', '2.0', 'NAN'],
        'column2': ['value1', 'value2', 'value3', 'value4', 'value5']
    }
    testDF = pd.DataFrame(testData)
    testDFColumn = testDF['column1']

    # Convert string 'nan' to np.nan
    from avaframe.in3Utils.cfgUtils import setStrnanToNan
    resultDF = setStrnanToNan(testDF, testDFColumn, 'column1')

    # Check results
    assert pd.isna(resultDF.at[0, 'column1'])  # 'nan'
    assert resultDF.at[1, 'column1'] == '1.0'  # number
    assert pd.isna(resultDF.at[2, 'column1'])  # 'NaN'
    assert resultDF.at[3, 'column1'] == '2.0'  # number
    assert pd.isna(resultDF.at[4, 'column1'])  # 'NAN'

    # Check other column remains unchanged
    assert resultDF['column2'].equals(testDF['column2'])


def test_setStrnanToNan_no_nans():
    """Test behavior when no 'nan' strings are present"""

    # Create test dataframe without 'nan' values
    testData = {
        'column1': ['1.0', '2.0', '3.0', '4.0', '5.0'],
        'column2': ['a', 'b', 'c', 'd', 'e']
    }
    testDF = pd.DataFrame(testData)
    testDFColumn = testDF['column1']

    # Process dataframe
    from avaframe.in3Utils.cfgUtils import setStrnanToNan
    resultDF = setStrnanToNan(testDF, testDFColumn, 'column1')

    # Verify no changes were made
    assert resultDF.equals(testDF)


def test_setStrnanToNan_empty_df():
    """Test handling of empty dataframe"""

    # Create empty dataframe
    testDF = pd.DataFrame(columns=['column1', 'column2'])
    testDFColumn = pd.Series(dtype='object')

    # Process empty dataframe
    from avaframe.in3Utils.cfgUtils import setStrnanToNan
    resultDF = setStrnanToNan(testDF, testDFColumn, 'column1')

    # Verify result is still empty
    assert resultDF.empty
    assert len(resultDF.columns) == 2


def test_setStrnanToNan_mixed_types():
    """Test handling of mixed data types"""

    # Create test dataframe with mixed types
    testData = {
        'column1': ['nan', 1.0, 'NaN', True, None],
        'column2': ['a', 'b', 'c', 'd', 'e']
    }
    testDF = pd.DataFrame(testData)
    testDFColumn = testDF['column1'].astype(str)

    # Process dataframe
    from avaframe.in3Utils.cfgUtils import setStrnanToNan
    resultDF = setStrnanToNan(testDF, testDFColumn, 'column1')

    # Check results
    assert pd.isna(resultDF.at[0, 'column1'])  # 'nan'
    assert resultDF.at[1, 'column1'] == 1.0    # number
    assert pd.isna(resultDF.at[2, 'column1'])  # 'NaN'
    assert resultDF.at[3, 'column1'] == True   # boolean
    assert pd.isna(resultDF.at[4, 'column1'])  # None


def test_setStrnanToNan_case_sensitivity():
    """Test case-insensitive matching of 'nan' strings"""

    # Create test dataframe with different cases of 'nan'
    testData = {
        'column1': ['nan', 'NaN', 'NAN', 'nAn', 'Nan'],
        'column2': ['a', 'b', 'c', 'd', 'e']
    }
    testDF = pd.DataFrame(testData)
    testDFColumn = testDF['column1']

    # Process dataframe
    from avaframe.in3Utils.cfgUtils import setStrnanToNan
    resultDF = setStrnanToNan(testDF, testDFColumn, 'column1')

    # Verify all variations of 'nan' were converted
    for i in range(5):
        assert pd.isna(resultDF.at[i, 'column1'])