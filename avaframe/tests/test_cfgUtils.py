"""Tests for module cfgUtils"""

from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import pathlib
import pytest
import configparser
import sys


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