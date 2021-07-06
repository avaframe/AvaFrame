"""Tests for module cfgUtils"""

from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import pathlib
import os
import configparser


def test_getModuleConfig():
    '''Test for module getModuleConfig'''

    avalancheDir = pathlib.Path(__file__).parents[0]
    # test with both a default and local .ini file
    cfg = cfgUtils.getModuleConfig(test_logUtils)
    sections = cfg.sections()

    assert sections == ['GENERAL', 'FLAGS', 'GOODSECTION1', 'GOODSECTION2', 'BADSECTION1']
    assert cfg['GENERAL']['inputDir'] == 'path/to/avalanche'
    assert cfg['GENERAL']['fullLog'] == 'True'
    assert cfg['GOODSECTION1']['goodKey1'] == '1'
    assert cfg['GOODSECTION1']['goodKey2'] == 'False'
    hasKey = cfg.has_option('GOODSECTION1', 'badKey1')
    assert hasKey == False

    # test reading a different file
    filename = avalancheDir / 'local_test_logUtilsCfg.ini'
    cfg = cfgUtils.getModuleConfig(test_logUtils, fileOverride=filename)
    sections = cfg.sections()
    assert sections == ['GENERAL', 'FLAGS', 'GOODSECTION1', 'GOODSECTION2', 'BADSECTION1']
    assert sections != ['GENERAL', 'FLAGS', 'GOODSECTION1', 'BADSECTION1']
    assert cfg['GENERAL']['inputDir'] == 'path/to/avalanche'
    assert cfg['GENERAL']['fullLog'] == 'True'
    assert cfg['GOODSECTION1']['goodKey1'] == '1'
    assert cfg['GOODSECTION1']['goodKey2'] == 'False'
    assert cfg['GOODSECTION2']['goodKey3'] == '0'
    assert cfg['GOODSECTION2']['goodKey4'] == 'False'


def test_cfgHash():
    '''Test for the uid hash generation '''

    avalancheDir = pathlib.Path(__file__).parents[0]

    filename = avalancheDir / 'local_test_logUtilsCfg.ini'
    cfg = cfgUtils.getModuleConfig(test_logUtils, fileOverride=filename)

    uid = cfgUtils.cfgHash(cfg)

    # test for the correct uid
    assert uid == 'bcc6c69699'

    # change and test again
    cfg['GOODSECTION1']['goodKey1'] = '1.5'
    uid = cfgUtils.cfgHash(cfg)
    # make sure it is not the same hash
    assert uid != 'bcc6c69699'


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


def test_orderSimFiles():
    """ test generating order of simulation results """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir
    inputDir = avaDir / 'Outputs' / 'com1DFA' / 'peakFiles'

    varParList = 'releaseScenario'

    simFilesDF = cfgUtils.orderSimFiles(avaDir, inputDir, varParList, True, specDir='')

    assert simFilesDF['simName'][0] == 'release1HS_ent_dfa_67dc2dc10a'


def test_createConfigurationInfo(tmp_path):
    """ test configuration info generation as DF """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir

    simDF = cfgUtils.createConfigurationInfo(avaDir, standardCfg='', writeCSV=False, specDir='')

    assert simDF.loc['67dc2dc10a']['releaseScenario'] == 'release1HS'
    assert simDF.loc['67dc2dc10a']['mu'] == 0.155
    assert simDF.loc['872f0101a4']['releaseScenario'] != 'release1HS'
    assert simDF.loc['872f0101a4']['relTh'] == 1.


def test_filterSims(tmp_path):
    """ test filtering of simulations using configuration files """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir

    parametersDict = {'releaseScenario': 'release1HS'}

    simNames = cfgUtils.filterSims(avaDir, parametersDict, specDir='')

    testRel = False
    testRel2 = False
    if 'release1HS' in simNames[0]:
        testRel = True
    if 'release2HS' in simNames[0]:
        testRel2 = True

    assert testRel
    assert testRel2 == False
    assert len(simNames) == 1
