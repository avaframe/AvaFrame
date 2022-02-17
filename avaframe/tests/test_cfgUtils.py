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

    simFilesDF = cfgUtils.orderSimFiles(avaDir, inputDir, varParList, True, specDir='', resFiles=True)

    assert simFilesDF['simName'][0] == 'release1HS_ent_dfa_d10bdc1e81'


def test_createConfigurationInfo(tmp_path):
    """ test configuration info generation as DF """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir

    simDF = cfgUtils.createConfigurationInfo(avaDir, standardCfg='', writeCSV=False, specDir='')

    assert simDF.loc['d10bdc1e81']['releaseScenario'] == 'release1HS'
    assert simDF.loc['d10bdc1e81']['mu'] == 0.155
    assert simDF.loc['e2145362b7']['releaseScenario'] != 'release1HS'
    assert simDF.loc['e2145362b7']['relTh0'] == 1.


def test_readAllConfigurationInfo():
    """ test readAllConfigurationInfo as DF """

    dirPath = pathlib.Path(__file__).parents[0]
    configDir = dirPath / 'data' / 'testCfgUtils'

    simDF, simDFName = cfgUtils.readAllConfigurationInfo(configDir, specDir=configDir)
    print(simDF)
    print(simDFName)
    for simHash, simDFrow in simDF.iterrows():
        assert simHash in ['88cfde04b5', '1adb9373f1']


def test_appendCgf2DF(tmp_path):
    """ test appendCgf2DF """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir

    simDF = ''
    cFile = avaDir / 'Outputs' / 'com1DFA' / 'configurationFiles' / 'release1HS_ent_dfa_d10bdc1e81.ini'
    simName = pathlib.Path(cFile).stem
    simHash = 'd10bdc1e81'
    cfgObject = cfgUtils.readCfgFile(avaDir, fileName=cFile)
    simDF = cfgUtils.appendCgf2DF(simHash, simName, cfgObject, simDF)

    assert simDF.loc['d10bdc1e81']['releaseScenario'] == 'release1HS'
    assert simDF.loc['d10bdc1e81']['mu'] == '0.15500'

    cFile = avaDir / 'Outputs' / 'com1DFA' / 'configurationFiles' / 'release2HS_ent_dfa_e2145362b7.ini'
    simName = pathlib.Path(cFile).stem
    simHash = 'e2145362b7'
    cfgObject = cfgUtils.readCfgFile(avaDir, fileName=cFile)
    simDF = cfgUtils.appendCgf2DF(simHash, simName, cfgObject, simDF)

    assert simDF.loc['d10bdc1e81']['releaseScenario'] == 'release1HS'
    assert simDF.loc['d10bdc1e81']['mu'] == '0.15500'
    assert simDF.loc['e2145362b7']['releaseScenario'] != 'release1HS'
    assert simDF.loc['e2145362b7']['relTh0'] == '1.0'
    simDF = cfgUtils.convertDF2numerics(simDF)

    assert simDF.loc['d10bdc1e81']['releaseScenario'] == 'release1HS'
    assert simDF.loc['d10bdc1e81']['mu'] == 0.155
    assert simDF.loc['e2145362b7']['releaseScenario'] != 'release1HS'
    assert simDF.loc['e2145362b7']['relTh0'] == 1.


def test_filterSims(tmp_path):
    """ test filtering of simulations using configuration files """

    avaTestDir = 'avaHockeyChannelPytest'
    dirPath = pathlib.Path(__file__).parents[0]
    avaDir = dirPath / '..' / '..' / 'benchmarks' / avaTestDir
    avaDir2 = dirPath / 'data' / 'avaFilterTest' / 'com1DFA'

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

    parametersDict = {'relTh': 1.}
    simNames2 = cfgUtils.filterSims(avaDir, parametersDict, specDir='')
    simNames2 = sorted(simNames2)
    assert len(simNames2) == 2
    assert simNames2[0] == 'release1HS_ent_dfa_d10bdc1e81'
    assert simNames2[1] == 'release2HS_ent_dfa_e2145362b7'

    parametersDict = {'~relTh': 1.}
    simNames3 = cfgUtils.filterSims(avaDir, parametersDict, specDir='')

    assert len(simNames3) == 0

    parametersDict = {'~releaseScenario': 'release1HS'}
    simNames4 = cfgUtils.filterSims(avaDir, parametersDict, specDir='')

    assert len(simNames4) == 1
    assert simNames4[0] == 'release2HS_ent_dfa_e2145362b7'

    parametersDict = {'relTh': 1.}

    simNames = cfgUtils.filterSims(avaDir2, parametersDict, specDir=avaDir2)
    simNames = sorted(simNames)

    print('simNames', simNames)

    assert len(simNames) == 2
    assert simNames == ['relGar_ent_dfa_b9b17dd019', 'relGar_null_dfa_6f35cbd808']

    parametersDict = {'relTh': [1., 1.5]}

    simNames = cfgUtils.filterSims(avaDir2, parametersDict, specDir=avaDir2)
    simNames =sorted(simNames)

    assert len(simNames) == 4
    assert simNames == ['relGar_ent_dfa_b9b17dd019', 'relGar_null_dfa_1022880a70',
        'relGar_null_dfa_6f35cbd808', 'relGar_null_dfa_c4337e50ac']
