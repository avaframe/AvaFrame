"""Tests for module logUtilis"""
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import os


def test_getModuleConfig(capfd):
    '''Simple test for module getModuleConfig'''
    dirname = os.path.dirname(__file__)
    avalancheDir = dirname
    # test with both a default and local .ini file
    cfg = cfgUtils.getModuleConfig(test_logUtils)
    sections = cfg.sections()

    assert sections == ['GENERAL', 'FLAGS', 'GOODSECTION1', 'GOODSECTION2']
    assert cfg['GENERAL']['inputDir'] == 'path/to/avalanche'
    assert cfg['GENERAL']['fullLog'] == 'True'
    assert cfg['GOODSECTION1']['goodKey1'] == '1'
    assert cfg['GOODSECTION1']['goodKey2'] == 'False'
    hasKey = cfg.has_option('GOODSECTION1', 'badKey1')
    assert hasKey == False

    # test reading a different file
    filemane = os.path.join(avalancheDir,'local_test_logUtilsCfg.ini')
    cfg = cfgUtils.getModuleConfig(test_logUtils, fileOverride=filemane)
    sections = cfg.sections()
    assert sections == ['GENERAL', 'FLAGS', 'GOODSECTION1', 'BADSECTION1']
    assert cfg['GENERAL']['inputDir'] == 'path/to/avalanche'
    assert cfg['GENERAL']['fullLog'] == 'True'
    assert cfg['GOODSECTION1']['goodKey1'] == '1'
    assert cfg['GOODSECTION1']['goodKey2'] == 'False'
    assert cfg['GOODSECTION1']['badKey1'] == 'Bonjour'
    assert cfg['BADSECTION1']['badKey2'] == '10'
    assert cfg['BADSECTION1']['badKey3'] == 'None'
