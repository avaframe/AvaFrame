"""Tests for module cfgUtils"""

from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import os


def test_getModuleConfig():
    '''Test for module getModuleConfig'''

    avalancheDir = os.path.dirname(__file__)
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
    filename = os.path.join(avalancheDir, 'local_test_logUtilsCfg.ini')
    cfg = cfgUtils.getModuleConfig(test_logUtils, fileOverride=filename)
    sections = cfg.sections()
    assert sections == ['GENERAL', 'FLAGS', 'GOODSECTION1', 'BADSECTION1']
    assert cfg['GENERAL']['inputDir'] == 'path/to/avalanche'
    assert cfg['GENERAL']['fullLog'] == 'True'
    assert cfg['GOODSECTION1']['goodKey1'] == '1'
    assert cfg['GOODSECTION1']['goodKey2'] == 'False'
    assert cfg['GOODSECTION1']['badKey1'] == 'Bonjour'
    assert cfg['BADSECTION1']['badKey2'] == '10'
    assert cfg['BADSECTION1']['badKey3'] == 'None'


def test_cfgHash():
    '''Test for the uid hash generation '''

    avalancheDir = os.path.dirname(__file__)

    filename = os.path.join(avalancheDir, 'local_test_logUtilsCfg.ini')
    cfg = cfgUtils.getModuleConfig(test_logUtils, fileOverride=filename)

    uid = cfgUtils.cfgHash(cfg)

    # test for the correct uid
    assert uid == 'aec6a993be'

    # change and test again
    cfg['GOODSECTION1']['goodKey1'] = '1.5'
    uid = cfgUtils.cfgHash(cfg)
    # make sure it is not the same hash
    assert uid != 'aec6a993be'
