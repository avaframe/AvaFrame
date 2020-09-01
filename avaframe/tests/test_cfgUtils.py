"""Tests for module logUtilis"""
from avaframe.in3Utils import cfgUtils
from avaframe.tests import test_logUtils
import os
import configparser
# import logging
# import logging.config


def test_getModuleConfig(capfd):
    '''Simple test for module getModuleConfig'''
    dirname = os.path.dirname(__file__)
    avalancheDir = dirname
    cfg = cfgUtils.getModuleConfig(test_logUtils)
    assert cfg['DICT1']['parameter1']=='600'
    assert cfg['DICT1']['parameter2']=='None'
    assert cfg['DICT2']['parameter5']=='False'
