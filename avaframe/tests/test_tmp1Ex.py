"""Example test for module template"""

from avaframe.tmp1Ex import tmp1Ex
from avaframe.in3Utils import cfgUtils


def test_tmp1ExMain(capfd):
    '''Simple test for module template main'''

    cfg = cfgUtils.getModuleConfig(tmp1Ex)

    tmp1Ex.tmp1ExMain(cfg)
    out, err = capfd.readouterr()
    assert out == 'In tmp1Example\n'
