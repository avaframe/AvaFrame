"""Example test for module template"""

from avaframe.tmp1Ex.tmp1Ex import tmp1ExMain

def test_tmp1ExMain(capfd):
    '''Simple test for module template main'''
    tmp1ExMain()
    out, err = capfd.readouterr()
    assert out == 'In tmp1Example\n'
