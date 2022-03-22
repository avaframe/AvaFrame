"""
    Pytest for module com1DFATools
"""

#  Load modules
import math
import pytest
import configparser


from avaframe.com1DFA import com1DFATools


def test_getPartInitMethod(capfd):
    cfg = configparser.ConfigParser()
    csz = 1
    relThForPart = 2
    cfg['GENERAL'] = {'rho': '3', 'massPerPart': '10', 'deltaTh': '0.25', 'sphKernelRadius': '1', 'nPPK0': '5',
                      'aPPK': '-1', 'sphKR0': '5', 'massPerParticleDeterminationMethod': 'MPPDIR'}
    massPerPart, nPPK = com1DFATools.getPartInitMethod(cfg['GENERAL'], csz, relThForPart)
    assert massPerPart == 10
    assert nPPK == 0

    cfg['GENERAL']['massPerParticleDeterminationMethod'] = 'MPPDH'
    massPerPart, nPPK = com1DFATools.getPartInitMethod(cfg['GENERAL'], csz, relThForPart)
    assert massPerPart == 0.75
    assert nPPK == 0

    cfg['GENERAL']['massPerParticleDeterminationMethod'] = 'MPPKR'
    massPerPart, nPPK = com1DFATools.getPartInitMethod(cfg['GENERAL'], csz, relThForPart)

    assert massPerPart == pytest.approx(math.pi * 6 / 25, abs=1e-6)
    assert nPPK == 25
