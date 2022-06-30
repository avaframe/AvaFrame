"""Tests for module com1DFAtools"""
import numpy as np
import configparser

# Local imports
import avaframe.com1DFA.timeDiscretizations as tD


def test_getSphKernelRadiusTimeStep(capfd):
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'cMax': '0.01'}
    dem = {'header': {'cellsize': 1}, 'headerNeighbourGrid': {'cellsize': 2}}

    dtStable = tD.getSphKernelRadiusTimeStep(dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.01

    cfg['GENERAL']['constrainCFL'] = 'True'
    dtStable = tD.getSphKernelRadiusTimeStep(dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.01

    dem['header']['cellsize'] = 2
    dtStable = tD.getSphKernelRadiusTimeStep(dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.02

    dem['header']['cellsize'] = 3
    dtStable = tD.getSphKernelRadiusTimeStep(dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.02

    cfg['GENERAL']['cMax'] = '0.02'
    dtStable = tD.getSphKernelRadiusTimeStep(dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.04
