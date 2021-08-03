"""Tests for module com1DFAtools"""
import numpy as np
import configparser

# Local imports
import avaframe.com1DFA.timeDiscretizations as tD


def test_getcflTimeStep(capfd):
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'cMax': '0.5', 'maxdT': '0.5', 'constrainCFL': 'False', 'mindT': '0.01'}
    dem = {'header': {'cellsize': 1}, 'headerNeighbourGrid': {'cellsize': 2}}
    particles = {}
    particles['ux'] = np.array([100.])
    particles['uy'] = np.array([0.])
    particles['uz'] = np.array([0.])

    dtStable = tD.getcflTimeStep(particles, dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.005

    cfg['GENERAL']['constrainCFL'] = 'True'
    dtStable = tD.getcflTimeStep(particles, dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.01

    particles['ux'] = np.array([0.5])
    dtStable = tD.getcflTimeStep(particles, dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.5

    particles['ux'] = np.array([5])
    dtStable = tD.getcflTimeStep(particles, dem, cfg['GENERAL'])
    print(dtStable)
    assert dtStable == 0.1
