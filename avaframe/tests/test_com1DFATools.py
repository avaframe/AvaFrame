"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
import logging
import pytest
import configparser
import pathlib
import os
import copy
import pickle
import pandas as pd
import shutil


from avaframe.com1DFA import com1DFATools
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


def test_extendTopBottomCom1DFAPath():
    """"""
    # setup required inputs
    cfg = configparser.ConfigParser()
    cfg['PATH'] = {'nCellsResample': '5', 'extTopOption': '1', 'nCellsMinExtend': '2',
                      'nCellsMaxExtend': '30'}
    avaProfile = {'x': np.array([10, 20, 30]), 'y': np.array([10, 20, 30]), 'z': np.array([40, 30, 20])}
    particlesIni = {'x': np.array([7, 8, 9, 10, 11, 12, 13]),
                    'y': np.array([7, 8, 9, 10, 11, 12, 13]),
                    'z': np.array([60, 42, 41, 40, 39, 38, 37])}
    dem = {'header': {'xllcenter': 0, 'yllcenter': 0, 'cellsize': 10, 'nrows': 10, 'ncols': 11},
           'rasterData': np.array([[50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0]])}

    avaProfileExt = com1DFATools.extendCom1DFAPath(cfg, dem, particlesIni, avaProfile)
    print(avaProfileExt)
    atol = 1e-10
    assert np.allclose(avaProfileExt['x'][0], 0.05050505, atol=atol)
    assert avaProfileExt['x'][-1] == 34
    assert np.allclose(avaProfileExt['y'][0], 0.05050505, atol=atol)
    assert avaProfileExt['y'][-1] == 34
    assert np.allclose(avaProfileExt['z'][0], 49.94949495, atol=atol)
    assert avaProfileExt['z'][-1] == 16
    assert np.allclose(avaProfileExt['s'], np.array([ 0., 14.0707107 , 28.21284632, 42.35498194, 48.01183619]), atol=atol)

    cfg = configparser.ConfigParser()
    cfg['PATH'] = {'nCellsResample': '5', 'extTopOption': '0', 'nCellsMinExtend': '2',
                      'nCellsMaxExtend': '30'}
    avaProfile = {'x': np.array([10, 20, 30]), 'y': np.array([10, 20, 30]), 'z': np.array([40, 30, 20])}
    particlesIni = {'x': np.array([7, 8, 9, 10, 11, 12, 13]),
                    'y': np.array([7, 8, 9, 10, 11, 12, 13]),
                    'z': np.array([43, 42, 41, 40, 39, 38, 37])}
    dem = {'header': {'xllcenter': 0, 'yllcenter': 0, 'cellsize': 10, 'nrows': 10, 'ncols': 11},
           'rasterData': np.array([[50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0],
                                   [50, 40, 30, 20, 10, 0, 0, 0, 0, 0, 0]])}

    avaProfileExt = com1DFATools.extendCom1DFAPath(cfg, dem, particlesIni, avaProfile)
    print(avaProfileExt)
    atol = 1e-10
    assert np.allclose(avaProfileExt['x'][0], 7, atol=atol)
    assert avaProfileExt['x'][-1] == 34
    assert np.allclose(avaProfileExt['y'][0], 7, atol=atol)
    assert avaProfileExt['y'][-1] == 34
    assert np.allclose(avaProfileExt['z'][0], 43, atol=atol)
    assert avaProfileExt['z'][-1] == 16
    assert np.allclose(avaProfileExt['s'], np.array([ 0., 4.24264069, 18.38477631, 32.52691193, 38.18376618]), atol=atol)
