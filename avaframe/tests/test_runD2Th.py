"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
import logging
import pytest
import configparser
import pathlib
import shutil

from avaframe.runScripts import runD2Th
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf


def test_runD2Th(tmp_path):
    """test check if sphOption and viscOption settings are consistent"""

    # setup required input data
    comMod = "com1DFA"
    resType = "pft"
    num = 6
    profileAxis = "x"
    profileIndex = 200

    dirName = pathlib.Path(__file__).parents[0]
    avaDir = pathlib.Path(tmp_path, "data", "avatestD")
    pDir = avaDir / "Outputs" / comMod / "peakFiles"
    inDir = avaDir / 'Inputs'
    fU.makeADir(pDir)
    fU.makeADir(inDir)
    demPath = dirName / ".." / "data" / "avaInclinedPlane" / "Inputs" / "DEM_IP_Topo.asc"
    demPath2 = avaDir / "Inputs" / "DEM_PF_Topo.asc"
    shutil.copy(demPath, demPath2)
#    print('DEMPATH', demPath2)

    # test depth
    testDepth = np.ones((401, 1001))
    outFile = pDir / 'release1IP_testID_C_L_null_dfa_pft.asc'
    testHeader ={'cellsize': 5., 'nrows': 401, 'ncols': 1001, 'xllcenter': 1000., 'yllcenter': -5000.,
                 'nodata_value': -9999}
    IOf.writeResultToAsc(testHeader, testDepth, outFile, flip=True)

    # call funct
    runD2Th.runD2Th(avaDir, comMod, resType, profileAxis, profileIndex)

    # read thickness
    thFile = pDir / 'transformed' / 'release1IP_testID_C_L_null_dfatransformed_pft.asc'
    thicknessField = IOf.readRaster(thFile)
    thickness = thicknessField['rasterData']

#    print('min ', np.amin(thickness), 'max ', np.amax(thickness), 'mean ', np.mean(thickness))
#    print('value', 1. * np.cos(np.deg2rad(34.)))
    assert np.isclose(np.mean(thickness), 1. * np.cos(np.deg2rad(34.)))

