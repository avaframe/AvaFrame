"""
    Pytest for module com1DFATools
"""

#  Load modules
import math
import pytest
import configparser
import pathlib
import shutil

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


def test_createSimDictFromCfgs(tmp_path):
    """ test creating a simDict from multiple cfg files """

    dirPath = pathlib.Path(__file__).parents[0]
    testPath = dirPath / 'data' / 'com1DFAConfigsTest'
    inputDir = dirPath / 'data' / 'testCom1DFA2'
    avaDir = pathlib.Path(tmp_path, 'testCom1DFA')
    shutil.copytree(inputDir, avaDir)
    cfgMain = configparser.ConfigParser()
    cfgMain['MAIN'] = {'avalancheDir': avaDir}

    with pytest.raises(FileNotFoundError) as e:
        # call function to be tested
        simDict, inputSimFiles, simDFExisting, outDir = com1DFATools.createSimDictFromCfgs(cfgMain, testPath)
    assert ("No configuration file found to create simulation runs in:") in str(e.value)


    testPath = dirPath / 'data' / 'com1DFAConfigs'
    avaDir = pathlib.Path(tmp_path, 'testCom1DFA2')
    shutil.copytree(inputDir, avaDir)
    cfgMain['MAIN'] = {'avalancheDir': avaDir}
    # call function to be tested
    simDict, inputSimFiles, simDFExisting, outDir = com1DFATools.createSimDictFromCfgs(cfgMain, testPath)

    for sim in simDict.keys():
        print('simDict: ', sim, simDict[sim])
    print('inputSimFiles', inputSimFiles)
    print('simDFExisting', simDFExisting)
    print('outDir', outDir)

    assert simDFExisting == None
    assert len(simDict) == 16
