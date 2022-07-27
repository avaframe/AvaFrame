"""
    Pytest for similarity solution test
    This file is part of Avaframe.
 """

#  Load modules
import pathlib
import configparser
import pytest
import numpy as np
import shutil

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in1Data import getInput as gI
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.outAna1Plots as outAna1Plots


def test_mainCompareSimSolCom1DFA(tmp_path):
    dirname = pathlib.Path(__file__).parents[0]
    simiSolCfgFile = dirname / '..' / 'tests' / 'data' / 'testSimiSol' / 'simiSol_com1DFACfg.ini'
    sourceDir = dirname / '..' / 'data' / 'avaSimilaritySol' / 'Inputs'
    destDir = tmp_path / 'avaSimilaritySol' / 'Inputs'
    avalancheDir = tmp_path / 'avaSimilaritySol'
    # setup work folder
    workPath = pathlib.Path(avalancheDir, 'Work', 'ana1Tests', 'simiSolTest')
    fU.makeADir(workPath)

    shutil.copytree(sourceDir, destDir)

    outDirTest = avalancheDir / 'Outputs' / 'ana1Tests'
    fU.makeADir(outDirTest)

    cfgMain = cfgUtils.getGeneralConfig()

    # Load configuration for similarity solution test
    simiSolCfg = cfgUtils.getModuleConfig(simiSolTest, fileOverride=simiSolCfgFile)
    # ++++++++++ set configurations for all the used modules and override ++++++++++++
    # get comDFA configuration and save to file
    com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
                                          onlyDefault=simiSolCfg['com1DFA_override']['defaultConfig'])
    com1DFACfg, simiSolCfg = cfgHandling.applyCfgOverride(com1DFACfg, simiSolCfg, com1DFA, addModValues=False)
    com1DFACfgFile = cfgUtils.writeCfgFile(avalancheDir, com1DFA, com1DFACfg, fileName='com1DFA_settings',
                                           filePath=workPath)
    # Define release thickness distribution
    demFile = gI.getDEMPath(avalancheDir)
    relDict = simiSolTest.getReleaseThickness(avalancheDir, simiSolCfg['SIMISOL'], demFile,
                                              simiSolCfg['SIMISOL'].getfloat('sphKernelRadius'))
    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=com1DFACfgFile)

    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)
    solSimi = simiSolTest.mainSimilaritySol(simiSolCfg['SIMISOL'], com1DFACfg)

    simDF = simiSolTest.postProcessSimiSol(avalancheDir, cfgMain, simiSolCfg, simDF, solSimi, outDirTest)

    # make convergence plot
    fig1, ax1 = outAna1Plots.plotErrorConvergence(simDF, outDirTest, simiSolCfg['SIMISOL'], 'nPart', 'vhErrorL2',
                                                  'aPPK', 'cMax', logScale=True)

    outAna1Plots.plotTimeCPULog(simDF, outDirTest, simiSolCfg['SIMISOL'], 'nPart', 'aPPK', 'nPPK0')

    simDF = simDF[simDF['nPPK0']==15]
    fig1, ax1 = outAna1Plots.plotPresentation(simDF, outDirTest, simiSolCfg['SIMISOL'], 'nPart', 'hErrorL2',
                              'aPPK', 'nPPK0', logScale=True, fit=True)
    simDF = simDF[simDF['aPPK']==-0.5]
    simiSolCfg['SIMISOL']['plotErrorTime'] = 'True'
    simiSolCfg['SIMISOL']['plotSequence'] = 'True'
    simiSolCfg['SIMISOL']['onlyLast'] = 'False'
    simDF = simiSolTest.postProcessSimiSol(avalancheDir, cfgMain, simiSolCfg, simDF, solSimi, outDirTest)
