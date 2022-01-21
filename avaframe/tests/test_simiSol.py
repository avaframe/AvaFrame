"""
    Pytest for module in1Data
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
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in1Data import getInput as gI
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.outAna1Plots as outAna1Plots


def test_mainCompareSimSolCom1DFA(tmp_path):
    dirname = pathlib.Path(__file__).parents[0]
    simiSolCfg = dirname / '..' / 'tests' / 'data' / 'testSimiSol' / 'simiSol_com1DFACfg.ini'
    sourceDir = dirname / '..' / 'data' / 'avaSimilaritySol'
    avalancheDir = tmp_path / 'avaSimilaritySol'

    shutil.copytree(sourceDir, avalancheDir)

    outDirTest = avalancheDir / 'Outputs' / 'ana1Tests'
    fU.makeADir(outDirTest)

    cfgMain = cfgUtils.getGeneralConfig()
    cfg = cfgUtils.getModuleConfig(com1DFA, simiSolCfg)
    # Define release thickness distribution
    demFile = gI.getDEMPath(avalancheDir)
    relDict = simiSolTest.getReleaseThickness(avalancheDir, cfg, demFile)
    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=simiSolCfg)

    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)
    solSimi = simiSolTest.mainSimilaritySol(simiSolCfg)

    simDF = simiSolTest.postProcessSimiSol(avalancheDir, cfgMain, cfg['SIMISOL'], simDF, solSimi, outDirTest)

    outAna1Plots.plotErrorRef(simDF, outDirTest, cfg['SIMISOL'], 'subgridMixingFactor', ['hErrorL2', 'vhErrorL2'],
                              'aPPK', 'cMax', logScale=False)

    # make convergence plot
    fig1, ax1, ax2, slopeU, slopeH = outAna1Plots.plotErrorConvergence(simDF, outDirTest, cfg['SIMISOL'], 'nPart', ['hErrorL2', 'vhErrorL2'],
                              'aPPK', 'cMax', logScale=True)

    outAna1Plots.plotTimeCPULog(simDF, outDirTest, cfg['SIMISOL'], 'nPart', 'aPPK', 'nPPK0')

    simDF = simDF[simDF['nPPK0']==15]
    simDF = simDF[simDF['aPPK']==-0.5]
    cfg['SIMISOL']['plotIntermediate'] = 'True'
    simDF = simiSolTest.postProcessSimiSol(avalancheDir, cfgMain, cfg['SIMISOL'], simDF, solSimi, outDirTest)
