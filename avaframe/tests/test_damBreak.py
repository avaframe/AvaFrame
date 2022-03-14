"""
    Pytest for dam break test
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
import avaframe.ana1Tests.damBreak as damBreak
import avaframe.out3Plot.outAna1Plots as outAna1Plots


def test_mainCompareSimSolCom1DFA(tmp_path):
    dirname = pathlib.Path(__file__).parents[0]
    damBreakCfg = dirname / '..' / 'tests' / 'data' / 'testDamBreak' / 'damBreak_com1DFACfg.ini'
    sourceDir = dirname / '..' / 'data' / 'avaDamBreak' / 'Inputs'
    destDir = tmp_path / 'avaDamBreak' / 'Inputs'
    avalancheDir = tmp_path / 'avaDamBreak'

    shutil.copytree(sourceDir, destDir)

    outDirTest = avalancheDir / 'Outputs' / 'ana1Tests'
    fU.makeADir(outDirTest)

    cfgMain = cfgUtils.getGeneralConfig()
    cfg = cfgUtils.getModuleConfig(com1DFA, damBreakCfg)
    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=damBreakCfg)

    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

    solDam = damBreak.damBreakSol(avalancheDir, cfgMain, cfg, outDirTest)

    simDF = damBreak.postProcessDamBreak(avalancheDir, cfgMain, cfg, simDF, solDam, outDirTest)

    # make convergence plot
    fig1, ax1 = outAna1Plots.plotErrorConvergence(simDF, outDirTest, cfg['DAMBREAK'], 'nPart', 'hErrorL2',
                              'aPPK', 'nPPK0', logScale=True)

    outAna1Plots.plotTimeCPULog(simDF, outDirTest, cfg['DAMBREAK'], 'nPart', 'aPPK', 'nPPK0')

    simDF = simDF[simDF['nPPK0']==15]
    fig1, ax1 = outAna1Plots.plotPresentation(simDF, outDirTest, cfg['DAMBREAK'], 'nPart', 'hErrorL2',
                              'aPPK', 'nPPK0', logScale=True, fit=True)
    simDF = simDF[simDF['aPPK']==-0.5]
    cfg['DAMBREAK']['plotErrorTime'] = 'True'
    cfg['DAMBREAK']['plotSequence'] = 'True'
    cfg['DAMBREAK']['onlyLast'] = 'False'
    simDF = damBreak.postProcessDamBreak(avalancheDir, cfgMain, cfg, simDF, solDam, outDirTest)
