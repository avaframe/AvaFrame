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
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.ana1Tests.damBreakTest as damBreak
import avaframe.out3Plot.outAna1Plots as outAna1Plots


def test_mainCompareSimSolCom1DFA(tmp_path):
    dirname = pathlib.Path(__file__).parents[0]
    damBreakCfgFile = dirname / '..' / 'tests' / 'data' / 'testDamBreak' / 'damBreak_com1DFACfg.ini'
    sourceDir = dirname / '..' / 'data' / 'avaDamBreak' / 'Inputs'
    destDir = tmp_path / 'avaDamBreak' / 'Inputs'
    avalancheDir = tmp_path / 'avaDamBreak'
    # setup work folder
    workPath = pathlib.Path(avalancheDir, 'Work', 'ana1Tests', 'damBreakTest')
    fU.makeADir(workPath)

    shutil.copytree(sourceDir, destDir)

    outDirTest = avalancheDir / 'Outputs' / 'ana1Tests'
    fU.makeADir(outDirTest)

    cfgMain = cfgUtils.getGeneralConfig()
    damBreakCfg = cfgUtils.getModuleConfig(damBreak, fileOverride=damBreakCfgFile)
    # ++++++++++ set configurations for all the used modules and override ++++++++++++
    # get comDFA configuration and save to file
    com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
                                          onlyDefault=damBreakCfg['com1DFA_override']['defaultConfig'])
    com1DFACfg, damBreakCfg = cfgHandling.applyCfgOverride(com1DFACfg, damBreakCfg, com1DFA, addModValues=False)
    com1DFACfgFile = cfgUtils.writeCfgFile(avalancheDir, com1DFA, com1DFACfg, fileName='com1DFA_settings',
                                           filePath=workPath)
    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=com1DFACfgFile)

    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

    solDam = damBreak.damBreakSol(avalancheDir, cfgMain, damBreakCfg['DAMBREAK'], com1DFACfg, outDirTest)

    simDF = damBreak.postProcessDamBreak(avalancheDir, cfgMain, damBreakCfg, simDF, solDam, outDirTest)

    # make convergence plot
    fig1, ax1 = outAna1Plots.plotErrorConvergence(simDF, outDirTest, damBreakCfg['DAMBREAK'], 'nPart', 'hErrorL2',
                              'aPPK', 'nPPK0', logScale=True)

    outAna1Plots.plotTimeCPULog(simDF, outDirTest, damBreakCfg['DAMBREAK'], 'nPart', 'aPPK', 'nPPK0')

    simDF = simDF[simDF['nPPK0']==15]
    fig1, ax1 = outAna1Plots.plotPresentation(simDF, outDirTest, damBreakCfg['DAMBREAK'], 'nPart', 'hErrorL2',
                              'aPPK', 'nPPK0', logScale=True, fit=True)
    simDF = simDF[simDF['aPPK']==-0.5]
    damBreakCfg['DAMBREAK']['plotErrorTime'] = 'True'
    damBreakCfg['DAMBREAK']['plotSequence'] = 'True'
    damBreakCfg['DAMBREAK']['onlyLast'] = 'False'
    simDF = damBreak.postProcessDamBreak(avalancheDir, cfgMain, damBreakCfg, simDF, solDam, outDirTest)
