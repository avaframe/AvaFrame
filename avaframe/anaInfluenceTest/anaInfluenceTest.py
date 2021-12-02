"""
Influence Test module

This module objetive is to test the influence of some parameters on the runout
area computed by the com1DFA module.
The chosen refence solution will be the one computed by com2AB module.
"""

# imports
import numpy as np
import math
import logging

# local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in1Data import getInput as gI
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.ana3AIMEC.ana3AIMEC as ana3AIMEC

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def mainAnaInfluenceTest(avalancheDir, cfgMain, influenceTestCfg, ABCfg):
    """ Compute com1DFA sol, similarity solution and compare
    Main script for the comparison of DFA simulations to the analytic similarity solution
    Parameters
    -----------
    avalancheDir: str or pathlib path
        avalanche directory
    cfgMain: confiparser
        avaframeCfg configuration
    DFACfg: pathlib path
        path to Dense-Flow-Avalanche configuration file
    ABCfg: pathlib path
        path to Alpha-Beta configuration file

    """

    cfgAIMEC = cfgUtils.getModuleConfig(ana3AIMEC)

    # Call to com2AB to compute the runout reference
    log.info('call to com2AB to compute the runout reference')
    resAB = com2AB.com2ABMain(ABCfg, avalancheDir)

    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    log.info('call com1DFA to perform simulations')
    _, _, _, _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=influenceTestCfg, relThField='', variationDict='')
    if isinstance(simDF, str):
        simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False)

    # compute the similartiy solution (this corresponds to our reference)
    log.info('Bite')
