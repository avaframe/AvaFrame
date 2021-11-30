"""
Similarity solution module

This module contains functions that compute the similarity solution
for a gliding avalanche on a inclined plane according to similarity solution from :
Hutter, K., Siegel, M., Savage, S.B. et al.
Two-dimensional spreading of a granular avalanche down an inclined plane
Part I. theory. Acta Mechanica 100, 37–68 (1993).
https://doi.org/10.1007/BF01176861
"""

# imports
import numpy as np
from scipy.integrate import ode
import math
import logging

# local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.in2Trans.ascUtils as IOf
import avaframe.ana1Tests.analysisTools as anaTools
import avaframe.out3Plot.outAna1Plots as outAna1Plots
import avaframe.in2Trans.shpConversion as shpConv


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def mainCompareSimSolCom1DFA(avalancheDir, cfgMain, simiSolCfg, outDirTest):
    """ Compute com1DFA sol, similarity solution and compare
    Main script for the comparison of DFA simulations to the analytic similarity solution
    Parameters
    -----------
    avalancheDir: str or pathlib path
        avalanche directory
    cfgMain: confiparser
        avaframeCfg configuration
    simiSolCfg: pathlib path
        path to simiSol configuration file
    outDirTest: pathlib path
        path to output directory

    """

    cfg = cfgUtils.getModuleConfig(com1DFA, simiSolCfg)

    # Define release thickness distribution
    demFile = gI.getDEMPath(avalancheDir)
    relDict = getReleaseThickness(avalancheDir, cfg, demFile)
    relTh = relDict['relTh']
    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    _, _, _, _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=simiSolCfg, relThField=relTh)

    if isinstance(simDF, str):
        simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False)

    # compute the similartiy solution (this corresponds to our reference)
    log.info('Computing similarity solution')
    solSimi = mainSimilaritySol(simiSolCfg)

    # now compare the simulations to the reference
    # first fetch info about all the simulations performed (and maybe order them)
    varParList = cfg['ANALYSIS']['varParList'].split('|')
    ascendingOrder = cfg['ANALYSIS']['ascendingOrder']
    # load info for all configurations and order them
    simDF = simDF.sort_values(by=varParList, ascending=ascendingOrder)
    simDF = postProcessSimiSol(avalancheDir, cfgMain, cfg['SIMISOL'], simDF, solSimi, outDirTest)
    outAna1Plots.plotError(simDF, outDirTest, cfg['SIMISOL'])
    outAna1Plots.plotErrorLog(simDF, outDirTest, cfg['SIMISOL'])
