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
from avaframe.out3Plot import outAB
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.com2AB.com2AB as com2AB
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

    # Analyse/ plot/ write results #
    reportDictList = []
    _, plotFile, writeFile = outAB.writeABpostOut(resAB, ABCfg, reportDictList)

    log.info('Plotted to: %s' % [str(plotFileName) for plotFileName in plotFile])
    log.info('Data written: %s' % [str(writeFileName) for writeFileName in writeFile])

    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    log.info('call com1DFA to perform simulations')
    _, _, _, _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=influenceTestCfg, relThField='', variationDict='')
    if isinstance(simDF, str):
        simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False)

    # Get the alpha-beta runout
#    abRunout = getABRunout()
    # Get each DFA runout
#    DFARunouts = getDFARunout()
    # Plot and compare and analyze the results
#    analyzeResults(abRunout, DFARunouts)

    # Check if last line is computed
    log.info('Main Function well computed')

def getABRunout(resAB):
    """ Get alpha-beta runout
    Parameters
    ----------
    resAB : dict
        dict with com2AB results
    Returns
    -------
    abRunout: float
    """
    AvaPath = resAB['AvaPath']
    NameAva = AvaPath['Name']
    print('STOOOOOOOOOOOOOOOOOOOOOOP')
    print(len(NameAva))
#    for i in range(len(NameAva)):
#        name = NameAva[i]
#        abRunout = resAB[name]['x']
#    return abRunout

def getDFARunout():
    """ Get DFA runout
    Parameters
    ----------
    ...
    Returns
    -------
    DFARunouts: dictionary
        dfaRunout: float array
    """
    abRunout = 0
    return DFARunouts

def analyzeResults(abRunout, DFARunout):
    """ Plot different runouts and store results in a data frame
    Parameters
    ----------
    abRunout
    DFARunout
    Returns
    -------
        dfaRunout: float array
    """
