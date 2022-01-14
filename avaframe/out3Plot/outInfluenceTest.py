"""
    Plotting and saving Influence Test results
"""

import logging
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Local imports
from avaframe.in3Utils import fileHandlerUtils as fU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)

def analyzeResults(simDF, avalancheDir, IDList, refID, resInfluenceTest):
    """ generate Lists to plot results
    Parameters
    ----------
    simDF, avalancheDir, samosIDList, ataIDList, refRunout, param2Study, visc2Study, otherParamValueList, otherParamNameList
        all described in previous functions
    Returns
    -------
    runoutArray: float List
        list of every Runnout
    MMPPRArray: float List
        list of every MMPPR
    MMPFDArray: float List
        list of every graph MMPFD
    """

    otherParamNameList = resInfluenceTest['otherParamNameList']
    otherParamValueList = resInfluenceTest['otherParamValueList']
    param2Study = resInfluenceTest['param2Study']
    xAxisArray = resInfluenceTest['xAxisArray']

    runoutArray = []
    MMPPRArray = []
    MMPFDArray = []

    for iN in range(len(otherParamNameList)):
        paramN = otherParamNameList[iN]
        for iV in range(len(otherParamValueList[iN])):
            paramV = otherParamValueList[iN][iV]
            runoutArray.append([])
            MMPPRArray.append([])
            MMPFDArray.append([])
            for xCoord in xAxisArray:
                for id in IDList:
                    if simDF.loc[id][paramN] == float(paramV) and simDF.loc[id][param2Study] == xCoord:
                        runoutArray[iV].append(simDF.loc[id]['runout'])
                        MMPPRArray[iV].append(simDF.loc[id]['MMPPR'])
                        MMPFDArray[iV].append(simDF.loc[id]['MMPFD'])

    refArray = [[], [], []]
    for xCoord in xAxisArray:
        refArray[0].append(simDF.loc[refID]['runout'])
        refArray[1].append(simDF.loc[refID]['MMPPR'])
        refArray[2].append(simDF.loc[refID]['MMPFD'])


    resInfluenceTest['runoutArray'] = runoutArray
    resInfluenceTest['MMPPRArray'] = MMPPRArray
    resInfluenceTest['MMPFDArray'] = MMPFDArray
    resInfluenceTest['refArray'] = refArray

    return resInfluenceTest



def plotResults(avalancheDir, resInfluenceTest):
    """ generate Plots contained in xAxisArray, runoutArray, MMPPRArray, MMPFDArray and refArray lists
    Parameters
    ----------
    param2Study, visc2Study, otherParamNameList, otherParamValueList, xAxisArray, runoutArray, MMPPRArray, MMPFDArray
        all described in previous functions
    """

    # If needed, create folders to store the simResults
    outDir = avalancheDir +'/Outputs/ana1Test/influenceTest'
    log.info(outDir)
    fU.makeADir(outDir)

    otherParamNameList = resInfluenceTest['otherParamNameList']
    otherParamValueList = resInfluenceTest['otherParamValueList']
    param2Study = resInfluenceTest['param2Study']
    visc2Study = resInfluenceTest['visc2Study']
    xAxisArray = resInfluenceTest['xAxisArray']
    runoutArray = resInfluenceTest['runoutArray']
    MMPPRArray = resInfluenceTest['MMPPRArray']
    MMPFDArray = resInfluenceTest['MMPFDArray']
    refArray = resInfluenceTest['refArray']

    for iN in range(len(otherParamNameList)):
        paramN = otherParamNameList[iN]
        plt.xlabel(param2Study)

        plt.ylabel('runout')
        for iV in range(len(otherParamValueList[iN])):
            paramV = otherParamValueList[iN][iV]
            plt.plot(xAxisArray,  runoutArray[iV], label=visc2Study + ' viscosity, ' + paramN + ' = ' + paramV)
        plt.plot(xAxisArray,  refArray[0], label='ref')
        plt.legend()
        plt.savefig(outDir + '/Runout_' + visc2Study + 'Viscosity_param2Study_' + param2Study + '_' + paramN + '_' + paramV + '.png')
        plt.show()

        plt.xlabel(param2Study)
        plt.ylabel('MMPPR')
        for iV in range(len(otherParamValueList[iN])):
            paramV = otherParamValueList[iN][iV]
            plt.plot(xAxisArray,  MMPPRArray[iV], label=visc2Study + ' viscosity, ' + paramN + ' = ' + paramV)
        plt.plot(xAxisArray,  refArray[1], label='ref')
        plt.legend()
        plt.savefig(outDir + '/MMPPR' + visc2Study + 'Viscosity_param2Study_' + param2Study + '_' + paramN + '_' + paramV + '.png')
        plt.show()

        plt.xlabel(param2Study)
        plt.ylabel('MMPFD')
        for iV in range(len(otherParamValueList[iN])):
            paramV = otherParamValueList[iN][iV]
            plt.plot(xAxisArray,  MMPFDArray[iV], label=visc2Study + ' viscosity, ' + paramN + ' = ' + paramV)
        plt.plot(xAxisArray,  refArray[2], label='ref')
        plt.legend()
        plt.savefig(outDir + '/MMPFD' + visc2Study + 'Viscosity_param2Study_' + param2Study + '_' + paramN + '_' + paramV + '.png')
        plt.show()
