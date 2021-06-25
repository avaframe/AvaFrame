"""

Function to determine statistics of datasets

"""

import os
import numpy as np
import logging
import pathlib
from matplotlib import pyplot as plt

from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def readAimecRunout(workingDir, avaName, cfg):
    """ Read runout length from aimec results

        Parameters
        ----------
        workingDir: str
            path to avalanche aimec directoy
        avaName: str
            name of avalanche directoy
        cfg : dict
            configuration read from ini file of aimec

        Returns
        --------
        Lrun: numpy array 1D
            runout length from aimec analaysis

    """

    # load configuration
    pLim = cfg['pressureLimit']
    dWidth = cfg['domainWidth']

    # set input file
    inputFileName = 'Results_%s__com1DFA__plim_%s_w_%s.txt' % (avaName, pLim, dWidth)
    inputFile = pathlib.Path(workingDir, inputFileName)
    dataset = np.loadtxt(inputFile, skiprows=7)
    Lrun = dataset[:, 3]

    return Lrun


def extractMaxValues(inputDir, cfgMain, avaDir, varPar, nameScenario='', parametersDict=''):
    """ Extract max values of result parameters and save to dictionary

        Parameters
        -----------
        inputDir: str
            path to directoy where peak files can be found
        cfgMain: dict
            configuration used to perform simulations
        avaDir: str
            path to avalanche directoy
        varPar: str
            parameter that has been varied when performing simulations (for example relTh)
        nameScenario: str
            parameter that shall be used for color coding of simulation results in plots (for example releaseScenario)
        parametersDict: dict
            dictionary with parameter and parameter values to filter simulations

        Returns
        --------
        peakValues: dict
            dictionary that contains max values for all result parameters for
            each simulation
        """

    # filter simulation results using parametersDict
    simNameList = cfgUtils.filterSims(avaDir, parametersDict)
    # load dataFrame of all simulation configurations
    simDF = cfgUtils.createConfigurationInfo(avaDir, standardCfg='')

    # load peakFiles of all simulations and generate dictionary
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)
    nSims = len(peakFilesDF['simName'])
    peakValues = {}
    for sName in simDF['simName'].tolist():
        peakValues[sName] = {}

    # Loop through peakFiles and compute probability
    for m in range(nSims):

        if peakFilesDF['simName'][m] in simNameList:

            # Load data
            fileName = peakFilesDF['files'][m]
            simName = peakFilesDF['simName'][m]
            data = np.loadtxt(fileName, skiprows=6)

            # compute max
            max = np.amax(data)

            # add statistical measures
            # fetch varPar value and nameScenario
            varParVal = simDF[simDF['simName'] == simName][varPar]
            if nameScenario != '':
                nameScenarioVal = simDF[simDF['simName'] == simName][nameScenario]
            log.info('Simulation parameter %s= %s for resType: %s and name %s' % (varPar, varParVal[0], peakFilesDF['resType'][m], nameScenarioVal[0]))
            peakValues[simName].update({peakFilesDF['resType'][m]: max})
            peakValues[simName].update({'varPar': float(varParVal)})
            peakValues[simName].update({'scenario': nameScenarioVal[0]})
        else:
            peakValues.pop(peakFilesDF['simName'][m], None)

    return peakValues
