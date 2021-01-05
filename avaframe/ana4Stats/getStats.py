"""

Function to determine statistics of datasets

"""

import os
import numpy as np
import logging
from matplotlib import pyplot as plt

import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf

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
    inputFile = workingDir+os.sep+ 'Results_%s__com1DFA__plim_%s_w_%s.txt' % (avaName, pLim, dWidth)
    dataset = np.loadtxt(inputFile, skiprows=7)
    Lrun = dataset[:,3]

    return Lrun


def extractMaxValues(inputDir, cfgMain, avaDir, nameScenario=''):
    """ Extract max values of result parameters and save to dictionary

        Parameters
        -----------
        inputDir: str
            path to directoy where peak files can be found
        cfgMain: dict
            configuration used to perform simulations
        avaDir: str
            path to avalanche directoy

        Returns
        --------
        peakValues: dict
            dictionary that contain max values for all result parameters for
            each simulation
        """

    # load simulation results and info
    varPar = cfgMain['PARAMETERVAR']['varPar']
    peakFiles = fU.makeSimDict(inputDir, varPar, avaDir)
    nSims = len(peakFiles['simName'])

    # load parameter variation values to check if they include default value
    itemsRaw = com1DFA.readVarPar(cfgMain)
    dVal = float(cfgMain['DEFVALUES'][varPar])
    flagValue = False
    if dVal in itemsRaw:
        flagValue = True

    # initialize dictionary to save values, if default value not in parameter variation, exclude this value
    peakValues = {}
    count = 0
    for simName in peakFiles['simName']:
        if flagValue == False:
            if peakFiles[varPar][count] != dVal:
                peakValues[simName] = {}
        else:
            peakValues[simName] = {}
        count = count + 1

    # Loop through peakFiles and compute probability
    for m in range(nSims):

        # Load peak fields
        # be aware of the standard simulation - so if default value should not be part of the analysis
        if flagValue == True:
            log.debug('Simulation parameter %s= %s' % (varPar, peakFiles[varPar][m]))

            # Load data
            fileName = peakFiles['files'][m]
            simName = peakFiles['simName'][m]
            data = np.loadtxt(fileName, skiprows=6)

            # compute max
            max = np.amax(data)

            # add statistical measures
            peakValues[simName].update({peakFiles['resType'][m]: max})
            peakValues[simName].update({'varPar': float(peakFiles[varPar][m])})
            if nameScenario != '':
                peakValues[simName].update({'scenario': nameScenario})
        else:
            if peakFiles[varPar][m] != dVal:
                log.debug('Simulation parameter %s= %s' % (varPar, peakFiles[varPar][m]))

                # Load data
                fileName = peakFiles['files'][m]
                simName = peakFiles['simName'][m]
                data = np.loadtxt(fileName, skiprows=6)

                # compute max
                max = np.amax(data)

                # add statistical measures
                peakValues[simName].update({peakFiles['resType'][m]: max})
                peakValues[simName].update({'varPar': float(peakFiles[varPar][m])})
                if nameScenario != '':
                    peakValues[simName].update({'scenario': nameScenario})

    return peakValues
