"""

Function to determine statistics of datasets

"""

import os
import numpy as np
import logging
from matplotlib import pyplot as plt

from avaframe.out3Plot.plotUtils import *
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


def extractMaxValues(inputDir, cfgMain, avaDir):
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

    # initialize dictionary to save values
    peakValues = {}
    count = 0
    for simName in peakFiles['simName']:
        if peakFiles[cfgMain['PARAMETERVAR']['varPar']][count] != cfgMain['DEFVALUES'][cfgMain['PARAMETERVAR']['varPar']]:
            peakValues[simName] = {}
        count = count + 1

    # Loop through peakFiles and compute probability
    for m in range(nSims):

        # Load peak fields
        # be aware of the standard simulation - so if default value should not be part of the analysis
        if peakFiles[cfgMain['PARAMETERVAR']['varPar']][m] != cfgMain['DEFVALUES'][cfgMain['PARAMETERVAR']['varPar']]:
            log.debug('Simulation parameter %s= %s' % (cfgMain['PARAMETERVAR']['varPar'], peakFiles[cfgMain['PARAMETERVAR']['varPar']][m]))

            # Load data
            fileName = peakFiles['files'][m]
            simName = peakFiles['simName'][m]
            data = np.loadtxt(fileName, skiprows=6)

            # compute max
            max = np.amax(data)

            # add statistical measures
            peakValues[simName].update({peakFiles['resType'][m]: max})
            peakValues[simName].update({'varPar': float(peakFiles[varPar][m])})

    return peakValues
