"""
    functions to perform sampling of the parameter space
"""

# Load modules
import logging
import configparser
import numpy as np
from scipy.stats import qmc
import pathlib

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import logUtils

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def createSample(cfgSA):
    """ create a sample for parameter variation

        Parameters
        -----------
        cfgSample: configparser object
            here used:
            nSample: int
                factor used to create the number of samples following, N*(2D+2) where D is the amount of parameters
            sampleMethod: str
                name of sample method available in SALib
                    - saltelli
            seed: int
                integer to initialise the random gnr.
            varParMin, varParMax, varParList

        Returns
        ---------
        paramValuesD: dict
            paramNames: list
                list of parameter names of paramValues
            paramValues: numpy array
                array of all sets of parameter values (N* (2D+2), one row per value set

    """

    # random generator object with seed
    seed = np.random.default_rng(cfgSA.getint('sampleSeed'))

    # create parameter variation info
    parameterNames = cfgSA['varParList'].split('|')
    lowerBoundsStr = cfgSA['varParMin'].split('|')
    upperBoundsStr = cfgSA['varParMax'].split('|')
    lowerBounds = [float(item) for item in lowerBoundsStr]
    upperBounds = [float(item) for item in upperBoundsStr]

    log.info('Parameter to be varied are: %s' % cfgSA['varParList'])
    log.info('Min values of these parameters are: %s' % cfgSA['varParMin'])
    log.info('Max values of these parameters are: %s' % cfgSA['varParMax'])

    if cfgSA['sampleMethod'].lower() == 'latin':
        # create sample using latin hypercube sampling from scipy
        sampler = qmc.LatinHypercube(d=len(parameterNames))
        sample = sampler.random(n=cfgSA.getint('nSample'))
        sample = qmc.scale(sample, lowerBounds, upperBounds)
        log.info('Parameter sample created using latin hypercube sampling')

    # create dictionary with all the info
    paramValuesD = {'names': parameterNames, 'values': sample, 'typeList': cfgSA['varParType'].split('|')}

    return paramValuesD


def createCfgFiles(paramValuesD, comMod, cfg, cfgPath='', fileOverride=''):
    """ create all config files required to run com Module from parameter variations using paramValues

        Parameters
        -----------
        paramValuesD: dict
            dictionary with parameter names and values (array of all sets of parameter values, one row per value set)
        comMod: com module
            computational module
        cfgPath: str
            path where cfg files should be saved to
        fileOverride: pathlib path
            if config of comMod should be overwritten provide file path

        Returns
        --------
        cfgFiles: list
            list of cfg file paths for comMod including the updated values of the parameters to vary

    """

    # get path of module
    modPath = pathlib.Path(comMod.__file__).resolve().parent

    # get filename of module
    modName = str(pathlib.Path(comMod.__file__).stem)

    # read initial configuration
    if fileOverride != '' and cfg['GENERAL'].getboolean('defaultConfig'):
        message = 'fileOverride provided AND defaultComModuleCfg set to True, only one is allowed'
        log.error(message)
        raise AssertionError(message)
    else:
        cfgStart = cfgUtils.getModuleConfig(comMod, fileOverride=fileOverride, toPrint=False,
                                            onlyDefault=cfg['GENERAL'].getboolean('defaultConfig'))

    # create one cfgFile with one line of the parameter values from the full parameter variation
    cfgFiles = []
    for count1, pVal in enumerate(paramValuesD['values']):
        for index, par in enumerate(paramValuesD['names']):
            cfgStart['GENERAL'][par] = str(pVal[index])
        cfgStart['VISUALISATION']['scenario'] = str(count1)
        cfgF = pathlib.Path(cfgPath, ('%d_%sCfg.ini' % (count1, modName)))
        with open(cfgF, 'w') as configfile:
            cfgStart.write(configfile)
        # append file path to list of cfg files
        cfgFiles.append(cfgF)

    return cfgFiles
