"""
    Create dictionary for parameter variations
"""

import logging
import numpy as np

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils


log = logging.getLogger(__name__)


def getVariationDict(avaDir, fullCfg, modDict):
    """ Create a dictionary with all the parameters that shall be varied from the standard configuration;
        either provide a cfgFile that contains info on parameters to be varied,
        or the local_ cfg file is used

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        fullCfg: configParser object
            full configuration potentially including variations of parameter
        modDict: dict
            info on modifications to standard configuration

        Returns
        -------
        variationDict: dict
            dictionary with the parameters that shall be varied and the chosen flags for the run

    """

    # look for parameters that are different than default in section GENERAL
    section = 'GENERAL'
    variations = {}
    for key, value in fullCfg.items(section):
        if key == 'resType':
            fullCfg = checkResType(fullCfg, section, key, value)
        # output saving options not relevant for parameter variation!
        if key not in ['resType', 'tSteps']:
            # if yes and if this value is different add this key to
            # the parameter variation dict
            if ':' in value or '|' in value:
                locValue = fU.splitIniValueToArraySteps(value)
                variations[key] = locValue
                defValue = modDict[section][key][1]
                log.info('%s: %s (default value was: %s)' % (key, locValue, defValue))

    # print modified parameters
    for sec in modDict:
        for value in modDict[sec]:
            if sec != section:
                log.info('%s: %s (default value was: %s)' % (value, modDict[sec][value][0], modDict[sec][value][1]))
            else:
                if value not in variations:
                    log.info('%s: %s (default value was: %s)' % (value, modDict[sec][value][0], modDict[sec][value][1]))

    return variations


def checkResType(fullCfg, section, key, value):
    """ Check if the resTypes asked for exist
    Warns the user if some do not, removes them from the resType list and
    updates the cfg

        Parameters
        -----------
        fullCfg: configParser object
            full configuration potentially including variations of parameter
        section: str
            section name
        key: str
            key name
        value: str
            corresponding value

        Returns
        --------
        fullCfg: configParser object
            full configuration updated with resType if this last one was modified

    """
    # check that the resType asked actually exists
    if value != '':
        resType = value.split('|')
        validResTypes = ['ppr', 'pfd', 'pfv', 'FD', 'FV', 'P', 'particles']
        message = (
            'The parameter % s is not a valid resType. It will not be saved')
        newResType = []
        for res in resType:
            if res not in (validResTypes or ['']):
                log.warning(message % res)
            else:
                newResType.append(res)
        fullCfg[section][key] = '|'.join(newResType)
    return fullCfg


def validateVarDict(variationDict, standardCfg):
    """ Check if all parameters in variationDict exist in default configuration and
        are provided in the correct format

        Parameters
        -----------
        variationDict: dict
            dictionary with parameters that shall be varied and a list for the parameter values for each parameter
        standardCfg: config Parser object
            default model configuration

        Returns
        --------
        variationDict: dict
            cleaned variation dict that meets the required structure

    """

    # check if values are provided as list or numpy array
    # check if parameters exist in model configuration
    ignoredParameters = []
    for parameter in variationDict:
        if parameter in standardCfg['GENERAL']:
            if not isinstance(variationDict[parameter], (list, np.ndarray)):
                variationDict[parameter] = [variationDict[parameter]]
        else:
            ignoredParameters.append(parameter)

    for ipar in ignoredParameters:
        log.warning('Parameter %s does not exist in model configuration - parameter is ignored' % ipar)
        del variationDict[ipar]

    return variationDict
