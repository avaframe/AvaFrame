"""
    Main functions for python DFA kernel
"""

import logging

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils


log = logging.getLogger(__name__)


def getVariationDict(avaDir, module, standardCfg, cfgFile=''):
    """ Create a dictionary with all the parameters that shall be varied from the standard configuration;
        either provide a cfgFile that contains info on parameters to be varied,
        or the local_ cfg file is used

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        module:
            computational module
        standardCfg: dict
            default configuration
        cfgFile: str
            path to configuration file that includes info on parameters to be varied - optional

        Returns
        -------
        variationDict: dict
            dictionary with the parameters that shall be varied

    """

    # default configurations for module
    defCfg = standardCfg

    # load modified configuration either from provided cfg file or from local_ cfg file
    if cfgFile != '':
        locCfg = cfgUtils.getModuleConfig(module, fileOverride=cfgFile)
    else:
        locCfg = cfgUtils.getModuleConfig(module)

    variations = {}
    # loop through all sections of the defCfg
    for section in defCfg.sections():
        # look for parameters that are different than default in section GENERAL
        if section == 'GENERAL':
            variations[section] = {}
            for key in defCfg.items(section):
                # output saving options not relevant for parameter variation!

                defValue = key[1]
                # check if key is also in the localCfg
                if locCfg.has_option(section, key[0]):
                    locValue = locCfg.get(section, key[0])
                    if locValue != defValue:
                        # if yes and if this value is different add this key to
                        # the parameter variation dict
                        if key[0] == 'resType' or key[0] == 'tSteps':
                            locValue = [locValue]
                        elif ':' in locValue or '|' in locValue:
                            locValue = fU.splitIniValueToArraySteps(locValue)
                        else:
                            # if only one item - create list
                            locValue = [locValue]
                        variations[section][key[0]] = locValue
                        log.info('\t\t%s : %s \t(default value was : %s)',
                                 key[0], locValue, defValue)
                    # remove the key from the localCfg
                    locCfg.remove_option(section, key[0])

    # Now check if there are some sections/ keys left in the local cfg and
    # that are not used
    for section in locCfg.sections():
        if section == 'GENERAL':
            for key in locCfg.items(section):
                log.warning('Key [\'%s\'] in section [\'%s\'] in the parameter variation Cfg file is not needed.' % (key[0], section))

    return variations['GENERAL']
