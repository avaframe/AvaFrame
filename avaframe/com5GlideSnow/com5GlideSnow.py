"""
    Glide snow tool implemented in com1DFA
"""

# Load modules
import logging
import numpy as np
import pathlib


# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import logUtils
from avaframe.com1DFA import com1DFA

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def runGlideSnow(cfgMain, glideSnowCfg):
    """ Run glide snow tool and adjust parameters to match glide snow settings for com1DFA run,
        result files, reports and plots are saved analog to a standard com1DFA model run

        Parameters
        -----------
        cfgMain: configparser object
            main AvaFrame settings
        glideSnowCfg: configparser object
            configuration settings for glide snow including com1DFA override parameters

    """

    avalancheDir = cfgMain['MAIN']['avalancheDir']

    # ++++++++++ set configurations for com1DFA and override ++++++++++++
    # get comDFA configuration and update with glide snow parameter set
    com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
                                          onlyDefault=glideSnowCfg['com1DFA_override'].getboolean('defaultConfig'))
    com1DFACfg, glideSnowCfg = cfgHandling.applyCfgOverride(com1DFACfg, glideSnowCfg, com1DFA, addModValues=False)

    # run the com1DFA module with glideSnow settings
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgInfo=com1DFACfg)

    # print info about simulation performed to log
    log.info('Com1DFA run performed with glideSnowTool settings')

    return dem, plotDict, reportDictList, simDF
