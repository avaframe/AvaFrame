"""
    Tools specific to the com1DFA computational kernel
"""

# Load modules
import logging
import math
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from deepdiff import DeepDiff

# local imports
import avaframe.com1DFA.deriveParameterSet as dP
from avaframe.in1Data import getInput as gI

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getPartInitMethod(cfg, csz, relThForPart):
    """ Get particle initialization parameters

    Get the massPerPart and nPPK corresponding to the desired initialization method

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    csz: float
        cell size
    relThForPart: float
        relTh value

    Returns
    -------
    massPerPart : float
        mass per particle desired
    nPPK : int
        number of particles per kernel radius desired
    """
    rho = cfg.getfloat('rho')
    massPerParticleDeterminationMethod = cfg['massPerParticleDeterminationMethod']
    nPPK = 0
    # derive mass per particle to define number of particles per cell:
    if massPerParticleDeterminationMethod == 'MPPDIR':
        massPerPart = cfg.getfloat('massPerPart')
        log.debug('Number of particles defined by: mass per particle %s' % cfg['massPerPart'])
    elif massPerParticleDeterminationMethod == 'MPPDH':
        deltaTh = cfg.getfloat('deltaTh')
        ds = min(csz, cfg.getfloat('sphKernelRadius'))
        massPerPart = rho * ds * ds * deltaTh
        log.debug('Number of particles defined by: release thickness per particle: %s' % cfg['deltaTh'])
        log.debug('mass per particle is %.2f' % massPerPart)
    elif massPerParticleDeterminationMethod == 'MPPKR':
        sphKernelRadius = cfg.getfloat('sphKernelRadius')
        cszMin = min(csz, sphKernelRadius)
        nPPK0 = cfg.getfloat('nPPK0')
        sphKR0 = cfg.getfloat('sphKR0')
        aPPK = cfg.getfloat('aPPK')
        nPPK = round(nPPK0 * (cszMin/sphKR0)**aPPK)
        massPerPart = rho * math.pi * cszMin * cszMin * relThForPart / nPPK

    return massPerPart, nPPK

def setRelThIni(avaDir, modName, cfgInitial):
    """ Add thickness values in configuration file according to thickness flags, ini settings or shapefile attributes
        and create inputSimFiles dictionary with file paths to input data

        Parameters
        -----------
        avaDir: str or pathlib path
            path to avalanche directory
        modName: module
            computational module
        cfgInitial: configparser object
            full configuration settings of com Module

        Returns
        --------
        inputSimFilesAll: dict
            dictionary with infos about input data file paths and flags for entrainment, resistance
    """

    # check if thickness settings in ini file are valid
    for thType in ['entTh', 'relTh', 'secondaryRelTh']:
        _ = dP.checkThicknessSettings(cfgInitial, thType)

    # fetch input data - dem, release-, entrainment- and resistance areas (and secondary release areas)
    inputSimFilesAll = gI.getInputDataCom1DFA(avaDir, cfgInitial)

    # get thickness of release and entrainment areas (and secondary release areas) -if thFromShp = True
    inputSimFilesAll, cfgInitial = gI.getThickness(inputSimFilesAll, avaDir, modName, cfgInitial)

    return inputSimFilesAll, cfgInitial


def compareSimCfgToDefaultCfgCom1DFA(simCfg):
    """ Compares the given simulation configuration (as dict) to the default 
        com1DFA configuration. Disregards values like avalancheDir that are expected to 
        change. Returns True if it is the default + an indentifier string: D = Default and 
        C = Changed

        Parameters
        -----------
        simCfg: dict
            simulation configuration

        Returns
        --------
        isDefault: bool
            True if it is the same as the default configuration 
        defaultIdentifierString: str
            D if default and C if changed

    """

    defaultIdentifierString = 'D'

    # Get default cfg and convert to dict for comparison
    defCfgObject = cfgUtils.getDefaultModuleConfig(com1DFA, toPrint=False)
    defCfg = cfgUtils.convertConfigParserToDict(defCfgObject)

    # Which changes to ignore (in the case of com1DFA). These are expected
    # to change...
    excludeItems = ["root['GENERAL']['avalancheDir']", 
                   "root['GENERAL']['secRelArea']",
                   "root['GENERAL']['simTypeList']",
                   "root['INPUT']['releaseScenario']",
                   ]

    # do the diff and analyse
    diff = DeepDiff(defCfg, simCfg, exclude_paths=excludeItems)


    if 'values_changed' in diff:
        log.info('Comparing to default cfg, values changed:')
        log.info(diff['values_changed'])
        log.info('Modified from default')
        isDefault = False 
        defaultIdentifierString = 'C'
    else: 
        diff['values_changed'] = None

    if 'dictionary_item_added' in diff:
        log.debug('Comparing to default cfg, added items:')
        log.debug(diff['dictionary_item_added'])

    return defaultIdentifierString, diff['values_changed']
