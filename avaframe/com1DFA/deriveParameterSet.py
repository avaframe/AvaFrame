"""
    Create dictionary for parameter variations
"""

import logging
import numpy as np

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.com1DFA import com1DFA

log = logging.getLogger(__name__)


def getVariationDict(avaDir, fullCfg, modDict):
    """ Create a dictionary with all the parameters that shall be varied from the standard configuration;
        ONLY takes care of variations in section GENERAL

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
        # percent variation info already used for updating thickness values
        if key not in ['resType', 'tSteps']:
            # if yes and if this value is different add this key to
            # the parameter variation dict
            if key in ['relThPercentVariation', 'entThPercentVariation', 'secondaryRelThPercentVariation'] and value != '':
                # here the factor for changing thValues is added to the variationDict instead of the
                # values directly
                locValue = splitVariationToArraySteps(value, key)
                variations[key] = locValue
                defValue = modDict[section][key][1]
                log.info('%s: %s (default value was: %s)' % (key, locValue, defValue))
            else:
                if ':' in value or '|' in value or '$' in value:
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
        validResTypes = ['ppr', 'pfd', 'pfv', 'pta', 'FD', 'FV', 'Vx', 'Vy', 'Vz', 'P', 'TA', 'particles']
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
        if isinstance(variationDict[parameter], str):
            if '|' in variationDict[parameter] or ':' in variationDict[parameter]:
                items = fU.splitIniValueToArraySteps(variationDict[parameter], returnList=False)
                variationDict[parameter] = items
        if parameter in standardCfg['GENERAL']:
            if not isinstance(variationDict[parameter], (list, np.ndarray)):
                variationDict[parameter] = [variationDict[parameter]]
        else:
            ignoredParameters.append(parameter)

    for ipar in ignoredParameters:
        log.warning('Parameter %s does not exist in model configuration - parameter is ignored' % ipar)
        del variationDict[ipar]

    return variationDict


def getParameterVariationInfo(avalancheDir, com1DFA, cfgFile):
    """ read info on which simulations shall be performed according to parameter variation

        Parameters
        -----------
        avalancheDir: str or pathlib Path
            path to avalanche directory
        cfgFile: str or pathlib Path
            path to override configuration file

        Returns
        --------
        modCfg: configparser object
            configuration of simulations to be performed
        variationDict: dict
            dictionary with information on parameter variations

    """

    # generate list of simulations from desired configuration
    # Load full configuration
    modCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=cfgFile, modInfo=True)
    variationDict = getVariationDict(avalancheDir, modCfg, modInfo)

    # add avalanche directory info to cfg
    modCfg['GENERAL']['avalancheDir'] = str(avalancheDir)

    return modCfg, variationDict


def checkRelEntThVariation(cfg, variationDict):
    """ check if release or entrainment thickness variation is working - due to where it is read from

        cfg: configparser object
            configuration settings
        variationDict: dict
            dictionary with info on parameter variation

    """

    # if parameter variation for release or entrainment thickness, print warning if thickness read from shape file
    thicknessTypes = {'relTh': 'release', 'entTh': 'entrainment'}
    for key, value in thicknessTypes.items():
        flag = 'use' + key[0].upper() + key[1:] + 'FromIni'
        if key in variationDict and cfg['GENERAL'].getboolean(flag) is False:
            log.warning('Parameter variation for %s thickness not working as %s read from shp file - \
                         consider setting %s to True' % (value, key, flag))


def getThicknessValue(cfg, inputSimFiles, fName, thType):
    """ set thickness values according to settings chosen and add info to cfg
        if thFromShp = True - add in INPUT section thickness and id info
        if thFromShp = False - check format of thickness value in GENERAL section

        Parameters
        -----------
        cfg: configparser object
            configuration settings
        inputSimFiles: dict
            dictionary with info on input files and attributes (id and thickness)
        fName: str
            name of scenario shp file (entrainment, release, ...)
        thType: str
            thickness parameter name (entTh, relTh, ...)

        Returns
        -------
        cfg: configparser object
            updated configuration settings with info on actual thickness values to be used for simulations

        """

    # fetch thickness values from shapefile
    thicknessList = inputSimFiles[fName]['thickness']
    idList = inputSimFiles[fName]['id']

    # create key names for flags and parameter variation info
    thFlag = thType + 'FromShp'
    thPercent = thType + 'PercentVariation'

    # if thickness should be read from shape file
    if cfg['GENERAL'].getboolean(thFlag):
        # if at least one but not all features in a shapefile have a thickness value - error
        if ('None' in thicknessList):
            message = 'Not all features in shape file have a thickness value - check shape file attributes: %s' % fName
            log.error(message)
            raise AssertionError(message)
        else:
            # set thickness value in ini file from info of shape file
            thId = idList[0]
            thThickness = thicknessList[0]
            for count, id in enumerate(idList[1:]):
                thId = thId +'|' + id
                thThickness = thThickness + '|'+ thicknessList[count+1]

            # add in INPUT section
            cfg['INPUT'][thType + 'Id'] = thId
            cfg['INPUT'][thType + 'Thickness'] = thThickness

    else:
        # if thickness should be read from ini file - check if format is correct
        if '$' in cfg['GENERAL'][thType] and len(cfg['GENERAL'][thType].split('$')) != 3:
            message = 'Format of relTh value in ini file is not correct - for variation from ini use refValue$percent$numberOfSteps'
            log.error(message)
            raise AssertionError(message)

    return cfg


def checkThicknessSettings(cfg, thName):
    """ check if thickness setting format is correct

        Parameters
        ------------
        cfg: configparser
            configuration settings
        thName: str
            thickness parameter name (entTh, ...)

    """

    # create key name for thickness flag
    thFlag = thName + 'FromShp'

    # check if flag is set correctly and thickness parameter has correct format
    if cfg['GENERAL'][thFlag] == 'True' or cfg['GENERAL'][thFlag] == 'False':
        if cfg['GENERAL'].getboolean(thFlag):
            if cfg['GENERAL'][thName] != '':
                message = 'If %s is set to True - it is not allowed to set a value for %s' % (thFlag, thName)
                log.error(message)
                raise AssertionError(message)
        else:
            if cfg['GENERAL'][thName] == '':
                message = 'If %s is set to False - it is required to set a value for %s' % (thFlag, thName)
                log.error(message)
                raise AssertionError(message)

    else:
        message = 'Check %s - needs to be True or False' % (thFlag)
        log.error(message)
        raise AssertionError(message)


def splitVariationToArraySteps(value, key):
    """ split variation in percent to create a list of factors to set paramerter value for variations
        or if a rangeVariation is given in absolute values

        Parameters
        -----------
        value: str
            value read from configuration
        key: str
            name of parameter

        Returns
        --------
        itemsArray: numpy array
            factor to change parameter values by multiplication (Percent) or addition (Range)
    """

    # check if positive or negative or both way variation
    itemsL = value.split('$')
    if 'Percent' in key:
        if '-' in itemsL[0]:
            itemsP = itemsL[0].split('-')[1]
            percentsStart = 1. - float(itemsP) / 100.
            percentsStop = 1. - float(itemsP) / 100.
        elif '+' in itemsL[0]:
            itemsP = itemsL[0].split('+')[1]
            percentsStart = 1. + float(itemsP) / 100.
            percentsStop = 1. + float(itemsP) / 100.
        else:
            percentsStart = 1. - float(itemsL[0]) / 100.
            percentsStop = 1. + float(itemsL[0]) / 100.
        # get number of steps
        steps = int(itemsL[1])
        # create array with variation factor
        itemsArray = np.linspace(percentsStart, percentsStop, steps)
    elif 'Range' in key:
        itemsArray = np.linspace(-1.*itemsL[0], itemsL[0], itemsL[1])

    return itemsArray


def setThicknessValueFromVariation(key, cfg, simType, row):
    """ set thickness value for thickness parameter for all features if multiple according to
        desired variation

        Parameters
        ------------
        key: str
            thickness variation info

        Returns
        --------
        cfg: dict
            updated dict with info on thickness

    """

    # only add entries to cfg if appropriate for chosen simType (e.g. entTh if ent or entres run)
    entCondition = (key == 'entThPercentVariation' and 'ent' in simType)
    secRelCondition = (key == 'secondaryRelThPercentVariation' and cfg['GENERAL']['secRelArea'] == 'True')
    relCondition = (key == 'relThPercentVariation')

    # fetch variation factor
    variationFactor = float(row._asdict()[key])

    # update thickness values accoridng to variation
    if entCondition or secRelCondition or relCondition:
        thType = key.split('Percent')[0]
        thFlag = thType + 'FromShp'

        # add thickness values for all features if thFromShape = True
        if cfg['GENERAL'][thFlag] == 'True':
            thId = thType + 'Id'
            thThickness = thType + 'Thickness'
            thicknessList = cfg['INPUT'][thThickness].split('|')
            idList = cfg['INPUT'][thId].split('|')
            for count, id in enumerate(idList):
                thNameId = thType + id
                cfg['GENERAL'][thNameId] = str(float(thicknessList[count]) * variationFactor)

            # set percentVaration parameter to actual variation in percent$steps
            if variationFactor <= 1:
                variationIni = '-' + str((1. - variationFactor) * 100) + '$1'
            else:
                variationIni = '+' + str((variationFactor - 1.) * 100) + '$1'

            # update parameter value
            cfg['GENERAL'][key] = variationIni

        else:
            # update ini thValue if thFromShape=False
            cfg['GENERAL'][thType] = str(float(cfg['GENERAL'][thType]) * variationFactor)
            # set parameter to '' as new thickness value is set for cfg['GENERAL'][thType] and read from here
            cfg['GENERAL'][key] = ''

    return cfg


def appendShpThickness(cfg, variationDict):
    """ append thickness values to GENERAL section if read from shp and not varied

        Parameters
        -----------
        cfg: dict
            configuration settings

        Returns
        --------
        cfg: dict
            updated configuartion settings

    """

    # first create which type of thickness settings are relevant for the current cfg dict
    thTypes = ['relTh']
    if cfg['GENERAL']['simTypeActual'] in ['ent', 'entres']:
        thTypes.append('entTh')
    if cfg['GENERAL']['secRelArea'] == 'True':
        thTypes.append('secondaryRelTh')

    # loop over all types and if thickness value read from shp file and no variation has been applied
    # (in this case already added to GENERAL section) - add to section GENERAL
    for thType in thTypes:
        thFlag = thType + 'FromShp'
        thPV = thType + 'PercentVariation'
        if cfg['GENERAL'][thFlag] == 'True' and cfg['GENERAL'][thPV] == '':
            thThickness = thType + 'Thickness'
            thId = thType + 'Id'
            thicknessList = cfg['INPUT'][thThickness].split('|')
            idList = cfg['INPUT'][thId].split('|')
            for count, id in enumerate(idList):
                thNameId = thType + id
                cfg['GENERAL'][thNameId] = str(float(thicknessList[count]))

    return cfg
