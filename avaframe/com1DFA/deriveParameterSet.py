"""
    Create dictionary for parameter variations
"""

import logging
import numpy as np
import sys
from datetime import datetime

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.com1DFA import com1DFA
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import geoTrans
import avaframe.in1Data.computeFromDistribution as cP
from avaframe.in1Data import getInput as gI


log = logging.getLogger(__name__)


def getVariationDict(avaDir, fullCfg, modDict):
    """ Create a dictionary with all the parameters that shall be varied from the standard configuration;
        ONLY accounts for variations in section GENERAL and INPUT/releaseScenario

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
            keyList = ['relThPercentVariation', 'entThPercentVariation',
                'secondaryRelThPercentVariation', 'relThRangeVariation',
                'entThRangeVariation', 'secondaryRelThRangeVariation',
                'relThDistVariation', 'entThDistVariation', 'relThRangeFromCiVariation',
                'entThRangeFromCiVariation', 'secondaryRelThRangeFromCiVariation',
                           'secondaryRelThDistVariation', ]
            if key in keyList and value != '':
                # here the factor for changing thValues is added to the variationDict instead of the
                # values directly
                locValue = splitVariationToArraySteps(value, key, fullCfg)
                variations[key] = locValue
                defValue = modDict[section][key][1]
                log.info('%s: %s (default value was: %s)' % (key, locValue, defValue))
            else:
                if any(c in value for c in [':', '|', '$']):
                    # look for default value. If it does not exist, it seems
                    # to be added, ignore it
                    try:
                        defValue = modDict[section][key][1]
                        locValue = fU.splitIniValueToArraySteps(value)
                        variations[key] = locValue
                        log.info('%s: %s (default value was: %s)' % (key, locValue, defValue))
                    except KeyError:
                        defValue = None
                        log.warning('Parameter %s: has a variation, seems to be added, is it acutally used? Ignored for now ' % (key))

    # add releaseScenario info to variations dict - also if this parameter is not varied
    variations['releaseScenario'] = fullCfg['INPUT']['releaseScenario'].split('|')

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
        validResTypes = ['ppr', 'pft', 'pfv', 'pta', 'pke', 'FT', 'FV', 'FM', 'Vx', 'Vy', 'Vz', 'P', 'TA', 'particles']
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


def getParameterVariationInfo(avalancheDir, module, cfgStart):
    """ create a variation dictionary according to parameter variation given in cfg object

        Parameters
        -----------
        avalancheDir: str or pathlib Path
            path to avalanche directory
        module: module

        cfgStart: configparser object
            full configuration object

        Returns
        --------
        modCfg: configparser object
            configuration of simulations to be performed
        variationDict: dict
            dictionary with information on parameter variations

    """

    # Load full configuration: default model settings to compare configuration with modifications
    defCfg = cfgUtils.getDefaultModuleConfig(module, toPrint=False)
    modInfo, modCfg = cfgUtils.compareTwoConfigs(defCfg, cfgStart, toPrint=True)
    # create variation dictionary from modification with respect to the default settings
    variationDict = getVariationDict(avalancheDir, modCfg, modInfo)

    # add avalanche directory info to cfg
    modCfg['GENERAL']['avalancheDir'] = str(avalancheDir)

    return modCfg, variationDict


def getThicknessValue(cfg, inputSimFiles, fName, thType):
    """ set thickness values according to settings chosen and add info to cfg
        if thFromShp = True - add in INPUT section thickness and id info and ci95
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
    ci95List = inputSimFiles[fName]['ci95']

    # create key name for flag
    thFlag = thType + 'FromShp'
    thDistVariation = thType + 'DistVariation'

    # create prefix for release area
    if thType == 'relTh':
        fNamePrefix = fName + '_'
    else:
        fNamePrefix = ''

    # if thickness should be read from shape file
    if cfg['GENERAL'].getboolean(thFlag):
        # if at least one but not all features in a shapefile have a thickness value - error
        if ('None' in thicknessList):
            message = 'Not all features in shape file have a thickness value - check shape file attributes: %s' % fName
            log.error(message)
            raise AssertionError(message)
        elif cfg['GENERAL'][thDistVariation] != '' and ('None' in ci95List):
            message = 'Not all features in shape file have a ci95 value - check shape file attributes: %s' % fName
            log.error(message)
            raise AssertionError(message)
        else:
            # set thickness value in ini file from info of shape file
            thId = idList[0]
            thThickness = thicknessList[0]
            thCi95 = ci95List[0]
            for count, id in enumerate(idList[1:]):
                thId = thId + '|' + id
                thThickness = thThickness + '|' + thicknessList[count+1]
                thCi95 = thCi95 + '|' + ci95List[count+1]

            # add in INPUT section for each feature thickness, id and ci95 values
            # append name of release Scenario in case multiple are provided
            cfg['INPUT'][fNamePrefix + thType + 'Id'] = thId
            cfg['INPUT'][fNamePrefix + thType + 'Thickness'] = thThickness
            cfg['INPUT'][fNamePrefix + thType + 'Ci95'] = thCi95

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

        Returns
        -------
        thicknessSettingsCorrect: bool
            True if settings are provided in the correct format

    """

    # create key name for thickness flag
    thFlag = thName + 'FromShp'
    thFile = thName + 'FromFile'

    # check if flag is set correctly and thickness parameter has correct format
    if cfg['GENERAL'][thFlag] == 'True' or cfg['GENERAL'][thFlag] == 'False':
        if cfg['GENERAL'].getboolean(thFlag):
            if cfg['GENERAL'][thName] != '':
                message = 'If %s is set to True - it is not allowed to set a value for %s' % (thFlag, thName)
                log.error(message)
                raise AssertionError(message)
        else:
            if cfg['GENERAL'][thName] == '' and cfg['GENERAL'].getboolean(thFile) is False:
                message = 'If %s is set to False - it is required to set a value for %s' % (thFlag, thName)
                log.error(message)
                raise AssertionError(message)

    else:
        message = 'Check %s - needs to be True or False' % (thFlag)
        log.error(message)
        raise AssertionError(message)

    # if release thickness should be read from file check other parameters
    if thName == 'relTh':
        if cfg['GENERAL'].getboolean(thFile) and (cfg['GENERAL'].getboolean(thFlag) != False or cfg['GENERAL'][thName] != ''):
            message = 'If %s is set to True - it is not allowed to set %s to True or provide a value in %s' % (thFile, thFlag, thName)
            log.error(message)
            raise AssertionError(message)

    thRV = thName + 'RangeVariation'
    thPV = thName + 'PercentVariation'
    thRCiV = thName + 'RangeFromCiVariation'

    if cfg['GENERAL'][thRV] != '' and cfg['GENERAL'][thPV] != '' and cfg['GENERAL'][thRCiV] != '':
        message = 'Only one variation type is allowed - check %s and %s' % (thRV, thPV)
        log.error(message)
        raise AssertionError(message)

    # if no error ocurred - thickness settings are correct
    thicknessSettingsCorrect = True

    return thicknessSettingsCorrect


def splitVariationToArraySteps(value, key, fullCfg):
    """ split variation in percent to create a list of factors to set parameter value for variations
        or if a rangeVariation is given in absolute values
        or if distVariation create an info string on how the distribution can be build
        (e.g. of format typeOfDistribution$numberOfSteps$ci95value$ci95$support and append the step
        of the current variation in front)

        Parameters
        -----------
        value: str
            value read from configuration
        key: str
            name of parameter
        fullCfg: configparser
            full configuration settings

        Returns
        --------
        itemsArray: numpy array
            factor to change parameter values by multiplication (Percent) or addition (Range)
            or info string on how to build the distribution and which step to draw from it
    """

    # check if positive or negative or both way variation
    itemsL = value.split('$')
    if 'Percent' in key:
        if '-' in itemsL[0]:
            itemsP = itemsL[0].split('-')[1]
            percentsStart = 1. - float(itemsP) / 100.
            percentsStop = 1.0
            # to be sure the right value is taken
            if int(itemsL[1]) == 1:
                percentsStop = percentsStart
        elif '+' in itemsL[0]:
            itemsP = itemsL[0].split('+')[1]
            percentsStart = 1.0
            percentsStop = 1. + float(itemsP) / 100.
            # this is required as the linspace would instead take 1 if only one step
            if int(itemsL[1]) == 1:
                percentsStart = percentsStop
        else:
            percentsStart = 1. - float(itemsL[0]) / 100.
            percentsStop = 1. + float(itemsL[0]) / 100.
        # get number of steps
        steps = int(itemsL[1])
        # create array with variation factor
        itemsArray = np.linspace(percentsStart, percentsStop, steps)
    elif 'RangeFromCi' in key:
        if len(itemsL) == 2:
            itemsArray = ['%d$' % i + value for i in range(int(itemsL[1]))]
        elif len(itemsL) == 3:
            itemsArray = [value]
        else:
            message = ('Format of %s is not correct' % value)
            log.error(message)
            raise AssertionError
    # TODO add standard value to range
    elif 'Range' in key:
        if '-' in itemsL[0]:
            itemsArray = np.linspace(float(itemsL[0]), 0.0, int(itemsL[1]))
        elif '+' in itemsL[0]:
            itemsArray = np.linspace(0.0, float(itemsL[0]), int(itemsL[1]))
        else:
            itemsArray = np.linspace(-1.*float(itemsL[0]), float(itemsL[0]), int(itemsL[1]))
    # if variaiton following normal distribution
    elif 'Dist' in key:
        itemsArray = []
        # if not already appended the step of the distribution values that shall be taken as value
        # add to string so it is clear which value shall be taken from distribution
        # first check format of string
        if len(itemsL) == 6:
            if 'distribution' in itemsL[0]:
                for i in range(int(itemsL[1])):
                    itemsArray.append('%d$' % i + value)
            else:
                message = ('Format of %s is not correct - required format: \
                    typeOfDistribution$numberOfSteps$ci95value$minMaxInterval$ci95$support, \
                    where the first item step is optional' % value)
                log.error(message)
                raise AssertionError
        elif len(itemsL) == 7:
            itemsArray = [value]

    return itemsArray


def setThicknessValueFromVariation(key, cfg, simType, row):
    """ set thickness value for thickness parameter for all features if multiple according to
        desired variation

        Parameters
        ------------
        key: str
            thickness variation info
        cfg: configparser object
            configuration settings of comModule
        simType: str
            simulation type (null, ent, entres, ..)
        row: pandas row
            info on variation of parameters

        Returns
        --------
        cfg: dict
            updated dict with info on thickness

    """

    # fetch info if variation is performed based on a given range or percentage
    if 'RangeFromCi' in key:
        varType = 'RangeFromCi'
    elif 'Range' in key:
        varType = 'Range'
    elif 'Percent' in key:
        varType = 'Percent'
    elif 'Dist' in key:
        varType = 'Dist'

    # only add entries to cfg if appropriate for chosen simType (e.g. entTh if ent or entres run)
    entCondition = (key == ('entTh%sVariation' % varType) and 'ent' in simType)
    secRelCondition = (key == ('secondaryRelTh%sVariation' % varType) and cfg['GENERAL']['secRelArea'] == 'True')
    relCondition = (key == ('relTh%sVariation' % varType))

    # fetch variation factor
    if varType == 'Dist' or varType == 'RangeFromCi':
        # if dist - this is still a string as distribution is only build once mean value is known -
        # if read from shp file - not available yet
        variationFactor = row._asdict()[key]
    else:
        variationFactor = float(row._asdict()[key])

    # update thickness values according to variation
    if entCondition or secRelCondition or relCondition:
        thType = key.split(varType)[0]
        thFlag = thType + 'FromShp'

        # add thickness values for all features if thFromShape = True
        if cfg['GENERAL'][thFlag] == 'True':
            cfg = setVariationForAllFeatures(cfg, key, thType, varType, variationFactor)

        else:
            # update ini thValue if thFromShape=False
            if varType == 'Range':
                cfg['GENERAL'][thType] = str(float(cfg['GENERAL'][thType]) + variationFactor)
            elif varType == 'RangeFromCi':
                message = ('Variation using RangeFromCi is only allowed if thFromShp is set to True')
                log.error(message)
                raise AssertionError(message)
            elif varType == 'Percent':
                cfg['GENERAL'][thType] = str(float(cfg['GENERAL'][thType]) * variationFactor)
            elif varType == 'Dist':
                distInfo = variationFactor.split('$')
                cfgDist = {'sampleSize': distInfo[2], 'mean': float(cfg['GENERAL'][thType]),
                    'buildValue': distInfo[3], 'minMaxInterval': distInfo[4], 'support': '10000',
                    'buildType': distInfo[5]}
                _, distValues, _, _ = cP.extractNormalDist(cfgDist)
                cfg['GENERAL'][thType] = distValues[int(distInfo[0])]
            # set parameter to '' as new thickness value is set for cfg['GENERAL'][thType] and read from here
            cfg['GENERAL'][key] = ''

    else:
        log.debug('%s set but simType is: %s so it is not relevant' % (key, simType))

    return cfg


def setVariationForAllFeatures(cfg, key, thType, varType, variationFactor):
    """ set thickness value for all features according to varType variation

        Parameters
        ----------
        cfg: configparser
            configuration settings of comModule for thickness
        key: str
            name of parameter
        thType: str
            thickness type (e.g. relTh, entTh, ...)
        varType: str
            type of variation (range or percent)
        variationFactor: float or str (if Dist or rangeFromCi)
            value used for variation

        Returns
        --------
        cfg: configparser
            updated configuration settings regarding thickness settings
    """

    # fetch thickness feature ids
    idList = cfg['INPUT'][thType + 'Id'].split('|')
    # fetch thickness list
    thicknessList = cfg['INPUT'][thType + 'Thickness'].split('|')
    # fetch ci95 list
    ci95List = cfg['INPUT'][thType + 'Ci95'].split('|')

    # do some preprocessing if varType is dist
    if varType == 'Dist':
        distInfo = variationFactor.split('$')
        if len(distInfo) != 7:
            message = 'Format of distVariation string is wrong'
            log.error(message)
            raise AssertionError
        else:
            cfgDist = {'sampleSize': distInfo[2], 'minMaxInterval': distInfo[4], 'support': distInfo[6],
                'buildType': distInfo[5]}

    # loop over all features
    for count, id in enumerate(idList):
        thNameId = thType + id

        # set thickness value per feature and fetch value for updating percent/range varation
        if varType == 'Percent':
            # set thickness value in in file for the feature with id Id
            cfg['GENERAL'][thNameId] = str(float(thicknessList[count]) * variationFactor)
            variationIni = setPercentVariation(variationFactor, variationFactor, thNameId)
        elif varType == 'RangeFromCi':
            fromCiVal = setRangeFromCiVariation(cfg, variationFactor,
                thicknessList[count], ci95List[count])
            cfg['GENERAL'][thNameId] = str(float(thicknessList[count]) + fromCiVal)
            variationIni = variationFactor
            strIni = cfg['INPUT']['thFromIni']
            cfg['INPUT']['thFromIni'] = (thNameId if strIni == '' else strIni + '|' + thNameId)
        elif varType == 'Range':
            # set thickness value in in file for the feature with id Id
            cfg['GENERAL'][thNameId] = str(float(thicknessList[count]) + variationFactor)
            variationIni = setRangeVariation(cfg, variationFactor, thNameId)
        elif varType == 'Dist':
            cfgDist['mean'] = str(float(thicknessList[count]))
            cfgDist['buildValue'] = str(float(ci95List[count]))
            _, distValues, _, _ = cP.extractNormalDist(cfgDist)
            cfg['GENERAL'][thNameId] = distValues[int(distInfo[0])]
            distInfo[3] = cfgDist['buildValue']
            variationIni = '$'.join(distInfo)

    # update variation parameter value in config file
    cfg['GENERAL'][key] = variationIni

    return cfg


def setPercentVariation(cfg, variationFactor, thNameId):
    """ determine thickness value if set from percentVariation and set from shp file and update
        percentVariation value that is used for exactely this sim

        this is required for reproducing this sim when using its configuration file - so that in the
        percentVariation parameter it is only one value e.g. +50$1 so +50% in 1 step

        Parameters
        -----------
        cfg: configparser object
            comModule configuration file with info on thickness settings
        variationFactor: float
            value of percent variation in terms of required multiplication of reference value
            (e.g. variationFactor= 0.5 - a variation of 50% performed by multiplication of reference
            value times variationFactor)
        thNameId: str
            name of thickness feature (e.g. relTh0 if release thickness and feature with id 0)

        Returns
        ---------
        variationIni: str
            percentVariation parameter value for this sim to be added in cfg file
    """

    # set percentVaration parameter to actual variation in percent$steps
    if variationFactor == 1.:
        # if variation is 0% set percentVaration = ''
        variationIni = ''
    elif variationFactor < 1:
        variationIni = '-' + str((1. - variationFactor) * 100) + '$1'
    elif variationFactor > 1:
        variationIni = '+' + str((variationFactor - 1.) * 100) + '$1'

    return variationIni


def setRangeVariation(cfg, variationFactor, thNameId):
    """ determine thickness value if set from rangeVariation and set from shp file and update
        rangeVariation value that is used for exactely this sim

        this is required for reproducing this sim when using its configuration file - so that in the
        rangeVariation parameter it is only one value e.g. +0.5$1 so +0.5m in 1 step

        Parameters
        -----------
        cfg: configparser object
            comModule configuration file with info on thickness settings
        variationFactor: float
            value of range variation in terms of required addition to the reference value
            (e.g. variationFactor= +0.5 - a variation of +0.5m performed by addition of reference
            value + variationFactor)
        thNameId: str
            name of thickness feature (e.g. relTh0 if release thickness and feature with id 0)

        Returns
        ---------
        variationIni: str
            rangeVariation parameter value for this sim to be added in cfg file
    """

    # set rangeVaration parameter to actual variation in range$steps
    if variationFactor == 0.:
        # if variation is 0 set RangeVaration = ''
        variationIni = ''
    else:
        variationIni = str(variationFactor) + '$1'

    return variationIni


def setRangeFromCiVariation(cfg, variationFactor, thValue, ciValue):
    """ determine thickness value if set from rangeFromCiVariation and set from shp file and update
        rangeFromCiVariation value that is used for exactely this sim

        this is required for reproducing this sim when using its configuration file - so that in the
        rangeFromCiVariation parameter it is only one value e.g. ci95$4$1 so -ci95m

        Parameters
        -----------
        cfg: configparser object
            comModule configuration file with info on thickness settings
        variationFactor: str
            value of range variation in terms of required addition to the reference value
            (e.g. variationFactor= chi95$4$0 - a variation of -ci95 value performed by addition of reference
            value + the variation value)
        thValue: str
            thickness value of thickness feature in meter
        ciValue: str
            ci value of thickness features in meter

        Returns
        ---------
        variationValue: float
            actual thickness value modified according to variation for feature
        variationIni: str
            rangeFromCiVariation parameter value for this sim to be added in cfg file
    """

    if ciValue == 'None':
        msg = ('ci95 values required in shape file for rangeFromCi variation but not provided')
        log.error(msg)
        raise AssertionError(msg)

    varValStep = int(variationFactor.split('$')[0])
    allSteps = int(variationFactor.split('$')[2])
    thicknessValues = np.linspace(-float(ciValue), float(ciValue), allSteps)
    variationValue = float(thicknessValues[varValStep])

    return variationValue


def appendShpThickness(cfg):
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

    cfgGen = cfg['GENERAL']

    # first create which type of thickness settings are relevant for the current cfg dict
    thTypes = ['relTh']
    if cfgGen['simTypeActual'] in ['ent', 'entres']:
        thTypes.append('entTh')
    if cfgGen['secRelArea'] == 'True':
        thTypes.append('secondaryRelTh')

    # loop over all types and if thickness value read from shp file and no variation has been applied
    # (in this case already added to GENERAL section) - add to section GENERAL
    for thType in thTypes:
        thFlag = thType + 'FromShp'
        thPV = thType + 'PercentVariation'
        thRV = thType + 'RangeVariation'
        thDV = thType + 'DistVariation'
        thRCiV = thType + 'RangeFromCiVariation'
        if cfgGen[thFlag] == 'True' and cfgGen[thPV] == '' and cfgGen[thRV] == '' and cfgGen[thDV] == '' and cfgGen[thRCiV] == '':
            thThickness = thType + 'Thickness'
            thId = thType + 'Id'
            thicknessList = cfg['INPUT'][thThickness].split('|')
            idList = cfg['INPUT'][thId].split('|')
            for count, id in enumerate(idList):
                thNameId = thType + id
                if thNameId in cfg['GENERAL'].keys():
                    log.info('Thickness value for %s already set in initial config file, \
                              read from there not from shp file' % thNameId)
                else:
                    cfgGen[thNameId] = str(float(thicknessList[count]))

    return cfg


def checkDEM(cfgSim, demFile, onlySearch=False):
    """ check if cell size of DEM in Inputs/ is same as desired meshCellSize
        if not - check for remeshed DEM or remesh the DEM

        Parameters
        -----------
        cfgSim: configparser object
            configuration settings of com module
        demFile: str or pathlib path
            path to dem in Inputs/
        onlySearch: bool
            if True - only searching for remeshed DEM but not remeshing if not found

        Returns
        --------
        pathToDem: str
            path to DEM with correct cellSize relative to Inputs/
    """

    # read header of DEM file
    headerDEM = IOf.readASCheader(demFile)

    # fetch info on desired meshCellSize
    meshCellSize = float(cfgSim['GENERAL']['meshCellSize'])
    meshCellSizeThreshold = float(cfgSim['GENERAL']['meshCellSizeThreshold'])

    # if cell size of DEM is different from desired meshCellSize - look for remeshed DEM or remesh
    if np.abs(meshCellSize - headerDEM['cellsize']) > meshCellSizeThreshold:
        pathToDem = geoTrans.remeshDEM(demFile, cfgSim, onlySearch=onlySearch)
    else:
        log.info('DEM taken from Inputs/')
        pathToDem = demFile.name

    log.info('path to DEM is: %s' % pathToDem)

    return pathToDem

def writeToCfgLine(values):
    """ write an array of values to a string of values separated by | for configuration

        Parameters
        -----------
        values: numpy array
            array of values

        Returns
        --------
        valString: str
            string of array values separated by |
    """

    valString = '%.12f' % values[0]
    for val in values[1:]:
        valString = valString + '|%.12f' % val

    return valString


def createSimDict(avalancheDir, com1DFA, cfgInitial, inputSimFiles, simNameExisting):
    """ Create a simDict with all the simulations that shall be performed

        Parameters
        -----------
        avalancheDir: pathlib path
            path to avalanche directory
        com1DFA: module
            computational module
        cfgStart: configparser object
            configuration settings for com1DFA
        inputSimFiles: dict
            dictionary with info in input files (release area, dem, ...)
        simNameExisting: list
            list with names of sims that already exist in outputs

        Returns
        --------
        simDict: dict
            dicionary with info on simHash, releaseScenario, release area file path,
            simType and contains full configuration configparser object for simulation run
    """

    # check if thickness settings in ini file are valid
    for thType in ['entTh', 'relTh', 'secondaryRelTh']:
        _ = checkThicknessSettings(cfgInitial, thType)
    # update thickness settings, e.g. fetch if th read from shp
    cfgInitial = gI.updateThicknessCfg(inputSimFiles, cfgInitial)

    # reset variationDict
    variationDict = ''

    # create a dictionary with information on which parameter shall be varied for individual simulations
    # compare cfgStart to default module config for this
    modCfg, variationDict = getParameterVariationInfo(avalancheDir, com1DFA, cfgInitial)

    # create a configuration object per simulation to run (from configuration) gathered in simDict
    # only new simulations are included in this simDict
    # key is simName and corresponds to one simulation
    simDict = {}
    simDict = com1DFA.prepareVarSimDict(modCfg, inputSimFiles, variationDict, simNameExisting=simNameExisting)

    # write full configuration (.ini file) to file
    date = datetime.today()
    fileName = 'sourceConfiguration_' + '{:%d_%m_%Y_%H_%M_%S}'.format(date)
    cfgUtils.writeCfgFile(avalancheDir, com1DFA, modCfg, fileName=fileName)

    return simDict
