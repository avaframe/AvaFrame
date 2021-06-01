'''
    Utilities for handling configuration files

'''

import configparser
import os
import logging
import pathlib
import hashlib
import json
import pandas as pd
import glob

# Local imports
import avaframe as avaf
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import fileHandlerUtils as fU


log = logging.getLogger(__name__)


def getGeneralConfig():
    ''' Returns the general configuration for avaframe
    returns a configParser object
    '''

    # get path of module
    modPath = pathlib.Path(avaf.__file__).resolve().parent

    localFile = modPath / 'local_avaframeCfg.ini'
    defaultFile = modPath / 'avaframeCfg.ini'

    if localFile.is_file():
        iniFile = localFile
        iniFile = [defaultFile, localFile]
        compare = True
    elif defaultFile.is_file():
        iniFile = defaultFile
        compare = False
    else:
        raise FileNotFoundError('None of the provided cfg files exist ')

    # Finally read it
    cfg = compareConfig(iniFile, 'General', compare)

    return cfg


def getModuleConfig(module, fileOverride='', modInfo=False):
    ''' Returns the configuration for a given module
    returns a configParser object

    module object: module : the calling function provides the already imported
           module eg.:
           from avaframe.com2AB import com2AB
           leads to getModuleConfig(com2AB)
           whereas
           from avaframe.com2AB import com2AB as c2
           leads to getModuleConfig(c2)

    Str: fileOverride : allows for a completely different file location

    modInfo: bool
        true if dictionary with info on differences to standard config

    Order is as follows:
    fileOverride -> local_MODULECfg.ini -> MODULECfg.ini

    '''

    # get path of module
    modPath = pathlib.Path(module.__file__).resolve().parent

    # get filename of module
    modName = str(pathlib.Path(module.__file__).stem)

    localFile = modPath / ('local_'+modName+'Cfg.ini')
    defaultFile = modPath / (modName+'Cfg.ini')

    log.debug('localFile: %s', localFile)
    log.debug('defaultFile: %s', defaultFile)

    # Decide which one to take
    if fileOverride:
        if os.path.isfile(fileOverride):
            iniFile = [defaultFile, fileOverride]
            compare = True
        else:
            raise FileNotFoundError('Provided fileOverride does not exist: '
                                    + fileOverride)

    elif localFile.is_file():
        iniFile = localFile
        iniFile = [defaultFile, localFile]
        compare = True
    elif defaultFile.is_file():
        iniFile = defaultFile
        compare = False
    else:
        raise FileNotFoundError('None of the provided cfg files exist ')

    # Finally read it
    cfg = compareConfig(iniFile, modName, compare, modInfo)

    return cfg


def getDefaultModuleConfig(module):
    ''' Returns the default configuration for a given module
    returns a configParser object

    module object: module : the calling function provides the already imported
           module eg.:
           from avaframe.com2AB import com2AB
           leads to getModuleConfig(com2AB)
           whereas
           from avaframe.com2AB import com2AB as c2
           leads to getModuleConfig(c2)

    '''

    # get path of module
    modPath = pathlib.Path(module.__file__).resolve().parent

    # get filename of module
    modName = str(pathlib.Path(module.__file__).stem)

    defaultFile = modPath / (modName+'Cfg.ini')

    log.debug('defaultFile: %s', defaultFile)

    # Finally read it
    cfg = compareConfig(defaultFile, modName, compare=False)

    return cfg


def compareConfig(iniFile, modName, compare, modInfo=False):
    ''' Compare configuration files (if a local and default are both provided)
    and inform user of the eventuel differences. Take the default as reference.

    Inputs:
            -iniFile: path to config file. Only one path if compare=False
            -compare: True if two paths are provided and a comparison is needed
            -modInfo: True if dictionary with modifications shall be returned

    Output: ConfigParser object
    '''

    modDict = {}
    if compare:
        log.info('Reading config from: %s and %s' % (iniFile[0], iniFile[1]))
        # initialize our final configparser object
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        # initialize configparser object to read
        defCfg = configparser.ConfigParser()
        defCfg.optionxform = str
        locCfg = configparser.ConfigParser()
        locCfg.optionxform = str
        # read default and local parser files
        defCfg.read(iniFile[0])
        locCfg.read(iniFile[1])
        # loop through all sections of the defCfg
        log.debug('Writing cfg for: %s', modName)
        for section in defCfg.sections():
            modDict[section] = {}
            cfg.add_section(section)
            log.info('\t%s', section)
            for key in defCfg.items(section):
                defValue = key[1]
                # check if key is also in the localCfg
                if locCfg.has_option(section, key[0]):
                    locValue = locCfg.get(section, key[0])
                    if locValue != defValue:
                        # if yes and if this value is different add this key to
                        # the cfg that will be returned
                        locValue = locCfg.get(section, key[0])
                        cfg.set(section, key[0], locValue)
                        log.info('\t\t%s : %s \t(default value was : %s)',
                                 key[0], locValue, defValue)
                        modString = [locValue, defValue]
                        modDict[section][key[0]] = modString
                    else:
                        cfg.set(section, key[0], defValue)
                        log.info('\t\t%s : %s', key[0], defValue)

                    # remove the key from the localCfg
                    locCfg.remove_option(section, key[0])
                else:
                    cfg.set(section, key[0], defValue)
                    log.info('\t\t%s : %s', key[0], defValue)

        # Now check if there are some sections/ keys left in the local cfg and
        # that are not used
        for section in locCfg.sections():
            for key in locCfg.items(section):
                log.warning('Key [\'%s\'] in section [\'%s\'] in the localCfg file is not needed.' % (key[0], section))

    else:
        log.info('Reading config from: %s', iniFile)
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        # Finally read it
        cfg.read(iniFile)
        # Write config to log file
        logUtils.writeCfg2Log(cfg, modName)

    if modInfo:
        return cfg, modDict
    else:
        return cfg


def writeCfgFile(avaDir, module, cfg, fileName=''):
    """ Save configuration used to text file in Outputs as moduleName_settings.ini
        or optional in Outputs/moduleName/configurationFiles/filenName.ini

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        module:
            module
        cfg: configparser object
            configuration settings
        fileName: str
            name of saved configuration file - optional

    """

    # get filename of module
    name = os.path.basename(module.__file__)
    modName = name.split('.')[0]

    # write to file
    if fileName != '':
        # set outputs
        modName = 'com1DFAPy'
        outDir = os.path.join(avaDir, 'Outputs', modName, 'configurationFiles')
        fU.makeADir(outDir)
        cfg.optionxform = str
        with open(os.path.join(outDir, '%s.ini' % (fileName)), 'w') as conf:
            cfg.write(conf)
    else:
        # set outputs
        outDir = os.path.join(avaDir, 'Outputs')
        cfg.optionxform = str
        with open(os.path.join(outDir, '%s_settings.ini' % (modName)), 'w') as conf:
            cfg.write(conf)


def readCfgFile(avaDir, module='', fileName=''):
    """ Read configuration from ini file, if module is provided, module configuration is read from Ouputs,
        if fileName is provided configuration is read from fileName

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        module:
            module
        fileName: str
            path to file that should be read - optional

        Returns
        --------
        cfg: configParser object
            configuration that is from file

    """

    # define file that should be read
    if fileName != '':
        inFile = fileName
    elif module != '':
        # get module name
        name = os.path.basename(module.__file__)
        modName = name.split('.')[0]
        # set input file
        inFile = os.path.join(avaDir, 'Outputs', '%s_settings.ini' % (modName))
    else:
        log.error('Please provide either a module or a fileName to read configuration from file')
        raise NameError

    # read configParser object from input file, case sensitive
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(inFile)
    cfg.optionxform = str

    return cfg


def cfgHash(cfg, typeDict=False):
    """ UID hash of a config. Given a configParser object cfg,
        or a dictionary - then typeDict=True, returns a uid
    hash

    Parameters
    ----------
    cfg: configParser object
    typeDict : dict
        dictionary

    Returns:
    --------
    uid: str
       uid hash
    """

    uidHash = hashlib.shake_256()

    if typeDict:
        cfgDict = cfg
    else:
        cfgDict = convertConfigParserToDict(cfg)

    jsonDict = json.dumps(cfgDict, sort_keys=True, ensure_ascii=True)
    encoded = jsonDict.encode()

    uidHash.update(encoded)
    uid = uidHash.hexdigest(5)

    return uid


def convertConfigParserToDict(cfg):
    """ create dictionary from configparser object """

    cfgDict = {}
    for section in cfg.sections():
        cfgDict[section] = {}
        for key, val in cfg.items(section):
            cfgDict[section][key] = val

    return cfgDict


def convertDictToConfigParser(cfgDict):
    """ create configParser object from dict """

    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    for section in cfgDict:
        cfg[section] = cfgDict[section]

    return cfg


def writeDictToJson(inDict, outFilePath):
    """ write a dictionary to a json file """

    jsonDict = json.dumps(inDict, sort_keys=True, ensure_ascii=True)
    f = open(outFilePath, "w")
    f.write(jsonDict)
    f.close()


def createConfigurationInfo(avaDir, standardCfg, writeCSV=False, specDir=''):
    """ Read configurations from all simulations configuration ini files from directory

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        standardCfg: dict
            standard configuration for module

        Returns
        --------
        simDF: pandas DataFrame
            DF with all the simulation configurations
    """

    # read default configuration of this module
    standardCfgDict = convertConfigParserToDict(standardCfg)
    # create pandas dataFrame
    simDF = pd.DataFrame(data=standardCfgDict['GENERAL'], index=['current standard'])

    # collect all configuration files for this module
    if specDir != '':
        inDir = os.path.join(specDir, 'configurationFiles')
    else:
        inDir = os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'configurationFiles')
    configFiles = glob.glob(inDir+os.sep+'*.ini')

    # create confiparser object, convert to json object, write to dataFrame
    # append all dataFrames
    for cFile in configFiles:
        if 'sourceConfiguration' not in cFile:
            simName = os.path.splitext(os.path.basename(cFile))[0]
            if '_AF_' in simName:
                nameParts = simName.split('_AF_')
                fNamePart = nameParts[0] + '_AF'
                infoParts = nameParts[1].split('_')

            else:
                nameParts = simName.split('_')
                fNamePart = nameParts[0]
                infoParts = nameParts[1:]
            simHash = infoParts[2]
            cfgObject = readCfgFile(avaDir, fileName=cFile)
            indexItem = [simHash]
            cfgDict = convertConfigParserToDict(cfgObject)
            simItemDF = pd.DataFrame(data=cfgDict['GENERAL'], index=indexItem)
            simItemDF = simItemDF.assign(simName=simName)
            simDF = pd.concat([simDF, simItemDF], axis=0)

    # if writeCSV, write dataFrame to csv file
    if writeCSV:
        outFile = os.path.join(inDir, 'allConfigurations.csv')
        simDF.to_csv(outFile)

    return simDF


def filterSims(avalancheDir, parametersDict, module, specDir=''):
    """ Filter simulations using a list of parameters and a pandas dataFrame of simulation configurations """

    standardCfg = getDefaultModuleConfig(module)
    simDF = createConfigurationInfo(avalancheDir, standardCfg, writeCSV=False, specDir=specDir)
    for key, value in parametersDict.items():
        if value != []:
            if isinstance(value, list) == False:
                value = [value]
            simDF = simDF[simDF[key].isin(value)]

    if 'current standard' in simDF.index:
        simDF = simDF.drop('current standard')

    simNameList = simDF['simName']

    return simNameList
