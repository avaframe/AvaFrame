'''
    Utilities for handling configuration files

'''

import configparser
import logging
import pathlib
import hashlib
import json
import pandas as pd
import numpy as np

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
    cfg, _ = compareConfig(iniFile, 'General', compare)

    return cfg


def getModuleConfig(module, fileOverride='', modInfo=False, toPrint=True, onlyDefault=False):
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
    onlyDefault: bool
        if True, only use the default configuration

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
        fileOverride = fU.checkPathlib(fileOverride)
        if fileOverride.is_file():
            iniFile = [defaultFile, fileOverride]
            compare = True
        else:
            raise FileNotFoundError('Provided fileOverride does not exist: ' +
                                    str(fileOverride))

    elif localFile.is_file() and not onlyDefault:
        iniFile = localFile
        iniFile = [defaultFile, localFile]
        compare = True
    elif defaultFile.is_file():
        iniFile = defaultFile
        compare = False
    else:
        raise FileNotFoundError('None of the provided cfg files exist ')

    # Finally read it
    cfg, modDict = compareConfig(iniFile, modName, compare, modInfo, toPrint)

    if modInfo:
        return cfg, modDict

    return cfg


def getDefaultModuleConfig(module, toPrint=True):
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
    cfg, _ = compareConfig(defaultFile, modName, compare=False, toPrint=toPrint)

    return cfg


def compareConfig(iniFile, modName, compare, modInfo=False, toPrint=True):
    ''' Compare configuration files (if a local and default are both provided)
    and inform user of the eventuel differences. Take the default as reference.

    Inputs:
            -iniFile: path to config file. Only one path if compare=False
            -compare: True if two paths are provided and a comparison is needed
            -modInfo: True if dictionary with modifications shall be returned
            -toPrint: True print configuration to terminal

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
            if defCfg.has_section(section):
                for key in locCfg.items(section):
                    # an exception is made for thickness values that are added for the features of a releaseScenario,
                    # entrainment Scenario or secondar. release scenario
                    # these are added to the configuration and also to the modDict if variation is applied
                    validItems = ['entrainmentScenario', 'DEM', 'secondaryReleaseScenario']
                    searchItems = ['relTh', 'entTh', 'secondaryRelTh']
                    if any(s in key[0] for s in searchItems) or key[0] in validItems:
                        locValue = locCfg.get(section, key[0])
                        cfg.set(section, key[0], locValue)
                        log.debug('\t\t%s : %s added to %s' % (key[0], locValue, section))
                        if '$' in locValue:
                            modString = [locValue, locValue.split('$')[0]]
                            modDict[section][key[0]] = modString
                    else:
                        log.warning('Additional Key [\'%s\'] in section [\'%s\'] is ignored.' % (key[0], section))
            else:
                cfg.add_section(section)
                log.info('Additional section [\'%s\'] is added to the configuration.' % (section))
                for key in locCfg.items(section):
                    log.info('Additional Key [\'%s\'] in section [\'%s\'] is added to the configuration.' %
                             (key[0], section))
                    cfg.set(section, key[0], key[1])
                    log.info('\t\t%s : %s', key[0], key[1])

    else:
        log.info('Reading config from: %s', iniFile)
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        # Finally read it
        cfg.read(iniFile)
        # Write config to log file
        if toPrint:
            logUtils.writeCfg2Log(cfg, modName)

    return cfg, modDict


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
    name = pathlib.Path(module.__file__).name
    modName = name.split('.')[0]

    # write to file
    if fileName != '':
        # set outputs
        outDir = pathlib.Path(avaDir, 'Outputs', modName, 'configurationFiles')
        fU.makeADir(outDir)
        cfg.optionxform = str
        with open(pathlib.Path(outDir, '%s.ini' % (fileName)), 'w') as conf:
            cfg.write(conf)
    else:
        # set outputs
        outDir = pathlib.Path(avaDir, 'Outputs')
        cfg.optionxform = str
        with open(pathlib.Path(outDir, '%s_settings.ini' % (modName)), 'w') as conf:
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
        name = pathlib.Path(module.__file__).name
        modName = name.split('.')[0]
        # set input file
        inFile = pathlib.Path(avaDir, 'Outputs', '%s_settings.ini' % (modName))
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


def createConfigurationInfo(avaDir, standardCfg='', writeCSV=False, specDir=''):
    """ Read configurations from all simulations configuration ini files from directory

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        standardCfg: dict
            standard configuration for module - option
        writeCSV: bool
            True if configuration dataFrame shall be written to csv file
        specDir: str
            path to a directory where simulation configuration files can be found - optional

        Returns
        --------
        simDF: pandas DataFrame
            DF with all the simulation configurations
    """

    # collect all configuration files for this module from directory
    if specDir != '':
        inDir = pathlib.Path(specDir, 'configurationFiles')
    else:
        inDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'configurationFiles')
    configFiles = inDir.glob('*.ini')

    if not inDir.is_dir():
        message = 'configuration file directory not found: %s' % (inDir)
        log.error(message)
        raise NotADirectoryError(message)
    elif configFiles == []:
        message = 'No configuration file found in: %s' % (inDir)
        log.error(message)
        raise FileNotFoundError(message)

    # create confiparser object, convert to json object, write to dataFrame
    # append all dataFrames
    simDF = ''
    for cFile in configFiles:
        if 'sourceConfiguration' not in str(cFile):
            simName = pathlib.Path(cFile).stem
            if '_AF_' in simName:
                nameParts = simName.split('_AF_')
                infoParts = nameParts[1].split('_')

            else:
                nameParts = simName.split('_')
                infoParts = nameParts[1:]
            simHash = infoParts[0]
            cfgObject = readCfgFile(avaDir, fileName=cFile)
            simDF = appendCgf2DF(simHash, simName, cfgObject, simDF)

    # convert numeric parameters to numerics
    simDF = convertDF2numerics(simDF)

    # add default configuration
    if standardCfg != '':
        # read default configuration of this module
        simDF = appendCgf2DF('current standard', 'current standard', standardCfg, simDF)

    # if writeCSV, write dataFrame to csv file
    if writeCSV:
        writeAllConfigurationInfo(avaDir, simDF, specDir=specDir)

    return simDF


def appendCgf2DF(simHash, simName, cfgObject, simDF):
    """ append simulation configuration to the simulation dataframe
        only account for sections GENERAL and INPUT

        Parameters
        -----------
        simHash: str
            hash of the simulation to append
        simName: str
            name of the simulation
        cfgObject: configParser
            configuration coresponding to the simulation
        simDF: pandas dataFrame
            configuration dataframe

        Returns
        --------
        simDF: pandas DataFrame
            DFappended with the new simulation configuration
    """
    indexItem = [simHash]
    cfgDict = convertConfigParserToDict(cfgObject)
    simItemDFGeneral = pd.DataFrame(data=cfgDict['GENERAL'], index=indexItem)
    simItemDFInput = pd.DataFrame(data=cfgDict['INPUT'], index=indexItem)
    simItemDF = pd.concat([simItemDFGeneral, simItemDFInput], axis=1)
    simItemDF = simItemDF.assign(simName=simName)
    if isinstance(simDF, str):
        simDF = simItemDF
    else:
        simDF = pd.concat([simDF, simItemDF], axis=0)
    return simDF


def appendTcpu2DF(simHash, tCPU, tCPUDF):
    """ append Tcpu dictionary to the dataframe

        Parameters
        -----------
        simHash: str
            hash of the simulation corresponding to the tCPU dict to append
        tCPU: dict
            cpu time dict of the simulation
        tCPUDF: pandas dataFrame
            tCPU dataframe

        Returns
        --------
        simDF: pandas DataFrame
            DFappended with the new simulation configuration
    """
    indexItem = [simHash]
    tCPUItemDF = pd.DataFrame(data=tCPU, index=indexItem)
    if isinstance(tCPUDF, str):
        tCPUDF = tCPUItemDF
    else:
        tCPUDF = pd.concat([tCPUDF, tCPUItemDF], axis=0)
    return tCPUDF


def convertDF2numerics(simDF):
    """ convert a string DF to a numerical one

        Parameters
        -----------
        simDF: pandas dataFrame
            dataframe

        Returns
        --------
        simDF: pandas DataFrame
    """

    for name, values in simDF.iteritems():
        simDFTest = simDF[name].str.replace('.', '', regex=True)
        # allow for - sign too
        simDFTest = simDFTest.replace('-', '', regex=True)
        # also include columns where nan is in first row - so check for any row
        if simDFTest.str.isdigit().any():
            # problem here is that it finds even if not present in | although not in ini
            simDFTest = simDF[name].str.replace('|', '§', regex=True)
            if simDFTest.str.contains('§').any() == False:
                simDF[name] = pd.to_numeric(simDF[name])
                log.debug('Converted to numeric %s' % name)
        else:
            log.debug('Not converted to numeric: %s' % name)

    return simDF


def readAllConfigurationInfo(avaDir, specDir=''):
    """ Read allConfigurations.csv file as dataFrame from directory

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        specDir: str
            path to a directory where simulation configuration files can be found - optional

        Returns
        --------
        simDF: pandas DataFrame
            DF with all the simulation configurations
        simDFName: array
            simName column of the dataframe
    """

    # collect all configuration files for this module from directory
    if specDir != '':
        inDir = pathlib.Path(specDir, 'configurationFiles')
    else:
        inDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'configurationFiles')
    configFiles = inDir / 'allConfigurations.csv'

    if configFiles.is_file():
        with open(configFiles, 'rb') as file:
            simDF = pd.read_csv(file, index_col=0, keep_default_na=False)
        simDFName = simDF['simName'].to_numpy()
    else:
        simDF = None
        simDFName = []

    return simDF, simDFName


def writeAllConfigurationInfo(avaDir, simDF, specDir=''):
    """ Write cfg configuration to allConfigurations.csv

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        simDF: pandas dataFrame
            daaframe of the configuration
        specDir: str
            path to a directory where simulation configuration shal be saved - optional

        Returns
        --------
        configFiles: pathlib Path
            path where the configuration dataframe was saved
    """

    # collect all configuration files for this module from directory
    if specDir != '':
        inDir = pathlib.Path(specDir, 'configurationFiles')
    else:
        inDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'configurationFiles')
    configFiles = inDir / 'allConfigurations.csv'

    simDF.to_csv(configFiles)

    return configFiles


def filterSims(avalancheDir, parametersDict, specDir=''):
    """ Filter simulations using a list of parameters and a pandas dataFrame of simulation configurations
        if ~ is used as a prefix for a parameter - it is filtered according to values that do NOT match the value
        provided with the ~Parameter

        Parameters
        -----------
        avalancheDir: str
            path to avalanche directory
        parametersDict: dict
            dictionary with parameter and parameter values for filtering
        specDir: str
            path to a directory where simulation configuration files can be found - optional

        Returns
        --------
        simNameList: list
            list of simNames that match filtering criteria
    """

    # load dataFrame for all configurations
    simDF = createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False, specDir=specDir)

    # filter simulations all conditions in the parametersDict have to be met
    if parametersDict != '':
        for key, value in parametersDict.items():
            # first check if values are valid
            if value != '' and value != []:
                # convert values to list
                if not isinstance(value, (list, np.ndarray)):
                    value = [value]

                # remove non matching simulations from simDF
                if key in ['relTh', 'entTh', 'secondaryRelTh', '~relTh', '~entTh', '~secondaryRelTh']:
                    simDF = filterCom1DFAThicknessValues(key, value, simDF)
                else:
                    simDF = removeSimsNotMatching(simDF, key, value)
            else:
                log.debug('Parameter %s is not used for filtering as no valid value is provided: %s' % (key, value))

    # list of simNames after filtering
    simNameList = simDF['simName'].tolist()
    return simNameList


def removeSimsNotMatching(simDF, key, value):
    """ remove simulations from simDF that do not match filtering critera

        Parameters
        -----------
        simDF: pandas dataframe
            dataframe with one row per simulation and info on its characteristics, parameters used,..
        key: str
            name of parameter that shall be used for filtering
        value: list
            list of parameter values used for filtering

        Returns
        ---------
        simDF: pandas dataframe
            updated dataframe with only those simulations that match filtering criteria
    """

    # check if negation in filtering criteria
    notIn = False
    if '~' in key:
        # only add simulations that do not match the value of ~key
        key = key.replace("~", "")
        notIn = True

    # only keep simulations in simDF that match filtering criteria
    if isinstance(value[0], str):
        if notIn:
            simDF = simDF[~simDF[key].isin(value)]
        else:
            simDF = simDF[simDF[key].isin(value)]
    else:
        # if float comparison allow for tolerance
        filterMask = np.isclose(simDF[key].values.reshape(-1, 1), value, atol=1.e-7, rtol=1.e-8).any(axis=1)
        if notIn:
            simDF = simDF[~filterMask]
        else:
            simDF = simDF[filterMask]

    return simDF


def orderSimFiles(avalancheDir, inputDir, varParList, ascendingOrder, specDir='', resFiles=False):
    """ Filter simulations results using a list of parameters and a flag if in ascending or descending order

        Parameters
        -----------
        avalancheDir: str
            path to avalanche directory
        inputDir: str
            path to simulation results
        varParList: str or list
            simulation configuration parameters for ordering simulations
        ascendingOrder: bool
            True if simulations shall be ordered in ascending order regarding varPar
        specDir: str
            path to a directory where simulation configuration files can be found - optional

        Returns
        --------
        dataDF: pandas dataFrame
            dataFrame of simulation results (fileName, ... and values for parameters in varParList)
    """

    # load dataFrame for all configurations
    simDF = createConfigurationInfo(avalancheDir, specDir=specDir)

    # make sure that parameters used for ordering are provided as list
    if isinstance(varParList, str):
        varParList = [varParList]

    if resFiles:
        # create dataframe for simulation results in inputDir
        dataDF = fU.makeSimDF(inputDir)
        # append 'simName' for merging of dataframes according to simNames
        columnNames = ['simName'] + varParList
        # merge varParList parameters as columns to dataDF for matching simNames
        dataDFNew = dataDF.merge(simDF[columnNames], left_on='simName',
                                 right_on='simName')
    else:
        dataDFNew = simDF

    # sort according to varParList and ascendingOrder flag
    dataDFNew = dataDFNew.sort_values(by=varParList, ascending=ascendingOrder)

    return dataDFNew


def filterCom1DFAThicknessValues(key, value, simDF):
    """ thickness settings different if read from shpfile - requires more complex filtering
        if read from shp - thickness values are provided per feature!!
        for example relTh = '' but relTh0 = 1 is appended for feature with id 0, relTh1 for feature
        with id 1, etc.

        Parameters
        -----------
        key: str
            name of parameter
        value: list
            list of values used for filtering
        simDF: pandas dataframe
            configuration info for each simulation

        Returns
        --------
        simDF: pandas data frame
            updated dataframe
    """

    # check if filter for values that do NOT match criteria
    notIn = False
    if '~' in key:
        key = key.split('~')[1]
        notIn = True

    # create required parameters for searching
    thFlag = key + 'FromShp'
    thId = key + 'Id'
    thThickness = key + 'Thickness'
    thPercentVariation = key + 'PercentVariation'

    # append identifier if simulation matches thickness filter criteria
    simDF['toBeAdded'] = False

    # initialize list for thickness parameter names (according to thickness configuration -
    # e.g. mutiple features)
    allThNames = []
    # loop over simDF and set identifier if filter criteria are matched
    for simHash, simDFrow in simDF.iterrows():
        if simDFrow[thFlag] == 'True':
            # inititialise thickness ids and thickness parameter names if thickness read from shp
            thIdList = str(simDFrow[thId]).split('|')
            thNames = [(key + id) for id in thIdList]
            allThNames = allThNames + thNames
            log.warning('Filtering applied for %s - multiple features found as %s was read \
                from shp file - only simulations where all features match %s will be added' %
                (key, key, value))
        else:
            # if thickness read from ini add thickness parameter name
            thIdList = [0]
            thNames = [key]
            allThNames = allThNames + [key]
        # check if filter criteria are met by thickness parameters for the sim in simDFrow
        for val in value:
            if (simDFrow[thNames].values == [val] * len(thIdList)).all():
                simDF.loc[simHash, 'toBeAdded'] = True

    # get a list with all thickness parameters included in search
    allThNames = list(set(allThNames))
    if notIn:
        # return all sims that do not match filter criteria
        simDF = simDF[simDF['toBeAdded'] == False]
    else:
        # return all sims that do match filter criteria
        simDF = simDF[simDF['toBeAdded'] == True]

    log.info('simulations for %s found with values: %s' % (key, simDF[allThNames]))

    return simDF
