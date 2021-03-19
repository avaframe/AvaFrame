'''
    Utilities for handling configuration files

    This file is part of Avaframe.

'''

import configparser
import os
import logging
import pathlib
# Local imports
import avaframe as avaf
from avaframe.in3Utils import logUtils


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


def getModuleConfig(module, fileOverride=''):
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
            iniFile = fileOverride
            compare = False
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
    cfg = compareConfig(iniFile, modName, compare)

    return cfg


def compareConfig(iniFile, modName, compare):
    ''' Compare configuration files (if a local and default are both provided)
    and inform user of the eventuel differences. Take the default as reference.

    Inputs:
            -iniFile: path to config file. Only one path if compare=False
            -compare: True if two paths are provided and a comparison is needed

    Output: ConfigParser object
    '''
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

    return cfg


def writeCfgFile(avaDir, module, cfg, suffix=''):
    """ Save configuration used to text file in Outputs """

    # write configuration to ini file
    modPath = os.path.dirname(module.__file__)
    # get filename of module
    name = os.path.basename(module.__file__)
    modName = name.split('.')[0]
    # set outputs
    outDir = os.path.join(avaDir, 'Outputs')

    # write to file
    with open(os.path.join(outDir, '%s%s_settings.ini' % (modName, suffix)), 'w') as conf:
        cfg.write(conf)


def readCfgFile(avaDir, module, suffix=''):
    """ Save configuration used to text file in Outputs """

    # write configuration to text file
    modPath = os.path.dirname(module.__file__)
    # get filename of module
    name = os.path.basename(module.__file__)
    modName = name.split('.')[0]
    # set outputs
    inFile = os.path.join(avaDir, 'Outputs', '%s%s_settings.ini' % (modName, suffix))

    # read configParser object
    cfg = configparser.ConfigParser()
    cfg.read(inFile)

    return cfg
