'''
    Utilities for handling configuration files

    This file is part of Avaframe.

'''

import configparser
import os
import logging
log = logging.getLogger(__name__)


def getGeneralConfig():
    ''' Returns the general configuration for avaframe
    returns a configParser object
    '''

    # initialize a configparser object
    cfg = configparser.ConfigParser()

    # get path of module
    modPath = os.path.dirname('avaframeCfg.ini')
    iniFile = os.path.join(modPath, 'avaframeCfg.ini')

    # Finally read it
    cfg.read(iniFile)

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
    modPath = os.path.dirname(module.__file__)

    # get filename of module
    modName = os.path.basename(module.__file__)
    modName = os.path.splitext(modName)[0]

    localFile = os.path.join(modPath, 'local_'+modName+'Cfg.ini')
    defaultFile = os.path.join(modPath, modName+'Cfg.ini')

    log.debug('localFile: %s', localFile)
    log.debug('defaultFile: %s', defaultFile)

    # Decide which one to take
    if fileOverride:
        if os.path.isfile(fileOverride):
            iniFile = fileOverride
            compare = False
        else:
            raise FileNotFoundError('Provided fileOverride does not exist')

    elif os.path.isfile(localFile):
        iniFile = localFile
        iniFile = [defaultFile, localFile]
        compare = True
    elif os.path.isfile(defaultFile):
        iniFile = defaultFile
        compare = False
    else:
        raise FileNotFoundError('None of the provided cfg files exist ')

    # Finally read it
    cfg = compareConfig(iniFile, compare)

    return cfg

def compareConfig(iniFile, compare):
    ''' Compare configuration files (if a local and default are both provided)
    and inform user of the eventuel differences. Take the default as reference.

    Inputs:
            -iniFile: path to config file. Only one path if compare=False
            -compare: True if two paths are provided and a comparison is needed

    Output: ConfigParser object
    '''
    if compare:
        log.info('Reading config from: %s and %s'% (iniFile[0], iniFile[1]))
        # initialize our final configparser object
        cfg = configparser.ConfigParser()
        # initialize configparser object to read
        defCfg = configparser.ConfigParser()
        locCfg = configparser.ConfigParser()
        # read default and local parser files
        defCfg.read(iniFile[0])
        locCfg.read(iniFile[1])
        # loop through all sections of the localCfg
        for section in locCfg.sections():
            # check if section is also in the default cfg
            if defCfg.has_section(section):
                # if yes add this section to the cfg that will be returned
                cfg.add_section(section)
                for key in locCfg.items(section):
                    # check if key is also in the default cfg
                    if defCfg.has_option(section, key[0]):
                        # if yes add this key to the cfg that will be returned
                        cfg.set(section, key[0], locCfg.get(section, key[0]))
                        # remove the key from the default cfg
                        defCfg.remove_option(section, key[0])
                    else :
                        # if not warn the user that this key no longer exists
                        log.warning('Key [\'%s\'] in section [\'%s\'] in config file %s no longer exists.' % (key[0], section, iniFile[1]))
                # remove the section from the default cfg if it is empty
                if not defCfg.items(section):
                    defCfg.remove_section(section)
            else :
                # if not warn the user that this section no longer exists
                log.warning('Section [\'%s\'] in config file %s no longer exists.' % (section, iniFile[1]))

        # Now check if there are some sections/ keys left in the default cfg and add them to the Cfg
        # loop through all sections of the defaultCfg
        for section in defCfg.sections():
            # check if section is already in the cfg
            if not cfg.has_section(section):
                # if no add this section to the cfg that will be returned and warn the user
                cfg.add_section(section)
                log.warning('Section [\'%s\'] added because it is not present in the config file %s but should be.' % (section, iniFile[1]))
            # loop on key in defaultCfg
            for key in defCfg.items(section):
                # add this key to the cfg that will be returned and warn the user
                cfg.set(section, key[0], defCfg.get(section, key[0]))
                log.warning('Key [\'%s\'] in section [\'%s\'] added because it is not present in the config file %s but should be.' % (key[0], section, iniFile[1]))

    else:
        log.info('Reading config from: %s', iniFile)
        cfg = configparser.ConfigParser()
        # Finally read it
        cfg.read(iniFile)

    return cfg
