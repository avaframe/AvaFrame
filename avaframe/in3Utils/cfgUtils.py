'''Utilities for handling configuration files'''

import configparser
import sys
import os
import logging
log = logging.getLogger(__name__)


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

    TODO:
    - pytest
    '''


    # initialize a configparser object
    cfg = configparser.ConfigParser()

    # get path of module
    modPath = os.path.dirname(module.__file__)

    # get filename of module
    modName = os.path.basename(module.__file__)
    modName = os.path.splitext(modName)[0]

    localFile = os.path.join(modPath,'local_'+modName+'Cfg.ini')
    defaultFile = os.path.join(modPath,modName+'Cfg.ini')

    # Decide which one to take
    if fileOverride:
        if os.path.isfile(fileOverride):
            iniFile = fileOverride
        else:
            raise FileNotFoundError('Provided fileOverride does not exist')

    elif os.path.isfile(localFile):
        iniFile = localFile
    elif os.path.isfile(defaultFile):
        iniFile = defaultFile
    else:
        raise FileNotFoundError('None of the provided cfg files exist ')

    log.info('Reading config from: %s', iniFile)

    # Finally read it
    cfg.read(iniFile)

    return cfg


