"""
    Defining a writable object to write the config to the log file

    This file is part of Avaframe.
"""
import logging
import logging.config
import os
import pathlib
from pathlib import Path, PureWindowsPath
from datetime import datetime
import avaframe.version as gv


def writeCfg2Log(cfg, cfgName='Unspecified'):
    """ write a configparser object to log file"""

    log = logging.getLogger('avaframe')

    log.debug('Writing cfg for: %s', cfgName)
    for section in cfg.sections():
        log.info('\t%s', section)
        for key in dict(cfg[section]):
            log.info('\t\t%s : %s', key, cfg[section][key])


def initiateLogger(targetDir, logName='runLog'):
    """ Initiates logger object based on setup in logging.conf

        Parameters
        ----------
        targetDir: str
            folder to save log file to
        logName : str
            filename of log file; optional; defaults to runLog.log

        Returns
        -------
        log : logging object

    """

    # datetime object containing current date and time
    now = datetime.now()
    dtString = now.strftime("%d.%m.%Y_%H:%M:%S")

    logFileName = pathlib.Path(targetDir, logName+'_'+dtString+'.log')

    # get path of module and generate logging.conf file path
    logConfPath = pathlib.Path(__file__).parents[0]
    logConfFile = logConfPath / 'logging.conf'

    logging.config.fileConfig(fname=logConfFile,
                              defaults={'logfilename': str(logFileName).replace('\\', '/')},
                              disable_existing_loggers=False)
    log = logging.getLogger('avaframe')

    log.info('Started logging at: %s', dtString)
    log.info('Also logging to: %s', logFileName)
    log.info('AvaFrame Version: %s', gv.getVersion())

    return log
