"""
    Defining a writable object to write the config to the log file

    This file is part of Avaframe.
"""
import logging
import logging.config
import os
from pathlib import Path, PureWindowsPath
from datetime import datetime


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
    logFileName = os.path.join(targetDir, logName+'.log')

    # get path of module and generate logging.conf file path
    logConfPath = os.path.dirname(__file__)
    logConfFile = os.path.join(logConfPath, 'logging.conf')

    logging.config.fileConfig(fname=logConfFile,
                              defaults={'logfilename': logFileName.replace('\\', '/')},
                              disable_existing_loggers=False)
    log = logging.getLogger('avaframe')

    # datetime object containing current date and time
    now = datetime.now()
    dtString = now.strftime("%d.%m.%Y %H:%M:%S")
    log.info('Started logging at: %s', dtString)
    log.info('Also logging to: %s', logFileName)

    return log
