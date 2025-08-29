"""
Defining a writable object to write the config to the log file

This file is part of Avaframe.
"""

import logging
import logging.config
import pathlib
from datetime import datetime
import avaframe.version as gv
from contextlib import contextmanager


def writeCfg2Log(cfg, cfgName="Unspecified"):
    """write a configparser object to log file"""

    log = logging.getLogger("avaframe")

    log.debug("Writing cfg for: %s", cfgName)
    for section in cfg.sections():
        log.info("\t%s", section)
        for key in dict(cfg[section]):
            log.info("\t\t%s : %s", key, cfg[section][key])


def initiateLogger(targetDir, logName="runLog", modelInfo=""):
    """Initiates logger object based on setup in logging.conf

    Parameters
    ----------
    targetDir: str
        folder to save log file to
    logName : str
        filename of log file; optional; defaults to runLog.log
    modelInfo: str
        if not emtpy add info on modelInfo to log

    Returns
    -------
    log : logging object

    """

    # datetime object containing current date and time
    now = datetime.now()
    dtString = now.strftime("%Y%m%d_%Hh%Mm%Ss")

    logFileName = pathlib.Path(targetDir, logName + "_" + dtString + ".log")

    # get path of module and generate logging.conf file path
    logConfPath = pathlib.Path(__file__).parents[0]
    logConfFile = logConfPath / "logging.conf"

    # clean logFileName for special characters
    # very hacky way, but no better way found working for both linux and windows
    # Add further special characters here as needed...
    cleanName = str(logFileName)
    cleanName = cleanName.replace("'", "\\'")
    cleanName = cleanName.replace("\\", "/")

    logging.config.fileConfig(
        fname=logConfFile,
        defaults={"logfilename": cleanName},
        disable_existing_loggers=False,
    )
    log = logging.getLogger("avaframe")

    log.info("Started logging at: %s", dtString)
    log.info("Also logging to: %s", logFileName)
    log.info("AvaFrame Version: %s", gv.getVersion())

    if modelInfo != "":
        log.info("Used by: %s" % modelInfo)

    return log


@contextmanager
def silentLogger(loggerName="avaframe"):
    """Context manager to temporarily silence a logger.

    Parameters
    ----------
    loggerName : str
        Name of the logger to silence. Defaults to 'avaframe'

    Example
    -------
    >>> with silentLogger():
    >>>     # This code block will not produce any log output
    >>>     function_with_noisy_logs()
    """
    logger = logging.getLogger(loggerName)
    # Store the original level
    originalLevel = logger.level
    # Temporarily set the level to ERROR (will suppress INFO and DEBUG messages)
    logger.setLevel(logging.ERROR)
    try:
        yield
    finally:
        # Restore the original level
        logger.setLevel(originalLevel)
