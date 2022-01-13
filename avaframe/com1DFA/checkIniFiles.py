"""
    Check if the .ini file is well parametered
"""

import logging
import sys

# create local logger
log = logging.getLogger(__name__)

def checkIniFiles(cfg):
    # Check if Ata Parameters all well seted
    sphOption = float(cfg['GENERAL']['sphOption'])
    viscOption = float(cfg['GENERAL']['viscOption'])
    if viscOption == 2 and sphOption != 2:
        # raise an error
        message = 'If viscOption is set to 2 (Ata viscosity), you have to use sphOption to 2, or implement the Ata voscosity for the other sphOption values'
        log.error(message)
        raise AssertionError(message)
    log.info('Ata Parameters are well seted')
