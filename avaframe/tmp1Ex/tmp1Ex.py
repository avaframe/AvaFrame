"""
This is the template for new modules, with the bare minimal required files
"""

import logging

# Local import
import avaframe.tmp1Ex.tmp1ExCfg as conf

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def tmp1ExMain():
    """ Main function for module tmp1Example"""

    print('In tmp1Example')
    log.debug('Here')
    log.debug('Input directory %s', conf.inputDir)
