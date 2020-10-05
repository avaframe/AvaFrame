"""
This is the template for new modules, with the bare minimal required files
"""

import logging

# create local logger
log = logging.getLogger(__name__)


def tmp1ExMain(cfg):
    """ Main function for module tmp1Example"""

    print('In tmp1Example')
    log.info('Input directory %s', cfg['GENERAL']['inputDir'])
