"""
    Check if the .ini file is consistent for com1DFA computations
"""

import logging

# create local logger
log = logging.getLogger(__name__)


def checkCfgConsistency(cfg):
    """ Check the provided configuration for necessary consistency/relation between
        parameters for com1DFA.

    Parameters
    -----------
    cfg : configuration object
        com1DFA configuration
    Returns
    --------
    True if checks are passed, otherwise error is thrown
    """

    # Check if Ata Parameters are consistent
    sphOption = float(cfg['GENERAL']['sphOption'])
    viscOption = float(cfg['GENERAL']['viscOption'])

    if viscOption == 2:
        if sphOption != 2:
            # raise an error
            message = ('If viscOption is set to 2 (Ata viscosity), sphOption = 2 is needed '
                       '(or implement the Ata viscosity for other sphOption values)')
            log.error(message)
            raise AssertionError(message)
        else:
            log.debug('Ata parameters are consistent')

    return True
