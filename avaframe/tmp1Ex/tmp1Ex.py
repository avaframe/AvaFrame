"""
This is the template for new modules, with the bare minimal required files
"""

import logging

# create local logger
log = logging.getLogger(__name__)


def tmp1ExMain(cfg):
    """Main function for module tmp1Example

    Parameters
    ----------
    foo : int, float, str, or tf.Tensor
      The foo to bar, which has a really really, reeeeeeeeeeeeeeeeally
      unnecessarily long multiline description.
    bar : str
      Bar to use on foo
    baz : float
      Baz to frobnicate

    Returns
    -------
    float
      The frobnicated baz
    """

    print('In tmp1Example')
    log.info('Input directory %s', cfg['GENERAL']['inputDir'])

