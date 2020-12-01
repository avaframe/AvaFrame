"""
This is the template for new modules, with the bare minimal required files
"""

import logging
import avaframe
from avaframe.out3Plot.plotUtils import *

# create local logger
log = logging.getLogger(__name__)


def tmp1ExMain(cfg):
  """Main function for module tmp1Example

  Args:
    foo (int): The foo to bar
    bar (str): Bar to use on foo
    baz (float): Baz to frobnicate

  Returns:
    float: The frobnicated baz
  """

  print('In tmp1Example')
  log.info('Input directory %s', cfg['GENERAL']['inputDir'])

def some_func(foo, bar, baz):
  """Does some stuff

  Args:
    foo (int): The foo to bar
    bar (str): Bar to use on foo
    baz (float): Baz to frobnicate

  Returns:
    float: The frobnicated baz
  """
  print('Hello')

