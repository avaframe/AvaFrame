"""
    Run script for generateTopo in module in3Utils
    This file is part of Avaframe.
"""

import sys
import logging

# Local imports
import avaframe as af
from avaframe.in3Utils import generateTopo as gT

# create logger, set to logging.DEBUG to see all messages
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(module)s:%(levelname)s - %(message)s')

gT.generateTopo()
