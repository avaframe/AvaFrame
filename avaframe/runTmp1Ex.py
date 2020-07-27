"""Run script for module template"""

import logging
import sys

# Local imports
import avaframe as af
from avaframe.tmp1Ex import tmp1Ex

# create logger, set to logging.DEBUG to see all messages
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(module)s:%(levelname)s - %(message)s')

# Different ways to call functions
tmp1Ex.tmp1ExMain()

# af.tmp1Ex.tmp1Ex.tmp1ExMain()
