"""
    Run script for generateTopo in module in3Utils
    This file is part of Avaframe.
"""

import sys
import logging

# Local imports
from avaframe.in3Utils import generateTopo as gT
from avaframe.out3SimpPlot import outGenerateTopo as oT

# create logger, set to logging.DEBUG to see all messages
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(module)s:%(levelname)s - %(message)s')

# Call main function to generate DEMs
[z, name_ext] = gT.generateTopo()

# Plot new topogrpahy
oT.plotDEM(z, name_ext)
