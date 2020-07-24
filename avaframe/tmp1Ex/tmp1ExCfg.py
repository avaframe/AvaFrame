"""
Config file
will be overridden by local_tmp1ExampleConf.py
So copy this file to local_tmp1ExmpleConf.py, adjust your variables there
This file serves as template file and general settings file
"""

# ##########################
# Required settings
# ##########################

# Input directory containing avalanche data
inputDir = 'Path to avalanche directory'

# path to work directory
workDir = 'Path to working directory'

# ##########################
# Optional settings
# ##########################

# Some example flag
fullLog = True

###########################
# Do not change anything below this line
###########################
try:
    from local_tmp1ExampleConf import *
except ImportError:
    pass
