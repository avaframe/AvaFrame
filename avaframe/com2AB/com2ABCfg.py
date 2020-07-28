"""
Config file
will be overridden by local_tmp1ExampleConf.py
So copy this file to local_tmp1ExmpleConf.py, adjust your variables there
This file serves as template file and general settings file
"""

# ##########################
# Required settings
# ##########################

ProfileLayer = 'YOUR PATH TO AVALANCHE PATH'

DGMSource = 'YOUR PATH TO DEM'

SplitPointSource = 'YOUR PATH TO SPLITPOINT'

saveOutPath = ''

# ##########################
# Optional settings
# ##########################


# if small avalanche set is wanted
smallAva = False
# if custom avalanche set is wanted
customParam = None
# customParam = {}
# customParam['k1'] = your value
# customParam['k2'] = your value
# customParam['k3'] = your value
# customParam['k4'] = your value
# customParam['SD'] = your value


# resampling step [m]
distance = 10

# Some example flag
fullLog = True

###########################
# Do not change anything below this line
###########################
try:
    from .local_com2ABCfg import *
except ImportError:
    pass
