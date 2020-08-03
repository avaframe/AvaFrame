"""
Config file
will be overridden by local_com2ABCfg.py
So copy this file to local_com2ABCfg.py, adjust your variables there
This file serves as template file and general settings file
"""
import os
dirname = os.path.dirname(__file__)
# ##########################
# Required settings
# ##########################

# 'YOUR PATH TO AVALANCHE PATH'
ProfileLayer = os.path.join(dirname, '../data/avaSlide/LINES/SlideProfile3')

DGMSource = os.path.join(dirname, '../data/avaSlide/slideTopo.asc')  # 'YOUR PATH TO DEM'

SplitPointSource = os.path.join(
    dirname, '../data/avaSlide/POINTS/slidePoints3')  # 'YOUR PATH TO SPLITPOINT'

saveOutPath = os.path.join(dirname, 'Outputs/')
outputName = 'SlideProfile3'

# ##########################
# Optional settings
# ##########################

fullLog = True
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

# plot save results
flags = {}
flags['PlotRes'] = True
flags['SavePlot'] = True
flags['WriteRes'] = True


###########################
# Do not change anything below this line
###########################
try:
    from local_com2ABCfg import *
except ImportError:
    pass
