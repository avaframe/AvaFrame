"""
    Run script for plotting a matrix dataset with one profile
"""

# Load modules
import glob
import os

# Local imports
from avaframe.out3Plot import outQuickPlot
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

"""
    Run script for creating directory structure for ISeeSnow models
"""
# Load modules
# importing general python modules
import pathlib

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.out3Plot import outQuickPlot as oQ


#+++++++++++ user input
level = 1.
resType = 'pfv'
secDir = ''
#+++++++++++++++++++++++

# fetch input directory
cfgMain = cfgUtils.getGeneralConfig()
avaDir = cfgMain['MAIN']['avalancheDir']

# call contour line plot
oQ.plotAllContours(avaDir, 'com1DFA', resType, level, specDir='')



