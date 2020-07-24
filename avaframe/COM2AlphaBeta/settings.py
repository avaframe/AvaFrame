# General setting file
# will be overridden by local_settings.py
# So copy this file to local_settings.py, adjust your variables and add that
# file to gitignore
# This file serves as template file and general settings file

# FSO--- Get full Samos output or not
fullOut = False

# FSO--- which path to use
ProfileLayer = 'YOUR PATH TO AVALANCHE PATH'
DGMSource = 'YOUR PATH TO DEM'
SplitPointSource = 'YOUR PATH TO SPLITPOINT'
try:
    from local_settings import *
except ImportError:
    pass
