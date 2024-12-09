#RunScarp.py

#from scarp import run_scarp_analysis
from com6RockAvalanche.scarp import runScarpAnalysis

# Specify the configuration file path
configFile = 'com6RockAvalanche/scarpCfg.ini'

# Call the scarp analysis function
runScarpAnalysis(configFile)
