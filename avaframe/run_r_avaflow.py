"""
    Run script for module r_avaflow
"""

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.r_avaflow import r_avaflow
#from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


""" run com1DFA module """

# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'run_r_avaflow'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Clean input directory of old work and output files from module
initProj.cleanModuleFiles(avalancheDir, r_avaflow, deleteOutput=True)
# Start logging
# log = logUtils.initiateLogger(avalancheDir, logName)
# log.info('MAIN SCRIPT')
# log.info('Current avalanche: %s', avalancheDir)

# call com1DFA and perform simulations
#dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile='')
r_avaflow.r_avaflow_Main(avalancheDir)