"""
    Run script for running the dam break problem on a sloping bed and compare to simulation results
"""

# Load modules
import pathlib
from configupdater import ConfigUpdater

# Local imports
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runDamBreakProblem'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaDamBreak'
cfgMain['MAIN']['avalancheDir'] = avalancheDir

# Clean input directory(ies) of old work and output files
# initProj.cleanModuleFiles(avalancheDir, com1DFA)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration
damBreakCfg = pathlib.Path(avalancheDir, 'Inputs', 'damBreak_com1DFACfg.ini')

# Overwriting the sphKernelRadius and meshCellSize from the ini file with the values specified in []
for sphKernelRadius in [5, 3]:
    updater = ConfigUpdater()
    updater.read(damBreakCfg)
    updater['GENERAL']['sphKernelRadius'].value = sphKernelRadius
    updater['GENERAL']['meshCellSize'].value = sphKernelRadius
    updater.update_file()

    # call com1DFA to perform simulation - provide configuration file and release thickness function
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=damBreakCfg)
