"""
    Run script for running python DFA kernel
"""

import pathlib
import shutil
from configupdater import ConfigUpdater

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


""" run com1DFA module """

# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runCom1DFA'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Clean input directory of old work and output files from module
initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=False)
# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)
localCfg = pathlib.Path('com1DFA', 'local_com1DFACfg.ini')


for sphKernelRadius in [10, 8, 6, 4, 5]:
    # update cfg
    updater = ConfigUpdater()
    updater.read(localCfg)
    updater['GENERAL']['sphKernelRadius'].value = sphKernelRadius
    updater['GENERAL']['meshCellSize'].value = sphKernelRadius
    updater.update_file()

    # copy dem with correct cell size to inputs
    demFile = pathlib.Path(avalancheDir, 'Inputs' , 'DEM' , 'DEM_PF_Topo' + str(round(sphKernelRadius)) + '.asc')
    demDestination = pathlib.Path(avalancheDir, 'Inputs' , 'DEM_PF_Topo.asc')
    if demFile.is_file():
        if demDestination.is_file():
            print('file aleready exist')
            # demDestination.unlink()
        shutil.copy(str(demFile), str(demDestination))
    # call com1DFA and perform simulations
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile='')

    # put back the original dem
    demFile = pathlib.Path(avalancheDir, 'Inputs' , 'DEM' , 'DEM_PF_Topo.asc')
    demDestination = pathlib.Path(avalancheDir, 'Inputs' , 'DEM_PF_Topo.asc')
    if demFile.is_file():
        if demDestination.is_file():
            print('file aleready exist')
            # demDestination.unlink()
        shutil.copy(str(demFile), str(demDestination))
