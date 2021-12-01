"""
    Run script for running module anaInfluenceTest
"""

import pathlib

# local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.com1DFA import com1DFA
from avaframe.com2AB import com2AB
from avaframe.anaInfluenceTest import anaInfluenceTest


# log file name; leave empty to use default runLog.log
logName = 'runAnaInfluenceTest'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Clean input directory of old work and output files from module
initProj.cleanModuleFiles(avalancheDir, anaInfluenceTest)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for similarity solution test
influenceTestCfg = pathlib.Path('anaInfluenceTest', 'anaInfluenceTestDFACfg.ini')
DFACfg = cfgUtils.getModuleConfig(com1DFA, influenceTestCfg)
ABCfg = cfgUtils.getModuleConfig(com2AB)
