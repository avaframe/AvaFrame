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
cfgDir = 'anaInfluenceTest'

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for similarity solution test
influenceTestCfg = pathlib.Path('anaInfluenceTest', 'anaInfluenceTestDFACfg.ini')
ABCfg = cfgUtils.getModuleConfig(com2AB)

# Run anaInfluenceTest -> have to decide if it returns a data frame or not
#simDF = anaInfluenceTest.mainAnaInfluenceTest(avalancheDir, cfgMain, DFACfg, ABCfg)
anaInfluenceTest.mainAnaInfluenceTest(avalancheDir, cfgMain, influenceTestCfg, ABCfg)
