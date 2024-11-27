import pathlib

import avaframe.in3Utils.initializeProject as initProj
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils, cfgHandling
from avaframe.in3Utils import logUtils

def processAvaDirCom1DFA(cfgMain, cfgCom7, avalancheDir):
    """Function to call com1DFA in each avalanche directory with regional override settings."""
    avalancheDir = pathlib.Path(avalancheDir)

    # Initialize log for each process
    log = logUtils.initiateLogger(avalancheDir, logName='runCom1DFA')
    log.info('PROCESS CALLED BY REGIONAL RUN')
    log.info('Current avalanche: %s', avalancheDir)

    # Update cfgMain to reflect the current avalancheDir
    cfgMain['MAIN']['avalancheDir'] = str(avalancheDir)

    # Clean input directory of old work and output files from module
    initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)

    # Create com1DFA configuration for the current avalanche directory and override with regional settings
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False,
                                          onlyDefault=cfgCom7['com1DFA_com1DFA_override'].getboolean('defaultConfig'))
    cfgCom1DFA, cfgCom7 = cfgHandling.applyCfgOverride(cfgCom1DFA, cfgCom7, com1DFA, addModValues=False)

    # Run com1DFA in the current avalanche directory
    com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgCom1DFA)

    return avalancheDir, "Success"