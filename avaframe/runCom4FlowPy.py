"""Run script for module com4FlowPy
"""

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.com4FlowPy import com4FlowPy
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def main():
    # log file name; leave empty to use default runLog.log
    logName = 'runcom4FlowPy'
    
    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']
    # Clean input directory of old work and output files from module
    initProj.cleanModuleFiles(avalancheDir, com4FlowPy, deleteOutput=False)
    
    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)
    
    # Load all input Parameters from config file
    # get the configuration of an already imported module
    # write config to log file
    cfg = cfgUtils.getModuleConfig(com4FlowPy)
    
    cfgSetup = cfg['SETUP']
    cfgFlags = cfg['FLAGS']
    
    # Extract input file locations
    cfgPath = com4FlowPy.readFlowPyinputs(avalancheDir, cfg)
    
    com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)

if __name__ == '__main__':
    main()

