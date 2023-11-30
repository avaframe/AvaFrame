"""
Run script for module com4FlowPy
"""

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.com4FlowPy import com4FlowPy
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def main():
    # log file name; leave empty to use default runLog.log
    logName = 'runcom4FlowPy'
    # Read main Config and com4FlowPy config from files
    cfgMain = cfgUtils.getGeneralConfig()
    cfg = cfgUtils.getModuleConfig(com4FlowPy)

    cfgSetup = cfg['SETUP']
    cfgFlags = cfg['FLAGS']
    cfgPaths = cfg['PATHS']

    if cfgPaths["useCustomPaths"]=='False':
        # Load avalanche directory from general configuration file
        avalancheDir = cfgMain['MAIN']['avalancheDir']
        # Clean input directory of old work and output files from module
        initProj.cleanModuleFiles(avalancheDir, com4FlowPy, deleteOutput=False)
        # Extract input file locations
        cfgPath = com4FlowPy.readFlowPyinputs(avalancheDir, cfg)
        # Start logging
        log = logUtils.initiateLogger(avalancheDir, logName)
        log.info('MAIN SCRIPT')
        log.info('Current avalanche: %s', avalancheDir)

        # IMPORTANT!! - this is a quick'n'dirty hack to set nSims to 99, so that
        # min(nCPU,nSims) in cfgUtils.getNumberOfProcesses returns nCPU
        cfgSetup['cpuCount'] = str(cfgUtils.getNumberOfProcesses(cfgMain,999))

        com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)

    elif cfgPaths["useCustomPaths"]=='True':
        print('True')
    else:
        pass

if __name__ == '__main__':
    main()

