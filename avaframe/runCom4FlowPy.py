"""
Run script for module com4FlowPy
"""
# Load modules
import pathlib
import os
from datetime import datetime

# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.com4FlowPy import com4FlowPy
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initialiseDirs as inDirs

def main():
    # log file name; leave empty to use default runLog.log
    logName = 'runcom4FlowPy'
    # Read main Config and com4FlowPy config from files
    cfgMain = cfgUtils.getGeneralConfig()
    cfg = cfgUtils.getModuleConfig(com4FlowPy)

    cfgSetup = cfg['SETUP']
    cfgFlags = cfg['FLAGS']
    cfgCustomPaths = cfg['PATHS']

    if cfgCustomPaths["useCustomPaths"]=='False':
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

    elif cfgCustomPaths["useCustomPaths"]=='True':
        cfgPath = {}

        # Handling Custom directory creation
        workDir = pathlib.Path(cfgCustomPaths['workDir'])
        
        time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        try:
            os.makedirs(workDir / 'res_{}'.format(time_string))
            res_dir = (workDir / 'res_{}'.format(time_string))
        except FileExistsError:
            res_dir = (workDir / 'res_{}'.format(time_string))

        try:
            os.makedirs(workDir / res_dir / 'temp')
            temp_dir = (workDir / res_dir / 'temp')
        except FileExistsError:
            temp_dir = (workDir / res_dir / 'temp')


        cfgPath["workDir"]=pathlib.Path(workDir)
        cfgPath["outDir"]=pathlib.Path(res_dir)
        cfgPath["tempDir"]=pathlib.Path(temp_dir)
        cfgPath["demPath"]=pathlib.Path(cfgCustomPaths["demPath"])
        cfgPath["releasePath"]=pathlib.Path(cfgCustomPaths["releasePath"])
        cfgPath["infraPath"]=pathlib.Path(cfgCustomPaths["infraPath"])

        log = logUtils.initiateLogger(cfgPath["outDir"], logName)

        cfgSetup["cpuCount"] = str(cfgUtils.getNumberOfProcesses(cfgMain,999))
        cfgPath["customDirs"] = cfgCustomPaths["useCustomPaths"]

        com4FlowPy.com4FlowPyMain(cfgPath, cfgSetup)

    else:
        log.info('INPUT SETTINGS incorrect - abort')
        sys.exit(1)
        pass

if __name__ == '__main__':
    main()

