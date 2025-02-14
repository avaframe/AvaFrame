import os
import subprocess
import platform
import logging

import avaframe.com1DFA.com1DFATools as com1DFATools
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in3Utils import cfgUtils

# create local logger
log = logging.getLogger(__name__)


# TODO move to whereever?
def _runAndCheck(command):
    if os.name == "nt":
        useShell = True
    elif platform.system() == "Darwin":
        useShell = False
    else:
        useShell = False

    # This starts the subprocess
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=useShell,
        encoding="utf-8",
        errors="replace",
        universal_newlines=True,
    )

    while True:
        realtimeOutput = process.stdout.readline()

        if realtimeOutput == "" and process.poll() is not None:
            break

        if realtimeOutput:
            line = realtimeOutput.strip()
            print(line)


# TODO move to utils
def cfgToRcf(cfg, fileName):
    with open(fileName, 'w') as f:
        for section in cfg.sections():
            if section in ("FOREST_EFFECTS", "ENTRAINMENT"):
                pass
            elif section in ("GENERAL"):
                continue
            else:
                f.write(f"# {section.replace('_', ' ')}\n")
                f.write("#\n")
            for key, value in cfg.items(section):
                # key = key.replace('_', ' ')
                key = key.strip()
                f.write(f"{key:<40}{value}\n")
            f.write("#\n")


def com8MoTPSAMain(cfgMain, cfgInfo=None):
    # Get all necessary information from the configuration files

    # Load avalanche directory from general configuration file
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    # fetch type of cfgInfo
    typeCfgInfo = com1DFATools.checkCfgInfoType(cfgInfo)

    # preprocessing to create configuration objects for all simulations to run
    simDict, outDir, inputSimFiles, simDFExisting = com1DFA.com1DFAPreprocess(cfgMain, typeCfgInfo, cfgInfo)

    log.info("The following simulations will be performed")
    for key in simDict:
        log.info("Simulation: %s" % key)
        exportFlag = simDict[key]["cfgSim"]["EXPORTS"].getboolean("exportData")

    # # Generate command and run via subprocess.run
    #
    # # Configuration that needs adjustment
    # cfgInfo['Run information']['Area of Interest'] = cfgMain['MAIN']['avalancheDir'].replace("/", "")
    # # cfgInfo['Run information']['UTM Zone']
    # # cfgInfo['Run information']['EPGS geodetic datum code']
    # cfgInfo['Run information']['Run name'] = cfgMain['MAIN']['avalancheDir'].replace("/", "")
    #
    #
    # rcfFile = 'outputTemplate.rcf'
    #
    # cfgToRcf(cfgInfo, rcfFile)
    #
    # command = ['./MoT-PSA', rcfFile]
    # _runAndCheck(command)
