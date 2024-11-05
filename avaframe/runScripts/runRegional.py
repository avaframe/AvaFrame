"""Run script for running simulations in parallel based on input from multiple avaFolders
and copying all results to one output folder"""

import logging
from avaframe.in3Utils.cfgUtils import getModuleConfig
from avaframe.in4Region import splitInputs

def runRegional(inFolders, splitInputs=True):

    getModuleConfig(splitInputs)
    for folder in inFolders:
        runCom1DFA
    return str(outFolder)