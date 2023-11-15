"""
    Pytest for module com1DFA
"""

#  Load modules
import numpy as np
import logging
import pytest
import configparser
import pathlib
import copy
import pickle
import shutil

from avaframe.com6RockAvalanche import com6RockAvalanche
# import avaframe.in2Trans.ascUtils as IOf
# import avaframe.in3Utils.fileHandlerUtils as fU
# import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils

def test_runCom6RockAvalanche(tmp_path, caplog):
    """Check that runCom1DFA produces the good outputs"""
    testDir = pathlib.Path(__file__).parents[0]
    inputDir = testDir / "data" / "testCom6RockAvalanche"
    avaDir = pathlib.Path(tmp_path, "testCom6RockAvalanche")
    shutil.copytree(inputDir, avaDir)
    cfgFile = avaDir / "test_com1DFACfg.ini"
    cfgMain = configparser.ConfigParser()
    cfgMain["MAIN"] = {"avalancheDir": str(avaDir), "nCPU": "auto", "CPUPercent": "90"}
    cfgMain["FLAGS"] = {
        "showPlot": "False",
        "savePlot": "True",
        "ReportDir": "True",
        "reportOneFile": "True",
        "debugPlot": "False",
    }
    # modCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=cfgFile, modInfo=True)
    modCfg, modInfo = cfgUtils.getModuleConfig(com6RockAvalanche, modInfo=True)
    print(modCfg)

    dem, plotDict, reportDictList, simDF = com6RockAvalanche.runRockAvalanche(cfgMain, modCfg)


    outDir = avaDir / "Outputs" / "com1DFA"
    for ext in ["FT", "pft", "pfv"]:
        assert (outDir / "peakFiles" / ("%s_%s.asc" % (simDF["simName"][0], ext))).is_file()

    assert (outDir / "configurationFiles" / ("%s.ini" % (simDF["simName"][0]))).is_file()
    assert (outDir / "configurationFiles" / ("allConfigurations.csv")).is_file()


