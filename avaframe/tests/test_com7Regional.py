"""Test functions for com7Regional module."""

import pathlib
import configparser
import shutil
from avaframe.com7Regional import splitInputs
from avaframe.com7Regional import com7Regional
from avaframe.in3Utils import generateTopo, getReleaseArea, cfgUtils, initializeProject


def setup_test_directory(tmp_path, avaName):
    """Set up test directory with required files."""
    avaDir = tmp_path / "RegionalAvalanches" / avaName
    initializeProject.createFolderStruct(avaDir)

    # Create topography
    cfgTopo = cfgUtils.getModuleConfig(generateTopo)
    cfgTopo["TOPO"]["dx"] = "5"
    cfgTopo["TOPO"]["xEnd"] = "500"
    cfgTopo["TOPO"]["yEnd"] = "500"
    cfgTopo["TOPO"]["meanAlpha"] = "30"
    cfgTopo["TOPO"]["demType"] = "FP"
    cfgTopo["TOPO"]["z0"] = "500"
    cfgTopo["DEMDATA"] = {
        "xl": "0",
        "yl": "0",
        "zEdit": "",
        "nodata_value": "-9999",
        "demName": "ava_topo",
    }

    generateTopo.generateTopo(cfgTopo, avaDir)

    # Create release area
    cfgRelease = cfgUtils.getModuleConfig(getReleaseArea)
    cfgRelease["FILE"]["relNo"] = "1"
    cfgRelease["FILE"]["relName"] = "release1"
    cfgRelease["GENERAL"]["dh"] = "1.0"
    cfgRelease["FILE"]["outputtxt"] = "False"
    cfgRelease["GEOMETRY"] = {"x0": "25", "y0": "25", "width": "10", "length": "10"}

    getReleaseArea.getReleaseArea(cfgTopo, cfgRelease, avaDir)

    return avaDir


def create_test_configs():
    """Create test configurations."""
    # Create main config
    cfgMainDict = {"MAIN": {"avalancheDir": "", "nCPU": "2", "flagDev": "True"}}

    # Create com7Regional config
    cfgDict = {
        "GENERAL": {
            "regionalDir": "RegionalAvalanches",
            "copyPeakFiles": "True",
            "moveInsteadOfCopy": "False",
            "mergeOutput": "True",
            "mergeTypes": "pfv",
            "mergeMethods": "max|count",
        },
        "com1DFA_com1DFA_override": {
            "defaultConfig": "False",
            "dt": "0.5",
            "simTypeList": "null",
            "resType": "ppr|pft|pfv",
            "tSteps": "1",
            "tEnd": "2",
            "meshCellSizeThreshold": "0.001",
        },
        "INPUT": {"releaseScenario": "release1"},
        "com1DFA": {"simTypeList": "null", "releaseScenario": "release1"},
        "FLAGS": {"debugPlot": "False"},
    }

    cfgMain = cfgUtils.convertDictToConfigParser(cfgMainDict)
    cfg = cfgUtils.convertDictToConfigParser(cfgDict)

    return cfgMain, cfg


def test_com7RegionalMain(tmp_path):
    """Test the main function of com7Regional module."""
    # Create avadirs with valid inputs
    avaNames = ["avaTest1", "avaTest2", "avaTest3"]
    avaDirs = []
    for avaName in avaNames:
        avaDir = setup_test_directory(tmp_path, avaName)
        avaDirs.append(avaDir)

    cfgMain, cfg = create_test_configs()
    cfgMain["MAIN"]["avalancheDir"] = str(tmp_path)

    allPeakFilesDir, mergedRastersDir = com7Regional.com7RegionalMain(cfgMain, cfg)

    # Check for expected output files and directories
    assert allPeakFilesDir is not None
    assert allPeakFilesDir.exists()
    assert mergedRastersDir is not None
    assert mergedRastersDir.exists()

    for avaDir in avaDirs:
        outputDir = avaDir / "Outputs" / "com1DFA"
        assert outputDir.exists()
        assert (outputDir / "peakFiles").exists()


def test_getTotalNumberOfSims(tmp_path):
    """Test the getTotalNumberOfSims function."""
    avaDir = setup_test_directory(tmp_path, "testAva")

    cfgMain, cfg = create_test_configs()
    cfgMain["MAIN"]["avalancheDir"] = str(tmp_path)

    totalSims = com7Regional.getTotalNumberOfSims([avaDir], cfgMain, cfg)

    assert totalSims > 0


def test_splitInputsMain(tmp_path):
    """Test splitInputsMain function using pre-generated test data."""
    # Set up test data
    test_data_dir = pathlib.Path(__file__).parent / "data" / "testIn4Region"
    inputDir = tmp_path / "avalancheDir"
    shutil.copytree(test_data_dir, inputDir)
    outputDir = tmp_path / "avalancheDir" / "com7Regional"

    # Configure test parameters
    cfg = configparser.ConfigParser()
    cfg["GENERAL"] = {"bufferSize": "2500"}
    cfgMain = configparser.ConfigParser()
    cfgMain["FLAGS"] = {"createReport": "True", "savePlot": "True"}

    # Run function
    splitInputs.splitInputsMain(inputDir, outputDir, cfg, cfgMain)

    # Verify outputs
    assert outputDir.exists()

    # Check group directories
    groupDirs = list(outputDir.glob("group*"))
    assert len(groupDirs) == 2

    for groupDir in groupDirs:
        # Check directory structure
        assert (groupDir / "Inputs").exists()
        assert (groupDir / "Inputs" / "REL").exists()
        assert len(list((groupDir / "Inputs").glob("*.tif"))) == 1

        # Check release areas were split by scenarios
        relDir = groupDir / "Inputs" / "REL"
        assert len(list(relDir.glob("*.shp"))) == 2

        # Check optional inputs
        assert (groupDir / "Inputs" / "ENT").exists()
        assert len(list((groupDir / "Inputs" / "ENT").glob("*.shp"))) == 1
        assert (groupDir / "Inputs" / "RES").exists()
        assert len(list((groupDir / "Inputs" / "RES").glob("*.shp"))) == 1

    # Check reports
    assert (outputDir / "splitInputs_scenarioReport.txt").exists()
    assert (outputDir / "splitInputs_visualReport_basic.png").exists()
    assert (outputDir / "splitInputs_visualReport_optional.png").exists()
