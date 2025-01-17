"""

    This file is part of Avaframe.

 """

#  Load modules
import os
import pathlib
from avaframe.in1Data import getInput
import configparser
import pytest
import shutil
import numpy as np
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import generateTopo
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.com1DFA.DFAtools as DFAtls
from avaframe.com1DFA import com1DFA
import logging


def test_readDEM():
    """test reading DEM"""

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = "avaHockeyChannel"
    avaDir = dirPath / ".." / "data" / avaName

    # call function to be tested
    dem = getInput.readDEM(avaDir)

    #    print('demHeader', dem['header'])

    assert dem["header"]["ncols"] == 1001
    assert dem["header"]["nrows"] == 401
    assert dem["header"]["xllcenter"] == 1000
    assert dem["header"]["yllcenter"] == -5000
    assert dem["header"]["cellsize"] == 5
    assert dem["rasterData"].shape == (401, 1001)


def test_getDEMFromConfig(tmp_path):
    """test get path of DEM"""

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = "avaHockeyChannel"
    avaDir = dirPath / ".." / "data" / avaName
    avaDirInputs = avaDir / "Inputs"

    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / "Inputs"
    shutil.copytree(avaDirInputs, avaTestDirInputs)

    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getDEMFromConfig(avaTestDir, fileName="Hutzly")
    assert str(e.value) == "Dem file: %s/Hutzly does not exist" % avaTestDirInputs

    demFile = getInput.getDEMFromConfig(avaTestDir, fileName="DEM_HS_Topo.asc")
    assert "DEM_HS_Topo.asc" in str(demFile)


def test_getDEMPath(tmp_path):
    """test get path of DEM"""

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = "avaHockeyChannel"
    avaDir = dirPath / ".." / "data" / avaName
    avaDirInputs = avaDir / "Inputs"
    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / "Inputs"
    shutil.copytree(avaDirInputs, avaTestDirInputs)

    # call function to be tested
    demPath = getInput.getDEMPath(avaTestDir)

    #    print('dem path', demPath)

    assert "DEM_HS_Topo.asc" in str(demPath)

    # call function to be tested
    inputFile = avaDirInputs / "DEM_HS_Topo.asc"
    testFile = avaTestDirInputs / "DEM_HS_Topo2.asc"
    shutil.copyfile(inputFile, testFile)

    with pytest.raises(AssertionError) as e:
        assert getInput.getDEMPath(avaTestDir)

    assert str(e.value) == "There should be exactly one topography .asc/.tif file in %s/Inputs/" % (
        avaTestDir
    )

    # call function to be tested
    avaDirTest3 = pathlib.Path(tmp_path, "avaTest2")
    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getDEMPath(avaDirTest3)
    assert str(e.value) == "No topography .asc / .tif file in %s/Inputs/" % (avaDirTest3)


def test_getInputData(tmp_path):
    """test check for input data"""

    # get input data
    dirPath = os.path.dirname(__file__)
    avaName = "avaHockeyChannel"
    avaDir = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, "Inputs")
    avaData = os.path.join(dirPath, "..", "data", avaName, "Inputs")
    shutil.copytree(avaData, avaInputs)

    # Initialise input in correct format
    cfg = configparser.ConfigParser()
    cfg["INPUT"] = {"releaseScenario": ""}

    # call function to be tested
    dem, rels, ent, res, wall, entResInfo = getInput.getInputData(avaDir, cfg["INPUT"])
    # second option
    cfg["INPUT"]["releaseScenario"] = "release1HS"
    dem2, rels2, ent2, res2, wall, entResInfo2 = getInput.getInputData(avaDir, cfg["INPUT"])
    # Test
    assert str(dem) == str(pathlib.Path(avaDir, "Inputs", "DEM_HS_Topo.asc"))
    assert rels == [
        pathlib.Path(avaDir, "Inputs", "REL", "release1HS.shp"),
        pathlib.Path(avaDir, "Inputs", "REL", "release2HS.shp"),
        pathlib.Path(avaDir, "Inputs", "REL", "release3HS.shp"),
    ]
    assert rels2 == [pathlib.Path(avaDir, "Inputs", "REL", "release1HS.shp")]
    assert res == ""
    assert str(ent) == str(pathlib.Path(avaDir, "Inputs", "ENT", "entrainment1HS.shp"))
    assert entResInfo["flagEnt"] == "Yes"
    assert entResInfo["flagRes"] == "No"
    assert entResInfo["flagWall"] == "No"
    assert wall is None

    # third option
    cfg["INPUT"]["releaseScenario"] = "release1HS.shp"
    dem3, rels3, ent3, res3, wall, entResInfo3 = getInput.getInputData(avaDir, cfg["INPUT"])
    assert str(ent3) == str(pathlib.Path(avaDir, "Inputs", "ENT", "entrainment1HS.shp"))
    assert entResInfo3["flagEnt"] == "Yes"
    assert entResInfo3["flagRes"] == "No"
    assert entResInfo["flagWall"] == "No"
    assert wall is None
    assert rels3 == [pathlib.Path(avaDir, "Inputs", "REL", "release1HS.shp")]

    # call function to be tested
    cfg["INPUT"]["releaseScenario"] = "release4HS"
    releaseF = pathlib.Path(avaDir, "Inputs", "REL", "release4HS.shp")
    with pytest.raises(FileNotFoundError) as e:
        assert getInput.getInputData(avaDir, cfg["INPUT"])
    assert str(e.value) == ("No release scenario called: %s" % releaseF)

    # fifth option
    cfg["INPUT"]["releaseScenario"] = "release1BL.shp"
    avaName = "avaBowl"
    avaDir = os.path.join(tmp_path, avaName)
    avaInputs = os.path.join(avaDir, "Inputs")
    avaData = os.path.join(dirPath, "..", "data", avaName, "Inputs")
    shutil.copytree(avaData, avaInputs)
    dem6, rels6, ent6, res6, wall, entResInfo6 = getInput.getInputData(avaDir, cfg["INPUT"])
    assert ent6 == ""
    assert res6 == ""
    assert str(wall) == os.path.join(avaDir, "Inputs", "DAM", "dam.shp")
    assert entResInfo6["flagWall"] == "Yes"
    assert entResInfo6["flagEnt"] == "No"
    assert entResInfo6["flagRes"] == "No"


def test_getInputDataCom1DFA(tmp_path):
    """test check for input data"""

    # get input data
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = "avaHockeyChannel"
    avaDir = pathlib.Path(tmp_path, avaName)
    avaInputs = avaDir / "Inputs"
    avaData = dirPath / ".." / "data" / avaName / "Inputs"
    shutil.copytree(avaData, avaInputs)

    # call function to be tested
    inputSimFiles = getInput.getInputDataCom1DFA(avaDir)

    # Test
    #    print(inputSimFiles['demFile'])
    #    print(avaDir / 'Inputs' / 'DEM_HS_Topo.asc')
    #    print(inputSimFiles['relFiles'])
    #    print([avaDir / 'Inputs' / 'REL' / 'release1HS.shp', avaDir / 'Inputs' / 'REL' / 'release2HS.shp', avaDir / 'Inputs' / 'REL' / 'release3HS.shp'])
    assert inputSimFiles["demFile"] == avaDir / "Inputs" / "DEM_HS_Topo.asc"
    assert inputSimFiles["relFiles"] == [
        avaDir / "Inputs" / "REL" / "release1HS.shp",
        avaDir / "Inputs" / "REL" / "release2HS.shp",
        avaDir / "Inputs" / "REL" / "release3HS.shp",
    ]
    assert inputSimFiles["resFile"] == None
    assert inputSimFiles["entFile"] == avaDir / "Inputs" / "ENT" / "entrainment1HS.shp"
    assert inputSimFiles["entResInfo"]["flagEnt"] == "Yes"
    assert inputSimFiles["entResInfo"]["flagRes"] == "No"
    assert inputSimFiles["relThFile"] == None


def test_getAndCheckInputFiles(tmp_path):
    """test fetching input files and checking if exist"""

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = "avaHockeyChannel"
    avaDir = dirPath / ".." / "data" / avaName
    avaDirInputs = avaDir / "Inputs"
    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / "Inputs"
    shutil.copytree(avaDirInputs, avaTestDirInputs)
    folder = "ENT"
    inputType = "entrainment"

    # call function to be tested
    outFile, available = getInput.getAndCheckInputFiles(avaTestDirInputs, folder, inputType, fileExt="shp")

    #    print('outfile', outFile)
    #    print('available', available)

    assert available == "Yes"
    assert "Inputs/ENT/entrainment1HS.shp" in str(outFile)

    # call function to be tested
    inputFile = avaDirInputs / "ENT" / "entrainment1HS.shp"
    testFile = avaTestDirInputs / "ENT" / "entrainment1HS2.shp"
    shutil.copyfile(inputFile, testFile)
    with pytest.raises(AssertionError) as e:
        assert getInput.getAndCheckInputFiles(avaTestDirInputs, folder, inputType, fileExt="shp")
    assert str(e.value) == (
            "More than one %s .shp file in %s/%s/ not allowed" % (inputType, avaTestDirInputs, folder)
    )


def test_getThicknessInputSimFiles(tmp_path):
    """test fetching thickness from shapefiles attributes"""

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = "avaHockeyChannel"
    avaDir = dirPath / ".." / "data" / avaName
    avaDirInputs = avaDir / "Inputs"
    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / "Inputs"
    shutil.copytree(avaDirInputs, avaTestDirInputs)
    outDir = avaTestDir / "Outputs" / "com1DFA"
    fU.makeADir(outDir)

    demFile = avaTestDirInputs / "DEM_HS_Topo.asc"
    relFile1 = avaTestDirInputs / "REL" / "release1HS.shp"
    relFile2 = avaTestDirInputs / "REL" / "release2HS.shp"
    entFile = avaTestDirInputs / "ENT" / "entrainment1HS.shp"
    inputSimFiles = {
        "demFile": demFile,
        "relFiles": [relFile1, relFile2],
        "entFile": entFile,
        "secondaryReleaseFile": None,
        "entResInfo": {"flagRes": "No", "flagEnt": "Yes", "flagSecondaryRelease": "No"},
        "relThFile": None,
    }

    inputSimFiles = getInput.getThicknessInputSimFiles(inputSimFiles)

    #    print('inputSimFiles', inputSimFiles)

    assert inputSimFiles["release1HS"]["thickness"] == ["1.0"]
    assert inputSimFiles["release2HS"]["thickness"] == ["1.0", "1.0"]
    assert inputSimFiles["release1HS"]["id"] == ["0"]
    assert inputSimFiles["release2HS"]["id"] == ["0", "1"]
    assert inputSimFiles["release1HS"]["ci95"] == ["None"]
    assert inputSimFiles["release2HS"]["ci95"] == ["None", "None"]
    assert inputSimFiles["entrainment1HS"]["thickness"] == ["0.3"]
    assert inputSimFiles["entrainment1HS"]["id"] == ["0"]


def test_updateThicknessCfg(tmp_path):
    """test fetching thickness from shapefiles attributes"""

    # setup required input
    dirPath = pathlib.Path(__file__).parents[0]
    avaName = "avaHockeyChannel"
    avaDir = dirPath / ".." / "data" / avaName
    avaDirInputs = avaDir / "Inputs"
    avaTestDir = pathlib.Path(tmp_path, avaName)
    avaTestDirInputs = avaTestDir / "Inputs"
    shutil.copytree(avaDirInputs, avaTestDirInputs)
    outDir = avaTestDir / "Outputs" / "com1DFA"
    fU.makeADir(outDir)

    cfg = configparser.ConfigParser()
    cfg = cfgUtils.getModuleConfig(com1DFA, toPrint=False, onlyDefault=True)

    cfg["GENERAL"]["relThFromFile"] = "False"
    cfg["GENERAL"]["simTypeList"] = "null|ent"
    cfg["GENERAL"]["secRelAra"] = "False"
    cfg["INPUT"] = {"releaseScenario": ""}

    demFile = avaTestDirInputs / "DEM_HS_Topo.asc"
    relFile1 = avaTestDirInputs / "REL" / "release1HS.shp"
    relFile2 = avaTestDirInputs / "REL" / "release2HS.shp"
    entFile = avaTestDirInputs / "ENT" / "entrainment1HS.shp"
    inputSimFiles = {
        "demFile": demFile,
        "relFiles": [relFile1, relFile2],
        "entFile": entFile,
        "secondaryReleaseFile": None,
        "entResInfo": {"flagRes": "No", "flagEnt": "Yes", "flagSecondaryRelease": "No"},
        "relThFile": None,
        "releaseScenarioList": ["release1HS", "release2HS"],
    }

    inputSimFiles["release1HS"] = {"thickness": ["1.0"], "id": ["0"], "ci95": ["None", "None"]}
    inputSimFiles["release2HS"] = {"thickness": ["1.0", "1.0"], "id": ["0", "1"], "ci95": ["None", "None"]}
    inputSimFiles["entrainment1HS"] = {"thickness": ["0.3"], "id": ["0"], "ci95": ["None"]}

    cfg = getInput.updateThicknessCfg(inputSimFiles, cfg)

    #    print('inputSimFiles', inputSimFiles)

    assert cfg["INPUT"]["releaseScenario"] == "release1HS|release2HS"
    assert cfg["INPUT"]["release1HS_relThId"] == "0"
    assert cfg["INPUT"]["release2HS_relThId"] == "0|1"
    assert cfg["INPUT"]["release1HS_relThThickness"] == "1.0"
    assert cfg["INPUT"]["release2HS_relThThickness"] == "1.0|1.0"
    assert cfg["INPUT"]["release1HS_relThCi95"] == "None"
    assert cfg["INPUT"]["release2HS_relThCi95"] == "None|None"

    assert cfg["GENERAL"]["relTh"] == ""
    assert cfg["GENERAL"].getboolean("relThFromShp") == True
    assert cfg["GENERAL"].getboolean("relThFromFile") == False
    assert cfg["GENERAL"]["entTh"] == ""
    assert cfg["GENERAL"].getboolean("entThFromShp") == True
    assert cfg["INPUT"]["entrainmentScenario"] == "entrainment1HS"
    assert cfg["INPUT"]["entThId"] == "0"
    assert cfg["INPUT"]["entThThickness"] == "0.3"


def test_selectReleaseFile(tmp_path):
    """testing selecting a release area scenario according to configuration settings"""

    # setup the required inputs
    testPath = pathlib.Path(tmp_path, "avaTest", "Inputs", "REL")
    rel1 = testPath / "rel1.shp"
    rel2 = testPath / "rel2.shp"

    inputSimFiles = {"relFiles": [rel1, rel2]}
    cfg = configparser.ConfigParser()
    cfg["INPUT"] = {"releaseScenario": "rel1"}

    # call function to be tested
    inputSimFiles = getInput.selectReleaseFile(inputSimFiles, cfg["INPUT"]["releaseScenario"])

    assert inputSimFiles["relFiles"][0].name == "rel1.shp"
    assert inputSimFiles["relFiles"][1].name == "rel2.shp"
    assert len(inputSimFiles["relFiles"]) == 2
    assert inputSimFiles["releaseScenario"] == rel1

    cfg = configparser.ConfigParser()
    cfg["INPUT"] = {"releaseScenario": "rel2"}
    inputSimFiles = {"relFiles": [rel1, rel2]}

    # call function to be tested
    inputSimFiles = getInput.selectReleaseFile(inputSimFiles, cfg["INPUT"]["releaseScenario"])

    assert inputSimFiles["relFiles"][0].name == "rel1.shp"
    assert inputSimFiles["relFiles"][1].name == "rel2.shp"
    assert len(inputSimFiles["relFiles"]) == 2
    assert inputSimFiles["releaseScenario"] == rel2


def test_fetchReleaseFile(tmp_path):
    """testing selecting a release area scenario according to configuration settings"""

    # setup the required inputs
    testPath = pathlib.Path(tmp_path, "avaTest", "Inputs", "REL")
    rel1 = testPath / "rel1.shp"
    rel2 = testPath / "rel2.shp"

    inputSimFiles = {"relFiles": [rel1, rel2]}
    cfg = configparser.ConfigParser()
    cfg["INPUT"] = {"releaseScenario": "rel1"}
    cfg["GENERAL"] = {"relThFromShp": False}
    releaseScenario = "rel1"
    releaseList = ["rel1", "rel2"]

    # call function to be tested
    releaseScenarioPath, cfg = getInput.fetchReleaseFile(inputSimFiles, releaseScenario, cfg, releaseList)

    assert releaseScenarioPath == rel1
    assert cfg["INPUT"]["releaseScenario"] == "rel1"

    cfg = configparser.ConfigParser()
    cfg["INPUT"] = {"releaseScenario": "rel2"}
    inputSimFiles = {"relFiles": [rel1, rel2]}
    cfg["GENERAL"] = {"relThFromShp": True}
    cfg["INPUT"] = {
        "rel2_relThId": "0",
        "rel2_relThThickness": "2.",
        "rel2_relThCi95": "",
        "rel1_relThId": "1",
        "rel1_relThThickness": "1.",
        "rel1_relThCi95": "",
    }
    # call function to be tested
    releaseScenarioPath, cfg = getInput.fetchReleaseFile(inputSimFiles, "rel2", cfg, releaseList)

    assert releaseScenarioPath == rel2
    assert cfg["INPUT"]["relThId"] == "0"
    assert cfg["INPUT"]["relThThickness"] == "2."
    assert cfg["INPUT"]["relThCi95"] == ""


def test_createReleaseStats(tmp_path):
    """test creating a release shp file info"""

    testPath = pathlib.Path(tmp_path, "avaTestGeo")
    testPathInputs = pathlib.Path(tmp_path, "avaTestGeo", "Inputs", "REL")
    fU.makeADir(testPathInputs)

    cfgGenTop = cfgUtils.getModuleConfig(generateTopo)
    cfgGenTop["TOPO"]["dx"] = "1."
    cfgGenTop["TOPO"]["demType"] = "IP"
    cfgGenTop["TOPO"]["meanAlpha"] = "27.5"
    cfgGenTop["TOPO"]["z0"] = "2000"
    cfgGenTop["TOPO"]["xEnd"] = "5000"
    cfgGenTop["TOPO"]["yEnd"] = "1000"
    cfgGenTop["TOPO"]["channel"] = "False"
    cfgGenTop["DEMDATA"]["xl"] = "0."
    cfgGenTop["DEMDATA"]["yl"] = "0."

    [z, name_ext, outDir] = generateTopo.generateTopo(cfgGenTop, testPath)

    print(testPath)
    # setup release line
    lineDict = {
        "x": np.asarray([100.0, 100.0, 150.0, 200.0, 200.0, 150.0, 100.0]),
        "y": np.asarray([100.0, 150.0, 150.0, 150.0, 100.0, 100.0, 100]),
    }
    fileName = pathlib.Path(testPath, "Inputs", "REL", "releaseIP.shp")
    lineName = "release1"
    shpConv.writeLine2SHPfile(lineDict, lineName, fileName, header="")

    # Load configuration file for probabilistic run and analysis
    cfg = cfgUtils.getModuleConfig(com1DFA)
    relDFDict = getInput.createReleaseStats(testPath, cfg)

    # compute parameter
    zMax = 2000.0 - np.tan(np.deg2rad(27.5)) * 100.0
    zMin = 2000.0 - np.tan(np.deg2rad(27.5)) * 200.0

    assert relDFDict["releaseIP"]["release feature"].iloc[0] == "release1"
    assert np.isclose(relDFDict["releaseIP"]["slope [deg]"].iloc[0], 27.5)
    print(20 * "-")
    print(relDFDict["releaseIP"]["MaxZ [m]"].iloc[0], zMax)
    print(20 * "-")
    assert np.isclose(relDFDict["releaseIP"]["MaxZ [m]"].iloc[0], zMax)

    print(20 * "-")
    print(relDFDict["releaseIP"]["MinZ [m]"].iloc[0], zMin)
    print(20 * "-")
    assert np.isclose(relDFDict["releaseIP"]["MinZ [m]"].iloc[0], zMin)
    assert np.isclose(relDFDict["releaseIP"]["projected area [ha]"].iloc[0], 0.5151)
    assert np.isclose(relDFDict["releaseIP"]["actual area [ha]"].iloc[0], 0.58071)


def test_computeAreasFromLines():
    """test computing areas with shapely from a lineDict"""

    # setup required inputs
    # setup release line
    lineDict = {
        "x": np.asarray([100.0, 100.0, 150.0, 200.0, 200.0, 150.0, 100.0]),
        "y": np.asarray([100.0, 150.0, 150.0, 150.0, 100.0, 100.0, 100]),
        "Start": np.asarray([0.0]),
        "Length": np.asarray([7]),
        "Name": [""],
    }

    # call function to be tested
    projectedAreas = getInput.computeAreasFromLines(lineDict)

    assert projectedAreas[0] == 5000.0
    assert len(projectedAreas) == 1


def test_computeAreasFromRasterAndLine(tmp_path):
    """test computation of areas using a lineDict and a dem raster"""

    # setup required inputs
    testPath = pathlib.Path(tmp_path, "avaTestGeo2")
    testPathInputs = pathlib.Path(tmp_path, "avaTestGeo2", "Inputs", "REL")
    fU.makeADir(testPathInputs)

    cfgGenTop = cfgUtils.getModuleConfig(generateTopo)
    cfgGenTop["TOPO"]["dx"] = "1."
    cfgGenTop["TOPO"]["demType"] = "IP"
    cfgGenTop["TOPO"]["meanAlpha"] = "27.5"
    cfgGenTop["TOPO"]["z0"] = "2000"
    cfgGenTop["TOPO"]["xEnd"] = "5000"
    cfgGenTop["TOPO"]["yEnd"] = "1000"
    cfgGenTop["TOPO"]["channel"] = "False"
    cfgGenTop["DEMDATA"]["xl"] = "0."
    cfgGenTop["DEMDATA"]["yl"] = "0."

    [z, name_ext, outDir] = generateTopo.generateTopo(cfgGenTop, testPath)
    dem = getInput.readDEM(testPath)

    dem["originalHeader"] = dem["header"].copy()
    methodMeshNormal = 1
    # get normal vector of the grid mesh
    dem = geoTrans.getNormalMesh(dem, num=methodMeshNormal)
    dem = DFAtls.getAreaMesh(dem, methodMeshNormal)

    lineDict = {
        "x": np.asarray([100.0, 100.0, 150.0, 200.0, 200.0, 150.0, 100.0]),
        "y": np.asarray([100.0, 150.0, 150.0, 150.0, 100.0, 100.0, 100]),
        "Start": np.asarray([0.0]),
        "Length": np.asarray([7]),
        "Name": [""],
    }

    # call function to be tested
    areaActualList, areaProjectedList, line = getInput.computeAreasFromRasterAndLine(lineDict, dem)

    assert np.isclose(areaActualList[0], 5807.14)
    assert areaProjectedList[0] == 5151.00


def test_computeRelStats(tmp_path):
    """test computing min, max eleavtions and other stats of a line and dem"""

    # setup required inputs
    testPath = pathlib.Path(tmp_path, "avaTestGeo3")
    testPathInputs = pathlib.Path(tmp_path, "avaTestGeo3", "Inputs", "REL")
    fU.makeADir(testPathInputs)

    cfgGenTop = cfgUtils.getModuleConfig(generateTopo)
    cfgGenTop["TOPO"]["dx"] = "1."
    cfgGenTop["TOPO"]["demType"] = "IP"
    cfgGenTop["TOPO"]["meanAlpha"] = "27.5"
    cfgGenTop["TOPO"]["z0"] = "2000"
    cfgGenTop["TOPO"]["xEnd"] = "5000"
    cfgGenTop["TOPO"]["yEnd"] = "1000"
    cfgGenTop["TOPO"]["channel"] = "False"
    cfgGenTop["DEMDATA"]["xl"] = "0."
    cfgGenTop["DEMDATA"]["yl"] = "0."

    [z, name_ext, outDir] = generateTopo.generateTopo(cfgGenTop, testPath)
    dem = getInput.readDEM(testPath)

    dem["originalHeader"] = dem["header"].copy()
    methodMeshNormal = 1
    # get normal vector of the grid mesh
    dem = geoTrans.getNormalMesh(dem, num=methodMeshNormal)
    dem = DFAtls.getAreaMesh(dem, methodMeshNormal)

    lineDict = {
        "x": np.asarray([100.0, 100.0, 150.0, 200.0, 200.0, 150.0, 100.0]),
        "y": np.asarray([100.0, 150.0, 150.0, 150.0, 100.0, 100.0, 100]),
        "Start": np.asarray([0.0]),
        "Length": np.asarray([7]),
        "Name": [""],
    }
    lineDict = geoTrans.prepareArea(lineDict, dem, 0.01, combine=False, checkOverlap=False)

    # call function to be tested
    lineDict = getInput.computeRelStats(lineDict, dem)

    # compute parameter
    zMax = 2000.0 - np.tan(np.deg2rad(27.5)) * 100.0
    zMin = 2000.0 - np.tan(np.deg2rad(27.5)) * 200.0

    assert np.isclose(lineDict["zMax"][0], zMax)
    assert np.isclose(lineDict["zMin"][0], zMin)
    assert np.isclose(lineDict["meanSlope"][0], 27.5)
    assert lineDict["featureNames"][0] == ""


# Based on AI suggestions:

def test_getAndCheckInputFiles_noFilesFound(mocker):
    mockInDir = mocker.MagicMock()
    mockInDir.glob.return_value = []
    mockPath = mocker.patch('pathlib.Path')
    mockPath.return_value = mockInDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    output_file, available = getInput.getAndCheckInputFiles(inputDir, folder, 'inputType', fileExt='shp')

    mockPath.assert_called_once_with(inputDir, folder)
    mockInDir.glob.assert_called_once_with('*.shp')
    assert output_file is None
    assert available == 'No'


def test_getAndCheckInputFilesone_valid_file(mocker):
    mock_file = mocker.MagicMock()
    mock_file.suffix = '.shp'
    mock_inDir = mocker.MagicMock()
    mock_inDir.glob.return_value = [mock_file]
    mock_path = mocker.patch('pathlib.Path')
    mock_path.return_value = mock_inDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    output_file, available = getInput.getAndCheckInputFiles(inputDir, folder, 'inputType', fileExt='shp')

    mock_inDir.glob.assert_called_once_with('*.shp')
    assert output_file == mock_file
    assert available == 'Yes'


def test_getAndCheckInputFilesmultiple_files_error(mocker):
    mock_file1 = mocker.MagicMock()
    mock_file1.suffix = '.shp'
    mock_file2 = mocker.MagicMock()
    mock_file2.suffix = '.shp'
    mock_inDir = mocker.MagicMock()
    mock_inDir.glob.return_value = [mock_file1, mock_file2]
    mock_path = mocker.patch('pathlib.Path')
    mock_path.return_value = mock_inDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    with pytest.raises(AssertionError) as excinfo:
        getInput.getAndCheckInputFiles(inputDir, folder, 'inputType', fileExt='shp')

    expected_msg = "More than one inputType .shp file in /fake/dir/fake_folder/ not allowed"
    assert expected_msg in str(excinfo.value)


def test_getAndCheckInputFilesunsupported_extension_error(mocker):
    mock_file = mocker.MagicMock()
    mock_file.suffix = '.txt'
    mock_inDir = mocker.MagicMock()
    mock_inDir.glob.return_value = [mock_file]
    mock_path = mocker.patch('pathlib.Path')
    mock_path.return_value = mock_inDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    with pytest.raises(AssertionError) as excinfo:
        getInput.getAndCheckInputFiles(inputDir, folder, 'inputType', fileExt='txt')

    assert "Unsupported file format found for OutputFile" in str(excinfo.value)


def test_getAndCheckInputFilesraster_extensions(mocker):
    mock_asc = mocker.MagicMock()
    mock_asc.suffix = '.asc'
    mock_tif = mocker.MagicMock()
    mock_tif.suffix = '.tif'
    mock_inDir = mocker.MagicMock()
    mock_inDir.glob.side_effect = lambda p: [mock_asc] if p == '*.asc' else [mock_tif]
    mock_path = mocker.patch('pathlib.Path')
    mock_path.return_value = mock_inDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    with pytest.raises(AssertionError) as excinfo:
        getInput.getAndCheckInputFiles(inputDir, folder, 'inputType', fileExt='raster')

    assert "More than one inputType .raster file" in str(excinfo.value)
    mock_inDir.glob.assert_any_call('*.asc')
    mock_inDir.glob.assert_any_call('*.tif')


def test_getAndCheckInputFilesfile_suffix_filter(mocker):
    mock_file = mocker.MagicMock()
    mock_file.suffix = '.shp'
    mock_inDir = mocker.MagicMock()
    mock_inDir.glob.return_value = [mock_file]
    mock_path = mocker.patch('pathlib.Path')
    mock_path.return_value = mock_inDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    fileSuffix = '_suffix'
    output_file, available = getInput.getAndCheckInputFiles(
        inputDir, folder, 'inputType', fileExt='shp', fileSuffix=fileSuffix
    )

    mock_inDir.glob.assert_called_once_with('*_suffix.shp')
    assert output_file == mock_file
    assert available == 'Yes'


def test_getAndCheckInputFilesempty_file_ext_with_suffix(mocker):
    mock_file = mocker.MagicMock()
    mock_file.suffix = ''
    mock_inDir = mocker.MagicMock()
    mock_inDir.glob.return_value = [mock_file]
    mock_path = mocker.patch('pathlib.Path')
    mock_path.return_value = mock_inDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    fileSuffix = '_suffix'
    with pytest.raises(AssertionError) as excinfo:
        getInput.getAndCheckInputFiles(inputDir, folder, 'inputType', fileExt='', fileSuffix=fileSuffix)

    mock_inDir.glob.assert_called_once_with('*_suffix.')
    assert "Unsupported file format" in str(excinfo.value)


def test_getAndCheckInputFilesempty_file_ext_and_suffix(mocker):
    mock_file = mocker.MagicMock()
    mock_file.suffix = ''
    mock_inDir = mocker.MagicMock()
    mock_inDir.glob.return_value = [mock_file]
    mock_path = mocker.patch('pathlib.Path')
    mock_path.return_value = mock_inDir

    inputDir = '/fake/dir'
    folder = 'fake_folder'
    with pytest.raises(AssertionError) as excinfo:
        getInput.getAndCheckInputFiles(inputDir, folder, 'inputType', fileExt='', fileSuffix='')

    mock_inDir.glob.assert_called_once_with('*.')
    assert "Unsupported file format" in str(excinfo.value)
