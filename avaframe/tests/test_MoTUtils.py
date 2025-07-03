import pytest
import numpy as np
import pathlib
from unittest.mock import patch, MagicMock, Mock, call
from avaframe.in3Utils import MoTUtils
import os
import platform


# Note 3.7.25: these tests are mostly AI-generated

def test_rewriteDEMtoZeroValues(tmp_path):
    # Create mock input data
    mockDemData = {
        "rasterData": np.array([[1.0, np.nan, 3.0],
                              [np.nan, 5.0, 6.0]]),
        "header": {
            "nodata_value": np.nan
        }
    }
    
    expectedRasterData = np.array([[1.0, 0.0, 3.0],
                                   [0.0, 5.0, 6.0]])
    
    # Create a mock Path object    avaTestDir = "avaHockeyChannelPytest"
    # mockPath = pathlib.Path(tmp_path)
    # outFile  = pathlib.Path(tmp_path, 'testDir', 'testDem.asc')
    # print(outFile)
    mockPath = MagicMock(
        spec=pathlib.Path,
        parent=MagicMock(
            spec=pathlib.Path,
            name="testDir"
        ),
        name="testDem.asc"
    )
    # mockPath = MagicMock(spec=pathlib.Path)
    # print(mockPath)
    # mockPath.stem            # = "testDem"
    
    # Mock the rasterUtils functions
    with patch('avaframe.in2Trans.rasterUtils.readRaster') as mockRead:
        with patch('avaframe.in2Trans.rasterUtils.writeResultToRaster') as mockWrite:
            # Set up the mock return value for readRaster
            mockRead.return_value = mockDemData
            
            # Call the function
            MoTUtils.rewriteDEMtoZeroValues(mockPath)
            # MoTUtils.rewriteDEMtoZeroValues(outFile)

            # Verify readRaster was called with correct argument
            mockRead.assert_called_once_with(mockPath)
            # mockRead.assert_called_once_with(outFile)

            # Verify the data was modified correctly
            np.testing.assert_array_equal(
                mockDemData["rasterData"],
                expectedRasterData
            )
            
            # Verify nodata value was updated
            assert mockDemData["header"]["nodata_value"] == 0.0
            
            # Verify writeResultToRaster was called with correct arguments
            mockWrite.assert_called_once()
            call_args = mockWrite.call_args[0]
            assert call_args[0] == mockDemData["header"]


def test_CopyMoTFiles(tmp_path):
    # Create temporary source and destination directories
    workDir = tmp_path / "work"
    outputDir = tmp_path / "output"
    workDir.mkdir()
    outputDir.mkdir()

    # Create some test files in work directory
    testFiles = [
        "test_p_max_001.txt",
        "test_p_max_002.txt",
        "another_p_max_file.txt"
    ]

    for fileName in testFiles:
        (workDir / fileName).touch()

    # Test parameters
    searchString = "p_max"
    replaceString = "ppr"

    # Create expected target file names
    expectedTargets = [
        outputDir / name.replace(searchString, replaceString)
        for name in testFiles
    ]

    # Mock shutil.copy2 to avoid actual file copying
    with patch('shutil.copy2') as mockCopy:
        # Run the function
        MoTUtils.copyMoTFiles(workDir, outputDir, searchString, replaceString)

        # Check that copy2 was called the correct number of times
        assert mockCopy.call_count == len(testFiles)

        # Verify each copy operation
        # Get actual calls and sort them by source path
        actualCopies = sorted(mockCopy.call_args_list,
                              key=lambda x: str(x[0][0]))
        sourceFilePaths = sorted(workDir.glob(f"*{searchString}*"))
        expectedTargets = sorted(expectedTargets)

        for i in range(len(testFiles)):
            callArgs = actualCopies[i][0]
            assert callArgs[0] == sourceFilePaths[i]
            assert callArgs[1] == expectedTargets[i]


def test_CopyMoTFiles_EmptyDirectory(tmp_path):
    # Test behavior with empty source directory
    workDir = tmp_path / "empty_work"
    outputDir = tmp_path / "empty_output"
    workDir.mkdir()
    outputDir.mkdir()

    with patch('shutil.copy2') as mockCopy:
        MoTUtils.copyMoTFiles(workDir, outputDir, "p_max", "ppr")
        mockCopy.assert_not_called()

def test_CopyMoTFiles_NoMatchingFiles(tmp_path):
    # Test behavior when no files match the search string
    workDir = tmp_path / "work"
    outputDir = tmp_path / "output"
    workDir.mkdir()
    outputDir.mkdir()

    # Create files that don't match the search pattern
    (workDir / "test1.txt").touch()
    (workDir / "test2.txt").touch()

    with patch('shutil.copy2') as mockCopy:
        MoTUtils.copyMoTFiles(workDir, outputDir, "p_max", "ppr")
        mockCopy.assert_not_called()

class MockProcess:
    def __init__(self, outputLines):
        self.outputLines = outputLines
        self.currentLine = 0
        self.stdout = self

    def readline(self):
        if self.currentLine >= len(self.outputLines):
            return ""
        line = self.outputLines[self.currentLine]
        self.currentLine += 1
        return line

    def poll(self):
        return None if self.currentLine < len(self.outputLines) else 0

@pytest.mark.parametrize("osName,platformSystem,expectedShell", [
    ("nt", "Windows", True),
    ("posix", "Darwin", False),
    ("posix", "Linux", False)
])
def test_RunAndCheckMoT_OSSpecific(osName, platformSystem, expectedShell):
    with patch('os.name', osName), \
            patch('platform.system', return_value=platformSystem), \
            patch('subprocess.Popen') as mockPopen:

        mockPopen.return_value = MockProcess([])
        MoTUtils.runAndCheckMoT("testCommand")

        mockPopen.assert_called_once()
        assert mockPopen.call_args[1]['shell'] == expectedShell



def test_RunAndCheckMoT_CommandTypes():
    testCases = [
        "singleCommand",
        ["command", "with", "arguments"],
        "command with spaces"
    ]

    for command in testCases:
        with patch('subprocess.Popen') as mockPopen:
            mockPopen.return_value = MockProcess([])
            MoTUtils.runAndCheckMoT(command)

            mockPopen.assert_called_once()
            assert mockPopen.call_args[0][0] == command

def test_RunAndCheckMoT_ProcessExit():
    # Test proper handling of process termination
    testOutput = ["Line 1\n", "Line 2\n"]

    with patch('subprocess.Popen') as mockPopen, \
            patch('avaframe.in3Utils.MoTUtils.log.info') as mockLog:

        mockPopen.return_value = MockProcess(testOutput)
        MoTUtils.runAndCheckMoT("testCommand")

        # Verify all output was processed
        assert mockLog.call_count == 2
        mockLog.assert_has_calls([
            call("Line 1"),
            call("Line 2")
        ])