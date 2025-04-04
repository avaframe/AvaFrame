"""
    Pytest for initialise directories

 """

#  Load modules
import pathlib
import pytest
from avaframe.in3Utils import initialiseDirs as iD
from avaframe.in3Utils import fileHandlerUtils as fU


def test_initialiseRunDirs(tmp_path):
    """ test initialising run directories """

    # set directory
    avaName = 'avaTest'
    avaDirtmp = pathlib.Path(tmp_path, avaName)
    modName = 'com1DFA'

    # define results
    workDirTest = avaDirtmp / 'Work' / modName
    outputDirTest = avaDirtmp / 'Outputs' / modName
    remeshedRastersTest = avaDirtmp / 'Inputs' / 'remeshedRasters'
    fU.makeADir(remeshedRastersTest)

    # call function to be tested
    workDir, outputDir = iD.initialiseRunDirs(avaDirtmp, modName, cleanRemeshedRasters=True)

    assert workDir == workDirTest
    assert outputDir == outputDirTest
    assert workDir.is_dir()
    assert outputDir.is_dir()
    assert remeshedRastersTest.is_dir() is False

    # call function to be tested
    with pytest.raises(AssertionError) as e:
        assert iD.initialiseRunDirs(avaDirtmp, modName, cleanRemeshedRasters=False)
    assert str(e.value) == 'Work directory %s already exists - delete first!' % (workDirTest)

    # set directory
    avaName = 'avaTest2'
    avaDirtmp = pathlib.Path(tmp_path, avaName)
    modName = 'com1DFA'

    # define results
    workDirTest = avaDirtmp / 'Work' / modName
    outputDirTest = avaDirtmp / 'Outputs' / modName
    remeshedRastersTest = avaDirtmp / 'Inputs' / 'remeshedRasters'
    fU.makeADir(remeshedRastersTest)

    # call function to be tested
    workDir, outputDir = iD.initialiseRunDirs(avaDirtmp, modName, cleanRemeshedRasters=False)

    assert workDir == workDirTest
    assert outputDir == outputDirTest
    assert workDir.is_dir()
    assert outputDir.is_dir()
    assert remeshedRastersTest.is_dir() is True
    assert (avaDirtmp / 'Outputs' / modName / 'configurationFiles' / 'configurationFilesDone').is_dir()
    assert (avaDirtmp / 'Outputs' / modName / 'configurationFiles' / 'configurationFilesLatest').is_dir()


