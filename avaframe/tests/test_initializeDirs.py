"""
    Pytest for initialise directories

 """

#  Load modules
import pathlib
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
    DEMremeshedTest = avaDirtmp / 'Inputs' / 'DEMremeshed'
    fU.makeADir(DEMremeshedTest)

    # call function to be tested
    workDir, outputDir = iD.initialiseRunDirs(avaDirtmp, modName, cleanDEMremeshed=True)

    assert workDir == workDirTest
    assert outputDir == outputDirTest
    assert workDir.is_dir()
    assert outputDir.is_dir()
    assert DEMremeshedTest.is_dir() is False

    # define results
    workDirTest = avaDirtmp / 'Work' / modName
    outputDirTest = avaDirtmp / 'Outputs' / modName
    DEMremeshedTest = avaDirtmp / 'Inputs' / 'DEMremeshed'
    fU.makeADir(DEMremeshedTest)

    # call function to be tested
    workDir, outputDir = iD.initialiseRunDirs(avaDirtmp, modName, cleanDEMremeshed=False)

    assert workDir == workDirTest
    assert outputDir == outputDirTest
    assert workDir.is_dir()
    assert outputDir.is_dir()
    assert DEMremeshedTest.is_dir() is True
