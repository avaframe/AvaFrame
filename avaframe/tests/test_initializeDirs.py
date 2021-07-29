"""
    Pytest for initialise directories

 """

#  Load modules
import pathlib
from avaframe.in3Utils import initialiseDirs as iD


def test_initialiseRunDirs(tmp_path):
    """ test initialising run directories """

    # set directory
    avaName = 'avaTest'
    avaDirtmp = pathlib.Path(tmp_path, avaName)
    modName = 'com1DFA'

    # define results
    workDirTest = avaDirtmp / 'Work' / modName
    outputDirTest = avaDirtmp / 'Outputs' / modName

    # call function to be tested
    workDir, outputDir = iD.initialiseRunDirs(avaDirtmp, modName)

    assert workDir == workDirTest
    assert outputDir == outputDirTest
    assert workDir.is_dir()
    assert outputDir.is_dir()
