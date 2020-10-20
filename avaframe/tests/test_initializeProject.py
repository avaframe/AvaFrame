'''
Tests for initialize and clean project routines
'''
import os
from avaframe.in3Utils import initializeProject as initProj


def test_cleanSingleAvaDir(tmp_path):

    # Make sure cleanSingleAvaDir catches:
    # - empty variable
    undefVar = initProj.cleanSingleAvaDir('')
    assert undefVar == 'AvaDir is empty'

    # - non string variable
    notAString = initProj.cleanSingleAvaDir(2)
    assert notAString == 'AvaDir is NOT a string'

    # create test directories and touch log files
    # then delete and check again
    avaDir = tmp_path / "avaTestDir"
    initProj.initializeFolderStruct(avaDir)

    # create and fill fake logs 
    fakeLog = avaDir / "fake.log"
    keepLog = avaDir / "keep.log"
    content = 'FakelogContent'
    fakeLog.write_text(content)
    keepLog.write_text(content)

    assert fakeLog.read_text() == content
    assert keepLog.read_text() == content

    # check with keep and deleteOutput false
    initProj.cleanSingleAvaDir(str(avaDir), keep='keep', deleteOutput=False)
    # 3 items should be left (Inputs, keep.log, Outputs)
    assert len(list(avaDir.iterdir())) == 3

    # check with keep
    initProj.cleanSingleAvaDir(str(avaDir), keep='keep')
    # 2 items should be left (inputs, keep.log)
    assert len(list(avaDir.iterdir())) == 2

    # check with no optional parameter
    initProj.cleanSingleAvaDir(str(avaDir))
    # 1 items should be left (inputs)
    assert len(list(avaDir.iterdir())) == 1

def test_initializeFolderStruct(tmp_path):


    avaDir = tmp_path / "avaTestDir"

    initProj.initializeFolderStruct(avaDir)

    assert os.path.exists(avaDir / 'Inputs')
    assert os.path.exists(avaDir / 'Inputs' / 'REL')
    assert os.path.exists(avaDir / 'Outputs')
