'''
Tests for initialize and clean project routines
'''
import os
from avaframe.in3Utils import initializeProject as initProj
from avaframe.com2AB import com2AB as c2


def test_cleanModuleFiles(tmp_path):

    # Make sure cleanModuleFile catches:
    # - empty variable
    undefVar = initProj.cleanModuleFiles('', 'Module')
    assert undefVar == 'avaDir is empty'

    # - non string variable
    notAString = initProj.cleanModuleFiles(2, 'Module')
    assert notAString == 'avaDir is NOT a string or PurePath'

    # create test directories
    # then delete and check again
    avaDir = tmp_path / "avaTestDir"
    initProj.initializeFolderStruct(avaDir)

    # make test dirs in outputs
    for dirName in ['com2AB', 'com1DFA', 'TestDir']:
        testDir = os.path.join(avaDir, 'Outputs', dirName)
        os.makedirs(testDir)

    # make test dirs in work
    for dirName in ['com2AB', 'com1DFA', 'TestDir']:
        testDir = os.path.join(avaDir, 'Work', dirName)
        os.makedirs(testDir)

    outputDir = os.path.join(avaDir, 'Outputs')
    workDir = os.path.join(avaDir, 'Work')

    # call the function
    initProj.cleanModuleFiles(str(avaDir), c2)

    # only 2 directories should remain in Work and Output
    assert len(os.listdir(workDir)) == 2
    assert len(os.listdir(outputDir)) == 2

    # call the function this time with the alternativeName
    initProj.cleanModuleFiles(str(avaDir), c2, alternativeName='TestDir')

    # only 1 directory should remain in Work and Output
    assert len(os.listdir(workDir)) == 1
    assert len(os.listdir(outputDir)) == 1


def test_cleanSingleAvaDir(tmp_path):

    # Make sure cleanSingleAvaDir catches:
    # - empty variable
    undefVar = initProj.cleanSingleAvaDir('')
    assert undefVar == 'avaDir is empty'

    # - non string variable
    notAString = initProj.cleanSingleAvaDir(2)
    print(notAString)
    assert notAString == 'avaDir is NOT a string or PurePath'

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

    initProj.initializeFolderStruct(avaDir, removeExisting=True)
    assert os.path.exists(avaDir / 'Inputs')
    assert os.path.exists(avaDir / 'Inputs' / 'REL')
    assert os.path.exists(avaDir / 'Outputs')
