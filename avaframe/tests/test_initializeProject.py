from avaframe.in3Utils import initializeProject as initProj

def test_cleanSingleAvaDir():

    # Make sure cleanSingleAvaDir catches:
    # - empty variable
    undefVar = initProj.cleanSingleAvaDir('')
    assert undefVar == 'AvaDir is empty'

    # - non string variable
    notAString = initProj.cleanSingleAvaDir(2)
    assert notAString == 'AvaDir is NOT a string'

    # TODO check correct deletion
    # create test directories and touch log files
    # then delete and check again
