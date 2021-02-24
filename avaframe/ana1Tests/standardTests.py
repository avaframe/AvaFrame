"""
    test
"""

# Load modules
import os
import glob
import logging
import json

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from benchmarks import simParameters


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def createDesDictTemplate():
    """ create a json file with all required test info """

    desDict = {'TAGS': '',
                'DESCRIPTION': '',
                'TYPE': '',
                'FILES': ''}

    return desDict


def writeDesDicttoJson(desDict, testName, outDir):
    """ create a json file with all required test info from desdict """

    jsonDict = json.dumps(desDict)
    fileName = os.path.join(outDir,"%s_desDict.json" % testName)
    f = open(fileName, "w")
    f.write(jsonDict)
    f.close()

    return fileName

def readDesDictFromJson(fileName):
    """ read dict from json file """

    with open(fileName) as f:
        desDict = json.load(f)

    return desDict
