"""Run script for module com2AB
"""

import logging
import logging.config
import sys
import os
import shutil
import configparser
from datetime import datetime

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3SimpPlot import outAB
from avaframe.configLogger import ConfigLogger

from avaframe.test_run.avaframetask import *
import avaframe.test_run.avaframedata


class AvaframeTaskInitialize(AvaframeTask):

    def __init__self(self):

        super(AvaframeTaskInitialize, self).__init__()

    def cliName(self):

        return '-init'

    def name(self):

        return 'Initializing sequence'

    def resultPath(self):

        return 'initialization'

    def description(self):

        return ('Running initialization sequence. Creating repository structure. Creating local copies of config files')

    def validateData(self, data):

        # TODO fixme
        return True

    def run(self, data, callBack=None):

        """Run initialization sequence
        """
        pathAvaName = data.pathAvalancheName
        logFileName = os.path.join(pathAvaName, 'initialization.log')
        logConfigLocName = os.path.join(pathAvaName, 'local_logging.conf')
        logConfigName = 'avaframe/logging.conf'
        if os.path.exists(pathAvaName):
            logConf = True
            if not os.path.exists(logConfigName):
                shutil.copyfile(logConfigName, logConfigLocName)
                logConf = False
            logging.config.fileConfig(fname=logConfigName, defaults={'logfilename': logFileName}, disable_existing_loggers=False)
            log = logging.getLogger(__name__)
            log.info(datetime.now().strftime("%H_%M_%d_%m_%Y"))
            log.info("Running initialization sequecnce for avalanche %s", str(pathAvaName))
            log.warning("Folder %s already exists. Checking folder structure and adding missing files and folders.", str(pathAvaName).split('/')[-1])
            createFolderStruct(data)
        else:
            try:
                os.makedirs(pathAvaName)
            except IOError as e:
                print('I/O error({0}): {1}'.format(e.errno, e.strerror))
                return
            shutil.copyfile(logConfigName, logConfigLocName)
            logging.config.fileConfig(fname=logConfigName, defaults={'logfilename': logFileName}, disable_existing_loggers=False)
            log = logging.getLogger(__name__)
            log.info(datetime.now().strftime("%H_%M_%d_%m_%Y"))
            log.info("Running initialization sequecnce for avalanche %s", str(pathAvaName))
            createFolderStruct(data)
            log.info("Done initializing avalanche %s", str(pathAvaName).split('/')[-1])



def createFolderStruct(data):
    log = logging.getLogger(__name__)
    pathAvaName = data.pathAvalancheName
    Inputs = os.path.join(pathAvaName, 'Inputs/')
    if not os.path.exists(Inputs):
        os.makedirs(Inputs)
    Res = os.path.join(Inputs, 'RES/')
    if not os.path.exists(Res):
        os.makedirs(Res)
    Rel = os.path.join(Inputs, 'REL/')
    if not os.path.exists(Rel):
        os.makedirs(Rel)
    Ent = os.path.join(Inputs, 'ENT/')
    if not os.path.exists(Ent):
        os.makedirs(Ent)
    Points = os.path.join(Inputs, 'POINTS/')
    if not os.path.exists(Points):
        os.makedirs(Points)
    Lines = os.path.join(Inputs, 'LINES/')
    if not os.path.exists(Lines):
        os.makedirs(Lines)

    Outputs = os.path.join(pathAvaName, 'Outputs/')
    if not os.path.exists(Outputs):
        os.makedirs(Outputs)

    Work = os.path.join(pathAvaName, 'Work/')
    if not os.path.exists(Work):
        os.makedirs(Work)

    if data.locCfgAB:
        ABCfgLocName = os.path.join(pathAvaName, 'local_com2ABCfg.ini')
        ABCfgName = 'avaframe/com2AB/com2ABCfg.ini'
        if not os.path.exists(ABCfgLocName):
            log.info("Creating local copy of %s to %s", str(ABCfgName), str(ABCfgLocName))
            shutil.copyfile(ABCfgName, ABCfgLocName)
        else:
            log.warning("%s already exists", str(ABCfgLocName))

    if data.locCfg1DFA:
        DFACfgLocName = os.path.join(pathAvaName, 'local_com1DFACfg.ini')
        DFACfgName = 'avaframe/com1DFA/com1DFACfg.ini'
        if not os.path.exists(DFACfgLocName):
            log.info("Creating local copy of %s to %s", str(DFACfgName), str(DFACfgLocName))
            shutil.copyfile(DFACfgName, DFACfgLocName)
        else:
            log.warning("%s already exists", str(DFACfgLocName))
