"""Run script for module com2AB
"""

import logging
import logging.config
import sys
import os
import configparser
import time
from datetime import datetime

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3SimpPlot import outAB
from avaframe.configLogger import ConfigLogger

import avaframe.test_run.avaframedata
from avaframe.test_run.avaframetask import *


class AvaframeTaskAlphaBeta(AvaframeTask):

    def __init__self(self):

        super(AvaframeTaskAlphaBeta, self).__init__()

    def cliName(self):

        return '-alphabeta'

    def name(self):

        return 'Alpha Beta computation module'

    def resultPath(self):

        return 'alphabeta'

    def description(self):

        return ('Runs the Alpha Beta model. Parameter \"pathAvalancheName=<value>\" is requiered.\n' +
                'If no local_com2ABCfg configuration file is provided, the standard setting will be used.\n' +
                'Run \"python Main_avaframe.py -init locCfgAB=1\" if you want to change the com2ABCfg file.')

    def validateData(self, data):

        # TODO fixme
        return True

    def run(self, data, callBack=None):
        """Run script for module com2AB
        """

        startTime = time.time()

        # ---------------------------------------------
        # ---------------------------------------------
        #############################################
        # Load all input Parameters from congig file #
        #############################################

        cfg = configparser.ConfigParser(allow_no_value=True)
        pathAvaName = data.pathAvalancheName
        ABCfgLocName = os.path.join(pathAvaName, 'local_com2ABCfg.ini')
        ABCfgName = 'avaframe/com2AB/com2ABCfg.ini'
        if os.path.exists(ABCfgLocName):
            cfg.read(ABCfgLocName)
        else:
            cfg.read(ABCfgName)

        cfgPath = cfg['INPATH']
        cfgsetup = cfg['ABSETUP']
        cfgFlags = cfg['FLAGS']

        #########################################
        # Load all input Parameters for logging #
        #########################################
        logFileName = os.path.join(cfgPath['saveOutPath'], "runCom2AB.log")
        if os.path.exists(logFileName):
            os.remove(logFileName)
        logConfigLocName = os.path.join(pathAvaName, 'local_logging.conf')
        logging.config.fileConfig(fname=logConfigLocName, defaults={
                                  'logfilename': logFileName}, disable_existing_loggers=False)
        log = logging.getLogger(__name__)
        log.info(datetime.now().strftime("%H_%M_%d_%m_%Y"))
        config_logger = ConfigLogger(log)
        config_logger(cfg)

        # ---------------------------------------------
        # Start ALPHABETA
        # ---------------------------------------------

        #################################
        # Read input data for ALPHABETA #
        # Preprocessing
        #################################
        log.info("Running com2ABMain model on test case DEM \n %s \n with profile \n %s ",
                 cfgPath['DGMSource'], cfgPath['ProfileLayer'])

        DGM = com2AB.readRaster(cfgPath['DGMSource'])
        Avapath = com2AB.readAvaPath(cfgPath['ProfileLayer'], cfgPath['outputName'], DGM['header'])
        SplitPoint = com2AB.readSplitPoint(cfgPath['SplitPointSource'], DGM['header'])

        #################################
        # Calculate ALPHABETA #
        # Processing
        #################################
        com2AB.com2ABMain(DGM, Avapath, SplitPoint, cfgPath['saveOutPath'],
                          cfgsetup)

        #################################
        # Analyse/ plot/ write results #
        # Postprocessing
        #################################
        fileNamePlot_ext, fileNameWrite_ext = outAB.writeABpostOut(DGM,
                                                                   Avapath, SplitPoint,
                                                                   cfgPath['saveOutPath'],
                                                                   cfgFlags)

        log.info('Plotted to: %s', fileNamePlot_ext)
        log.info('Data written: %s', fileNameWrite_ext)

        endTime = time.time()
        log.info('Done with AlphaBeta module')
        log.info('Took %.2f seconds to run AlphaBeta module and process results.' %
                 (endTime - startTime))
