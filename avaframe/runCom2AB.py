"""Run script for module com2AB
"""

import logging
import logging.config
import sys
import os
import configparser

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3SimpPlot import outAB

#---------------------------------------------
#---------------------------------------------
#############################################
# Load all input Parameters from congig file #
#############################################
cfg = configparser.ConfigParser(allow_no_value=True)
if os.path.isfile('avaframe/com2AB/local_com2ABCfg.ini'):
    cfg.read('avaframe/com2AB/local_com2ABCfg.ini')
else:
    cfg.read('avaframe/com2AB/com2ABCfg.ini')

cfgPath = cfg['INPATH']
cfgsetup = cfg['ABSETUP']
cfgFlags = cfg['FLAGS']

#########################################
# Load all input Parameters for logging #
#########################################
logFileName = os.path.join(cfgPath['saveOutPath'], "runCom2AB.log")
logging.config.fileConfig(fname='avaframe/logging.conf', defaults={'logfilename': logFileName}, disable_existing_loggers=False)
log = logging.getLogger(__name__)


#---------------------------------------------
# Start ALPHABETA
#---------------------------------------------

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
