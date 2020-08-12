"""Run script for module com2AB
"""

import logging
import logging.config
import sys
import os

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3SimpPlot import outAB
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils.configLogger import ConfigLogger

#-----------Required settings-----------------

# TODO: move this to a main cfg file! (alongside flags for intermediate data/results/plots)
# and have function cfgUtils.readGlobalCfg that allows for a local_MainCfg

avalancheDir = 'data/avaSlide'

#---------------------------------------------

#############################################
# Load all input Parameters from config file #
#############################################

# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(com2AB)

cfgSetup = cfg['ABSETUP']
cfgFlags = cfg['FLAGS']
print(cfg.sections())

#########################################
# Load all input Parameters for logging #
#########################################
logFileName = os.path.join(cfgPath['saveOutPath'], "runCom2AB.log")
logging.config.fileConfig(fname='avaframe/in3Utils/logging.conf', defaults={'logfilename': logFileName})
log = logging.getLogger(__name__)
config_logger = ConfigLogger(log)
config_logger(cfg)

#---------------------------------------------
# Start ALPHABETA
#---------------------------------------------

#################################
# Read input data for ALPHABETA #
# Preprocessing
#################################

# extract input file locations
cfgPath = com2AB.readABinputs(avalancheDir)

log.info("Running com2ABMain model on test case DEM \n %s \n with profile \n %s ",
         cfgPath['DGMSource'], cfgPath['ProfileLayer'])

DGM = com2AB.readRaster(cfgPath['DGMSource'])
Avapath = com2AB.readAvaPath(cfgPath['ProfileLayer'], cfgPath['DefaltName'], DGM['header'])
SplitPoint = com2AB.readSplitPoint(cfgPath['SplitPointSource'], DGM['header'])

#################################
# Calculate ALPHABETA #
# Processing
#################################
com2AB.com2ABMain(DGM, Avapath, SplitPoint, cfgPath['saveOutPath'],
                  cfgSetup)


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
