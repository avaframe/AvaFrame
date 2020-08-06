"""Run script for module com2AB
"""

import logging
import sys
import os

# Local imports
from avaframe.com2AB import com2AB
import avaframe.com2AB.com2ABCfg as conf
from avaframe.out3SimpPlot import outAB

# create logger, set to logging.DEBUG to see all messages
# logging.basicConfig(filename='com2AB.log', filemode='w', level=logging.INFO,
#                     format='%(module)s:%(levelname)s - %(message)s')
# logging.basicConfig(stream=sys.stdout, level=logging.INFO,
#                     format='%(module)s:%(levelname)s - %(message)s')
logFileName = os.path.join(conf.saveOutPath, "runCom2AB.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(module)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(logFileName, "w"),
        logging.StreamHandler()
    ]
)
log = logging.getLogger(__name__)

log.info("Running com2ABMain model on test case DEM \n %s \n with profile \n %s ",
         conf.DGMSource, conf.ProfileLayer)

DGM = com2AB.readRaster(conf.DGMSource)
Avapath = com2AB.readAvaPath(conf.ProfileLayer, conf.outputName, DGM['header'])
SplitPoint = com2AB.readSplitPoint(conf.SplitPointSource, DGM['header'])

com2AB.com2ABMain(DGM, Avapath, SplitPoint, conf.saveOutPath,
                  conf.smallAva, conf.customParam, conf.distance)

fileNamePlot_ext, fileNameWrite_ext = outAB.writeABpostOut(DGM,
                                                           Avapath, SplitPoint,
                                                           conf.saveOutPath,
                                                           conf.flags)

log.info('Plotted to: %s', fileNamePlot_ext)
log.info('Data written: %s', fileNameWrite_ext)
