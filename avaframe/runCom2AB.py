"""Run script for module com2AB
"""

import logging
import sys

# Local imports
from avaframe.com2AB import com2AB
import avaframe.com2AB.com2ABCfg as conf
from avaframe.out3SimpPlot import outAB

# create logger, set to logging.DEBUG to see all messages
# logging.basicConfig(filename='com2AB.log', filemode='w', level=logging.INFO,
#                     format='%(module)s:%(levelname)s - %(message)s')
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(module)s:%(levelname)s - %(message)s')
log = logging.getLogger(__name__)

log.info("Running com2ABMain model on test case DEM %s with profile %s ",
         conf.DGMSource, conf.ProfileLayer)

[header, rasterdata] = com2AB.readRaster(conf.DGMSource)
avapath = com2AB.readAvaPath(conf.ProfileLayer, conf.outputName, header)
splitPoint = com2AB.readSplitPoint(conf.SplitPointSource, header)

com2AB.com2ABMain(header, rasterdata, avapath, splitPoint,
                    conf.saveOutPath, conf.smallAva, conf.distance)

fileNamePlot_ext, fileNameWrite_ext = outAB.writeABpostOut(header, rasterdata,
                                                           avapath, splitPoint,
                                                           conf.saveOutPath,
                                                           conf.flags)

log.info('Plotted to: %s', fileNamePlot_ext)
log.info('Data written: %s', fileNameWrite_ext)
