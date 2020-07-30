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

log.info("[com2AB] Running com2ABMain model on test case DEM : ",
      conf.DGMSource, 'with profile:', conf.ProfileLayer)

[header, rasterdata] = com2AB.readRaster(conf.DGMSource)
avapath = com2AB.readAvaPath(conf.ProfileLayer, header)
splitPoint = com2AB.readSplitPoint(conf.SplitPointSource, header)

com2AB.com2ABMain(header, rasterdata, avapath, splitPoint,
                  conf.saveOutPath, conf.smallAva, conf.distance)
FileNamePlot_ext, FileNameWrite_ext = outAB.writeABpostOut(header, rasterdata, avapath, splitPoint, conf.saveOutPath, conf.flags)

print(FileNamePlot_ext)
print(FileNameWrite_ext)
