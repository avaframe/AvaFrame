"""Run script for module com2AB
"""

# Local imports
from avaframe.com2AB import com2AB
import avaframe.com2AB.com2ABCfg as conf

print("[com2AB] Running com2ABMain model on test case DEM : ",
      conf.DGMSource, 'with profile:', conf.ProfileLayer)

[header, rasterdata] = com2AB.readRaster(conf.DGMSource)
avapath = com2AB.readAvaPath(conf.ProfileLayer, header)
splitPoint = com2AB.readSplitPoint(conf.SplitPointSource, header)

# com2AB.com2ABMain(header, rasterdata, avapath, splitPoint, conf.saveOutPath, conf.smallAva)
