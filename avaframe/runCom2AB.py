"""Run script for module com2AB
"""

# Local imports
import avaframe as af
from avaframe.com2AB import com2AB
import avaframe.com2AB.com2ABCfg as conf

print("[com2AB] Running com2ABMain model on test case DEM : ",
      conf.DGMSource, 'with profile:', conf.ProfileLayer)

[header, rasterdata] = com2AB.ReadRaster(conf.DGMSource)
avapath = com2AB.ReadAvaPath(conf.ProfileLayer, header)
splitPoint = com2AB.ReadSplitPoint(conf.SplitPointSource, header)

# com2AB.com2ABMain(header, rasterdata, avapath, splitPoint, conf.saveOutPath, conf.smallAva)
