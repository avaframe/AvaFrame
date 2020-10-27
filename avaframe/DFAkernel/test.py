import os
import glob
import logging
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.DFAkernel.tools as tools
# import avaframe.DFAkernel.test as test
import avaframe.in3Utils.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

if __name__ == "__main__":
    # log file name; leave empty to use default runLog.log
    logName = 'testKernel'

    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    inputDir = os.path.join(avalancheDir, 'Inputs')
    relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.shp')
    demFile = glob.glob(inputDir+os.sep+'*.asc')
    dem = IOf.readRaster(demFile[0])
    releaseLine = shpConv.readLine(relFiles[0], 'release1', dem)
    print(releaseLine['Name'])
    print(releaseLine['x'])
    print(releaseLine['Start'])
    print(releaseLine['y'])
    print(dem['header'])
    relRaster = tools.getRelease(dem, releaseLine)
    relTh = 1
    relRasterD = relRaster * relTh
    relRasterDSparse = sparse.csr_matrix(relRaster)
    indX, indY = np.nonzero(relRaster)
    print(indX, indY)
    print(relRasterDSparse)
