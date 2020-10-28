import os
import glob
import logging
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

    # fetch input data
    inputDir = os.path.join(avalancheDir, 'Inputs')
    relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.shp')
    demFile = glob.glob(inputDir+os.sep+'*.asc')
    dem = IOf.readRaster(demFile[0])
    releaseLine = shpConv.readLine(relFiles[0], 'release1', dem)

    # process release to get it as a raster
    relRaster = tools.polygon2Raster(dem['header'], releaseLine)
    relTh = 1
    relRasterD = relRaster * relTh

    # initialize simulation : create particles
    particles = tools.initializeRelease(relRaster, dem)

    # get normal vector of the grid mesh
    Nx, Ny, Nz = tools.getNormalVect(dem['rasterData'], dem['header'].cellsize)
    mesh = {}
    # mesh['header'] = dem['header']
    mesh['Nx'] = Nx
    mesh['Ny'] = Ny
    mesh['Nz'] = Nz

    Npart = particles['Npart']
    csz = dem['header'].cellsize
    # loop on particles
    for i in range(Npart):
        m = particles['m'][i]
        x = particles['x'][i]
        y = particles['y'][i]
        z = particles['z'][i]
        d = particles['d'][i]
        ux = particles['ux'][i]
        uy = particles['uy'][i]
        uz = particles['uz'][i]
        # get normal at the particle location
        Point = {}
        Point['x'] = x
        Point['y'] = y

        nx = geoTrans.projectOnRasterRoot(x, y, Nx, csz=csz)
        ny = geoTrans.projectOnRasterRoot(x, y, Ny, csz=csz)
        nz = geoTrans.projectOnRasterRoot(x, y, Nz, csz=csz)
        nx, ny, nz = tools.normalize(nx, ny, nz)

        xEnd = x + dt * ux
        yEnd = y + dt * uy
        nxEnd = geoTrans.projectOnRasterRoot(xEnd, yEnd, Nx, csz=csz)
        nyEnd = geoTrans.projectOnRasterRoot(xEnd, yEnd, Ny, csz=csz)
        nzEnd = geoTrans.projectOnRasterRoot(xEnd, yEnd, Nz, csz=csz)
        nxEnd, nyEnd, nzEnd = tools.normalize(nxEnd, nyEnd, nzEnd)

        nxAvg = nx + nxEnd
        nyAvg = ny + nyEnd
        nzAvg = nz + nzEnd
        nxAvg, nyAvg, nzAvg = tools.normalize(nxAvg, nyAvg, nzAvg)

        accxNormCurv = - ux * (nxEnd-nx) / dt
        accyNormCurv = - uy * (nyEnd-ny) / dt
        acczNormCurv = - uz * (nzEnd-nz) / dt
