import os
import glob
import logging
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.DFAkernel.tools as tools
from avaframe.DFAkernel.setParam import *
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

    # ------------------------
    # fetch input data
    inputDir = os.path.join(avalancheDir, 'Inputs')
    relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.shp')
    demFile = glob.glob(inputDir+os.sep+'*.asc')
    dem = IOf.readRaster(demFile[0])
    releaseLine = shpConv.readLine(relFiles[0], 'release1', dem)

    # ------------------------
    # process release to get it as a raster
    relRaster = tools.polygon2Raster(dem['header'], releaseLine)
    relTh = 1
    # could do something more advanced if we want varying release depth
    relRasterD = relRaster * relTh

    # ------------------------
    # initialize simulation : create particles, create resistance and
    # entrqinment matrix
    particles, Cres, Ment = tools.initializeSimulation(relRaster, dem)

    # get normal vector of the grid mesh
    Nx, Ny, Nz = tools.getNormalVect(dem['rasterData'], dem['header'].cellsize)
    dem['Nx'] = Nx
    dem['Ny'] = Ny
    dem['Nz'] = Nz

    # ------------------------
    # Start time step computation
    iterate = True
    t = 0
    X = particles['x'][0]
    Y = particles['y'][0]
    Z = particles['z'][0]
    UX = particles['ux'][0]
    UY = particles['uy'][0]
    UZ = particles['uz'][0]
    while t < Tend and iterate:
        log.info('Computing time step t = %f s', t)
        # get particles location (neighbours for sph)
        particles = tools.getNeighbours(particles, dem)
        # get forces
        force = tools.computeForce(particles, dem, Ment, Cres)
        # get forces sph
        # forceSPH = tools.computeForceSPH(particles, dem)
        # update velocity and particle position
        particles = tools.updatePosition(particles, dem, force)
        X = np.append(X, particles['x'][0])
        Y = np.append(Y, particles['y'][0])
        Z = np.append(Z, particles['z'][0])
        UX = np.append(UX, particles['ux'][0])
        UY = np.append(UY, particles['uy'][0])
        UZ = np.append(UZ, particles['uz'][0])
        # tools.plotPosition(particles, dem)
        t = t + dt
    print(tools.norm(UX, UY, UZ))
    print(gravAcc * math.sin(34*math.pi/180) * np.arange(451) * dt)
    # print(gravAcc * math.cos(34*math.pi/180) * (math.tan(34*math.pi/180) - mu) * np.arange(201) * dt)
