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

    # ------------------------
    # Start time step computation
    Npart = particles['Npart']
    csz = dem['header'].cellsize
    # initialize
    Fnormal = np.zeros(Npart)
    forceX = np.zeros(Npart)
    forceY = np.zeros(Npart)
    forceZ = np.zeros(Npart)
    # loop on particles
    for i in range(Npart):
        mass = particles['m'][i]
        x = particles['x'][i]
        y = particles['y'][i]
        z = particles['z'][i]
        h = particles['h'][i]
        ux = particles['ux'][i]
        uy = particles['uy'][i]
        uz = particles['uz'][i]
        indCellX, indCellY = particles['InCell'][i]
        # deduce area
        A = mass / (h * rho)
        # get velocity verctor direction
        uxDir, uyDir, uzDir = tools.normalize(ux, uy, uz)
        # get normal at the particle location
        nx, ny, nz = tools.getNormal(x, y, Nx, Ny, Nz, csz)

        xEnd = x + dt * ux
        yEnd = y + dt * uy
        nxEnd, nyEnd, nzEnd = tools.getNormal(xEnd, yEnd, Nx, Ny, Nz, csz)

        nxAvg = nx + nxEnd
        nyAvg = ny + nyEnd
        nzAvg = nz + nzEnd
        nxAvg, nyAvg, nzAvg = tools.normalize(nxAvg, nyAvg, nzAvg)

        # acceleration due to curvature
        accNormCurv = -(ux*(nxEnd-nx) + uy*(nyEnd-ny) + uz*(nzEnd-nz))/dt

        gravAccNorm = - gravAcc * nzAvg
        effAccNorm = gravAccNorm + accNormCurv
        if(effAccNorm > 0.0):
            Fnormal[i] = mass * effAccNorm

        # body forces
        gravAccTangX =         + gravAccNorm * nxAvg
        gravAccTangY =         + gravAccNorm * nyAvg
        gravAccTangZ = gravAcc + gravAccNorm * nzAvg
        # adding gravity force contribution
        forceX[i] = forceX[i] + gravAccTangX * mass
        forceY[i] = forceY[i] + gravAccTangY * mass
        forceZ[i] = forceZ[i] + gravAccTangZ * mass

        uMag = tools.norm(ux, uy, uz)

        # Calculating bottom sheer and normal stress
        if(effAccNorm < 0.0):
            # if fluid detatched
            tau = 0.0
        else:
            # bottom normal stress sigmaB
            sigmaB = effAccNorm * rho * h
            # SamosAT friction type (bottom shear stress)
            tau = tools.SamosATfric(uMag, sigmaB, h)

        # adding bottom shear resistance contribution
        forceBotTang = -A * tau
        forceX[i] = forceX[i] + forceBotTang * uxDir
        forceY[i] = forceY[i] + forceBotTang * uyDir
        forceZ[i] = forceZ[i] + forceBotTang * uzDir

        # compute entrained mass
        dm = 0
        if Ment[indCellY][indCellX] > 0:
            # either erosion or ploughing but not both
            # width of the particle
            width = math.sqrt(A)
            # bottom area covered by the particle during dt
            ABotSwiped = width * uMag * dt
            if(entEroEnergy > 0):
                # erosion: erode according to shear and erosion energy
                dm = A * tau * uMag * dt / entEroEnergy
                Aent = A
            else:
                # ploughing in at avalanche front: erode full area weight
                # mass available in the cell [kg/m²]
                rhoHent = Ment[indCellY][indCellX]
                dm = rhoHent * ABotSwiped
                Aent = rhoHent / rhoEnt
            # adding mass balance contribution
            forceX[i] = forceX[i] - dm / dt * ux
            forceY[i] = forceY[i] - dm / dt * uy
            forceZ[i] = forceZ[i] - dm / dt * uz

            # adding force du to entrained mass
            Fent = width * (entShearResistance + dm / Aent * entDefResistance)
            forceX[i] = forceX[i] + Fent * uxDir
            forceY[i] = forceY[i] + Fent * uyDir
            forceZ[i] = forceZ[i] + Fent * uzDir

        # adding resistance force du to obstacles
        if Cres[indCellY][indCellX] > 0:
            if(h < hRes):
                hResEff = h
            cres = - rho * A * hResEff * Cres * uMag
            forceX[i] = forceX[i] + cres * ux
            forceY[i] = forceY[i] + cres * uy
            forceZ[i] = forceZ[i] + cres * uz
