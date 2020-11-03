import os
import glob
import time
import logging
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.DFAkernel.DFAtools as DFAtools
# from avaframe.DFAkernel.setParam import *
from avaframe.out3Plot.plotUtils import *
import avaframe.in3Utils.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'testKernel'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

cfg = cfgUtils.getModuleConfig(DFAtools)['GENERAL']

# ------------------------
# fetch input data
inputDir = os.path.join(avalancheDir, 'Inputs')
relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.shp')
demFile = glob.glob(inputDir+os.sep+'*.asc')
demOri = IOf.readRaster(demFile[0])
releaseLine = shpConv.readLine(relFiles[0], 'release1', demOri)

# ------------------------
# process release to get it as a raster
relRaster = DFAtools.polygon2Raster(demOri['header'], releaseLine)
relTh = 1
# could do something more advanced if we want varying release depth
relRasterD = relRaster * relTh

TForce = []
TPos = []
TNeigh = []
TField = []
MassPart = [5000, 2500, 1000, 750, 500, 250, 100, 75, 50]
NP = []
for massPart in MassPart:
    cfg['massPerPart'] = str(massPart)
    # ------------------------
    # initialize simulation : create particles, create resistance and
    # entrqinment matrix
    dem = demOri.copy()
    dem['header'].xllcenter = 0
    dem['header'].yllcenter = 0
    dem['header'].xllcorner = 0
    dem['header'].yllcorner = 0
    particles, fields, Cres, Ment = DFAtools.initializeSimulation(cfg, relRaster, dem)
    log.info('Initializted simulation. M = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
    # get particles location (neighbours for sph)
    particles = DFAtools.getNeighbours(particles, dem)
    # update fields (compute grid values)
    t = 0
    particles['t'] = t
    fields = DFAtools.updateFields(cfg, particles, dem, fields)
    # get normal vector of the grid mesh
    Nx, Ny, Nz = DFAtools.getNormalVect(dem['rasterData'], dem['header'].cellsize)
    dem['Nx'] = Nx
    dem['Ny'] = Ny
    dem['Nz'] = Nz

    # ------------------------
    # Start time step computation
    iterate = True
    results = [particles.copy()]
    Fields = [fields.copy()]
    TcpuForce = 0.
    TcpuPos = 0.
    TcpuNeigh = 0.
    TcpuField = 0.
    Tend = cfg.getfloat('Tend')
    dtSave = cfg.getfloat('dtSave')
    dt = cfg.getfloat('dt')
    nSave = 1
    niter = 0
    while t < Tend and iterate:
        t = t + dt
        niter = niter + 1
        log.debug('Computing time step t = %f s', t)
        particles['t'] = t
        NP.append(particles['Npart'])
        # get forces
        startTime = time.time()
        force = DFAtools.computeForce(cfg, particles, dem, Ment, Cres)
        tcpuForce = time.time() - startTime
        TcpuForce = TcpuForce + tcpuForce
        # log.info(('cpu time Force = %s s' % (tcpuForce)))
        # get forces sph
        # forceSPH = tools.computeForceSPH(particles, dem)
        # update velocity and particle position
        startTime = time.time()
        particles = DFAtools.updatePosition(cfg, particles, dem, force)
        tcpuPos = time.time() - startTime
        TcpuPos = TcpuPos + tcpuPos
        # log.info(('cpu time Position = %s s' % (tcpuPos)))
        # get particles location (neighbours for sph)
        startTime = time.time()
        particles = DFAtools.getNeighbours(particles, dem)
        tcpuNeigh = time.time() - startTime
        TcpuNeigh = TcpuNeigh + tcpuNeigh
        # log.info(('cpu time Neighbour = %s s' % (tcpuNeigh)))
        # update fields (compute grid values)
        startTime = time.time()
        fields = DFAtools.updateFields(cfg, particles, dem, fields)
        tcpuField = time.time() - startTime
        TcpuField = TcpuField + tcpuField
        # log.info(('cpu time Fields = %s s' % (tcpuField)))
        if t >= nSave * dtSave:
            log.info('Saving results for time step t = %f s', t)
            results.append(particles.copy())
            Fields.append(fields.copy())
            nSave = nSave + 1

    log.info(('cpu time Force = %s s' % (TcpuForce / niter)))
    TForce.append(TcpuForce)
    log.info(('cpu time Position = %s s' % (TcpuPos / niter)))
    TPos.append(TcpuPos)
    log.info(('cpu time Neighbour = %s s' % (TcpuNeigh / niter)))
    TNeigh.append(TcpuNeigh)
    log.info(('cpu time Fields = %s s' % (TcpuField / niter)))
    TField.append(TcpuField)
    # tools.plotPosition(particles, dem)
    partRef = results[0]
    print(partRef['z'][0])
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    mu = cfg.getfloat('mu')
    for part, field in zip(results, Fields):
        print(part['t'])
        print(DFAtools.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
        # exact solution for inclined plane with no friction
        # print(gravAcc * math.sin(34*math.pi/180) * part['t'])
        # print(math.sqrt(2 * gravAcc * abs(partRef['z'][0] - part['z'][0])))

        # exact solution for inclined plane with friction
        print(gravAcc * math.cos(34*math.pi/180) * (math.tan(34*math.pi/180) - mu) * part['t'])

    #     DFAtools.plotPosition(part, dem, dem['rasterData'])
    #     DFAtools.plotPosition(part, dem, field['PFD'])
    # plt.show()
    # fieldRef = Fields[0]
    # DFAtools.plotPosition(part, dem, fieldRef['PFD'])

print(TForce)
print(TPos)
print(TNeigh)
print(TField)
fig, ax = plt.subplots(figsize=(figW, figH))
ax.plot(NP, TForce, 'ok', linestyle='-')
ax.plot(NP, TPos, 'sk', linestyle='-')
ax.plot(NP, TNeigh, 'dk', linestyle='-')
ax.plot(NP, TField, '+k', linestyle='-')
plt.show()
