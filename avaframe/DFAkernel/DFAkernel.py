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
import avaframe.DFAkernel.DFAtools as tools
from avaframe.DFAkernel.setParam import *
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

# ------------------------
# fetch input data
inputDir = os.path.join(avalancheDir, 'Inputs')
relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.shp')
demFile = glob.glob(inputDir+os.sep+'*.asc')
demOri = IOf.readRaster(demFile[0])
releaseLine = shpConv.readLine(relFiles[0], 'release1', demOri)

# ------------------------
# process release to get it as a raster
relRaster = tools.polygon2Raster(demOri['header'], releaseLine)
relTh = 1
# could do something more advanced if we want varying release depth
relRasterD = relRaster * relTh

# ------------------------
# initialize simulation : create particles, create resistance and
# entrqinment matrix
dem = demOri.copy()
dem['header'].xllcenter = 0
dem['header'].yllcenter = 0
dem['header'].xllcorner = 0
dem['header'].yllcorner = 0
particles, fields, Cres, Ment = tools.initializeSimulation(relRaster, dem)
log.info('Initializted simulation. M = %f kg, %s particles' % (particles['mTot'], particles['Npart']))

# get normal vector of the grid mesh
Nx, Ny, Nz = tools.getNormalVect(dem['rasterData'], dem['header'].cellsize)
dem['Nx'] = Nx
dem['Ny'] = Ny
dem['Nz'] = Nz

# ------------------------
# Start time step computation
iterate = True
t = 0
particles['t'] = t
results = [particles.copy()]
Fields = [fields.copy()]
nSave = 1
while t < Tend and iterate:
    t = t + dt
    log.info('Computing time step t = %f s', t)
    particles['t'] = t
    # get particles location (neighbours for sph)
    particles = tools.getNeighbours(particles, dem)
    # get forces
    force = tools.computeForce(particles, dem, Ment, Cres)
    # get forces sph
    # forceSPH = tools.computeForceSPH(particles, dem)
    # update velocity and particle position
    particles = tools.updatePosition(particles, dem, force)
    # update fields (compute grid values)
    fields = tools.updateFields(particles, dem, fields)
    if t >= nSave * dtSave:
        log.info('Saving results for time step t = %f s', t)
        results.append(particles.copy())
        Fields.append(fields.copy())
        nSave = nSave + 1

# tools.plotPosition(particles, dem)
partRef = results[0]
print(partRef['z'][0])
fig, ax = plt.subplots(figsize=(figW, figH))
fig1, ax1 = plt.subplots(figsize=(figW, figH))
for part, field in zip(results, Fields):
    print(part['t'])
    print(tools.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
    # exact solution for inclined plane with no friction
    # print(gravAcc * math.sin(34*math.pi/180) * part['t'])
    # print(math.sqrt(2 * gravAcc * abs(partRef['z'][0] - part['z'][0])))

    # exact solution for inclined plane with friction
    print(gravAcc * math.cos(34*math.pi/180) * (math.tan(34*math.pi/180) - mu) * part['t'])
    fig, ax = tools.plotPosition(part, dem, dem['rasterData'], fig, ax)
    fig1, ax1 = tools.plotPosition(part, dem, field['FD'], fig1, ax1)
    plt.pause(1)
plt.show()
