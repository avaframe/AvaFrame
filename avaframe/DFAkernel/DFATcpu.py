import os
import glob
import copy
import logging
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.DFAkernel.DFAtools as DFAtools
# from avaframe.DFAkernel.setParam import *
from avaframe.out3Plot.plotUtils import *
import avaframe.in2Trans.ascUtils as IOf
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
dem = copy.deepcopy(demOri)
dem['header'].xllcenter = 0
dem['header'].yllcenter = 0
dem['header'].xllcorner = 0
dem['header'].yllcorner = 0

# ------------------------
# process release to get it as a raster
relRaster = DFAtools.prepareArea(releaseLine, demOri)
relTh = 1
# could do something more advanced if we want varying release depth
relRasterD = relRaster * relTh

TForce = []
TForceVect = []
TPos = []
TNeigh = []
TField = []
MassPart = [5000, 1000, 500, 200, 100, 50]
cfg['dt'] = str(0.1)
cfg['Tend'] = str(0.1)
NP = []
for massPart in MassPart:
    cfg['massPerPart'] = str(massPart)
    # ------------------------
    # initialize simulation : create particles, create resistance and
    # entrainment matrix, initialize fields, get normals and neighbours
    dem, particles, fields, Cres, Ment = DFAtools.initializeSimulation(cfg, relRaster, dem)
    NP.append(particles['Npart'])
    log.info('Initializted simulation. M = %f kg, %s particles' % (particles['mTot'], particles['Npart']))

    # ------------------------
    #  Start time step computation
    Tcpu = {}
    Tcpu['Force'] = 0.
    Tcpu['ForceVect'] = 0.
    Tcpu['Pos'] = 0.
    Tcpu['Neigh'] = 0.
    Tcpu['Field'] = 0.

    Particles, Fields, Tcpu = DFAtools.DFAIterate(cfg, particles, fields, dem, Ment, Cres, Tcpu)

    log.info(('cpu time Force = %s s' % (Tcpu['Force'] / Tcpu['niter'])))
    TForce.append(Tcpu['Force'])
    log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / Tcpu['niter'])))
    TForceVect.append(Tcpu['ForceVect'])
    log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / Tcpu['niter'])))
    TPos.append(Tcpu['Pos'])
    log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / Tcpu['niter'])))
    TNeigh.append(Tcpu['Neigh'])
    log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / Tcpu['niter'])))
    TField.append(Tcpu['Field'])

# -------------------------------
fig, ax = plt.subplots(figsize=(figW, figH))
m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForce))
cm1lab = "TForce : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'b--', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TForce), 'ok', linestyle='-', label='Tcpu Force')
# ---------------------------------
m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForceVect))
cm1lab = "TForceVect : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'b-.', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TForceVect), '*k', linestyle='-', label='Tcpu Force Vect')
# -----------------------------------
m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TPos[2:]))
cm1lab = "TPos : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'g--', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TPos), 'sk', linestyle='-', label='Tcpu Position')
m, c, r, p, se1 = stats.linregress(np.log(NP[0:3]), np.log(TNeigh[0:3]))
cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
m, c, r, p, se1 = stats.linregress(np.log(NP[4:]), np.log(TNeigh[4:]))
cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'r-.', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TNeigh), 'dk', linestyle='-', label='Tcpu Neighbours')
ax.plot(np.log(NP), np.log(TField), '+k', linestyle='-', label='Tcpu Fields')
plt.legend(loc='upper left')
# plt.show()

fig1, ax1 = plt.subplots(figsize=(figW, figH))
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForce))
# cm1lab = "TForce : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax1.loglog(NP, m*np.log(NP)+c, 'b--', linewidth=2, label=cm1lab)
ax1.loglog(NP, TForce, 'ok', linestyle='-', label='Tcpu Force')
ax1.plot(NP, TForceVect, '*k', linestyle='-', label='Tcpu Force Vect')
ax1.loglog(NP, TPos, 'sk', linestyle='-', label='Tcpu Position')
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TNeigh))
# cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax1.loglog(NP, m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
ax1.loglog(NP, TNeigh, 'dk', linestyle='-', label='Tcpu Neighbours')
ax1.loglog(NP, TField, '+k', linestyle='-', label='Tcpu Fields')
plt.legend(loc='upper left')
plt.show()
