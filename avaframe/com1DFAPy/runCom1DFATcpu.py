import os
import glob
import copy
import logging
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Local imports
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA
# from avaframe.DFAkernel.setParam import *
import avaframe.out3Plot.plotUtils as pU
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

# Load configuration
cfg = cfgUtils.getModuleConfig(com1DFA)
cfgGen = cfg['GENERAL']
flagDev = cfg['FLAGS'].getboolean('flagDev')


# ------------------------
# fetch input data
demFile, relFiles, entFiles, resFile, flagEntRes = gI.getInputData(
    avalancheDir, cfg['FLAGS'], flagDev=False)
demOri = IOf.readRaster(demFile)
releaseLine = shpConv.readLine(relFiles[0], 'release1', demOri)
dem = copy.deepcopy(demOri)
dem['header'].xllcenter = 0
dem['header'].yllcenter = 0
dem['header'].xllcorner = 0
dem['header'].yllcorner = 0

# -----------------------
# Initialize mesh
dem = com1DFA.initializeMesh(dem)
# ------------------------
# process release to get it as a raster
relRaster = com1DFA.prepareArea(releaseLine, demOri)
relTh = 1
# could do something more advanced if we want varying release depth
relRasterD = relRaster * relTh

TForce = []
TForceVect = []
TForceSPH = []
TPos = []
TNeigh = []
TField = []
MassPart = [5000, 1000, 500, 250, 150]
cfgGen['dt'] = str(0.1)
cfgGen['Tend'] = str(0.1)
NP = []
for massPart in MassPart:
    cfgGen['massPerPart'] = str(massPart)
    # ------------------------
    # initialize simulation : create particles, create resistance and
    # entrainment matrix, initialize fields, get normals and neighbours
    particles, fields, Cres, Ment = com1DFA.initializeSimulation(cfgGen, relRaster, dem)
    NP.append(particles['Npart'])
    log.info('Initializted simulation. M = %f kg, %s particles' %
             (particles['mTot'], particles['Npart']))

    # ------------------------
    #  Start time step computation
    Tcpu = {}
    Tcpu['Force'] = 0.
    Tcpu['ForceVect'] = 0.
    Tcpu['ForceSPH'] = 0.
    Tcpu['Pos'] = 0.
    Tcpu['Neigh'] = 0.
    Tcpu['Field'] = 0.

    T, U, Z, S, Particles, Fields, Tcpu = com1DFA.DFAIterate(
        cfgGen, particles, fields, dem, Ment, Cres, Tcpu)

    log.info(('cpu time Force = %s s' % (Tcpu['Force'] / Tcpu['nIter'])))
    TForce.append(Tcpu['Force'])
    log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / Tcpu['nIter'])))
    TForceVect.append(Tcpu['ForceVect'])
    log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / Tcpu['nIter'])))
    TForceSPH.append(Tcpu['ForceSPH'])
    log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / Tcpu['nIter'])))
    TPos.append(Tcpu['Pos'])
    log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / Tcpu['nIter'])))
    TNeigh.append(Tcpu['Neigh'])
    log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / Tcpu['nIter'])))
    TField.append(Tcpu['Field'])

# -------------------------------
fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
# -------------------------------
m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForce))
cm1lab = "TForce : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'b--', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TForce), 'ok', linestyle='-', label='Tcpu Force')
# ---------------------------------
m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TForceVect[2:]))
cm1lab = "TForceVect : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'b-.', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TForceVect), '*k', linestyle='-', label='Tcpu Force Vect')
# ---------------------------------
m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TForceSPH[2:]))
cm1lab = "TForceSPH : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'g--', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TForceSPH), '^k', linestyle='-', label='Tcpu Force SPH')
# -----------------------------------
# m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TPos[2:]))
# cm1lab = "TPos : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax.plot(np.log(NP), m*np.log(NP)+c, 'g--', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TPos), 'sk', linestyle='-', label='Tcpu Position')
# -----------------------------------
m, c, r, p, se1 = stats.linregress(np.log(NP[3:]), np.log(TNeigh[3:]))
cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
ax.plot(np.log(NP), m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
# m, c, r, p, se1 = stats.linregress(np.log(NP[4:]), np.log(TNeigh[4:]))
# cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax.plot(np.log(NP), m*np.log(NP)+c, 'r-.', linewidth=2, label=cm1lab)
ax.plot(np.log(NP), np.log(TNeigh), 'dk', linestyle='-', label='Tcpu Neighbours')
# -----------------------------------
ax.plot(np.log(NP), np.log(TField), '+k', linestyle='-', label='Tcpu Fields')
plt.legend(loc='upper left')
# plt.show()

fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForce))
# cm1lab = "TForce : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax1.loglog(NP, m*np.log(NP)+c, 'b--', linewidth=2, label=cm1lab)
ax1.loglog(NP, TForce, 'ok', linestyle='-', label='Tcpu Force')
ax1.plot(NP, TForceVect, '*k', linestyle='-', label='Tcpu Force Vect')
ax1.plot(NP, TForceSPH, '^k', linestyle='-', label='Tcpu Force SPH')
ax1.loglog(NP, TPos, 'sk', linestyle='-', label='Tcpu Position')
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TNeigh))
# cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax1.loglog(NP, m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
ax1.loglog(NP, TNeigh, 'dk', linestyle='-', label='Tcpu Neighbours')
ax1.loglog(NP, TField, '+k', linestyle='-', label='Tcpu Fields')
plt.legend(loc='upper left')
plt.show()
