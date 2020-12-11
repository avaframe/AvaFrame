import time
import glob
import copy
import logging
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
from avaframe.out3Plot.plotUtils import *
import avaframe.com1DFAPy.timeDiscretizations as tD
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.SPHfunctions as SPH
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.in2Trans.ascUtils as IOf
import avaframe.com1DFAPy.frictionLaws as fricLaws
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# import pyximport
# pyximport.install()
import SPHfunctionsCython as SPHC
import Tcpu as Tcpucore

# log file name; leave empty to use default runLog.log
logName = 'TcpuKernel'

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
rho = cfgGen.getfloat('rho')

# ------------------------
Lx = 50
Ly = 100
dt = 0.1
Lam = [1, 2, 3, 4, 5]
Nlam = len(Lam)
CSZ = [1, 2.5, 5]
Ncsz = len(CSZ)
TForce = np.zeros((Nlam, Ncsz))
TForceVect = np.zeros((Nlam, Ncsz))
TForceSPH = np.zeros((Nlam, Ncsz))
TForceSPHC = np.zeros((Nlam, Ncsz))
TFDSPH = np.zeros((Nlam, Ncsz))
TFDSPHC = np.zeros((Nlam, Ncsz))
TPos = np.zeros((Nlam, Ncsz))
TNeigh = np.zeros((Nlam, Ncsz))
TNeighC = np.zeros((Nlam, Ncsz))
TField = np.zeros((Nlam, Ncsz))
NP = np.zeros((Nlam, Ncsz))
count1 = 0
for csz in CSZ:
    count2 = 0
    dem = Tcpucore.createDEM(Lx, Ly, csz)
    # -----------------------
    # Initialize mesh
    dem = com1DFA.initializeMesh(dem)
    Ment = np.zeros(np.shape(dem['rasterData']))
    Cres = np.zeros(np.shape(dem['rasterData']))
    for lam in Lam:
        particles, fields = Tcpucore.createParticles(dem, lam, rho)
        NP[count2, count1] = particles['Npart']
        log.info('Initializted simulation. M = %f kg, %s particles' %
                 (particles['mTot'], particles['Npart']))
        # Xpart = particles['x']
        # Ypart = particles['y']
        # fig, ax = plt.subplots(figsize=(figW, figH))
        # ax.plot(Xpart, Ypart, color='r', marker='.', linestyle='None')
        # plt.show()
        # ------------------------
        #  Start time step computation
        Tcpu = {}
        Tcpu['Force'] = 0.
        Tcpu['ForceVect'] = 0.
        Tcpu['ForceSPH'] = 0.
        Tcpu['ForceSPHC'] = 0.
        Tcpu['FDSPH'] = 0.
        Tcpu['FDSPHC'] = 0.
        Tcpu['Pos'] = 0.
        Tcpu['Neigh'] = 0.
        Tcpu['NeighC'] = 0.
        Tcpu['Field'] = 0.

        # get forces
        # loop version of the compute force
        startTime = time.time()
        # forceLoop = computeForce(cfg, particles, dem, Ment, Cres)
        tcpuForce = time.time() - startTime
        Tcpu['Force'] = Tcpu['Force'] + tcpuForce
        # vectorized version of the compute force
        startTime = time.time()
        force = com1DFA.computeForceVect(cfgGen, particles, dem, Ment, Cres, dt)
        tcpuForceVect = time.time() - startTime
        Tcpu['ForceVect'] = Tcpu['ForceVect'] + tcpuForceVect
        # compute lateral force (SPH component of the calculation)
        startTime = time.time()
        # particles, force = com1DFA.computeForceSPH(cfgGen, particles, force, dem)
        tcpuForceSPH = time.time() - startTime
        Tcpu['ForceSPH'] = Tcpu['ForceSPH'] + tcpuForceSPH

        startTime = time.time()
        particles, force = SPHC.computeForceSPHC(cfgGen, particles, force, dem)
        tcpuForceSPHC = time.time() - startTime
        Tcpu['ForceSPHC'] = Tcpu['ForceSPHC'] + tcpuForceSPHC

        particles2 = copy.deepcopy(particles)
        # update velocity and particle position
        startTime = time.time()
        particles2 = com1DFA.updatePosition(cfgGen, particles2, dem, force)
        tcpuPos = time.time() - startTime
        Tcpu['Pos'] = Tcpu['Pos'] + tcpuPos

        # get particles location (neighbours for sph)
        startTime = time.time()
        # particles = getNeighbours(particles, dem)
        particles = SPH.getNeighboursVect(particles, dem)
        tcpuNeigh = time.time() - startTime
        Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

        particles3 = copy.deepcopy(particles)
        startTime = time.time()
        particles3 = SPHC.getNeighboursC(particles3, dem)
        tcpuNeighCython = time.time() - startTime
        Tcpu['NeighC'] = Tcpu['NeighC'] + tcpuNeighCython

        # get SPH flow depth
        startTime = time.time()
        # particles = SPH.computeFlowDepth(cfgGen, particles, dem)
        tcpuFD = time.time() - startTime
        Tcpu['FDSPH'] = Tcpu['FDSPH'] + tcpuFD
        # particles = SPH.computeFlowDepth(cfg, particles, dem)
        startTime = time.time()
        header = dem['header']
        Nx = dem['Nx']
        Ny = dem['Ny']
        Nz = dem['Nz']
        indX = (particles['InCell'][:, 0]).astype('int')
        indY = (particles['InCell'][:, 1]).astype('int')
        nx, ny, nz = DFAtls.getNormalArray(particles['x'], particles['y'], Nx, Ny, Nz, csz)
        H = SPHC.computeFDC(particles, header, nx, ny, nz, indX, indY)
        H = np.asarray(H)
        particles['hSPH'] = H
        tcpuFDC = time.time() - startTime
        Tcpu['FDSPHC'] = Tcpu['FDSPHC'] + tcpuFDC


        # update fields (compute grid values)
        startTime = time.time()
        particles, fields = com1DFA.updateFields(cfgGen, particles, force, dem, fields)
        tcpuField = time.time() - startTime
        Tcpu['Field'] = Tcpu['Field'] + tcpuField

        log.info(('cpu time Force = %s s' % (Tcpu['Force'])))
        TForce[count2, count1] = Tcpu['Force']
        log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'])))
        TForceVect[count2, count1] = Tcpu['ForceVect']
        log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'])))
        TForceSPH[count2, count1] = Tcpu['ForceSPH']
        log.info(('cpu time ForceSPHC = %s s' % (Tcpu['ForceSPHC'])))
        TForceSPHC[count2, count1] = Tcpu['ForceSPHC']
        log.info(('cpu time Position = %s s' % (Tcpu['Pos'])))
        TPos[count2, count1] = Tcpu['Pos']
        log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'])))
        TNeigh[count2, count1] = Tcpu['Neigh']
        log.info(('cpu time NeighbourC = %s s' % (Tcpu['NeighC'])))
        TNeighC[count2, count1] = Tcpu['NeighC']
        log.info(('cpu time FD = %s s' % (Tcpu['FDSPH'])))
        TFDSPH[count2, count1] = Tcpu['FDSPH']
        log.info(('cpu time FDC = %s s' % (Tcpu['FDSPHC'])))
        TFDSPHC[count2, count1] = Tcpu['FDSPHC']
        log.info(('cpu time Fields = %s s' % (Tcpu['Field'])))
        TField[count2, count1] = Tcpu['Field']
        count2 = count2 + 1

    count1 = count1 + 1

fig, ax = plt.subplots(figsize=(figW, figH))
fig1, ax1 = plt.subplots(figsize=(figW, figH))
for ncell in range(Ncsz):
    # -------------------------------

    # ---------------------------------
    # m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TForceVect[2:]))
    # cm1lab = "TForceVect : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax.plot(np.log(NP), m*np.log(NP)+c, 'b-.', linewidth=2, label=cm1lab)
    ax.plot(np.log(NP[:, ncell]), np.log(TForceVect[:, ncell]), '*k', linestyle='-', label='Tcpu Force Vect')
    # ---------------------------------
    # m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TForceSPH[2:]))
    # cm1lab = "TForceSPH : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax.plot(np.log(NP), m*np.log(NP)+c, 'g--', linewidth=2, label=cm1lab)
    ax.plot(np.log(NP[:, ncell]), np.log(TForceSPH[:, ncell]), '<k', linestyle='-', label='Tcpu Force SPH')
    # ---------------------------------
    # m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TForceSPHC[2:]))
    # cm1lab = "TForceSPHC : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax.plot(np.log(NP), m*np.log(NP)+c, 'k--', linewidth=2, label=cm1lab)
    ax.plot(np.log(NP[:, ncell]), np.log(TForceSPHC[:, ncell]), '>k', linestyle='-', label='Tcpu Force SPH C')

    # -----------------------------------
    # m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TPos[2:]))
    # cm1lab = "TPos : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax.plot(np.log(NP), m*np.log(NP)+c, 'g--', linewidth=2, label=cm1lab)
    ax.plot(np.log(NP[:, ncell]), np.log(TPos[:, ncell]), 'sk', linestyle='-', label='Tcpu Position')
    # -----------------------------------
    # m, c, r, p, se1 = stats.linregress(np.log(NP[3:]), np.log(TNeigh[3:]))
    # cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax.plot(np.log(NP), m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
    ax.plot(np.log(NP[:, ncell]), np.log(TNeigh[:, ncell]), 'dk', linestyle='-', label='Tcpu Neighbours')
    # -----------------------------------
    # m, c, r, p, se1 = stats.linregress(np.log(NP[3:]), np.log(TNeighC[3:]))
    # cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax.plot(np.log(NP), m*np.log(NP)+c, 'y--', linewidth=2, label=cm1lab)
    ax.plot(np.log(NP[:, ncell]), np.log(TNeighC[:, ncell]), '^k', linestyle='-', label='Tcpu NeighboursC')
    # -----------------------------------
    ax.plot(np.log(NP[:, ncell]), np.log(TField[:, ncell]), '+k', linestyle='-', label='Tcpu Fields')
    plt.legend(loc='upper left')
    # plt.show()


    # m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForce))
    # cm1lab = "TForce : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax1.loglog(NP, m*np.log(NP)+c, 'b--', linewidth=2, label=cm1lab)
    ax1.loglog(NP[:, ncell], TForce[:, ncell], 'ok', linestyle='-', label='Tcpu Force')
    ax1.plot(NP[:, ncell], TForceVect[:, ncell], '*k', linestyle='-', label='Tcpu Force Vect')
    ax1.plot(NP[:, ncell], TForceSPH[:, ncell], '^k', linestyle='-', label='Tcpu Force SPH')
    ax1.loglog(NP[:, ncell], TPos[:, ncell], 'sk', linestyle='-', label='Tcpu Position')
    # m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TNeigh))
    # cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
    # ax1.loglog(NP, m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
    ax1.loglog(NP[:, ncell], TNeigh[:, ncell], 'dk', linestyle='-', label='Tcpu Neighbours')
    ax1.loglog(NP[:, ncell], TField[:, ncell], '+k', linestyle='-', label='Tcpu Fields')
    plt.legend(loc='upper left')
plt.show()
