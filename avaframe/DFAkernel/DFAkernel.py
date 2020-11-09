import os
import glob
import time
import copy
import logging
import numpy as np
import math
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib as mpl

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
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

startTime = time.time()
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

massPart = 1250  # [200, 100, 50, 25, 10, 7.5, 5]
cfg['massPerPart'] = str(massPart)
# ------------------------
# initialize simulation : create particles, create resistance and
# entrainment matrix, initialize fields, get normals and neighbours
dem, particles, fields, Cres, Ment = DFAtools.initializeSimulation(cfg, relRaster, dem)

# ------------------------
#  Start time step computation
Tcpu = {}
Tcpu['Force'] = 0.
Tcpu['ForceVect'] = 0.
Tcpu['ForceSPH'] = 0.
Tcpu['Pos'] = 0.
Tcpu['Neigh'] = 0.
Tcpu['Field'] = 0.

T, U, Z, S, Particles, Fields, Tcpu = DFAtools.DFAIterate(cfg, particles, fields, dem, Ment, Cres, Tcpu)

log.info(('cpu time Force = %s s' % (Tcpu['Force'] / Tcpu['niter'])))
log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / Tcpu['niter'])))
log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / Tcpu['niter'])))
log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / Tcpu['niter'])))
log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / Tcpu['niter'])))
log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / Tcpu['niter'])))

tcpuDFA = time.time() - startTime
log.info(('cpu time DFA = %s s' % (tcpuDFA)))
# -------------------------------
# Analyse resutls
# tools.plotPosition(particles, dem)
partRef = Particles[0]
Z0 = partRef['z'][0]
rho = cfg.getfloat('rho')
gravAcc = cfg.getfloat('gravAcc')
mu = cfg.getfloat('mu')
repeat = True
while repeat == True:
    fig, ax = plt.subplots(figsize=(figW, figH))
    T = np.array([0])
    Z = np.array([0])
    U = np.array([0])
    S = np.array([0])
    for part, field in zip(Particles, Fields):
        T = np.append(T, part['t'])
        S = np.append(S, part['s'][0])
        Z = np.append(Z, part['z'][0])
        U = np.append(U, DFAtools.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
        # print(part['t'])
        # print(DFAtools.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
        # exact solution for inclined plane with no friction
        # print(gravAcc * math.sin(34*math.pi/180) * part['t'])
        # exact solution with no friction
        # print(math.sqrt(2 * gravAcc * abs(partRef['z'][0] - part['z'][0])))

        # exact solution for inclined plane with friction
        # print(gravAcc * math.cos(34*math.pi/180) * (math.tan(34*math.pi/180) - mu) * part['t'])
        # exact solution with friction
        # print(math.sqrt(2 * gravAcc * ((partRef['z'][0] - part['z'][0]) - mu * part['s'][0])))
        #
        fig, ax = DFAtools.plotPosition(part, demOri, field['PFD'], cmapPres, fig, ax, plotPart=True)
        # fig1, ax1 = DFAtools.plotPosition(part, dem, dem['rasterData'], fig1, ax1)
    # plt.show()
    # repeat = False
    value = input("[y] to repeat:\n")
    if value != 'y':
        repeat = False
fieldRef = Fields[-1]
fig1, ax1 = plt.subplots(figsize=(figW, figH))
fig2, ax2 = plt.subplots(figsize=(figW, figH))
fig2, ax2 = DFAtools.plotPosition(particles, demOri, fields['FD'], cmapPres, fig2, ax2, plotPart=False)
fig1, ax1 = DFAtools.plotPosition(particles, demOri, fields['PFD'], cmapPres, fig1, ax1, plotPart=False)
plt.show()


# fig, ax = plt.subplots(figsize=(figW, figH))
# ax.plot(T, U, 'k', linestyle='-', linewidth=2)
# ax.plot(T, Z, 'b', linestyle='-')
# ax.plot(T, np.sqrt(2 * gravAcc * ((Z0-Z) - mu * S)), 'r', linestyle='-', linewidth=1)
# # ax.plot(T, np.sqrt(2 * gravAcc * ((Z0-Z))), 'r', linestyle='-', linewidth=1)
# plt.show()
