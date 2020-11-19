"""

    This file is part of Avaframe.
"""


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
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.com1DFAPy.DFAtools as DFAtls

# from avaframe.DFAkernel.setParam import *
from avaframe.out3Plot.plotUtils import *
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'testKernel'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']
modName = 'com1DFAPy'

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

cfg = cfgUtils.getModuleConfig(com1DFA)['GENERAL']
cfgFull = cfgUtils.getModuleConfig(com1DFA)

startTime = time.time()
# ------------------------
# fetch input data
demFile, relFiles, entFiles, resFile, flagEntRes = gI.getInputData(avalancheDir, cfgFull['FLAGS'])
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

massPart = 1000  # [200, 100, 50, 25, 10, 7.5, 5]
cfg['massPerPart'] = str(massPart)
# ------------------------
# initialize simulation : create directories
workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName)
# create particles, create resistance and
# entrainment matrix, initialize fields, get normals and neighbours
particles, fields, Cres, Ment = com1DFA.initializeSimulation(cfg, relRaster, dem)

# ------------------------
#  Start time step computation
Tcpu = {}
Tcpu['Force'] = 0.
Tcpu['ForceVect'] = 0.
Tcpu['ForceSPH'] = 0.
Tcpu['Pos'] = 0.
Tcpu['Neigh'] = 0.
Tcpu['Field'] = 0.

T, U, Z, S, Particles, Fields, Tcpu = com1DFA.DFAIterate(cfg, particles, fields, dem, Ment, Cres, Tcpu)

log.info(('cpu time Force = %s s' % (Tcpu['Force'] / Tcpu['nIter'])))
log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / Tcpu['nIter'])))
log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / Tcpu['nIter'])))
log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / Tcpu['nIter'])))
log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / Tcpu['nIter'])))
log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / Tcpu['nIter'])))

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
        U = np.append(U, DFAtls.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
        # print(part['t'])
        # print(DFAtls.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
        # exact solution for inclined plane with no friction
        # print(gravAcc * math.sin(34*math.pi/180) * part['t'])
        # exact solution with no friction
        # print(math.sqrt(2 * gravAcc * abs(partRef['z'][0] - part['z'][0])))

        # exact solution for inclined plane with friction
        # print(gravAcc * math.cos(34*math.pi/180) * (math.tan(34*math.pi/180) - mu) * part['t'])
        # exact solution with friction
        # print(math.sqrt(2 * gravAcc * ((partRef['z'][0] - part['z'][0]) - mu * part['s'][0])))
        #
        fig, ax = com1DFA.plotPosition(part, demOri, field['pfd'], cmapPres, 'm', fig, ax, plotPart=True)
        # fig1, ax1 = com1DFA.plotPosition(part, dem, dem['rasterData'], fig1, ax1)
    # plt.show()
    # repeat = False
    value = input("[y] to repeat:\n")
    if value != 'y':
        repeat = False
fieldRef = Fields[-1]
fig1, ax1 = plt.subplots(figsize=(figW, figH))
fig2, ax2 = plt.subplots(figsize=(figW, figH))
fig3, ax3 = plt.subplots(figsize=(figW, figH))
fig1, ax1 = com1DFA.plotPosition(particles, demOri, fields['pfd'], cmapPres, 'm', fig1, ax1, plotPart=False)
fig2, ax2 = com1DFA.plotPosition(particles, demOri, fields['pv'], cmapPres, 'm/s', fig2, ax2, plotPart=False)
fig3, ax3 = com1DFA.plotPosition(particles, demOri, fields['ppr']/1000, cmapPres, 'kPa', fig3, ax3, plotPart=False)
plt.show()


# Result parameters to be exported
resTypesString = cfg['resType']
resTypes = resTypesString.split('_')
finalFields = Fields[-1]
for resType in resTypes:
    resField = finalFields[resType]
    if resType == 'ppr':
        resField = resField * 0.001
    relName = os.path.splitext(os.path.basename(relFiles[0]))[0]
    dataName = relName + '_' + 'null' + '_'  + 'dfa' + '_'  + '0.155' + '_'  + resType
    fU.writePeakResult(outDir, resField, demOri['header'], dataName)
    log.info('Results parameter: %s has been exported to Outputs/peakFiles' % resType)

# Generata plots for all peakFiles
plotDict = oP.plotAllPeakFields(avalancheDir, cfgFull, cfgMain['FLAGS'], modName)


# fig, ax = plt.subplots(figsize=(figW, figH))
# ax.plot(T, U, 'k', linestyle='-', linewidth=2)
# ax.plot(T, Z, 'b', linestyle='-')
# ax.plot(T, np.sqrt(2 * gravAcc * ((Z0-Z) - mu * S)), 'r', linestyle='-', linewidth=1)
# # ax.plot(T, np.sqrt(2 * gravAcc * ((Z0-Z))), 'r', linestyle='-', linewidth=1)
# plt.show()
