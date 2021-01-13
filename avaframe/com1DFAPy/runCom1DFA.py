"""
    Run script for running python DFA kernel
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
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.DFAfunctionsCython as DFAfunC

# from avaframe.DFAkernel.setParam import *
import avaframe.out3Plot.plotUtils as pU
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

debugPlot = True

# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'testKernel'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']
# set module name, reqiured as long we are in dev phase
# - because need to create e.g. Output folder for com1DFAPy to distinguish from
# current com1DFA
modName = 'com1DFAPy'

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration
cfg = cfgUtils.getModuleConfig(com1DFA)
cfgGen = cfg['GENERAL']
flagDev = cfg['FLAGS'].getboolean('flagDev')

# for timing the sims
startTime = time.time()


# +++++++++Inputs++++++++++++++++++++++++
# ------------------------
# fetch input data - dem, release-, entrainment- and resistance areas
demFile, relFiles, entFiles, resFile, flagEntRes = gI.getInputData(
    avalancheDir, cfg['FLAGS'], flagDev)
demOri = IOf.readRaster(demFile)
# derive line from release area polygon
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
# process release info to get it as a raster
relRaster = com1DFA.prepareArea(releaseLine, demOri)
relTh = cfgGen.getfloat('relTh')
relRaster = relRaster * relTh

# ------------------------
# initialize simulation : create directories
workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName)
# create particles, create resistance and
# entrainment matrix, initialize fields, get normals and neighbours
particles, fields, Cres, Ment = com1DFA.initializeSimulation(cfgGen, relRaster, dem)

# +++++++++PERFORM SIMULAITON++++++++++++++++++++++
# ------------------------
#  Start time step computation
Tcpu = {}
Tcpu['Force'] = 0.
Tcpu['ForceVect'] = 0.
Tcpu['ForceSPH'] = 0.
Tcpu['Pos'] = 0.
Tcpu['Neigh'] = 0.
Tcpu['Field'] = 0.

Tsave, T, U, Z, S, Particles, Fields, Tcpu = com1DFA.DFAIterate(
    cfgGen, particles, fields, dem, Ment, Cres, Tcpu)

log.info(('cpu time Force = %s s' % (Tcpu['Force'] / Tcpu['nIter'])))
log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / Tcpu['nIter'])))
log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / Tcpu['nIter'])))
log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / Tcpu['nIter'])))
log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / Tcpu['nIter'])))
log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / Tcpu['nIter'])))

tcpuDFA = time.time() - startTime
log.info(('cpu time DFA = %s s' % (tcpuDFA)))


# +++++++++POSTPROCESS++++++++++++++++++++++++
# -------------------------------
# Analyse resutls
# tools.plotPosition(particles, dem)
partRef = Particles[0]
Z0 = partRef['z'][0]
rho = cfgGen.getfloat('rho')
gravAcc = cfgGen.getfloat('gravAcc')
mu = cfgGen.getfloat('mu')
repeat = debugPlot
while repeat:
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    T = np.array([0])
    Z = np.array([0])
    U = np.array([0])
    S = np.array([0])
    for part, field in zip(Particles, Fields):
        T = np.append(T, part['t'])
        S = np.append(S, part['s'][0])
        Z = np.append(Z, part['z'][0])
        U = np.append(U, DFAtls.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
        fig, ax = com1DFA.plotPosition(
            part, demOri, dem['Nz'], pU.cmapDEM2, 'm', fig, ax, plotPart=True)
    fig, ax = com1DFA.plotPosition(
            part, demOri, dem['Nz'], pU.cmapDEM2, 'm', fig, ax, plotPart=True, last=True)
    value = input("[y] to repeat:\n")
    if value != 'y':
        repeat = False
fieldRef = Fields[-1]
fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
fig3, ax3 = plt.subplots(figsize=(pU.figW, pU.figH))
fig1, ax1 = com1DFA.plotPosition(
    particles, demOri, fields['FD'], pU.cmapPres, 'm', fig1, ax1, plotPart=False)
fig2, ax2 = com1DFA.plotPosition(
    particles, demOri, fields['V'], pU.cmapPres, 'm/s', fig2, ax2, plotPart=False)
fig3, ax3 = com1DFA.plotPosition(
    particles, demOri, fields['P']/1000, pU.cmapPres, 'kPa', fig3, ax3, plotPart=False)
fig4, ax4 = plt.subplots(figsize=(pU.figW, pU.figH))
ax4.plot(np.linspace(0, 400, 81), fields['FD'][40,:])
plt.show()


# +++++++++EXPORT RESULTS AND PLOTS++++++++++++++++++++++++
# Result parameters to be exported
resTypesString = cfgGen['resType']
resTypes = resTypesString.split('_')
tSteps = fU.getTimeIndex(cfgGen, Fields)
print(tSteps)
for tStep in tSteps:
    finalFields = Fields[tStep]
    for resType in resTypes:
        resField = finalFields[resType]
        if resType == 'ppr':
            resField = resField * 0.001
        relName = os.path.splitext(os.path.basename(relFiles[0]))[0]
        dataName = relName + '_' + 'null' + '_' + 'dfa' + '_' + '0.155' + '_' + resType + '_'  + 't%.2f' % (Tsave[tStep]) +'.asc'
        # create directory
        outDirPeak = os.path.join(outDir, 'peakFiles')
        fU.makeADir(outDirPeak)
        outFile = os.path.join(outDirPeak, dataName)
        IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
        if tStep == -1:
            log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f - FINAL time step ' % (resType,Tsave[tStep]))
        else:
            log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f ' % (resType,Tsave[tStep]))


# Generata plots for all peakFiles
plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'], modName)


# fig, ax = plt.subplots(figsize=(figW, figH))
# ax.plot(T, U, 'k', linestyle='-', linewidth=2)
# ax.plot(T, Z, 'b', linestyle='-')
# ax.plot(T, np.sqrt(2 * gravAcc * ((Z0-Z) - mu * S)), 'r', linestyle='-', linewidth=1)
# # ax.plot(T, np.sqrt(2 * gravAcc * ((Z0-Z))), 'r', linestyle='-', linewidth=1)
# plt.show()
