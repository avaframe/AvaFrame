import os
import glob
import time
import copy
import logging
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA

# from avaframe.DFAkernel.setParam import *
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# log file name; leave empty to use default runLog.log
logName = 'convergenceDFA'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']
modName = 'com1DFAPy'

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)
# initialize simulation : create directories
workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration
cfg = cfgUtils.getModuleConfig(com1DFA)
cfgGen = cfg['GENERAL']
flagDev = cfg['FLAGS'].getboolean('flagDev')

startTime = time.time()
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
MassPart = [1000, 500, 250]  # , 200]
cfgGen['dt'] = str(0.1)
cfgGen['Tend'] = str(10)
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
    TForce.append(Tcpu['Force'] / Tcpu['nIter'])
    log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / Tcpu['nIter'])))
    TForceVect.append(Tcpu['ForceVect'] / Tcpu['nIter'])
    log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / Tcpu['nIter'])))
    TForceSPH.append(Tcpu['ForceSPH'] / Tcpu['nIter'])
    log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / Tcpu['nIter'])))
    TPos.append(Tcpu['Pos'] / Tcpu['nIter'])
    log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / Tcpu['nIter'])))
    TNeigh.append(Tcpu['Neigh'] / Tcpu['nIter'])
    log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / Tcpu['nIter'])))
    TField.append(Tcpu['Field'] / Tcpu['nIter'])

    # Result parameters to be exported
    resTypesString = cfgGen['resType']
    resTypes = resTypesString.split('_')
    finalFields = Fields[-1]
    for resType in resTypes:
        resField = finalFields[resType]
        if resType == 'ppr':
            resField = resField * 0.001
        relName = os.path.splitext(os.path.basename(relFiles[0]))[0]
        dataName = relName + '_' + 'null' + '_' + 'dfa' + '_' + '0.155' + '_' + resType + \
            '_' + 'massPart' + '_' + str(massPart) + '_' + 'CFL' + '_' + cfgGen['cMax']
        # dataName = relName + '_' + 'massPart' + '_'  + str(massPart) + '_'  + 'CFL' + '_'  + cfg['cMax']
        # fU.writePeakResult(outDir, resField, demOri['header'], dataName)
        log.info('Results parameter: %s has been exported to Outputs/peakFiles' % resType)

    # Generata plots for all peakFiles
    # plotDict = oP.plotAllPeakFields(avalancheDir, cfgFull, cfgMain['FLAGS'], modName)

# -------------------------------
log.info('CPU times')
# log.info(NP)
# log.info(TForce)
# log.info(TForceVect)
# log.info(TForceSPH)
# log.info(TPos)
# log.info(TNeigh)
# log.info(TField)
log.info(('{: <17}' + '{:<10d}'*int(np.size(NP))).format(
    'Mass per Part ', *[np for np in NP]))
log.info(('{: <17}' + '{:<10.4f}'*int(np.size(NP))).format(
    'T force [s] ', *[tf for tf in TForce]))
log.info(('{: <17}' + '{:<10.4f}'*int(np.size(NP))).format(
    'T forceVect [s] ', *[tf for tf in TForceVect]))
log.info(('{: <17}' + '{:<10.4f}'*int(np.size(NP))).format(
    'T forceSPH [s] ', *[tf for tf in TForceSPH]))
log.info(('{: <17}' + '{:<10.4f}'*int(np.size(NP))).format(
    'T Pos [s] ', *[tf for tf in TPos]))
log.info(('{: <17}' + '{:<10.4f}'*int(np.size(NP))).format(
    'Time Neigh [s] ', *[tf for tf in TNeigh]))
log.info(('{: <17}' + '{:<10.4f}'*int(np.size(NP))).format(
    'Time Fields [s] ', *[tf for tf in TField]))

FileName_ext = os.path.join(avalancheDir, 'Outputs', modName, 'CPUTime.txt')
with open(FileName_ext, 'w') as outfile:
    outfile.write(('{: <17}' + '{:<10d}'*int(np.size(NP)) + '\n').format(
        'Mass per Part ', *[np for np in NP]))
    outfile.write(('{: <17}' + '{:<10.4f}'*int(np.size(NP)) + '\n').format(
        'T force [s] ', *[tf for tf in TForce]))
    outfile.write(('{: <17}' + '{:<10.4f}'*int(np.size(NP)) + '\n').format(
        'T forceVect [s] ', *[tf for tf in TForceVect]))
    outfile.write(('{: <17}' + '{:<10.4f}'*int(np.size(NP)) + '\n').format(
        'T forceSPH [s] ', *[tf for tf in TForceSPH]))
    outfile.write(('{: <17}' + '{:<10.4f}'*int(np.size(NP)) + '\n').format(
        'T Pos [s] ', *[tf for tf in TPos]))
    outfile.write(('{: <17}' + '{:<10.4f}'*int(np.size(NP)) + '\n').format(
        'Time Neigh [s] ', *[tf for tf in TNeigh]))
    outfile.write(('{: <17}' + '{:<10.4f}'*int(np.size(NP)) + '\n').format(
        'Time Fields [s] ', *[tf for tf in TField]))
# fig, ax = plt.subplots(figsize=(figW, figH))
# # -------------------------------
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForce))
# cm1lab = "TForce : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax.plot(np.log(NP), m*np.log(NP)+c, 'b--', linewidth=2, label=cm1lab)
# ax.plot(np.log(NP), np.log(TForce), 'ok', linestyle='-', label='Tcpu Force')
# # ---------------------------------
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForceVect))
# cm1lab = "TForceVect : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax.plot(np.log(NP), m*np.log(NP)+c, 'b-.', linewidth=2, label=cm1lab)
# ax.plot(np.log(NP), np.log(TForceVect), '*k', linestyle='-', label='Tcpu Force Vect')
# # ---------------------------------
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForceSPH))
# cm1lab = "TForceSPH : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax.plot(np.log(NP), m*np.log(NP)+c, 'g--', linewidth=2, label=cm1lab)
# ax.plot(np.log(NP), np.log(TForceSPH), '^k', linestyle='-', label='Tcpu Force SPH')
# # -----------------------------------
# # m, c, r, p, se1 = stats.linregress(np.log(NP[2:]), np.log(TPos[2:]))
# # cm1lab = "TPos : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# # ax.plot(np.log(NP), m*np.log(NP)+c, 'g--', linewidth=2, label=cm1lab)
# ax.plot(np.log(NP), np.log(TPos), 'sk', linestyle='-', label='Tcpu Position')
# # -----------------------------------
# m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TNeigh))
# cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# ax.plot(np.log(NP), m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
# # m, c, r, p, se1 = stats.linregress(np.log(NP[4:]), np.log(TNeigh[4:]))
# # cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# # ax.plot(np.log(NP), m*np.log(NP)+c, 'r-.', linewidth=2, label=cm1lab)
# ax.plot(np.log(NP), np.log(TNeigh), 'dk', linestyle='-', label='Tcpu Neighbours')
# # -----------------------------------
# ax.plot(np.log(NP), np.log(TField), '+k', linestyle='-', label='Tcpu Fields')
# plt.legend(loc='upper left')
# # plt.show()
#
# fig1, ax1 = plt.subplots(figsize=(figW, figH))
# # m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TForce))
# # cm1lab = "TForce : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# # ax1.loglog(NP, m*np.log(NP)+c, 'b--', linewidth=2, label=cm1lab)
# ax1.loglog(NP, TForce, 'ok', linestyle='-', label='Tcpu Force')
# ax1.plot(NP, TForceVect, '*k', linestyle='-', label='Tcpu Force Vect')
# ax1.plot(NP, TForceSPH, '^k', linestyle='-', label='Tcpu Force SPH')
# ax1.loglog(NP, TPos, 'sk', linestyle='-', label='Tcpu Position')
# # m, c, r, p, se1 = stats.linregress(np.log(NP), np.log(TNeigh))
# # cm1lab = "TNeigh : $" + ('y=%2.2fx+%2.2f, r^2=%1.2f' % (m, c, r**2)) + "$"
# # ax1.loglog(NP, m*np.log(NP)+c, 'r--', linewidth=2, label=cm1lab)
# ax1.loglog(NP, TNeigh, 'dk', linestyle='-', label='Tcpu Neighbours')
# ax1.loglog(NP, TField, '+k', linestyle='-', label='Tcpu Fields')
# plt.legend(loc='upper left')
# plt.show()
