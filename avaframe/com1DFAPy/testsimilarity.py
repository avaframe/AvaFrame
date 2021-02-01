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
import avaframe.in3Utils.geoTrans as geoTrans
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
import avaframe.out3Plot.makePalette as makePalette
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFAPy.simiSol as simiSol


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
cfgGenDict = cfgUtils.parser2dict(cfgGen)

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
########## custom release for tests
nrows = demOri['header'].nrows
ncols = demOri['header'].ncols
xllc = demOri['header'].xllcenter
yllc = demOri['header'].yllcenter
csz = demOri['header'].cellsize
x = np.linspace(0, ncols-1, ncols)*csz+xllc
y = np.linspace(0, nrows-1, nrows)*csz+yllc
X, Y = np.meshgrid(x, y)
cos = math.cos(math.pi*35/180)
sin = math.sin(math.pi*35/180)
X1 = X/cos
Y1 = Y
r = np.sqrt((X*X)/(cos*cos)+(Y*Y))
H0 = 8
relTh = H0 * (1 - (r/80) * (r/80))
relTh = np.where(relTh < 0, 0, relTh)
# relTh1 = -(np.sqrt((X-100)*(X-100)+0*(Y-100)*(Y-100))*5-102.5)*0.04765
# relTh2 = -(np.sqrt(0*(X-100)*(X-100)+(Y-100)*(Y-100))*5-102.5)*0.04
# relTh = np.minimum(relTh1, relTh2)
# could do something more advanced if we want varying release depth
########################
# relTh = 1
relRaster = relRaster * relTh

# ------------------------
# initialize simulation : create directories
workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName)
# create particles, create resistance and
# entrainment matrix, initialize fields, get normals and neighbours
particles, fields, Cres, Ment = com1DFA.initializeSimulation(cfgGenDict, relRaster, dem)

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
    cfgGenDict, particles, fields, dem, Ment, Cres, Tcpu)

log.info(('cpu time Force = %s s' % (Tcpu['Force'] / Tcpu['nIter'])))
log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / Tcpu['nIter'])))
log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / Tcpu['nIter'])))
log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / Tcpu['nIter'])))
log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / Tcpu['nIter'])))
log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / Tcpu['nIter'])))

tcpuDFA = time.time() - startTime
log.info(('cpu time DFA = %s s' % (tcpuDFA)))

log.info('Computing similarity solution')
solSimi = simiSol.runSimilarity()

# +++++++++POSTPROCESS++++++++++++++++++++++++
# -------------------------------
# Analyse resutls
# tools.plotPosition(particles, dem)
partRef = Particles[0]
Z0 = partRef['z'][0]
rho = cfgGen.getfloat('rho')
gravAcc = cfgGen.getfloat('gravAcc')
mu = cfgGen.getfloat('mu')
repeat = True
while repeat:
    fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
    for part, field in zip(Particles, Fields):
        t = part['t']
        ind_time = np.searchsorted(solSimi['Time'], t)
        hSimi = simiSol.h(solSimi, X1, Y1, ind_time)
        hSimi = np.where(hSimi <= 0, 0, hSimi)
        fig, ax, cmap, lev = com1DFA.plotContours(
            fig, ax, part, demOri, field['FD'], pU.cmapDepth, 'm')
        CS = ax.contour(X, Y, hSimi, levels=lev, origin='lower', cmap=cmap,
                        linewidths=2, linestyles='dashed')
        plt.pause(1)

    fig, ax, cmap, lev = com1DFA.plotContours(
        fig, ax, part, demOri, field['FD'], pU.cmapDepth, 'm', last=True)
    CS = ax.contour(X, Y, hSimi, levels=lev, origin='lower', cmap=cmap,
                    linewidths=2, linestyles='dashed')
    ax.clabel(CS, inline=1, fontsize=8)
    value = input("[y] to repeat:\n")
    if value != 'y':
        repeat = False

value = input("give time step to plot (float in s):\n")
try:
    value = float(value)
except ValueError:
    value = 'n'
while isinstance(value, float):
    ind_t = np.searchsorted(Tsave, value)
    fields = Fields[ind_t]
    particles = Particles[ind_t]
    ind_time = np.searchsorted(solSimi['Time'], value)
    # get simi sol
    hSimi = simiSol.h(solSimi, X1, Y1, ind_time)
    hSimi = np.where(hSimi <= 0, 0, hSimi)
    print(np.sum(hSimi))
    uxSimi = simiSol.u(solSimi, X1, Y1, ind_time)
    uxSimi = np.where(hSimi <= 0, 0, uxSimi)
    uySimi = simiSol.v(solSimi, X1, Y1, ind_time)
    uySimi = np.where(hSimi <= 0, 0, uySimi)
    vSimi = np.sqrt(uxSimi*uxSimi + uySimi*uySimi)
    xCenter = simiSol.xc(solSimi, X1, Y1, ind_time)*cos
    x = particles['x']
    y = particles['y']
    m = particles['m']
    ind = np.where(((particles['y']+yllc > -2.5) & (particles['y']+yllc < 2.5)))
    x = particles['x'][ind]+xllc
    x1 = x/cos
    h = particles['h'][ind]
    hsph = particles['hSPH'][ind]
    ux = particles['ux'][ind]
    uy = particles['uy'][ind]
    uz = particles['uz'][ind]
    Ux = DFAtls.scalProd(ux, uy, uz, cos, 0, -sin)
    Uy = DFAtls.scalProd(ux, uy, uz, 0, 1, 0)
    v = np.sqrt(ux*ux + uy*uy + uz*uz)
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    ax1.set_title('X cut of the solution at t=%.2f, %.2f s' % (Tsave[ind_t], solSimi['Time'][ind_time]))
    ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['FD'][50,:], 'k')
    ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['V'][50,:], 'g')
    ax1.plot(x, h, '.k', linestyle='None')
    ax1.plot(x, hsph, '*k', linestyle='None')
    # ax1.plot(x, Ux, '.b', linestyle='None')
    # ax1.plot(x, Uy, '.r', linestyle='None')
    ax1.plot(x, v, '.g', linestyle='None')
    ax1.plot(X[51,:], hSimi[50,:], '--k')
    # ax1.plot(X[51,:], uxSimi[51,:], '--b')
    # ax1.plot(X[51,:], uySimi[51,:], '--r')
    ax1.plot(X[51,:], vSimi[50,:], '--g')
    ax1.axvline(x=xCenter, linestyle=':')
    # ax1.set_xlim([x.min(), x.max()])
    # ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel('x in [m]')
    # ax1.set_ylabel('x in [m]')

    ind = np.where(((particles['x']+xllc > xCenter-2.5) & (particles['x']+xllc < xCenter+2.5)))
    x = particles['x'][ind]+xllc
    y = particles['y'][ind]+yllc
    h = particles['h'][ind]
    hsph = particles['hSPH'][ind]
    ux = particles['ux'][ind]
    uy = particles['uy'][ind]
    uz = particles['uz'][ind]
    Ux = DFAtls.scalProd(ux, uy, uz, cos, 0, -sin)
    Uy = DFAtls.scalProd(ux, uy, uz, 0, 1, 0)
    v = np.sqrt(ux*ux + uy*uy + uz*uz)
    indc = int(np.floor((xCenter - xllc)/csz))
    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    ax2.set_title('Y cut of the solution at t=%.2f, %.2f s' % (T[ind_t], solSimi['Time'][ind_time]))
    ax2.plot(np.linspace(yllc, yllc+(nrows-1)*csz, nrows), fields['FD'][:, indc], 'k')
    ax2.plot(np.linspace(yllc, yllc+(nrows-1)*csz, nrows), fields['V'][:, indc], 'g')
    ax2.plot(y, h, '.k', linestyle='None')
    ax2.plot(y, hsph, '*k', linestyle='None')
    # ax2.plot(y, Ux, '.b', linestyle='None')
    # ax2.plot(y, Uy, '.r', linestyle='None')
    ax2.plot(y, v, '.g', linestyle='None')
    ax2.plot(Y[:, indc], hSimi[:, indc], '--k')
    # ax2.plot(Y[:, indc], uxSimi[:, indc], '--b')
    # ax2.plot(Y[:, indc], uySimi[:, indc], '--r')
    ax2.plot(Y[:, indc], vSimi[:, indc], '--g')
    # ax2.set_xlim([x.min(), x.max()])
    # ax2.set_ylim([y.min(), y.max()])
    ax2.set_xlabel('y in [m]')
    # ax2.set_ylabel('x in [m]')
    plt.show()

    value = input("give time step to plot (float in s):\n")
    try:
        value = float(value)
    except ValueError:
        value = 'n'
