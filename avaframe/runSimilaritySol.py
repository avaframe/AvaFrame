"""
    Run script for running python DFA kernel and similarity solution test
"""
import os
import time
import copy
import numpy as np
import math
import matplotlib.pyplot as plt

# Local imports
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA
from avaframe.com1DFAPy import runCom1DFA
import avaframe.com1DFAPy.DFAtools as DFAtls

# from avaframe.DFAkernel.setParam import *
import avaframe.out3Plot.plotUtils as pU
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana1Tests.simiSol as simiSol


# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaSimilaritySol'
# set module name, reqiured as long we are in dev phase
# - because need to create e.g. Output folder for com1DFAPy to distinguish from
# current com1DFA
modName = 'com1DFAPy'

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration
simiSolCfg = os.path.join(avalancheDir, 'Inputs', 'simiSol_com1DFACfg.ini')
cfg = cfgUtils.getModuleConfig(com1DFA, simiSolCfg)
cfgGen = cfg['GENERAL']
# Load configuration
cfgSimi = cfg['SIMISOL']
planeinclinationAngleDeg = float(cfgSimi['planeinclinationAngle'])
# Dimensioning parameters L
L_x = cfgSimi.getfloat('L_x')
L_y = cfgSimi.getfloat('L_y')
Hini = cfgGen.getfloat('relTh')

# set development flag
flagDev = cfg['FLAGS'].getboolean('flagDev')

# for timing the sims
startTime = time.time()

# create output directory for test result plots
outDirTest = os.path.join(avalancheDir, 'Outputs', 'ana1Tests')
fU.makeADir(outDirTest)

# Define release thickness distribution
demFile, relFiles, entFiles, resFile, flagEntRes = gI.getInputData(
    avalancheDir, cfg['FLAGS'], flagDev)
demOri = IOf.readRaster(demFile)
nrows = demOri['header'].nrows
ncols = demOri['header'].ncols
xllc = demOri['header'].xllcenter
yllc = demOri['header'].yllcenter
csz = demOri['header'].cellsize
x = np.linspace(0, ncols-1, ncols)*csz+xllc
y = np.linspace(0, nrows-1, nrows)*csz+yllc
X, Y = np.meshgrid(x, y)
cos = math.cos(math.pi*planeinclinationAngleDeg/180)
sin = math.sin(math.pi*planeinclinationAngleDeg/180)
X1 = X/cos
Y1 = Y
r = np.sqrt((X*X)/(cos*cos)+(Y*Y))
H0 = Hini
relTh = H0 * (1 - (r/L_x) * (r/L_y))
relTh = np.where(relTh < 0, 0, relTh)

# call com1DFAPy to perform simulation - provide configuration file and release thickness function
Particles, Fields, Tsave = runCom1DFA.runCom1DFAPy(avaDir=avalancheDir, cfgFile=simiSolCfg, relTh=relTh, flagAnalysis=False)

# compute similartiy solution
log.info('Computing similarity solution')
solSimi = simiSol.runSimilarity()

# +++++++++POSTPROCESS++++++++++++++++++++++++
# -------------------------------
# Analyse resutls
partRef = Particles[0]
Z0 = partRef['z'][0]
rho = cfgGen.getfloat('rho')
gravAcc = cfgGen.getfloat('gravAcc')
mu = cfgGen.getfloat('mu')
repeat = True

if cfgMain['FLAGS'].getboolean('showPlot'):
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

        # option for user interaction
        if cfgSimi.getboolean('flagInteraction'):
            value = input("[y] to repeat:\n")
            if value != 'y':
                repeat = False
        else:
            repeat = False

# option for user interaction
if cfgSimi.getboolean('flagInteraction'):
    value = input("give time step to plot (float in s):\n")
else:
    value = cfgSimi.getfloat('dtSol')

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
    indy = int(nrows *0.5) -1
    log.info('indy %d' % indy)
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    ax1.set_title('Profile along flow at t=%.2f (com1DFAPy), %.2f s (simiSol)' % (Tsave[ind_t], solSimi['Time'][ind_time]))
    ax2 = ax1.twinx()
    ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['FD'][indy,:], 'k', label='Field flow depth')
    ax2.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['V'][indy,:], 'g', label='Field flow velocity')
    ax1.plot(x, h, '.k', linestyle='None', label='Part flow depth')
    # ax1.plot(x, hsph, '*k', linestyle='None', label='Part flow depth SPH')
    # ax1.plot(x, Ux, '.b', linestyle='None')
    # ax1.plot(x, Uy, '.r', linestyle='None')
    ax2.plot(x, v, '.g', linestyle='None', label='Part flow velocity')
    ax1.plot(X[indy,:], hSimi[indy,:], '--k', label='SimiSol flow depth')
    # ax1.plot(X[51,:], uxSimi[51,:], '--b')
    # ax1.plot(X[51,:], uySimi[51,:], '--r')
    ax2.plot(X[indy,:], vSimi[indy,:], '--g', label='SimiSol flow velocity')
    ax1.axvline(x=xCenter, linestyle=':')
    # ax1.set_xlim([x.min(), x.max()])
    # ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel('x in [m]')
    ax1.set_ylabel('flow depth [m]')
    color = 'tab:green'
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel('flow velocity [ms-1]', color=color)
    ax2.legend(loc='upper right')
    ax1.legend(loc='upper left')

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
    log.info('indc %d' % indc)
    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    ax3 = ax2.twinx()
    ax2.set_title('Profile accross flow at t=%.2f (com1DFAPy), %.2f (simiSol) s' % (Tsave[ind_t], solSimi['Time'][ind_time]))
    ax2.plot(np.linspace(yllc, yllc+(nrows-1)*csz, nrows), fields['FD'][:, indc], 'k', label='Field flow depth')
    ax3.plot(np.linspace(yllc, yllc+(nrows-1)*csz, nrows), fields['V'][:, indc], 'g', label='Field flow velocity')
    ax2.plot(y, h, '.k', linestyle='None', label='Part flow depth')
    # ax2.plot(y, hsph, '*k', linestyle='None', label='Part flow depth SPH')
    # ax2.plot(y, Ux, '.b', linestyle='None')
    # ax2.plot(y, Uy, '.r', linestyle='None')
    ax3.plot(y, v, '.g', linestyle='None',  label='Part flow velocity')
    ax2.plot(Y[:, indc], hSimi[:, indc], '--k', label='SimiSol flow depth')
    # ax2.plot(Y[:, indc], uxSimi[:, indc], '--b')
    # ax2.plot(Y[:, indc], uySimi[:, indc], '--r')
    ax3.plot(Y[:, indc], vSimi[:, indc], '--g', label='SimiSol flow velocity')
    # ax2.set_xlim([x.min(), x.max()])
    # ax2.set_ylim([y.min(), y.max()])
    ax2.set_xlabel('y in [m]')
    ax2.set_ylabel('flow depth [m]')
    color = 'tab:green'
    ax3.tick_params(axis='y', labelcolor=color)
    ax3.set_ylabel('flow velocity [ms-1]', color=color)
    ax3.legend(loc='upper right')
    ax2.legend(loc='upper left')

    if cfgMain['FLAGS'].getboolean('showPlot'):
        plt.show()

    fig1.savefig(os.path.join(outDirTest, 'xCutSol.%s' % (pU.outputFormat)))
    fig2.savefig(os.path.join(outDirTest, 'yCutSol.%s' % (pU.outputFormat)))

    # option for user interaction
    if cfgSimi.getboolean('flagInteraction'):
        value = input("give time step to plot (float in s):\n")
        try:
            value = float(value)
        except ValueError:
            value = 'n'
    else:
        value = 'n'
