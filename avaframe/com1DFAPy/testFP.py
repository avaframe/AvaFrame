"""
    Run script for running python DFA kernel
"""
import os
import time
import numpy as np
import math
import matplotlib.pyplot as plt

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA
from avaframe.com1DFAPy import runCom1DFA
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.DFAfunctionsCython as DFAfunC

# from avaframe.DFAkernel.setParam import *
import avaframe.out3Plot.plotUtils as pU
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.out3Plot.makePalette as makePalette

debugPlot = True

# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runFlatPlaneTest'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaFPtest'
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
FPCfg = os.path.join(avalancheDir, 'Inputs', 'FlatPlane_com1DFACfg.ini')
cfg = cfgUtils.getModuleConfig(com1DFA, FPCfg)
cfgGen = cfg['GENERAL']
# Load configuration
cfgFP = cfg['FPSOL']
H0 = float(cfgFP['H0'])
deltaX = float(cfgFP['deltaX'])
slope = float(cfgFP['slope'])
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
r = np.sqrt((X*X)+(Y*Y))
relTh = H0 - (r-deltaX)*slope
relTh = np.where(relTh < 0, 0, relTh)
relTh = np.where(relTh > H0, H0, relTh)

Particles, Fields, Tsave, dem = runCom1DFA.runCom1DFAPy(avaDir=avalancheDir, cfgFile=FPCfg, relTh=relTh, flagAnalysis=False)

# +++++++++POSTPROCESS++++++++++++++++++++++++
# -------------------------------
# Analyse resutls
# tools.plotPosition(particles, dem)
partRef = Particles[0]
Z0 = partRef['z'][0]
rho = cfgGen.getfloat('rho')
gravAcc = cfgGen.getfloat('gravAcc')
mu = cfgGen.getfloat('mu')
sphOption = cfgGen.getint('sphOption')
repeat = False
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
            fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True)
    fig, ax = com1DFA.plotPosition(
            fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True, last=True)
    value = input("[y] to repeat:\n")
    if value != 'y':
        repeat = False


# option for user interaction
if cfgFP.getboolean('flagInteraction'):
    value = input("give time step to plot (float in s):\n")
else:
    value = cfgFP.getfloat('dtSol')
try:
    value = float(value)
except ValueError:
    value = 'n'
while isinstance(value, float):
    ind_t = min(np.searchsorted(Tsave, value), len(Tsave)-1)
    fields = Fields[ind_t]
    particles = Particles[ind_t]
    x = particles['x']
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    m = particles['m']
    force2 = {}
    particles, force2 = DFAfunC.computeForceSPHC(cfgGen, particles, force2, dem, SPHOption=sphOption, gradient=1)
    force3 = {}
    particles, force3 = DFAfunC.computeForceSPHC(cfgGen, particles, force3, dem, SPHOption=4, gradient=1)
    gradNorm = DFAtls.norm(force2['forceSPHX'], force2['forceSPHY'], force2['forceSPHZ'])
    gradNorm1 = DFAtls.norm(force3['forceSPHX'], force3['forceSPHY'], force3['forceSPHZ'])
    # indNorm = np.where(np.abs(gradNorm) > 0.1)
    x1, y1, z1, = DFAtls.normalize(x+xllc, y+yllc, 0)
    uMag = DFAtls.norm(ux, uy, 0)
    v = DFAtls.scalProd(ux, uy, 0, x1, y1, z1)
    grad = DFAtls.scalProd(force2['forceSPHX'], force2['forceSPHY'], force2['forceSPHZ'], x1, y1, z1)
    Grad = np.zeros((nrows, ncols))
    MassBilinear = np.zeros((nrows, ncols))
    MassBilinear = DFAfunC.pointsToRasterC(x, y, m, MassBilinear, csz=5)
    Grad = DFAfunC.pointsToRasterC(x, y, m*gradNorm, Grad, csz=5)
    indMass = np.where(MassBilinear > 0)
    Grad[indMass] = Grad[indMass]/MassBilinear[indMass]
    # Grad = np.where(np.abs(Grad) > 0.01, 1, 0)
    x = particles['x']+xllc
    y = particles['y']+yllc
    m = particles['m']
    Theta = [0, 30, 45, 60]
    Col = ['k', 'b', 'r', 'g']
    d = 2.5

    print(np.max(np.abs(grad)), np.max(gradNorm), np.max(uMag))
    fig3, ax3 = plt.subplots(figsize=(pU.figW, pU.figH))
    com1DFA.plotPosition(fig3, ax3, particles, dem, fields['FD'], pU.cmapDepth, '', plotPart=False, last=False)
    variable = gradNorm
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        pU.cmapDepth, np.amin(variable), np.amax(variable), continuous=True)
    # set range and steps of colormap
    cc = variable
    sc = ax3.scatter(particles['x'], particles['y'], c=cc, cmap=cmap, marker='.')
    pU.addColorBar(sc, ax3, ticks, 'm', 'gradient')

    fig = plt.figure(figsize=(pU.figW, pU.figH))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    # for theta, col in zip(Theta, Col):
    #     theta = theta*math.pi/180
    #     cos = math.cos(theta)
    #     sin = math.sin(theta)
    #     tan = math.tan(theta)
    #     xx = r*cos
    #     yy = r*sin
    #     ind = np.where(((y > tan*x-d/cos) & (y < tan*x+d/cos)))
    #     r = np.sqrt(x*x + y*y)
    #     r = r[ind]
    #     h = particles['h'][ind]
    #     hsph = particles['hSPH'][ind]
    #     ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['FD'][100,:], 'k')
    #     ax1.plot(r, h, color=col, marker='.', linestyle='None')
    #     # ax1.plot(r, hsph, color=col, marker='*', linestyle='None')
    #
    #     ax2.plot(r, grad[ind], color=col, marker='*', linestyle='None')
    #     ax2.plot(r, v[ind], color=col, marker='o', linestyle='None')
    #     # ax2.plot(r, gradNorm[ind], color=col, marker='s', linestyle='None')

    r = np.sqrt(x*x + y*y)
    r = r
    h = particles['h']
    hsph = particles['hSPH']
    ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['FD'][100,:], '--b', label='field flow depth')
    ax1.plot(r, h, color='b', marker='.', linestyle='None', label='particle flow depth')
    # ax1.plot(r, hsph, color=col, marker='*', linestyle='None')

    # ax2.plot(r, grad, color='b', marker='.', linestyle='None')
    ax2.plot(r, gradNorm, color='k', marker='o', linestyle='None', label='SPH gradient used')
    ax2.plot(r, gradNorm1, color='r', marker='.', linestyle='None', label='SPH accurate gradient')
    # ax2.plot(r, v, color='b', marker='.', linestyle='None')

    ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), relTh[100, :], '--k')
    ax1.plot(r, H0-mu*(r-deltaX), '-k', label='initial expected flow depth')
    ax1.set_xlabel('r in [m]')
    ax1.set_title('flow depth, t=%.2f s' % (Tsave[ind_t]))

    ax2.plot(r, mu*np.ones(np.shape(r)), '-k', label='friction threashold')
    ax2.set_xlabel('r in [m]')
    ax2.set_title('Gradient of the flow depth')
    ax1.legend()
    ax2.legend()
    if cfgMain['FLAGS'].getboolean('showPlot'):
        plt.show()

    fig.savefig(os.path.join(outDirTest, 'radialCutSol.%s' % (pU.outputFormat)))

    # option for user interaction
    if cfgFP.getboolean('flagInteraction'):
        value = input("give time step to plot (float in s):\n")
        try:
            value = float(value)
        except ValueError:
            value = 'n'
    else:
        value = 'n'
