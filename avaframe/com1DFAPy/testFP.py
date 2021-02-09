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
########## custom release for tests
nrows = demOri['header'].nrows
ncols = demOri['header'].ncols
xllc = demOri['header'].xllcenter
yllc = demOri['header'].yllcenter
csz = demOri['header'].cellsize
x = np.linspace(0, ncols-1, ncols)*csz+xllc
y = np.linspace(0, nrows-1, nrows)*csz+yllc
X, Y = np.meshgrid(x, y)
cos = math.cos(math.pi*0/180)
sin = math.sin(math.pi*0/180)
X1 = X/cos
Y1 = Y
r = np.sqrt((X*X)/(cos*cos)+(Y*Y))
H0 = 5
deltaX = 20
relTh = H0 * (1 - (r-deltaX)/50)
relTh = np.where(relTh < 0, 0, relTh)
relTh = np.where(relTh > H0, H0, relTh)
# relTh = r
# relTh = np.where(relTh < 30, 1, 0)
# relTh1 = -(np.sqrt((X)*(X)+0*(Y)*(Y))*5-100)*0.1
# relTh1 = np.where(relTh1 < 0, 0, relTh1)
# relTh2 = -(np.sqrt(0*(X)*(X)+(Y)*(Y))*5-100)*0.1
# relTh2 = np.where(relTh2 < 0, 0, relTh2)
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
particles, fields, Cres, Ment = com1DFA.initializeSimulation(cfgGen, relRaster, dem)

# +++++++++PERFORM SIMULAITON++++++++++++++++++++++
# ------------------------
#  Start time step computation
Tsave, T, U, Z, S, Particles, Fields, Tcpu = com1DFA.DFAIterate(
    cfgGen, particles, fields, dem, Ment, Cres)
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


value = input("give time step to plot (float in s):\n")
try:
    value = float(value)
except ValueError:
    value = 'n'
while isinstance(value, float):
    ind_t = np.searchsorted(Tsave, value)
    fields = Fields[ind_t]
    particles = Particles[ind_t]
    x = particles['x']
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    m = particles['m']
    force2 = {}
    particles, force2 = DFAfunC.computeForceSPHC(cfgGen, particles, force2, dem, SPHOption=sphOption, gradient=1)
    gradNorm = DFAtls.norm(force2['forceSPHX'], force2['forceSPHY'], force2['forceSPHZ'])
    indNorm = np.where(np.abs(gradNorm) > 0.1)
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
    ax3.plot(x[indNorm], y[indNorm], color='k', marker='.', linestyle='None')
    com1DFA.plotPosition(fig3, ax3, particles, dem, Grad, pU.cmapDepth, '', plotPart=False, last=False)
    fig = plt.figure(figsize=(pU.figW, pU.figH))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    for theta, col in zip(Theta, Col):
        theta = theta*math.pi/180
        cos = math.cos(theta)
        sin = math.sin(theta)
        tan = math.tan(theta)
        xx = r*cos
        yy = r*sin
        ind = np.where(((y > tan*x-d/cos) & (y < tan*x+d/cos)))
        r = np.sqrt(x*x + y*y)
        r = r[ind]
        h = particles['h'][ind]
        hsph = particles['hSPH'][ind]
        ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['FD'][100,:], 'k')
        ax1.plot(r, h, color=col, marker='.', linestyle='None')
        # ax1.plot(r, hsph, color=col, marker='*', linestyle='None')

        ax2.plot(r, grad[ind], color=col, marker='*', linestyle='None')
        ax2.plot(r, v[ind], color=col, marker='o', linestyle='None')
        # ax2.plot(r, gradNorm[ind], color=col, marker='s', linestyle='None')

    ax1.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), relTh[100,:], '--b')
    ax1.plot(r, H0-mu*(r-deltaX), '-k')
    ax1.set_xlabel('r in [m]')
    ax1.set_title('flow depth')

    ax2.plot(r, -mu*np.ones(np.shape(r)), '-k')
    ax2.set_xlabel('r in [m]')
    ax2.set_title('Gradient of the flow depth')
    plt.show()

    value = input("give time step to plot (float in s):\n")
    try:
        value = float(value)
    except ValueError:
        value = 'n'



fieldRef = Fields[-1]
fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
fig3, ax3 = plt.subplots(figsize=(pU.figW, pU.figH))
fig1, ax1 = com1DFA.plotPosition(
    fig1, ax1, particles, demOri, fields['FD'], pU.cmapPres, 'm', plotPart=False)
fig2, ax2 = com1DFA.plotPosition(
    fig2, ax2, particles, demOri, fields['V'], pU.cmapPres, 'm/s', plotPart=False)
fig3, ax3 = com1DFA.plotPosition(
    fig3, ax3, particles, demOri, fields['P']/1000, pU.cmapPres, 'kPa', plotPart=False)

# fig4, ax4 = plt.subplots(figsize=(pU.figW, pU.figH))
# ax4.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols-1), (fields['FD'][101, 1:ncols]-fields['FD'][101, 0:(ncols-1)])/5*math.cos(math.pi*3/180))
# force2 = {}
# particles = Particles[-1]
# particles, force2 = DFAfunC.computeForceSPHC(cfgGen, particles, force2, dem, SPHOption=2, gradient=1)
# grad = DFAtls.norm(force2['forceSPHX'], force2['forceSPHY'], force2['forceSPHZ'])
# # grad = DFAtls.norm(force['forceSPHX'], force['forceSPHY'], force['forceSPHZ'])
# x = particles['x']
# y = particles['y']
# m = particles['m']
# ind = np.where(((particles['y']+yllc > -2.5) & (particles['y']+yllc < 2.5)))
# Grad = np.zeros((nrows, ncols))
# MassBilinear = np.zeros((nrows, ncols))
# MassBilinear = DFAfunC.pointsToRasterC(x, y, m, MassBilinear, csz=5)
# Grad = DFAfunC.pointsToRasterC(x, y, m*grad, Grad, csz=5)
# indMass = np.where(MassBilinear > 0)
# Grad[indMass] = Grad[indMass]/MassBilinear[indMass]
# fig5, ax5 = plt.subplots(figsize=(pU.figW, pU.figH))
# fig5, ax5 = com1DFA.plotPosition(particles, dem, Grad, pU.cmapPres, '', fig5, ax5, plotPart=False)
# fig6, ax6 = plt.subplots(figsize=(pU.figW, pU.figH))
# ax6.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), Grad[101,:], 'b')
# ax6.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['FD'][101,:], 'r')
# ax6.plot(np.linspace(xllc, xllc+(ncols-1)*csz, ncols), fields['V'][101,:], 'g')
# ax6.plot(particles['x'][ind]+xllc, grad[ind], '.r', linestyle='None')
# #     plt.show()
# plt.show()


# # +++++++++EXPORT RESULTS AND PLOTS++++++++++++++++++++++++
# # Result parameters to be exported
# resTypesString = cfgGen['resType']
# resTypes = resTypesString.split('_')
# tSteps = fU.getTimeIndex(cfgGen, Fields)
# for tStep in tSteps:
#     finalFields = Fields[tStep]
#     for resType in resTypes:
#         resField = finalFields[resType]
#         if resType == 'ppr':
#             resField = resField * 0.001
#         relName = os.path.splitext(os.path.basename(relFiles[0]))[0]
#         dataName = relName + '_' + 'null' + '_' + 'dfa' + '_' + '0.155' + '_' + resType + '_'  + 't%.2f' % (Tsave[tStep]) +'.asc'
#         # create directory
#         outDirPeak = os.path.join(outDir, 'peakFiles')
#         fU.makeADir(outDirPeak)
#         outFile = os.path.join(outDirPeak, dataName)
#         IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
#         if tStep == -1:
#             log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f - FINAL time step ' % (resType,Tsave[tStep]))
#         else:
#             log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f ' % (resType,Tsave[tStep]))
#
#
# # Generata plots for all peakFiles
# plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'], modName)
