"""
    Run script for running the dam break problem on a sloping bed and compare to simulation results
"""

# Load modules
import os
import numpy as np
import matplotlib.pyplot as plt
import glob


# Local imports
import avaframe.com1DFAPy.com1DFA as com1DFAPy
from avaframe.com1DFAPy import runCom1DFA
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.ana1Tests import damBreak
import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runDamBreakProblem'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avaDir = 'data/avaDamBreak'

# Start logging
log = logUtils.initiateLogger(avaDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avaDir)


# Load configuration
damBreakCfg = os.path.join(avaDir, 'Inputs', 'damBreak_com1DFACfg.ini')
cfg = cfgUtils.getModuleConfig(com1DFAPy, damBreakCfg)
cfgGen = cfg['GENERAL']

# Set configuration for outputs from com1DFAPy
resTypeFD = 'FD'
unitFD = pU.cfgPlotUtils['unit%s' % resTypeFD]
nameFD = pU.cfgPlotUtils['name%s' % resTypeFD]
resTypeV = 'V'
unitV = pU.cfgPlotUtils['unit%s' % resTypeV]
nameV = pU.cfgPlotUtils['name%s' % resTypeV]
dt2 = cfg['DAMBREAK'].getfloat('dtStep')

# Load flow depth from analytical solution
hL, hR, uR, phi = damBreak.damBreakSol(avaDir, cfgMain, cfg)
xR = np.linspace(-200, 200, 1000)*np.cos(phi)  # projected on the horizontal plane
tR = int(dt2 * 100.0)

# call com1DFAPy to perform simulation - provide configuration file and release thickness function
Particles, Fields, Tsave = runCom1DFA.runCom1DFAPy(avaDir=avaDir, cfgFile=damBreakCfg, flagAnalysis=True)

# grab results from simulations
inputDir = 'data/avaDamBreak/Outputs/com1DFAPy/peakFiles/'
name1FD = glob.glob(inputDir+os.sep + '*%s*t%d.*.asc' % (resTypeFD, int(dt2)))[0]
data1FD = np.loadtxt(name1FD, skiprows=6)
name1V = glob.glob(inputDir+os.sep + '*%s*t%d.*.asc' % (resTypeV, int(dt2)))[0]
data1V = np.loadtxt(name1V, skiprows=6)

log.info('File for flow depth: %s' % name1FD)
log.info('File for flow velocity: %s' % name1V)

# also load initial state
nameIniFD = 'data/avaDamBreak/Outputs/com1DFAPy/peakFiles/rel1_null_dfa_0.155_FD_t0.00.asc'
dataIniFD = np.loadtxt(nameIniFD, skiprows=6)
nameIniV = 'data/avaDamBreak/Outputs/com1DFAPy/peakFiles/rel1_null_dfa_0.155_V_t0.00.asc'
dataIniV = np.loadtxt(nameIniV, skiprows=6)
header = IOf.readASCheader(name1FD)
cellSize = header.cellsize
ny = data1FD.shape[0]
nx = data1FD.shape[1]
xllc = header.xllcenter
yllc = header.yllcenter

# Location of Profiles
nx_loc = int(ny *0.5)

# set x Vector
x = np.arange(xllc, xllc + nx*cellSize, cellSize)
y = np.zeros(len(x))
y[x<0] = hL
y[x>=0] = 0.0

fig1, ax = plt.subplots(nrows=1, sharex=True)
ax.plot(x, y, 'grey', linestyle='--')
ax.plot(x, dataIniFD[nx_loc, :], 'k--', label='init')
ax.plot(x, data1FD[nx_loc, :], 'b', label='com1DFAPy')
ax.plot(xR, hR[:,tR], 'r-', label='analyt')
ax.set_xlabel('Along track [ncols]')
ax.set_ylabel('Flow depth [m]')
plt.legend()
ax.set_title('%s at time step %.02f s' % (nameFD, dt2))

outDir = os.path.join(avaDir, 'Outputs', 'ana1Tests')
fU.makeADir(outDir)
fig1.savefig(os.path.join(outDir, 'CompareDamBreakH.%s' % (pU.outputFormat)))

y = np.zeros(len(x))
fig2, ax = plt.subplots(nrows=1, sharex=True)
ax.plot(x, y, 'grey', linestyle='--')
ax.plot(x, dataIniV[nx_loc, :], 'k--', label='init')
ax.plot(x, data1V[nx_loc, :], 'b', label='com1DFAPy')
ax.plot(xR, uR[:,tR], 'r-', label='analyt')
ax.set_xlabel('Along track [ncols]')
ax.set_ylabel('Flow velocity [ms-1]')
plt.legend()
ax.set_title('%s at time step %.02f s' % (nameV, dt2))

fig2.savefig(os.path.join(outDir, 'CompareDamBreakVel.%s' % (pU.outputFormat)))

if cfgMain['FLAGS'].getboolean('showPlot'):
    plt.show()
