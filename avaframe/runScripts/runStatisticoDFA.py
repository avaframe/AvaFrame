"""
Run a combination of the DFA kernel to get an alphaBeta path
to run alphaBeta model to get the alpha angle
and run the DFA kernel again
"""

import pathlib
import time
from configupdater import ConfigUpdater
import matplotlib.pyplot as plt
import numpy as np

# Local imports
# import config and init tools
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in1Data import getInput
from avaframe.ana5Hybrid import hybridTools
import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import geoTrans

# import computation modules
from avaframe.com1DFA import com1DFA
from avaframe.com2AB import com2AB

# import analysis tools
from avaframe.out3Plot import outAB
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools

# import plotting tools
from avaframe.out3Plot import outCom1DFA


# Time the whole routine
startTime = time.time()

# log file name; leave empty to use default runLog.log
logName = 'runSuperSmartANDFastModel'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# ----------------
# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

# Load configuration for hybrid model
hybridModelCfg = pathlib.Path('ana5Hybrid', 'hybridModel_com1DFACfg.ini')
# ----------------
# Run dense flow with coulomb friction
particlesList, fieldsList, tSaveList, dem, plotDict, reportDictList, _ = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=hybridModelCfg)

# postprocess to extract path and energy line# read DEM
demOri = getInput.readDEM(avalancheDir)
headerOri = demOri['header']
xllcOri = headerOri['xllcenter']
yllcOri = headerOri['yllcenter']
avaProfilePart, avaProfileMass, avaProfileKE, _, _ = hybridTools.getCom1DFAPath(particlesList, tSaveList, demOri)
avaProfilePart, avaProfileMass, avaProfileKE = hybridTools.elongateCom1DFAPath(demOri, particlesList[0], avaProfilePart, avaProfileMass, avaProfileKE)
# save profile as AB profile in Inputs

pathAB = pathlib.Path(avalancheDir, 'Inputs', 'LINES', 'pathAB_aimec')
name = 'massAvaPath'
hybridTools.writeLine2SHPfile(avaProfileMass, name, pathAB)

# Run Alpha Beta
cfgAB = cfgUtils.getModuleConfig(com2AB)
# take the path extracted from the DFA model as input
resAB = com2AB.com2ABMain(cfgAB, avalancheDir)


reportDictList = []
_, plotFile, writeFile = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

log.info('Plotted to: %s', plotFile)

# ----------------
# Run dense flow with new coulomb friction parameter
alpha = resAB[name]['alpha']
updater = ConfigUpdater()
updater.read(hybridModelCfg)
updater['GENERAL']['mu'].value = str(np.tan(np.radians(alpha)))
updater.update_file()

particlesList, fieldsList, tSaveList, _, plotDict, reportDictList, _ = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=hybridModelCfg)

avaProfilePartNew, avaProfileMassNew, avaProfileKENew, V2Path, EkinPath = hybridTools.getCom1DFAPath(particlesList, tSaveList, demOri)

# Load all input Parameters from config file
# get the configuration of an already imported module
cfgAIMEC = cfgUtils.getModuleConfig(ana3AIMEC)

initProj.cleanModuleFiles(avalancheDir, ana3AIMEC)

# write configuration to file
cfgUtils.writeCfgFile(avalancheDir, ana3AIMEC, cfgAIMEC)

cfgSetup = cfgAIMEC['AIMECSETUP']
anaMod = cfgSetup['anaMod']

# Setup input from com1DFA
pathDict = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfgAIMEC)

# TODO: define referenceFile
pathDict['numSim'] = len(pathDict['pfd'])

# define reference simulation
pathDict = aimecTools.fetchReferenceSimNo(pathDict, cfgSetup)

pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)

log.info("Running ana3AIMEC model on test case DEM: \n %s \n with profile: \n %s ",
         pathDict['demSource'], pathDict['profileLayer'])
# Run AIMEC postprocessing
rasterTransfo, newRasters, resAnalysis = ana3AIMEC.mainAIMEC(pathDict, cfgAIMEC)

# Plots
sAIMEC = resAnalysis['runout'][0]
xAIMEC = resAnalysis['runout'][1]
yAIMEC = resAnalysis['runout'][2]

ids_alpha = resAB[name]['ids_alpha']
indSplit = resAB[name]['indSplit']
xAB = resAB[name]['x'][ids_alpha]
yAB = resAB[name]['y'][ids_alpha]
sAB = resAB[name]['s'][ids_alpha]
zAB = resAB[name]['z'][ids_alpha]
fAB = resAB[name]['f']
# splitPoint = resAB[name]['splitPoint']

g = 9.81
fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
ax1 = outCom1DFA.addDem2Plot(ax1, dem, what='slope')
ax1 = outCom1DFA.addResult2Plot(ax1, dem, fieldsList[-1]['pta'], 'pta')
ax1.plot(xAIMEC - xllcOri, yAIMEC - yllcOri, 'X', markersize=8, color='0.7', label='Auslauf (A) com1DFA (AIMEC pfd=0m)')
ax1.plot(xAB - xllcOri, yAB - yllcOri, 'x', markersize=8, label='A com2AB')
ax1.plot(avaProfileMass['x'] - xllcOri, avaProfileMass['y'] - yllcOri, 'b-', label='AB path first')
ax1.plot(avaProfileMassNew['x'] - xllcOri, avaProfileMassNew['y'] - yllcOri, 'b-', label='AB path last')

# set titel of output png
title = ('myTitle1')

path = pathlib.Path(avalancheDir, 'Outputs')
pU.saveAndOrPlot({'pathResult': path}, title, fig1)

fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
unit = 'm'
cmap = pU.cmapPlasma
cmap.set_under(color='w')
projectedZ = 'yes'
if projectedZ == 'yes':
    avaProfileMassNew, _ = geoTrans.projectOnRaster(demOri, avaProfileMassNew, interp='bilinear')
Zene = avaProfileMassNew['z'] + V2Path/(2*g)
# Zene = z_m[0] + V2Path/(2*g)
# Colorbar: kinietische Energie [J]
scat = ax2.scatter(avaProfileMassNew['sCor'], Zene, marker='s', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Gesamtenergie(s_mod)')
scat = ax2.scatter(avaProfileMassNew['s'], Zene, marker='o', cmap=cmap, s=2*pU.ms, c=EkinPath, label='Gesamtenergie(s_real)')
cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
cbar2.ax.set_ylabel('kinetische Energie [J]')

ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'], 'b-', label='Z_av(s_real)')
ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'], 'b-.', label='Z_av(s_mod)')
ax2.plot(avaProfileMassNew['s'], avaProfileMassNew['z'], 'k-', label='Z_true(s_real)')
ax2.plot(avaProfileMassNew['sCor'], avaProfileMassNew['z'], 'k-.', label='Z_true(s_mod)')
# ax2.plot(part_aimec[:, 3], z_aimec, 'k--', label='lin. extrapol. Lp_m_projZ')
# ax2.plot(part_aimec[:, 3], part_aimec[:, 2], 'k-.', label='lin. extrapol. Lp_m')
GK = avaProfileMassNew['sCor'][-1] * np.tan(alpha*np.pi/180)
z_enda = avaProfileMassNew['z'][0] - GK
s_geomL = [avaProfileMassNew['sCor'][0], avaProfileMassNew['sCor'][-1]]
z_geomL = [avaProfileMassNew['z'][0], z_enda]
ax2.plot(s_geomL, z_geomL, 'r-', linewidth=0.3, label='alpha line from Z_av func s_mod')
GK = avaProfileMassNew['s'][-1] * np.tan(alpha*np.pi/180)
z_enda = avaProfileMassNew['z'][0] - GK
s_geomL = [avaProfileMassNew['s'][0], avaProfileMassNew['s'][-1]]
z_geomL = [avaProfileMassNew['z'][0], z_enda]
ax2.plot(s_geomL, z_geomL, 'g-', linewidth=0.3, label='alpha line from Z_true func s_mod')
# ax2.plot(sAB, fAB, 'b-', linewidth=0.3, label='Alphalinie com2AB')
ax2.plot(sAB, zAB, 'x', markersize=8, label='A com2AB   = %1.2f' % zAB + 'm')

ax2.set_xlabel('s [m]', fontsize=pU.fs)
ax2.set_ylabel('Höhe [m]', fontsize=pU.fs)
ax2.legend(loc='lower left')

# set titel of output png
title = ('myTitle')

path = pathlib.Path(avalancheDir, 'Outputs')
pU.saveAndOrPlot({'pathResult': path}, title, fig2)
