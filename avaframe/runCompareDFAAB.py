"""
    trackParticle
    get information from Outputs/com1DFA/particles/particlesxxxx.xxxx.p
    create paths as .shp
    create plots, compare runout
"""

import math
import shapefile
import time
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from configupdater import ConfigUpdater
import pathlib
import numpy.ma as ma
import shutil
from time import sleep

# Local imports
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in2Trans.ascUtils as IOf
import avaframe.com1DFA.DFAtools as DFAtls
from avaframe.out3Plot.plotUtils import *
from avaframe.com2AB import com2AB
from avaframe.out3Plot import outAB
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
from avaframe.in3Utils import initializeProject as iP
import avaframe.in3Utils.geoTrans as geoTrans
from avaframe.in1Data import getInput
import avaframe.out3Plot.plotUtils as pU

# changes local_....ini files
# INPUT
# path to AvaFrame
pathAvaFrame = '/home/marie/ava0/AvaFrame'
# name of avalanche
avaName = 'avaSei'
avalancheDir = 'data/avaSei'
logName = 'runCompareDFAAB'
# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)
# alpha [°] (avaParabola: 15, avaHelix: 30, avaSei: 26.5)
alpha = 26.5
# one Particle simulation: oneParticle = 'yes'
#  -> in local_com1DFA: deltaTh >> relTh
oneParticle = 'ys'
# clear Inputs/LINES pathelongation method changed
# pathelongation via linear extrapolation (betapoint): pathelongation = 'linearextrapolationbeta'
pathelongation = 'linearextrapolationbeta'
# pathelongation via small alpha/Reibungswinkel: pathelongation = 'smallalpha'
# Delete Output directory after first com1DFA computation!
# pathelongation = 'smallalpha'
# z-coordinates projected on the dem: projectedZ = 'yes'
projectedZ = 'yes'
# diffrence between modified s (s_mod) and real s (s_real) with linearextrapolationbeta, projectedZ
sreal = 'yes'
#########################################################


def writeLine2SHPfile(part, lineName, fileName):
    w = shapefile.Writer(fileName)
    w.field('name', 'C')
    w.line([part])
    w.record(lineName)
    w.close()


# small alpha for elongated path
if avaName == 'avaSei':
    kl = 0.90
else:
    kl = 0.80
sa = kl * alpha
print('Kleines Alpha = sa = %1.1f' % sa + '° = %1.1f' % kl + ' * alpha, alpha = %1.1f' % alpha + '°')
print('oneParticle = ' + oneParticle + ', pathelongation = ' + pathelongation + ', projectedZ = ' + projectedZ + ', reals = ' + sreal)

# update local_....ini files
avalancheDir = os.path.join('data', avaName)

k4 = str(alpha)

file_path_avaframe = os.path.join(pathAvaFrame, 'avaframe', 'local_avaframeCfg.ini')
updater = ConfigUpdater()
updater.read(file_path_avaframe)
updater['MAIN']['avalancheDir'].value = avalancheDir
updater.update_file()

file_path_AB = os.path.join(pathAvaFrame, 'avaframe', 'com2AB', 'local_com2ABCfg.ini')
updater = ConfigUpdater()
updater.read(file_path_AB)
updater['ABSETUP']['k4'].value = k4
updater.update_file()

# calculate mu_small=tan(sa)
mu_smalla_exact = np.tan(sa*np.pi/180)
mu_smalla_float = round(mu_smalla_exact, 5)
mu_smalla = str(mu_smalla_float)

file_path_DFA = os.path.join(pathAvaFrame, 'avaframe', 'com1DFA', 'local_com1DFACfg.ini')
updater = ConfigUpdater()
updater.read(file_path_DFA)
updater['GENERAL']['mu'].value = mu_smalla
updater.update_file()

# in case local files were changed
# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']
avalancheTitle = avalancheDir[5:]

cfgCom2AB = cfgUtils.getModuleConfig(com2AB)
alpha_str = cfgCom2AB['ABSETUP']['k4']
alpha = float(alpha_str)

cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
mu_smalla_str = cfgCom1DFA['GENERAL']['mu']
mu_smalla = float(mu_smalla_str)

avaDir = os.path.join(pathAvaFrame, 'avaframe', avalancheDir)

for file in os.listdir(os.path.join(avaDir, 'Inputs')):
    if file.endswith('.asc'):
        nameDEM = file

inDir = os.path.join(avaDir)
inDir = pathlib.Path(inDir, 'Outputs', 'com1DFA', 'particles')
pathAB_smalla = os.path.join(avaDir, 'Inputs', 'LINES', 'pathAB_aimec_smalla')
fileName_smalla = os.path.join(avaDir, 'Inputs', 'LINES', 'pathAB_aimec_smalla')

# read DEM
dem = getInput.readDEM(avalancheDir)

# get DEM Path
header = dem['header']
ncols = header['ncols']
nrows = header['nrows']
xllc = header['xllcenter']
yllc = header['yllcenter']
csz = header['cellsize']
relField = np.zeros((nrows, ncols))
# z = dem['rasterData']

if oneParticle == 'yes':
    if avaName == 'avaParabola':
        relField[round((-4000-yllc)/csz), round((2700-xllc)/csz)] = 1
    if avaName == 'avaHelix':
        relField[round((-4809-yllc)/csz), round((2263-xllc)/csz)] = 1
    if avaName == 'avaSei':
        relField[round((380955-yllc)/csz), round((252638-xllc)/csz)] = 1
else:
    relField = ''

"""elongate path with small alpha
"""
if pathelongation == 'smallalpha':

    particlesList, fieldsList, tSave, dem, plotDict, reportDictList = com1DFA.com1DFAMain(
        avalancheDir, cfgMain, cfgFile='', relThField=relField, variationDict='')

    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    XX = PointsX[0, :]
    YY = PointsY[:, 0]
    ZZ = dem['rasterData']
    Particles, TimeStepInfo = com1DFA.readPartFromPickle(inDir, flagAvaDir=False)

    part_smalla = np.empty((0, 4))
    count = 0
    for t in TimeStepInfo[0:]:
        particles = Particles[count]
        m = particles['m']
        x = particles['x'] + xllc
        y = particles['y'] + yllc
        z = particles['z']
        s = particles['s']
        Npart = particles['Npart']
        pond = m
        pondSum = np.sum(m)
        # mass-averaged path
        temp_smalla = np.array([[np.sum(m*x)/np.sum(m), np.sum(m*y)/np.sum(m), np.sum(m*z)/np.sum(m), np.sum(m*s)/np.sum(m)]])
        part_smalla = np.append(part_smalla, temp_smalla, axis=0)
        count = count + 1

    lineName_smalla = 'myLine_m_smalla'
    if projectedZ == 'yes':
        z_proj_tuble = geoTrans.projectOnGrid(part_smalla[:, 0], part_smalla[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
        z_proj = np.array(z_proj_tuble[0])
        z_proj.flatten()
        projZ_smalla = part_smalla
        projZ_smalla[:, 2] = z_proj
        writeLine2SHPfile(projZ_smalla, lineName_smalla, fileName_smalla)
    else:
        writeLine2SHPfile(part_smalla, lineName_smalla, fileName_smalla)

    # delete Outputs (not working automaticly so far)
    # Outputs = os.path.join(avaDir, 'Outputs')
    # shutil.rmtree('Outputs', ignore_errors=True)
    print('5 seconds to delete Outputs directory')
    sleep(5)
    print('timeover')

"""paths
"""

# calculate mu=tan(alpha)
mu_exact = np.tan(alpha*np.pi/180)
mu_float = round(mu_exact, 5)
mu = str(mu_float)

file_path_DFA = os.path.join(pathAvaFrame, 'avaframe', 'com1DFA', 'local_com1DFACfg.ini')
updater = ConfigUpdater()
updater.read(file_path_DFA)
updater['GENERAL']['mu'].value = mu
updater.update_file()

cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
mu_str = cfgCom1DFA['GENERAL']['mu']
mu = float(mu_str)

pathAB_px = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_px')
pathAB_p = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_p')
pathAB_m = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_m')
pathAB_kE = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_kE')
fileName_px = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_px')
fileName_p = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_p')
fileName_m = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_m')
fileName_kE = os.path.join(avaDir, 'Outputs', 'pathAB_aimec_kE')
fileName_long_m = os.path.join(avaDir, 'Inputs', 'LINES', 'pathAB_aimec_m')

particlesList, fieldsList, tSave, dem, plotDict, reportDictList = com1DFA.com1DFAMain(
    avalancheDir, cfgMain, cfgFile='', relThField=relField, variationDict='')

xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
PointsX, PointsY = np.meshgrid(xgrid, ygrid)
XX = PointsX[0, :]
YY = PointsY[:, 0]
ZZ = dem['rasterData']
Particles, TimeStepInfo = com1DFA.readPartFromPickle(inDir, flagAvaDir=False)

V2Path = np.empty((0, 1))
EkinPath = np.empty((0, 1))
EpotPath = np.empty((0, 1))
ax = plt.subplot(111)
part_px = np.empty((0, 4))
part_p = np.empty((0, 4))
part_m = np.empty((0, 4))
part_kE = np.empty((0, 4))
count = 0
for t in TimeStepInfo[0:]:
    particles = Particles[count]
    m = particles['m']
    x = particles['x'] + xllc
    y = particles['y'] + yllc
    z = particles['z']
    ux = particles['ux']
    uy = particles['uy']
    uz = particles['uz']
    s = particles['s']
    u = DFAtls.norm(ux, uy, uz)
    U2 = u*u
    Npart = particles['Npart']
    kineticEne = 0.5*m*u*u
    kineticEneSum = np.sum(kineticEne)
    # mass-averaged path
    if kineticEneSum <= 100:
        pond = np.ones(np.shape(kineticEne))
        pondSum = Npart
    else:
        pond = kineticEne
        pondSum = kineticEneSum

    pond = m
    pondSum = np.sum(m)
    v2coE = np.sum(pond*U2)/pondSum
    V2Path = np.append(V2Path, v2coE)

    # update energy
    EkinPath = np.append(EkinPath, kineticEneSum)
    # kinetic energy-averaged path
    if kineticEneSum > 0:
        temp_kE = np.array([[np.sum(kineticEne*x)/kineticEneSum, np.sum(kineticEne*y)/kineticEneSum,
                np.sum(kineticEne*z)/kineticEneSum, np.sum(kineticEne*s)/kineticEneSum]])
    else:
        temp_kE = np.array([[np.sum(m*x)/np.sum(m), np.sum(m*y)/np.sum(m),
                np.sum(m*z)/np.sum(m), np.sum(m*s)/np.sum(m)]])
    # mass-averaged path
    temp_m = np.array([[np.sum(m*x)/np.sum(m), np.sum(m*y)/np.sum(m),
            np.sum(m*z)/np.sum(m), np.sum(m*s)/np.sum(m)]])
    # particle-averaged path
    temp_p = np.array([[np.sum(x)/Npart, np.sum(y)/Npart, np.sum(z)/Npart,
            np.sum(s)/Npart]])
    # one single particle [number of particle]
    temp_px = np.array([[x[0], y[0], z[0], s[0]]])
    # add a line to part
    part_px = np.append(part_px, temp_px, axis=0)
    part_p = np.append(part_p, temp_p, axis=0)
    part_m = np.append(part_m, temp_m, axis=0)
    part_kE = np.append(part_kE, temp_kE, axis=0)
    count = count + 1

    variable = particles['h']
    cc = variable

lineName_px = 'myLine_px'
lineName_p = 'myLine_p'
lineName_m = 'myLine_m'
lineName_kE = 'myLine_kE'
writeLine2SHPfile(part_px, lineName_px, fileName_px)
writeLine2SHPfile(part_p, lineName_p, fileName_p)
writeLine2SHPfile(part_m, lineName_m, fileName_m)
writeLine2SHPfile(part_kE, lineName_kE, fileName_kE)

"""elongate path with linear linextrapolation
"""
if pathelongation == 'linearextrapolationbeta':
    lineName_long = 'myLine_m_ext'

    # input path com2AB
    length = s[-1]/3
    print(length)
    start_extend3 = len(part_m)/3
    s_e3 = int(round(-start_extend3))
    len_ext3 = math.sqrt(pow(part_m[s_e3, 0] - part_m[-1, 0], 2.0) + pow(
        part_m[s_e3, 1] - part_m[-1, 1], 2.0))
    x_ext3 = np.array([[part_m[-1, 0] + (part_m[-1, 0] - part_m[s_e3, 0]) / len_ext3 * length]])
    y_ext3 = np.array([[part_m[-1, 1] + (part_m[-1, 1] - part_m[s_e3, 1]) / len_ext3 * length]])
    # z_ext = part_m[-1, 2] + (part_m[-1, 2] - part_m[s_e, 2]) / len_ext * length
    s_ext3 = np.array([[part_m[-1, 3] + (part_m[-1, 3] - part_m[s_e3, 3]) / len_ext3 * length]])
    z_ext3 = geoTrans.projectOnGrid(x_ext3[:, 0], y_ext3[:, 0], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
    x_ext3_f = float(x_ext3)
    y_ext3_f = float(y_ext3)
    z_ext_f = float(z_ext3[0])
    s_ext3_f = float(s_ext3)
    part_ext = np.array([[x_ext3_f, y_ext3_f, z_ext_f, s_ext3_f]])
    part_long = np.append(part_m, part_ext, axis=0)
    writeLine2SHPfile(part_long, lineName_long, fileName_long_m)

    # copied from runCom2AB.py
    # log file name; leave empty to use default runLog.log
    logName = 'runCom2AB'

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    # Load all input Parameters from config file
    # get the configuration of an already imported module
    # write config to log file
    cfg = cfgUtils.getModuleConfig(com2AB)

    # Calculate ALPHABETA
    resAB = com2AB.com2ABMain(cfg, avalancheDir)

    # Analyse/ plot/ write results #
    reportDictList = []
    _, plotFile, writeFile = outAB.writeABpostOut(resAB, cfg, reportDictList)

    log.info('Plotted to: %s', plotFile)
    log.info('Data written: %s', writeFile)

    # extrapolate path starting at beta Point
    name = 'myLine_m_ext'

    ids10Point = resAB[name]['ids10Point']
    x_b10 = resAB[name]['x'][resAB[name]['ids10Point']]
    y_b10 = resAB[name]['y'][resAB[name]['ids10Point']]
    s_b10 = resAB[name]['s'][resAB[name]['ids10Point']]
    z_b10 = resAB[name]['z'][resAB[name]['ids10Point']]

    # input path com2AB
    len_ext = math.sqrt(pow(x_b10 - part_m[-1, 0], 2.0) + pow(y_b10 - part_m[-1, 1], 2.0))
    x_ext = np.array([[part_m[-1, 0] + (part_m[-1, 0] - x_b10) / len_ext * length]])
    y_ext = np.array([[part_m[-1, 1] + (part_m[-1, 1] - y_b10) / len_ext * length]])
    # z_ext = np.array([[part_m[-1, 2] + (part_m[-1, 2] - z_b10) / len_ext * length]])
    s_ext = np.array([[part_m[-1, 3] + (part_m[-1, 3] - s_b10) / len_ext * length]])
    z_ext = geoTrans.projectOnGrid(x_ext[:, 0], y_ext[:, 0], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
    x_ext_f = float(x_ext)
    y_ext_f = float(y_ext)
    z_ext_f = float(z_ext[0])
    s_ext_f = float(s_ext)
    part_ext = np.array([[x_ext_f, y_ext_f, z_ext_f, s_ext_f]])
    part_aimec = np.append(part_m, part_ext, axis=0)
    writeLine2SHPfile(part_long, lineName_long, fileName_long_m)

"""copied from runCom2AB.py
"""
# log file name; leave empty to use default runLog.log
logName = 'runCom2AB'

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
# write config to log file
cfg = cfgUtils.getModuleConfig(com2AB)

# Calculate ALPHABETA
resAB = com2AB.com2ABMain(cfg, avalancheDir)

# Analyse/ plot/ write results #
reportDictList = []
_, plotFile, writeFile = outAB.writeABpostOut(resAB, cfg, reportDictList)

log.info('Plotted to: %s', plotFile)
log.info('Data written: %s', writeFile)

"""copied from runAna3AIMEC.py
"""
# -----------Required settings-----------------
# log file name; leave empty to use default runLog.log
logName = 'runAna3AIMEC'

# ---------------------------------------------

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
# write config to log file
cfg = cfgUtils.getModuleConfig(ana3AIMEC)

iP.cleanModuleFiles(avalancheDir, ana3AIMEC)

# write configuration to file
cfgUtils.writeCfgFile(avalancheDir, ana3AIMEC, cfg)

cfgSetup = cfg['AIMECSETUP']
anaMod = cfgSetup['anaMod']

# Setup input from com1DFA
pathDict = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfg)

# TODO: define referenceFile
pathDict['numSim'] = len(pathDict['ppr'])

# define reference simulation
pathDict = aimecTools.fetchReferenceSimNo(pathDict, cfgSetup)

pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)

startTime = time.time()

log.info("Running ana3AIMEC model on test case DEM: \n %s \n with profile: \n %s ",
         pathDict['demSource'], pathDict['profileLayer'])
# Run AIMEC postprocessing
ana3AIMEC.mainAIMEC(pathDict, cfg)

endTime = time.time()

log.info(('Took %s seconds to calculate.' % (endTime - startTime)))

rasterTransfo, newRasters, resAnalysis = ana3AIMEC.mainAIMEC(pathDict, cfg)
# create local logger
log = logging.getLogger(__name__)

"""plots
"""

if pathelongation == 'linearextrapolationbeta':
    name = 'myLine_m_ext'
else:
    name = 'myLine_m_smalla'

#  Plot results
s_aimec = resAnalysis['runout'][0]
x_aimec = resAnalysis['runout'][1]
y_aimec = resAnalysis['runout'][2]
deltaH = resAnalysis['deltaH']
s = resAB[name]['s']
z = resAB[name]['z']
f = resAB[name]['f']

if projectedZ == 'yes':
    z_aimec_tuble = geoTrans.projectOnGrid(x_aimec, y_aimec, ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
    z_aimec = z_aimec_tuble[0]
    z_aimec.flatten()
else:
    z_aimec = f[0] - deltaH

print(x_aimec, y_aimec, s_aimec, z_aimec, resAB[name]['x'][resAB[name]['ids_alpha']],
     resAB[name]['y'][resAB[name]['ids_alpha']], resAB[name]['s'][resAB[name]['ids_alpha']],
     resAB[name]['z'][resAB[name]['ids_alpha']])

ids_alpha = resAB[name]['ids_alpha']
indSplit = resAB[name]['indSplit']
# splitPoint = resAB[name]['splitPoint']

g = 9.81
fig = plt.figure(figsize=(2*figW, figH))
ax1 = plt.subplot(211)
cmap = cmapReds
unit = 'm'
cmap.set_under(color='w')
# ax1.plot(part_p[:, 0], part_p[:, 1], 'g--', linewidth=1,
#         label='Lawinenpfad ungewichtet (Lp_p)')
# ax1.plot(part_m[:, 0], part_m[:, 1], 'b-', linewidth=1,
#         label='Lp massengemittelt (Lp_m)')
# ax1.plot(part_kE[:, 0], part_kE[:, 1], 'r-.', linewidth=1,
#         label='Lp kinetischer Energie gewichtet (Lp_kE)')
if pathelongation == 'smallalpha':
    ax1.plot(part_m[:, 0], part_m[:, 1], 'b-', label='Lawinenpfad massengemittelt (Lp_m) α=%1.2f' % alpha + '°')
    ax1.plot(part_smalla[:, 0], part_smalla[:, 1], 'k--', label='Lp_m α=%1.2f' % sa + '°')

if pathelongation == 'linearextrapolationbeta':
    ax1.plot(part_m[:, 0], part_m[:, 1], 'b-', label='Lawinenpfad massengemittelt (Lp_m)')
    # ax1.plot(part_m[-1, 0], part_m[-1, 1], 'x', color='b', markersize=8, label='Pfadende')
    ax1.plot(part_aimec[:, 0], part_aimec[:, 1], 'k--', label='verlängerter Lp_m')
    # ax1.plot(part_m[s_e3, 0], part_m[s_e3, 1], 'x', color='r', markersize=8, label='Start Extrapol. 1/3')
    # ax1.plot(x_b10, y_b10, 'x', color='g', markersize=8, label='Start Extrapol. b10')
    # ax1.plot(resAB[name]['x'][resAB[name]['splitPoint']], resAB[name]['y'][resAB[name]['splitPoint']], 'x', color='c', markersize=8, label='splitPoint')
    ax1.plot(resAB[name]['x'][resAB[name]['ids10Point']], resAB[name]['y'][resAB[name]['ids10Point']],
        'x', color='g', markersize=8, label='Beta')

ax1.plot(x_aimec[-1], y_aimec[-1], 'X', markersize=8, color='0.7',
     label='Auslauf (A) com1DFA (AIMEC pfd=0m)')
ax1.plot(resAB[name]['x'][resAB[name]['ids_alpha']], resAB[name]['y'][resAB[name]['ids_alpha']],
     'x', markersize=8, label='A com2AB')

pfdDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'peakFiles')
pfdFiles = list(pfdDir.glob('*_pfd.asc'))
pfdFile = str(pfdFiles[0])

# Fließhöhe
raster = IOf.readRaster(pfdFile, noDataToNan=True)
data1 = raster['rasterData']
header = IOf.readASCheader(pfdFile)
cellSize = header['cellsize']
ny = data1.shape[0]
nx = data1.shape[1]
Ly = ny * cellSize
Lx = nx * cellSize

plt.imshow(data1, extent=[1000, Lx+1000, -5000, Ly-5000])
cmap, _, ticks, norm = pU.makeColorMap(cmap, np.nanmin(data1), np.nanmax(data1), continuous=pU.contCmap)
cmap.set_bad('w')
data1P = ma.masked_where(data1 == 0.0, data1)
im1 = plt.imshow(data1P, cmap=cmap, extent=[1000, Lx+1000, -5000, Ly-5000], origin='lower', aspect='auto', norm=norm)
cbarpfd = ax1.figure.colorbar(im1, ax=ax1, use_gridspec=True)
cbarpfd.ax.set_ylabel('max. Fließhöhe (pfd) [m]')

ax1.legend()

if avaName == 'avaSei':
    plt.imshow(data1, extent=[251855.110199999996, Lx+251855.110199999996,  376544.749233085488, Ly+376544.749233085488])
    im1 = plt.imshow(data1P, cmap=cmap, extent=[251855.110199999996, Lx+251855.110199999996,  376544.749233085488, Ly+376544.749233085488], origin='lower', aspect='auto', norm=norm)
    ax1.set_xlim([252400, 253400])
    ax1.set_ylim([378000, 381500])
    ax1.legend(loc='lower left')

if avaName == 'avaParabola':
    ax1.set_ylim([-4150, -3850])
    ax1.set_xlim([2500, 5200])
    ax1.legend(loc='lower left')

if avaName == 'avaHelix':
    ax1.set_xlim([2000, 3200])
    ax1.set_ylim([-5000, -3200])
    ax1.legend(loc='upper left')

# Partikelhöhe
# cmap = cmapBlues
# scat = ax1.scatter(x, y, c=cc, cmap=cmap, marker='.', label='Partikelhöhe letzter Zeitschritt')
# cbar1 = ax1.figure.colorbar(scat, ax=ax1, use_gridspec=True)
# cbar1.ax.set_ylabel('Partikelhöhe [m]')

# givenPath
# pathAB_given = os.path.join(avaDir, 'Inputs', 'LINES', 'givenPaths', 'pathAB.shp')
# listx = []
# listy = []
# test = shapefile.Reader(pathAB_given)
# for sr in test.shapeRecords():
#     for xNew, yNew in sr.shape.points:
#         listx.append(xNew)
#         listy.append(yNew)
# ax1.plot(listx, listy, 'g--', linewidth=1.0, label='givenPath')

ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')

cmap = cmapPlasma
cmap.set_under(color='w')
ax2 = plt.subplot(212)

if sreal == 'yes':
    dx = part_m[1:, 0] - part_m[:-1, 0]
    dy = part_m[1:, 1] - part_m[:-1, 1]
    dx = np.hstack((np.array([0]), dx))
    dy = np.hstack((np.array([0]), dy))
    sStar = np.cumsum(np.sqrt(dx*dx + dy*dy))
    z_m, _ = geoTrans.projectOnGrid(part_m[:, 0], part_m[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
    z_aimec_tuble = geoTrans.projectOnGrid(part_aimec[:, 0], part_aimec[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
    z_aimec = np.array(z_aimec_tuble[0])
    z_aimec.flatten()
    ax2.plot(sStar, z_m, 'b-', label='Z_av(s_real)')
    ax2.plot(part_m[:, 3], z_m, 'b-.', label='Z_av(s_mod)')
    ax2.plot(sStar, part_m[:, 2], 'k-', label='Z_true(s_real)')
    ax2.plot(part_m[:, 3], part_m[:, 2], 'k-.', label='Z_true(s_mod)')
    GK = part_m[-1, 3] * np.tan(alpha*np.pi/180)
    z_ende = part_m[0, 2] - GK
    s_geomL = [part_m[0, 3], part_m[-1, 3]]
    z_geomL = [part_m[0, 2], z_ende]
    ax2.plot(s_geomL, z_geomL, 'r-', linewidth=0.3, label='alpha line from Z_av func s_mod')
    GK = sStar[-1] * np.tan(alpha*np.pi/180)
    z_ende = z_m[0] - GK
    s_geomL = [sStar[0], sStar[-1]]
    z_geomL = [z_m[0], z_ende]
    ax2.plot(s_geomL, z_geomL, 'g-', linewidth=0.3, label='alpha line from Z_true func s_mod')
    plt.axvline(x=s[ids10Point], color='g', linewidth=1, linestyle='-.', label='Beta')

    Zene = part_m[:, 2] + V2Path/(2*g)
    # Colorbar: kinietische Energie [J]
    scat = ax2.scatter(part_m[:, 3], Zene, marker='s', cmap=cmap, s=2*ms, c=EkinPath, label='Gesamtenergie(s_mod)')
    scat = ax2.scatter(sStar, Zene, marker='o', cmap=cmap, s=2*ms, c=EkinPath, label='Gesamtenergie(s_real)')
    cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('kinetische Energie [J]')

else:
    if projectedZ == 'yes':
        # z_px = geoTrans.projectOnGrid(part_px[:, 0], part_px[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
        z_p = geoTrans.projectOnGrid(part_p[:, 0], part_p[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
        z_kE = geoTrans.projectOnGrid(part_kE[:, 0], part_kE[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
        z_m, _ = geoTrans.projectOnGrid(part_m[:, 0], part_m[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
        # ax2.plot(part_px[:, 3], z_px[0], 'm:', label='Lp_px_projZ')
        # ax2.plot(part_p[:, 3], z_p[0], 'g-', label='Lp_p_projZ')
        # ax2.plot(part_m[:, 3], z_m[0], 'b--', label='Lp_m_projZ')
        # ax2.plot(part_kE[:, 3], z_kE[0], 'r--', label='Lp_kE_projZ')

        if pathelongation == 'smallalpha':
            z_smalla_tuble = geoTrans.projectOnGrid(part_smalla[:, 0], part_smalla[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
            z_smalla = np.array(z_smalla_tuble[0])
            z_smalla.flatten()
            ax2.plot(part_m[:, 3], z_m, 'b-', label='Lp_m_projZ α=%1.2f' % alpha + '°')
            ax2.plot(part_smalla[:, 3], z_smalla, 'k--', label='Lp_m_projZ α=%1.2f' % sa + '°')
            GK = part_smalla[-1, 3] * np.tan(alpha*np.pi/180)
            z_ende = z_smalla[0] - GK
            s_geomL = [part_smalla[0, 3], part_smalla[-1, 3]]
            z_geomL = [z_smalla[0], z_ende]

        if pathelongation == 'linearextrapolationbeta':
            z_aimec_tuble = geoTrans.projectOnGrid(part_aimec[:, 0], part_aimec[:, 1], ZZ, csz=5, xllc=xllc, yllc=yllc, interp='bilinear')
            z_aimec = np.array(z_aimec_tuble[0])
            z_aimec.flatten()
            ax2.plot(part_m[:, 3], z_m, 'b-', label='Lp_m_projZ')
            ax2.plot(part_aimec[:, 3], z_aimec, 'k--', label='lin. extrapol. Lp_m_projZ')
            GK = part_aimec[-1, 3] * np.tan(alpha*np.pi/180)
            z_ende = part_aimec[0, 2] - GK
            s_geomL = [part_aimec[0, 3], part_aimec[-1, 3]]
            z_geomL = [part_aimec[0, 2], z_ende]
            # plt.axvline(x=s[splitPoint], color='0.8', linewidth=1, linestyle='--', label='Split point')
            plt.axvline(x=s[ids10Point], color='g', linewidth=1, linestyle='-.', label='Beta')

    else:
        # ax2.plot(part_p[:, 3], part_p[:, 2], 'g--', label='Lp_p')
        # ax2.plot(part_m[:, 3], part_m[:, 2], 'b-', label='Lp_m')
        # ax2.plot(part_kE[:, 3], part_kE[:, 2], 'r-.', label='Lp_kE')

        if pathelongation == 'smallalpha':
            ax2.plot(part_m[:, 3], part_m[:, 2], 'b-', label='Lp_m α=%1.2f' % alpha + '°')
            ax2.plot(part_smalla[:, 3], part_smalla[:, 2], 'k--', label='Lp_m α=%1.2f' % sa + '°')
            GK = part_smalla[-1, 3] * np.tan(alpha*np.pi/180)
            z_ende = part_smalla[0, 2] - GK
            s_geomL = [part_smalla[0, 3], part_smalla[-1, 3]]
            z_geomL = [part_smalla[0, 2], z_ende]

        if pathelongation == 'linearextrapolationbeta':
            ax2.plot(part_m[:, 3], part_m[:, 2], 'b-', label='Lp_m')
            ax2.plot(part_aimec[:, 3], part_aimec[:, 2], 'k--', label='lin. extrapol. Lp_m')
            GK = part_aimec[-1, 3] * np.tan(alpha*np.pi/180)
            z_ende = part_aimec[0, 2] - GK
            s_geomL = [part_aimec[0, 3], part_aimec[-1, 3]]
            z_geomL = [part_aimec[0, 2], z_ende]
            # plt.axvline(x=s[splitPoint], color='0.8', linewidth=1, linestyle='--', label='Split point')
            plt.axvline(x=s[ids10Point], color='g', linewidth=1, linestyle='-.', label='Beta')

    ax2.plot(s_geomL, z_geomL, 'k-', linewidth=0.3, label='geometrische Lösung')
    ax2.plot(s, f, 'g-', linewidth=0.3, label='Alphalinie com2AB')

    Zene = part_m[:, 2] + V2Path/(2*g)
    # Colorbar: kinietische Energie [J]
    # scat = ax2.scatter(part_m[:, 3], Zene, marker='s', cmap=cmap, s=2*ms, c= EkinPath,
    #     label='Gesamtenergie com1DFA')
    # cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    # cbar2.ax.set_ylabel('kinetische Energie [J]')

    # Colorbar: Energiehöhe [m]
    Zene_EH = V2Path/(2*g)
    scat_EH = ax2.scatter(part_m[:, 3], Zene, marker='s', cmap=cmap, s=2*ms, c= Zene_EH,
        label='Gesamtenergie com1DFA')
    cbar2 = ax2.figure.colorbar(scat_EH, ax=ax2, use_gridspec=True)
    # cbar2.ax.set_ylabel('Energiehöhe [m] aus Ekin')
    cbar2.ax.set_ylabel(r'$\frac{Ekin}{2g}$ [m]')

ax2.plot(s_aimec[-1], z_aimec[-1], 'X', markersize=8, color='0.7', label='A com1DFA = %1.2f' % s_aimec[-1] + 'm')
ax2.plot(resAB[name]['s'][resAB[name]['ids_alpha']], resAB[name]['z'][resAB[name]['ids_alpha']],
    'x', markersize=8, label='A com2AB   = %1.2f' % resAB[name]['s'][resAB[name]['ids_alpha']] + 'm')

ax2.set_xlabel('s [m]', fontsize=fs)
ax2.set_ylabel('Höhe [m]', fontsize=fs)
ax2.legend(loc='lower left')
fig.tight_layout()

# set titel of output png
if sreal == 'yes':
    if oneParticle == 'yes':
        title = (avalancheTitle + '_%1.1f' % alpha + '_1p_linext_projZ_s_real.png')
    else:
        title = (avalancheTitle + '_%1.1f' % alpha + '_linext_projZ_s_real.png')
else:
    if oneParticle == 'yes':
        if pathelongation == 'linearextrapolationbeta':
            if projectedZ == 'yes':
                title = (avalancheTitle + '_%1.1f' % alpha + '_1p_linext_projZ.png')
            else:
                title = (avalancheTitle + '_%1.1f' % alpha + '_1p_linext.png')
        else:
            if projectedZ == 'yes':
                title = (avalancheTitle + '_%1.1f' % alpha + '_1p_projZ.png')
            else:
                title = (avalancheTitle + '_%1.1f' % alpha + '_1p.png')
    else:
        if pathelongation == 'linearextrapolationbeta':
            if projectedZ == 'yes':
                title = (avalancheTitle + '_%1.1f' % alpha + '_linext_projZ.png')
            else:
                title = (avalancheTitle + '_%1.1f' % alpha + '_linext.png')
        else:
            if projectedZ == 'yes':
                title = (avalancheTitle + '_%1.1f' % alpha + '_projZ.png')
            else:
                title = (avalancheTitle + '_%1.1f' % alpha + '.png')

path = os.path.join(avaDir, title)
plt.savefig(path)
plt.show()
