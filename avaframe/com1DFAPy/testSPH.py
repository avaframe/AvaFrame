import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.ascUtils as IOf
from avaframe.out3Plot.plotUtils import *
import avaframe.com1DFAPy.timeDiscretizations as tD
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.com1DFAPy.SPHfunctions as SPH
from avaframe.in3Utils import cfgUtils
cfg = cfgUtils.getModuleConfig(com1DFA)['GENERAL']
cfgFull = cfgUtils.getModuleConfig(com1DFA)


Lx = 30
Ly = 35

dx = 0.5
dy = 0.5

DX = 5
DY = 5
csz = DX

coef = 1
rho = 200

# define particles
nx = np.int(Lx/dx)
ny = np.int(Ly/dy)
Npart = nx*ny
x = np.linspace(0, Lx-dx, nx)
y = np.linspace(0, Ly-dy, ny)
xx, yy = np.meshgrid(x, y)
xx = xx.flatten()
yy = yy.flatten()
Xpart = xx + (np.random.rand(Npart) - 0.5) * dx * coef
Ypart = yy + (np.random.rand(Npart) - 0.5) * dy * coef
Hpart = (np.sin(math.pi*Xpart/Lx))*np.sin(math.pi*Xpart/Lx) + Ypart/Ly
Mpart = Hpart * dx * dy * rho
# create dictionnary to store particles properties
particles = {}
particles['Npart'] = Npart
particles['mTot'] = np.sum(Mpart)
particles['x'] = Xpart
particles['y'] = Ypart
# adding z component
Zpart = np.zeros(np.shape(Ypart))
particles['z'] = Zpart
particles['s'] = Zpart
particles['ux'] = Zpart
particles['uy'] = Zpart
particles['uz'] = Zpart

particles['m'] = Mpart
particles['h'] = Hpart
# define grid
NX = np.int(Lx/DX + 1)
NY = np.int(Ly/DY + 1)
header = IOf.cASCheader()
header.ncols = NX
header.nrows = NY
header.cellsize = csz
dem = {}
dem['header'] = header
X = np.linspace(0, Lx, NX)
Y = np.linspace(0, Ly, NY)
XX, YY = np.meshgrid(X, Y)
# find neighbours
particles = SPH.getNeighboursVect(particles, dem)

indOut = np.where(DFAtls.norm(Xpart-Lx/2, Ypart-Ly/2, Zpart) > 4*Lx/2.5)
mask = np.ones(len(Xpart), dtype=bool)
mask[indOut] = False

nRemove = len(mask)-np.sum(mask)
if nRemove > 0:
    particles = com1DFA.removePart(particles, mask, nRemove)
    print('removed %s particles' % (nRemove))


fig, ax = plt.subplots(figsize=(figW, figH))
Xpart = particles['x']
Ypart = particles['y']
Mpart = particles['m']
ax.plot(Xpart, Ypart, color='r', marker='.', linestyle='None')
ax.plot(XX, YY, color='k', marker='.', linestyle='None')
# plt.show()

# find neighbours
particles = SPH.getNeighboursVect(particles, dem)

# Compute sph FD
particles = SPH.computeFD(cfg, particles, dem)

# Compute SPH gradient

Npart = particles['Npart']
# initialize
GHX = np.zeros(Npart)
GHY = np.zeros(Npart)
GHZ = np.zeros(Npart)
# loop on particles
# TcpuSPH = 0
# Tcpuadd = 0
for j in range(Npart):
    mass = particles['m'][j]
    # adding lateral force (SPH component)
    # startTime = time.time()
    gradhX, gradhY,  gradhZ, _ = SPH.calcGradHSPHVect(particles, j, NX, NY, csz)
    # tcpuSPH = time.time() - startTime
    # TcpuSPH = TcpuSPH + tcpuSPH
    # startTime = time.time()
    GHX[j] = GHX[j] - gradhX / rho
    GHY[j] = GHY[j] - gradhY / rho
    GHZ[j] = GHZ[j] - gradhZ / rho
    # tcpuadd = time.time() - startTime
    # Tcpuadd = Tcpuadd + tcpuadd

# log.info(('cpu time SPH = %s s' % (TcpuSPH / Npart)))
# log.info(('cpu time SPH add = %s s' % (Tcpuadd / Npart)))


# Update fields using a nearest interpolation
MassNearest = np.zeros((NY, NX))
# startTime = time.time()
iC = particles['InCell'][:, 2]
MassNearest = MassNearest.flatten()
np.add.at(MassNearest, iC, Mpart)
MassNearest = np.reshape(MassNearest, (NY, NX))
FDNearest = MassNearest / (DX*DY * rho)
hNN, _ = geoTrans.projectOnRasterVectRoot(Xpart, Ypart, FDNearest, csz=csz, interp='nearest')
hNB, _ = geoTrans.projectOnRasterVectRoot(Xpart, Ypart, FDNearest, csz=csz, interp='bilinear')

# Update fields using a bilinear interpolation
MassBilinear = np.zeros((NY, NX))
#
# # startTime = time.time()
MassBilinear = geoTrans.pointsToRaster(Xpart, Ypart, Mpart, MassBilinear, csz=csz, interp='bilinear')
FDBilinear = MassBilinear / (DX*DY * rho)
hBB, _ = geoTrans.projectOnRasterVectRoot(Xpart, Ypart, FDBilinear, csz=csz, interp='bilinear')

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111, projection='3d')
# ax2.scatter(particles['x'], particles['y'], GHX, 'r')
# ax2.scatter(particles['x'], particles['y'], np.ones(np.shape(particles['x']))/Lx, 'b')
#
fig3 = plt.figure()
ax3 = fig3.add_subplot(111, projection='3d')
ax3.scatter(particles['x'], particles['y'], hBB, 'k')
ax3.scatter(particles['x'], particles['y'], particles['hSPH'], 'r')
ax3.scatter(particles['x'], particles['y'], ((np.sin(math.pi*Xpart/Lx))*np.sin(math.pi*Xpart/Lx)) + Ypart/Ly, 'b')



ind = np.where(((particles['y']>Ly/2-dy/2) & (particles['y']<Ly/2+dy/2)))
fig1, ax1 = plt.subplots(figsize=(2*figW, figH))
ax1.plot(particles['x'][ind], ((np.sin(math.pi*Xpart/Lx))*np.sin(math.pi*Xpart/Lx))[ind] + Ypart[ind]/Ly, color='b', marker='None')
ax1.plot(particles['x'][ind], particles['hSPH'][ind], color='r', marker='.', linestyle = 'None')
ax1.plot(particles['x'][ind], hBB[ind], color='k', marker='.', linestyle = 'None')
ax1.plot(particles['x'][ind], hNN[ind], color='y', marker='.', linestyle = 'None')
ax1.plot(particles['x'][ind], hNB[ind], color='g', marker='.', linestyle = 'None')

ind = np.where(((particles['x']>Lx/2-dx/2) & (particles['x']<Lx/2+dx/2)))
fig2, ax2 = plt.subplots(figsize=(2*figW, figH))
ax2.plot(particles['y'][ind], ((np.sin(math.pi*Xpart/Lx))*np.sin(math.pi*Xpart/Lx))[ind] + particles['y'][ind]/Ly, color='b', marker='None')
ax2.plot(particles['y'][ind], particles['hSPH'][ind], color='r', marker='.', linestyle = 'None')
ax2.plot(particles['y'][ind], hBB[ind], color='k', marker='.', linestyle = 'None')
ax1.plot(particles['y'][ind], hNN[ind], color='y', marker='.', linestyle = 'None')
ax1.plot(particles['y'][ind], hNB[ind], color='g', marker='.', linestyle = 'None')

ind = np.where(((particles['y']>Ly/2-dy/2) & (particles['y']<Ly/2+dy/2)))
fig3, ax3 = plt.subplots(figsize=(2*figW, figH))
ax3.plot(particles['x'][ind], Xpart[ind]/(Lx*Lx), color='b', marker='None')
ax3.plot(particles['x'][ind], GHX[ind], color='r', marker='.', linestyle = 'None')

ind = np.where(((particles['x']>Lx/2-dx/2) & (particles['x']<Lx/2+dx/2)))
fig4, ax4 = plt.subplots(figsize=(2*figW, figH))
ax4.plot(particles['y'][ind], np.ones(np.shape(particles['y'][ind]))/Ly, color='b', marker='None')
ax4.plot(particles['y'][ind], GHY[ind], color='r', marker='.', linestyle = 'None')

plt.show()
