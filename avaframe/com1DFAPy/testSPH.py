import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

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

########################################################################
## CHOOSE YOUR SETUP
##########################################################################
# Choose the snow depth you want to use (h function)
def Hfunction(x, y, z):
    h = np.ones(np.shape(x))
    GHx = np.zeros(np.shape(x))
    GHy = np.zeros(np.shape(x))
    GHz = np.zeros(np.shape(x))
    return h, GHx, GHy, GHz


# Choose the snow depth you want to use (h function)
def Sfunction(x, y, Lx, Ly):
    # horizontal plane
    Z = np.ones(np.shape(x))
    area = 1

    # plane
    # Z = x*np.tan(slopeAngle) + y*np.tan(slopeAngle)
    # area = np.sqrt(1 + (np.tan(slopeAngle))*(np.tan(slopeAngle)) + (np.tan(slopeAngle))*(np.tan(slopeAngle)))
    # quadratic surface
    Z = 1/Lx * x*x - 2*x + Lx + (1/Ly * y*y - 2*y + Ly)
    area = np.sqrt(1 + (2/Lx*x - 2)*(2/Lx*x - 2) + (2/Ly*y - 2)*(2/Ly*y - 2))
    return Z, area


slopeAngle = 30*math.pi/180

# set the size of the mesh grid [m]
NDX = [5]
# NDX = [7.5, 5, 3, 2, 1]
# choose the number of particles per DX and DY
# if you choose 3, you will have 3*3 = 9 particles per grid cell
# NPartPerD = [40]
NPartPerD = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40]

# choose if the particles should be randomly distributed.
# 0 no random part, up to 1, random fluctuation of dx/2 and dy/2
coef = 0.5
rho = 200
##############################################################################
## END CHOOSE SETUP
###############################################################################


fig1, ax1 = plt.subplots(figsize=(2*figW, figH))
fig2, ax2 = plt.subplots(figsize=(2*figW, figH))
fig3, ax3 = plt.subplots(figsize=(2*figW, figH))
fig4, ax4 = plt.subplots(figsize=(2*figW, figH))


# color = cmapAimec(np.linspace(1, 0, 2*len(NDX) + 3, dtype=float))
color = cmapAimec(np.linspace(1, 0, 2*len(NPartPerD) + 3, dtype=float))
markers = ['o', 's', 'd', '*', 'p', 'P', '^', '>', '<', 'X', 'h']
count = 0
for DX in NDX:
    DY = DX
    csz = DX

    # set the extend of your mesh
    Lx = 3*DX
    Ly = 3*DY

    for nPartPerD in NPartPerD:
        dx = DX/nPartPerD
        dy = DY/nPartPerD

        # ------------------------------------------
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
        # adding z component
        Zpart, area = Sfunction(Xpart, Ypart, Lx, Ly)
        Hpart, _, _, _ = Hfunction(Xpart, Ypart, Zpart)
        Mpart = Hpart * dx * dy * area * rho
        # create dictionnary to store particles properties
        particles = {}
        particles['Npart'] = Npart
        particles['mTot'] = np.sum(Mpart)
        particles['x'] = Xpart
        particles['y'] = Ypart
        particles['z'] = Zpart
        particles['s'] = np.zeros(np.shape(Ypart))
        particles['ux'] = np.zeros(np.shape(Ypart))
        particles['uy'] = np.zeros(np.shape(Ypart))
        particles['uz'] = np.zeros(np.shape(Ypart))

        particles['m'] = Mpart
        particles['h'] = Hpart

        # ------------------------------------------
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
        ZZ, _ = Sfunction(XX, YY, Lx, Ly)
        dem['rasterData'] = ZZ
        # Initialize mesh
        dem = com1DFA.initializeMesh(dem, num=4)

        # ------------------------------------------
        # find neighbours
        particles = SPH.getNeighboursVect(particles, dem)

        # ------------------------------------------
        # indOut = np.where(DFAtls.norm(Xpart-Lx/2, Ypart-Ly/2, Zpart) >6* Lx/2)
        # mask = np.ones(len(Xpart), dtype=bool)
        # mask[indOut] = False
        #
        # nRemove = len(mask)-np.sum(mask)
        # if nRemove > 0:
        #     particles = com1DFA.removePart(particles, mask, nRemove)
        #     print('removed %s particles' % (nRemove))
        #
        # # find neighbours
        # particles = SPH.getNeighboursVect(particles, dem)

        Xpart = particles['x']
        Ypart = particles['y']
        Mpart = particles['m']

        indPartInCell = particles['indPartInCell']
        partInCell = particles['partInCell']


        Nx = dem['Nx']
        Ny = dem['Ny']
        Nz = dem['Nz']
        for ic in range(NX*NY):
            part = partInCell[indPartInCell[ic]: indPartInCell[ic+1]]
            fig, ax = plt.subplots(figsize=(figW, figH))
            ax.plot(Xpart, Ypart, color='r', marker='.', linestyle='None')
            ax.plot(Xpart[part], Ypart[part], color='b', marker='.', linestyle='None')
            ax.plot(XX, YY, color='k', marker='.', linestyle='None')
            # fig2 = plt.figure()
            # ax2 = fig2.add_subplot(111, projection='3d')
            # ax2.scatter(Xpart, Ypart, Zpart, color='r', marker='.')
            # ax2.scatter(Xpart[part], Ypart[part], Zpart[part], color='b', marker='.')
            # ax2.scatter(XX, YY, ZZ, color='k', marker='.')
            plt.close()
            # plt.show()

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

            x = particles['x'][j]
            y = particles['y'][j]
            nx, ny, nz = DFAtls.getNormal(x, y, Nx, Ny, Nz, csz)
            gradhX, gradhY,  gradhZ, _ = SPH.calcGradHSPHVect(particles, j, NX, NY, csz, nx, ny, nz)
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

        Area = dem['Area']
        # Update fields using a nearest interpolation
        MassNearest = np.zeros((NY, NX))
        # startTime = time.time()
        iC = particles['InCell'][:, 2]
        MassNearest = MassNearest.flatten()
        np.add.at(MassNearest, iC, Mpart)
        MassNearest = np.reshape(MassNearest, (NY, NX))
        FDNearest = MassNearest / (Area * rho)
        hNN, _ = geoTrans.projectOnRasterVectRoot(Xpart, Ypart, FDNearest, csz=csz, interp='nearest')
        hNB, _ = geoTrans.projectOnRasterVectRoot(Xpart, Ypart, FDNearest, csz=csz, interp='bilinear')

        # Update fields using a bilinear interpolation
        MassBilinear = np.zeros((NY, NX))
        #
        # # startTime = time.time()
        MassBilinear = geoTrans.pointsToRaster(Xpart, Ypart, Mpart, MassBilinear, csz=csz, interp='bilinear')
        FDBilinear = MassBilinear / (Area * rho)
        hBB, _ = geoTrans.projectOnRasterVectRoot(Xpart, Ypart, FDBilinear, csz=csz, interp='bilinear')

        # fig2 = plt.figure()
        # ax2 = fig2.add_subplot(111, projection='3d')
        # ax2.scatter(particles['x'], particles['y'], GHX, 'r')
        # ax2.scatter(particles['x'], particles['y'], np.ones(np.shape(particles['x']))/Lx, 'b')
        #
        # fig3 = plt.figure()
        # ax3 = fig3.add_subplot(111, projection='3d')
        # ax3.scatter(particles['x'], particles['y'], hBB, 'k')
        # ax3.scatter(particles['x'], particles['y'], particles['hSPH'], 'r')
        h, Ghx, Ghy, Ghz = Hfunction(particles['x'], particles['y'], particles['z'])
        # ax3.scatter(particles['x'], particles['y'], h, 'b')

        count = count + 1
        col = color[2*count]
        mark = markers[count-1]
        ind = np.where(((particles['y']>Ly/2-dy) & (particles['y']<Ly/2+dy)))
        if count == 1:
            ax1.plot(particles['x'][ind], hBB[ind], color='r', marker=mark, linestyle = 'None', label='HBB flow depth')
            ax1.plot(particles['x'][ind], hNN[ind], color='y', marker=mark, linestyle = 'None', label='HNN flow depth')
            ax1.plot(particles['x'][ind], hNB[ind], color='g', marker=mark, linestyle = 'None', label='HNB flow depth')
        else:
            ax1.plot(particles['x'][ind], hBB[ind], color='r', marker=mark, linestyle = 'None')
            ax1.plot(particles['x'][ind], hNN[ind], color='y', marker=mark, linestyle = 'None')
            ax1.plot(particles['x'][ind], hNB[ind], color='g', marker=mark, linestyle = 'None')


        ax1.plot(particles['x'][ind], h[ind], color='b', marker='None')#, label='real flow depth')
        ax1.plot(particles['x'][ind], particles['hSPH'][ind], color=col, marker=mark, linestyle = 'None', label='SPH N=' + str(nPartPerD))

        ind = np.where(((particles['x']>Lx/2-dx) & (particles['x']<Lx/2+dx)))
        if count == 1:
            ax2.plot(particles['y'][ind], hBB[ind], color='r', marker=mark, linestyle = 'None', label='HBB flow depth')
            ax2.plot(particles['y'][ind], hNN[ind], color='y', marker=mark, linestyle = 'None', label='HNN flow depth')
            ax2.plot(particles['y'][ind], hNB[ind], color='g', marker=mark, linestyle = 'None', label='HNB flow depth')
        else:
            ax2.plot(particles['y'][ind], hBB[ind], color='r', marker=mark, linestyle = 'None')
            ax2.plot(particles['y'][ind], hNN[ind], color='y', marker=mark, linestyle = 'None')
            ax2.plot(particles['y'][ind], hNB[ind], color='g', marker=mark, linestyle = 'None')


        ax2.plot(particles['y'][ind], h[ind], color='r', marker='None')#, label='real flow depth')
        ax2.plot(particles['y'][ind], particles['hSPH'][ind], color=col, marker=mark, linestyle = 'None', label='SPH N=' + str(nPartPerD))

        ind = np.where(((particles['y']>Ly/2-dy) & (particles['y']<Ly/2+dy)))
        ax3.plot(particles['x'][ind], Ghx[ind], color='r', marker='None')#, label='real gradH')
        ax3.plot(particles['x'][ind], GHX[ind], color=col, marker=mark, linestyle = 'None', label='SPH N=' + str(nPartPerD))


        ind = np.where(((particles['x']>Lx/2-dx) & (particles['x']<Lx/2+dx)))
        ax4.plot(particles['y'][ind], Ghy[ind], color='r', marker='None')#, label='real gradH')
        ax4.plot(particles['y'][ind], GHY[ind], color=col, marker=mark, linestyle = 'None', label='SPH N=' + str(nPartPerD))
        plt.show(block=False)
        # plt.draw()
        # plt.pause(2)

fig1.legend()
fig2.legend()
fig3.legend()
fig4.legend()
plt.show()
