import pyximport
pyximport.install()
import numpy as np
import time
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
import avaframe.out3Plot.plotUtils as pU
import avaframe.com1DFAPy.timeDiscretizations as tD
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.com1DFAPy.SPHfunctions as SPH
from avaframe.in3Utils import cfgUtils
from SPHfunctionsCython import *
cfg = cfgUtils.getModuleConfig(com1DFA)['GENERAL']
cfgFull = cfgUtils.getModuleConfig(com1DFA)

########################################################################
# CHOOSE YOUR SETUP
##########################################################################
# Choose the snow depth you want to use (h function)


def Hfunction(x, y, z):
    # h = x*x/1000 + 1
    GHx = 2*x*y/5000
    GHy = x*x/5000
    h = x*x*y/5000 + 1
    # GHx = np.ones(np.shape(x))/10
    # h = np.ones(np.shape(x))
    # GHx = np.zeros(np.shape(x))
    # h = y/50 + 1
    # GHy = np.ones(np.shape(x))/50
    # GHy = np.zeros(np.shape(x))
    GHz = np.zeros(np.shape(x))
    return h, GHx, GHy, GHz


# Choose the surface shape you want to use
def Sfunction(x, y, Lx, Ly):
    # horizontal plane
    # Z = np.ones(np.shape(x))
    # area = 1

    # plane
    Z = x*np.tan(slopeAnglex) + y*np.tan(slopeAngley)
    sx = np.tan(slopeAnglex)*np.ones(np.shape(x))
    sy = np.tan(slopeAngley)*np.ones(np.shape(x))
    area = np.sqrt(1 + (np.tan(slopeAnglex))*(np.tan(slopeAnglex)) + (np.tan(slopeAngley))*(np.tan(slopeAngley)))*np.ones(np.shape(x))
    # quadratic surface
    # Z = 1/(Lx*Lx) * x*x + 2/(Ly*Ly) * y*y
    # area = np.sqrt(1 + (2/(Lx*Lx)*x)*(2/(Lx*Lx)*x) + (2*2/(Ly*Ly)*y)*(2*2/(Ly*Ly)*y))
    # sx = 2/(Lx*Lx)*x
    # sy = 2*2/(Ly*Ly)*y
    # other function
    # Z = Lx*np.exp(1/(Lx*Lx) * x*x) + 1/(Ly) * y*y
    # sx = (Lx*2/(Lx*Lx)*x)*np.exp(1/(Lx*Lx) * x*x)
    # sy = (2/(Ly)*y)
    # area = np.sqrt(1 + (Lx*2/(Lx*Lx)*x)*np.exp(1/(Lx*Lx) * x*x)*(Lx*2/(Lx*Lx)*x)*np.exp(1/(Lx*Lx) * x*x) + (2/(Ly)*y)*(2/(Ly)*y))
    return Z, sx, sy, area


slopeAnglex = 40*math.pi/180
slopeAngley = 20*math.pi/180
# set the size of the mesh grid [m]
NDX = [5]
# NDX = [7.5, 5, 3, 2, 1]
# choose the number of particles per DX and DY
# if you choose 3, you will have 3*3 = 9 particles per grid cell
# NPartPerD = [40]
NPartPerD = [39] #[2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40]

# choose if the particles should be randomly distributed.
# 0 no random part, up to 1, random fluctuation of dx/2 and dy/2
coef = 0.
rho = 200
##############################################################################
# END CHOOSE SETUP
###############################################################################


fig1, ax1 = plt.subplots(figsize=(2*pU.figW, pU.figH))
ax1.set_title('h(x)')
ax1.set_xlabel('x [m]')
ax1.set_ylabel('h [m]')
fig2, ax2 = plt.subplots(figsize=(2*pU.figW, pU.figH))
ax2.set_title('h(y)')
ax2.set_xlabel('y [m]')
ax2.set_ylabel('h [m]')
fig3, ax3 = plt.subplots(figsize=(2*pU.figW, pU.figH))
ax3.set_title('Gradh(x)')
ax3.set_xlabel('x [m]')
ax3.set_ylabel('Gradh []')
fig4, ax4 = plt.subplots(figsize=(2*pU.figW, pU.figH))
ax4.set_title('Gradh(y)')
ax4.set_xlabel('y [m]')
ax4.set_ylabel('Gradh []')


# color = cmapAimec(np.linspace(1, 0, 2*len(NDX) + 3, dtype=float))
color = pU.cmapAimec(np.linspace(1, 0, 2*len(NPartPerD) + 3, dtype=float))
color = pU.cmapAimec(np.linspace(1, 0, 4*len(NPartPerD) + 3, dtype=float))
markers = ['o', 's', 'd', '*', 'p', 'P', '^', '>', '<', 'X', 'h']
count = 0
for DX in NDX:
    DY = DX
    csz = DX

    # set the extend of your mesh
    Lx = 5*DX
    Ly = 5*DY

    for nPartPerD in NPartPerD:
        dx = DX/nPartPerD
        dy = DY/nPartPerD

        # ------------------------------------------
        # define particles
        nx = np.int(Lx/dx)-1
        ny = np.int(Ly/dy)-1
        Npart = nx*ny
        x = np.linspace(dx, Lx-dx, nx)
        y = np.linspace(dy, Ly-dy, ny)
        xx, yy = np.meshgrid(x, y)
        xx = xx.flatten()
        yy = yy.flatten()
        Xpart = xx + (np.random.rand(Npart) - 0.5) * dx * coef
        Ypart = yy + (np.random.rand(Npart) - 0.5) * dy * coef
        # adding z component
        Zpart, sx, sy, area = Sfunction(Xpart, Ypart, Lx, Ly)
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
        ZZ, _, _, _ = Sfunction(XX, YY, Lx, Ly)
        dem['rasterData'] = ZZ
        # Initialize mesh
        dem = com1DFA.initializeMesh(dem, num=4)

        # ------------------------------------------
        # find neighbours
        particles = SPH.getNeighboursVect(particles, dem)

        # ------------------------------------------
        indOut = np.where(DFAtls.norm(Xpart-Lx/2, Ypart-Ly/2, Zpart) >6*Lx/2)
        mask = np.ones(len(Xpart), dtype=bool)
        mask[indOut] = False

        nRemove = len(mask)-np.sum(mask)
        if nRemove > 0:
            particles = com1DFA.removePart(particles, mask, nRemove)
            print('removed %s particles' % (nRemove))

        # find neighbours
        particles = SPH.getNeighboursVect(particles, dem)

        Xpart = particles['x']
        Ypart = particles['y']
        Mpart = particles['m']

        indPartInCell = particles['indPartInCell']
        partInCell = particles['partInCell']

        Nx = dem['Nx']
        Ny = dem['Ny']
        Nz = dem['Nz']
        for ic in range(1): #range(NX*NY):
            part = partInCell[indPartInCell[ic]: indPartInCell[ic+1]]
            # part = np.where(((particles['x'] > Lx/2-dx) & (particles['x'] < Lx/2+dx)))
            fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
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


        # Compute SPH gradient
        m = particles['m']
        x = particles['x']
        y = particles['y']
        z = particles['z']
        ux = particles['ux']
        uy = particles['uy']
        uz = particles['uz']

        indPartInCell = (particles['indPartInCell']).astype('int')
        partInCell = (particles['partInCell']).astype('int')
        indX = particles['indX'].astype('int')
        indY = particles['indY'].astype('int')
        startTime = time.time()
        GHX, GHY, GHZ = computeGradC(cfg, particles, header, Nx, Ny, Nz, indX, indY, SPHOption=1, gradient=1)
        GHX = np.asarray(GHX)
        GHY = np.asarray(GHY)
        GHZ = np.asarray(GHZ)
        tottime = time.time() - startTime
        print('Time SPHOption 1: ', tottime)
        startTime = time.time()
        GHX2, GHY2, GHZ2 = computeGradC(cfg, particles, header, Nx, Ny, Nz, indX, indY, SPHOption=2, gradient=1)
        GHX2 = np.asarray(GHX2)
        GHY2 = np.asarray(GHY2)
        GHZ2 = np.asarray(GHZ2)
        tottime = time.time() - startTime
        print('Time SPHOption 2: ', tottime)
        startTime = time.time()
        startTime = time.time()
        GHX3, GHY3, GHZ3 = computeGradC(cfg, particles, header, Nx, Ny, Nz, indX, indY, SPHOption=3, gradient=1)
        GHX3 = np.asarray(GHX3)
        GHY3 = np.asarray(GHY3)
        GHZ3 = np.asarray(GHZ3)
        tottime = time.time() - startTime
        print('Time SPHOption 3: ', tottime)
        startTime = time.time()
        # Compute sph FD
        H = computeFDC(cfg, particles, header, Nx, Ny, Nz, indX, indY)
        H = np.asarray(H)
        particles['hSPH'] = H
        tottime = time.time() - startTime
        print('Time FD: ',tottime)

        Area = dem['Area']
        # Update fields using a nearest interpolation
        MassNearest = np.zeros((NY, NX))
        # startTime = time.time()
        iC = particles['InCell']
        MassNearest = MassNearest.flatten()
        np.add.at(MassNearest, iC, Mpart)
        MassNearest = np.reshape(MassNearest, (NY, NX))
        FDNearest = MassNearest / (Area * rho)
        hNN, _ = geoTrans.projectOnRasterVectRoot(
            Xpart, Ypart, FDNearest, csz=csz, interp='nearest')
        hNB, _ = geoTrans.projectOnRasterVectRoot(
            Xpart, Ypart, FDNearest, csz=csz, interp='bilinear')

        # Update fields using a bilinear interpolation
        MassBilinear = np.zeros((NY, NX))
        #
        # # startTime = time.time()
        MassBilinear = geoTrans.pointsToRaster(
            Xpart, Ypart, Mpart, MassBilinear, csz=csz, interp='bilinear')
        FDBilinear = MassBilinear / (Area * rho)
        hBB, _ = geoTrans.projectOnRasterVectRoot(
            Xpart, Ypart, FDBilinear, csz=csz, interp='bilinear')

        h, Ghx, Ghy, Ghz = Hfunction(particles['x'], particles['y'], particles['z'])

        count = count + 1
        col = color[2*count]
        col3 = color[-1]
        col2 = color[round(len(color)/2)]
        mark = markers[count-1]
        ind = np.where(((particles['y'] > Ly/2-0.5*dy) & (particles['y'] < Ly/2+0.5*dy)))
        # if count == 1:
        #     ax1.plot(particles['x'][ind], hBB[ind], color='r',
        #              marker=mark, linestyle='None', label='HBB flow depth')
        #     ax1.plot(particles['x'][ind], hNN[ind], color='y',
        #              marker=mark, linestyle='None', label='HNN flow depth')
        #     ax1.plot(particles['x'][ind], hNB[ind], color='g',
        #              marker=mark, linestyle='None', label='HNB flow depth')
        # else:
        #     ax1.plot(particles['x'][ind], hBB[ind], color='r', marker=mark, linestyle='None')
        #     ax1.plot(particles['x'][ind], hNN[ind], color='y', marker=mark, linestyle='None')
        #     ax1.plot(particles['x'][ind], hNB[ind], color='g', marker=mark, linestyle='None')

        # , label='real flow depth')
        ax1.plot(particles['x'][ind], h[ind], color='r', marker='.', linestyle='None')
        ax1.plot(particles['x'][ind], particles['hSPH'][ind], color=col,
                 marker=mark, linestyle='None', label='SPH N=' + str(nPartPerD))

        ind = np.where(((particles['x'] > Lx/2-0.5*dx) & (particles['x'] < Lx/2+0.5*dx)))
        # if count == 1:
        #     ax2.plot(particles['y'][ind], hBB[ind], color='r',
        #              marker=mark, linestyle='None', label='HBB flow depth')
        #     ax2.plot(particles['y'][ind], hNN[ind], color='y',
        #              marker=mark, linestyle='None', label='HNN flow depth')
        #     ax2.plot(particles['y'][ind], hNB[ind], color='g',
        #              marker=mark, linestyle='None', label='HNB flow depth')
        # else:
        #     ax2.plot(particles['y'][ind], hBB[ind], color='r', marker=mark, linestyle='None')
        #     ax2.plot(particles['y'][ind], hNN[ind], color='y', marker=mark, linestyle='None')
        #     ax2.plot(particles['y'][ind], hNB[ind], color='g', marker=mark, linestyle='None')

        # , label='real flow depth')
        ax2.plot(particles['y'][ind], h[ind], color='r', marker='.', linestyle='None')
        ax2.plot(particles['y'][ind], particles['hSPH'][ind], color=col,
                 marker=mark, linestyle='None', label='SPH N=' + str(nPartPerD))

        ind = np.where(((particles['y'] > Ly/2-0.5*dy) & (particles['y'] < Ly/2+0.5*dy)))

        if count == 1:
            ax3.plot(particles['x'][ind], Ghx[ind], color='r', marker='.', linestyle='None' , label='real gradHX')
            ax3.plot(particles['x'][ind], Ghy[ind], color='g', marker='.', linestyle='None' , label='real gradHY')
            ax3.plot(particles['x'][ind], sx[ind]*Ghx[ind] + sy[ind]*Ghy[ind], color='k', marker='.', linestyle='None' , label='real gradHZ')
        else:
            ax3.plot(particles['x'][ind], Ghx[ind], color='r', marker='.', linestyle='None')  # , label='real gradH')
            ax3.plot(particles['x'][ind], Ghy[ind], color='g', marker='.', linestyle='None')  # , label='real gradH')
            ax3.plot(particles['x'][ind], sx[ind]*Ghx[ind] + sy[ind]*Ghy[ind], color='k', marker='.', linestyle='None')  # , label='real gradH')

        ax3.plot(particles['x'][ind], GHX[ind], color=col, marker=mark, markersize=5,
                 linestyle='None', label='SPH1 N=' + str(nPartPerD))
        ax3.plot(particles['x'][ind], GHY[ind], color=col, marker=mark, markersize=5,
                 linestyle='None')
        ax3.plot(particles['x'][ind], GHZ[ind], color=col, marker=mark, markersize=5,
                 linestyle='None')
        ax3.plot(particles['x'][ind], GHX2[ind], color=col2, marker=mark, markersize=3,
                 linestyle='None', label='SPH2 N=' + str(nPartPerD))
        ax3.plot(particles['x'][ind], GHY2[ind], color=col2, marker=mark, markersize=3,
                 linestyle='None')
        ax3.plot(particles['x'][ind], GHZ2[ind], color=col2, marker=mark, markersize=3,
                 linestyle='None')
        ax3.plot(particles['x'][ind], GHX3[ind], color=col3, marker=mark, markersize=1,
                 linestyle='None', label='SPH3 N=' + str(nPartPerD))
        ax3.plot(particles['x'][ind], GHY3[ind], color=col3, marker=mark, markersize=1,
                 linestyle='None')
        ax3.plot(particles['x'][ind], GHZ3[ind], color=col3, marker=mark, markersize=1,
                 linestyle='None')

        ind = np.where(((particles['x'] > Lx/2-0.5*dx) & (particles['x'] < Lx/2+0.5*dx)))
        if count == 1:
            ax4.plot(particles['y'][ind], Ghx[ind], color='r', marker='.', linestyle='None' , label='real gradHX')
            ax4.plot(particles['y'][ind], Ghy[ind], color='g', marker='.', linestyle='None' , label='real gradHY')
            ax4.plot(particles['y'][ind], sx[ind]*Ghx[ind] + sy[ind]*Ghy[ind], color='k', marker='.', linestyle='None' , label='real gradHZ')
        else:
            ax4.plot(particles['y'][ind], Ghx[ind], color='r', marker='.', linestyle='None')  # , label='real gradH')
            ax4.plot(particles['y'][ind], Ghy[ind], color='g', marker='.', linestyle='None')  # , label='real gradH')
            ax4.plot(particles['y'][ind], sx[ind]*Ghx[ind] + sy[ind]*Ghy[ind], color='k', marker='.', linestyle='None')  # , label='real gradH')
        ax4.plot(particles['y'][ind], GHX[ind], color=col, marker=mark, markersize=5,
                 linestyle='None', label='SPH1 N=' + str(nPartPerD))
        ax4.plot(particles['y'][ind], GHY[ind], color=col, marker=mark, markersize=5,
                 linestyle='None')
        ax4.plot(particles['y'][ind], GHZ[ind], color=col, marker=mark, markersize=5,
                 linestyle='None')
        ax4.plot(particles['y'][ind], GHX2[ind], color=col2, marker=mark, markersize=3,
                 linestyle='None', label='SPH2 N=' + str(nPartPerD))
        ax4.plot(particles['y'][ind], GHY2[ind], color=col2, marker=mark, markersize=3,
                 linestyle='None')
        ax4.plot(particles['y'][ind], GHZ2[ind], color=col2, marker=mark, markersize=3,
                 linestyle='None')
        plt.show(block=False)
        ax4.plot(particles['y'][ind], GHX3[ind], color=col3, marker=mark, markersize=1,
                 linestyle='None', label='SPH3 N=' + str(nPartPerD))
        ax4.plot(particles['y'][ind], GHY3[ind], color=col3, marker=mark, markersize=1,
                 linestyle='None')
        ax4.plot(particles['y'][ind], GHZ3[ind], color=col3, marker=mark, markersize=1,
                 linestyle='None')
        plt.show(block=False)
        # plt.draw()
        # plt.pause(2)

fig1.legend()
fig2.legend()
fig3.legend()
fig4.legend()
plt.show()
