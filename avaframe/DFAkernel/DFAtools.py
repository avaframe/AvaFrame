import logging
import numpy as np
from scipy.interpolate import griddata
import math
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.ascUtils as IOf
from avaframe.out3Plot.plotUtils import *
from avaframe.DFAkernel.setParam import *

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def polygon2Raster(demHeader, Line):
    # adim and center dem and polygon
    ncols = demHeader.ncols
    nrows = demHeader.nrows
    xllc = demHeader.xllcenter
    yllc = demHeader.yllcenter
    csz = demHeader.cellsize
    xCoord = (Line['x'] - xllc) / csz
    yCoord = (Line['y'] - yllc) / csz
    # get the raster corresponding to the polygon
    mask = geoTrans.poly2maskSimple(xCoord, yCoord, ncols, nrows)

    if debugPlot:
        fig, ax = plt.subplots(figsize=(figW, figH))
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        im = plt.imshow(mask, cmap, origin='lower')
        ax.plot(xCoord, yCoord, 'k')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        plt.show()

    return mask


def initializeSimulation(relRaster, dem):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    S = csz * csz
    # initialize
    partPerCell = np.zeros(np.shape(relRaster), dtype=np.int64)
    Npart = 0
    Xpart = np.empty(0)
    Ypart = np.empty(0)
    Mpart = np.empty(0)
    Hpart = np.empty(0)
    InCell = np.empty((0, 3), int)
    # find all non empty cells
    indY, indX = np.nonzero(relRaster)
    # loop on non empty cells
    for indx, indy in zip(indX, indY):
        # number of particles for this cell
        h = relRaster[indy][indx]
        V = S * h
        mass = V * rho
        nPart = np.ceil(mass / massPerPart).astype('int')
        Npart = Npart + nPart
        partPerCell[indy][indx] = nPart
        mPart = mass / nPart
        xpart = csz * (np.random.rand(nPart) - 0.5 + indx)
        ypart = csz * (np.random.rand(nPart) - 0.5 + indy)
        # if one particle in center of cell
        # xpart = csz * indx
        # ypart = csz * indy
        Xpart = np.append(Xpart, xpart)
        Ypart = np.append(Ypart, ypart)
        Mpart = np.append(Mpart, mPart * np.ones(nPart))
        Hpart = np.append(Hpart, h * np.ones(nPart))
        ic = indx + ncols * indy
        InCell = np.append(InCell, np.tile(np.array([indx, indy, ic]), (nPart, 1)), axis=0)

    # create dictionnary to store particles properties
    particles = {}
    particles['Npart'] = Npart
    particles['mTot'] = np.sum(Mpart)
    particles['x'] = Xpart
    particles['y'] = Ypart
    # adding z component
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')

    particles['m'] = Mpart
    particles['h'] = Hpart
    particles['InCell'] = InCell
    particles['ux'] = np.zeros(np.shape(Xpart))
    particles['uy'] = np.zeros(np.shape(Xpart))
    particles['uz'] = np.zeros(np.shape(Xpart))

    Cres = np.zeros(np.shape(dem['rasterData']))
    Ment = np.zeros(np.shape(dem['rasterData']))
    PV = np.zeros(np.shape(dem['rasterData']))
    PP = np.zeros(np.shape(dem['rasterData']))
    PFD = np.zeros(np.shape(dem['rasterData']))
    fields = {}
    fields['PV'] = PV
    fields['PP'] = PP
    fields['PFD'] = PFD
    fields['V'] = PV
    fields['P'] = PP
    fields['FD'] = PFD

    if debugPlot:
        x = np.arange(ncols) * csz
        y = np.arange(nrows) * csz
        fig, ax = plt.subplots(figsize=(figW, figH))
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        ref0, im = NonUnifIm(ax, x, y, relRaster, 'x [m]', 'y [m]',
                             extent=[x.min(), x.max(), y.min(), y.max()],
                             cmap=cmap, norm=None)
        ax.plot(Xpart, Ypart, 'or', linestyle='None')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        plt.show()

    return particles, fields, Cres, Ment


def computeForce(particles, dem, Ment, Cres):
    Npart = particles['Npart']
    csz = dem['header'].cellsize
    ncols = dem['header'].ncols
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    # initialize
    Fnormal = np.zeros(Npart)
    forceX = np.zeros(Npart)
    forceY = np.zeros(Npart)
    forceZ = np.zeros(Npart)
    dM = np.zeros(Npart)
    # loop on particles
    for j in range(Npart):
        mass = particles['m'][j]
        x = particles['x'][j]
        y = particles['y'][j]
        z = particles['z'][j]
        h = particles['h'][j]
        ux = particles['ux'][j]
        uy = particles['uy'][j]
        uz = particles['uz'][j]
        indCellX, indCellY, Cell = particles['InCell'][j]
        # deduce area
        A = mass / (h * rho)
        # get velocity magnitude and direction
        uMag = norm(ux, uy, uz)
        uxDir, uyDir, uzDir = normalize(ux, uy, uz)
        # get normal at the particle location
        nx, ny, nz = getNormal(x, y, Nx, Ny, Nz, csz)
        # get normal at the particle estimated end location
        xEnd = x + dt * ux
        yEnd = y + dt * uy
        nxEnd, nyEnd, nzEnd = getNormal(xEnd, yEnd, Nx, Ny, Nz, csz)
        # get average of those normals
        nxAvg = nx + nxEnd
        nyAvg = ny + nyEnd
        nzAvg = nz + nzEnd
        nxAvg, nyAvg, nzAvg = normalize(nxAvg, nyAvg, nzAvg)

        # acceleration due to curvature
        accNormCurv = (ux*(nxEnd-nx) + uy*(nyEnd-ny) + uz*(nzEnd-nz)) / dt
        # normal component of the acceleration of gravity
        gravAccNorm = - gravAcc * nzAvg
        effAccNorm = gravAccNorm + accNormCurv
        if(effAccNorm < 0.0):
            Fnormal[j] = mass * effAccNorm

        # body forces (tangential component of acceleration of gravity)
        gravAccTangX =          - gravAccNorm * nxAvg
        gravAccTangY =          - gravAccNorm * nyAvg
        gravAccTangZ = -gravAcc - gravAccNorm * nzAvg
        # adding gravity force contribution
        forceX[j] = forceX[j] + gravAccTangX * mass
        forceY[j] = forceY[j] + gravAccTangY * mass
        forceZ[j] = forceZ[j] + gravAccTangZ * mass

        # Calculating bottom sheer and normal stress
        if(effAccNorm > 0.0):
            # if fluid detatched
            log.info('fluid detatched for particle %s', j)
            tau = 0.0
        else:
            # bottom normal stress sigmaB
            sigmaB = - effAccNorm * rho * h
            # SamosAT friction type (bottom shear stress)
            tau = SamosATfric(uMag, sigmaB, h)
            # coulomb friction type (bottom shear stress)
            tau = mu * sigmaB

        # adding bottom shear resistance contribution
        forceBotTang = - A * tau
        forceX[j] = forceX[j] + forceBotTang * uxDir
        forceY[j] = forceY[j] + forceBotTang * uyDir
        forceZ[j] = forceZ[j] + forceBotTang * uzDir

        # compute entrained mass
        dm = 0
        if Ment[indCellY][indCellX] > 0:
            # either erosion or ploughing but not both
            # width of the particle
            width = math.sqrt(A)
            # bottom area covered by the particle during dt
            ABotSwiped = width * uMag * dt
            if(entEroEnergy > 0):
                # erosion: erode according to shear and erosion energy
                dm = A * tau * uMag * dt / entEroEnergy
                Aent = A
            else:
                # ploughing in at avalanche front: erode full area weight
                # mass available in the cell [kg/m²]
                rhoHent = Ment[indCellY][indCellX]
                dm = rhoHent * ABotSwiped
                Aent = rhoHent / rhoEnt
            dM[j] = dm
            # adding mass balance contribution
            forceX[j] = forceX[j] - dm / dt * ux
            forceY[j] = forceY[j] - dm / dt * uy
            forceZ[j] = forceZ[j] - dm / dt * uz

            # adding force du to entrained mass
            Fent = width * (entShearResistance + dm / Aent * entDefResistance)
            forceX[j] = forceX[j] + Fent * uxDir
            forceY[j] = forceY[j] + Fent * uyDir
            forceZ[j] = forceZ[j] + Fent * uzDir

        # adding resistance force du to obstacles
        if Cres[indCellY][indCellX] > 0:
            if(h < hRes):
                hResEff = h
            cres = - rho * A * hResEff * Cres * uMag
            forceX[j] = forceX[j] + cres * ux
            forceY[j] = forceY[j] + cres * uy
            forceZ[j] = forceZ[j] + cres * uz

        # adding lateral force (SPH component)
        gradhX, gradhY,  gradhZ = calcGradHSPH(particles, j, ncols)
        forceX[j] = forceX[j] - gradhX * mass * gravAcc / rho
        forceY[j] = forceY[j] - gradhY * mass * gravAcc / rho
        forceZ[j] = forceZ[j] - gradhZ * mass * gravAcc / rho

    force = {}
    force['dM'] = dM
    force['forceX'] = forceX
    force['forceY'] = forceY
    force['forceZ'] = forceZ

    return force


def updatePosition(particles, dem, force):
    csz = dem['header'].cellsize
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    dM = force['dM']
    forceX = force['forceX']
    forceY = force['forceY']
    forceZ = force['forceZ']
    mass = particles['m']
    x = particles['x']
    y = particles['y']
    z = particles['z']
    h = particles['h']
    ux = particles['ux']
    uy = particles['uy']
    uz = particles['uz']

    # procede to time integration
    # update velocity
    uxNew = ux + forceX * dt / mass
    uyNew = uy + forceY * dt / mass
    uzNew = uz + forceZ * dt / mass
    # update mass
    massNew = mass + dM
    # update position
    xNew = x + dt * 0.5 * (ux + uxNew)
    yNew = y + dt * 0.5 * (uy + uyNew)
    zNew = z + dt * 0.5 * (uz + uzNew)

    particles['mTot'] = np.sum(massNew)
    particles['x'] = xNew
    particles['y'] = yNew
    # make sure particle is on the mesh (recompute the z component)
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')

    nx, ny, nz = getNormalArray(xNew, yNew, Nx, Ny, Nz, csz)
    particles['m'] = massNew
    # normal component of the velocity
    uN = uxNew*nx + uyNew*ny + uzNew*nz
    # print(norm(ux, uy, uz), uN)
    # remove normal component of the velocity
    particles['ux'] = uxNew - uN * nx
    particles['uy'] = uyNew - uN * ny
    particles['uz'] = uzNew - uN * nz
    return particles


def updateFields(particles, dem, fields):
    header = dem['header']
    csz = dem['header'].cellsize
    S = csz * csz
    ncols = header.ncols
    nrows = header.nrows
    Npart = particles['Npart']
    m = particles['m']
    ux = particles['ux']
    uy = particles['uy']
    uz = particles['uz']
    PV = fields['PV']
    PP = fields['PP']
    PFD = fields['PFD']
    Mass = np.zeros((nrows, ncols))
    MomX = np.zeros((nrows, ncols))
    MomY = np.zeros((nrows, ncols))
    MomZ = np.zeros((nrows, ncols))
    # loop on particles
    for j in range(Npart):
        indx, indy, ic = particles['InCell'][j]
        Mass[indy, indx] = Mass[indy, indx] + m[j]
        MomX[indy, indx] = MomX[indy, indx] + m[j] * ux[j]
        MomY[indy, indx] = MomY[indy, indx] + m[j] * uy[j]
        MomZ[indy, indx] = MomZ[indy, indx] + m[j] * uz[j]

    VX = np.where(Mass > 0, MomX/Mass, MomX)
    VY = np.where(Mass > 0, MomY/Mass, MomY)
    VZ = np.where(Mass > 0, MomZ/Mass, MomZ)
    V = norm(VX, VY, VZ)
    FD = Mass / (S * rho)
    P = V * V * rho
    PV = np.where(V > PV, V, PV)
    PP = np.where(V > PP, P, PP)
    PFD = np.where(V > PFD, FD, PFD)

    fields['V'] = V
    fields['P'] = P
    fields['FD'] = FD

    return fields


def plotPosition(particles, dem, data, fig, ax):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    x = particles['x']
    y = particles['y']
    xx = np.arange(ncols) * csz
    yy = np.arange(nrows) * csz
    ax.clear()
    ax.set_title('Particles on dem at t=%.2f s' % particles['t'])
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    ref0, im = NonUnifIm(ax, xx, yy, data, 'x [m]', 'y [m]',
                         extent=[x.min(), x.max(), y.min(), y.max()],
                         cmap=cmap, norm=None)
    ax.plot(x, y, 'or', linestyle='None')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax)
    # ax.set_ylim([510, 530])
    # ax.set_xlim([260, 300])

    return fig, ax


def getNeighbours(particles, dem):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    Npart = particles['Npart']
    x = particles['x']
    y = particles['y']
    check = np.zeros((nrows, ncols))
    indPartInCell = np.zeros(ncols*nrows + 1)
    partInCell = np.zeros(Npart)
    InCell = np.empty((0, 3), int)
    # Count number of particles in each cell
    indx = ((x + csz/2) / csz).astype(int)
    indy = ((y + csz/2) / csz).astype(int)
    check[indy, indx] = check[indy, indx] + 1
    ic = indx + ncols * indy
    ##################################
    # TODO: test speed between add.at and bincount
    # indPartInCell = np.bincount(ic, minlength=len(indPartInCell))
    # or
    np.add.at(indPartInCell, ic, 1)
    ##################################
    indPartInCell = np.cumsum(indPartInCell)
    # make the list of which particles are in which cell
    indPartInCell2 = np.copy(indPartInCell)
    indPartInCell[-1] = 0
    for j in range(Npart):
        indx = int((x[j] + csz/2) / csz)
        indy = int((y[j] + csz/2) / csz)
        ic = indx + ncols * indy
        partInCell[int(indPartInCell2[ic])-1] = j
        indPartInCell2[ic] = indPartInCell2[ic] - 1
        InCell = np.append(InCell, np.tile(np.array([indx, indy, ic]), (1, 1)), axis=0)
    particles['InCell'] = InCell
    particles['indPartInCell'] = indPartInCell
    particles['partInCell'] = partInCell

    # xx = np.arange(ncols) * csz
    # yy = np.arange(nrows) * csz
    # fig, ax = plt.subplots(figsize=(figW, figH))
    # cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    # ref0, im = NonUnifIm(ax, xx, yy, check, 'x [m]', 'y [m]',
    #                      extent=[x.min(), x.max(), y.min(), y.max()],
    #                      cmap=cmap, norm=None)
    # ax.plot(x, y, 'or', linestyle='None')
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="5%", pad=0.1)
    # fig.colorbar(im, cax=cax)
    # plt.show()

    return particles


def calcGradHSPH(particles, j, ncols):
    indx, indy, _ = particles['InCell'][j]
    indPartInCell = particles['indPartInCell']
    partInCell = particles['partInCell']
    x = particles['x'][j]
    y = particles['y'][j]
    z = particles['z'][j]
    gradhX = 0
    gradhY = 0
    gradhZ = 0
    # find all the neighbour boxes
    for n in range(-1, 2):
        ic = (indx - 1) + ncols * (indy + n)
        iPstart = int(indPartInCell[ic-1])
        iPend = int(indPartInCell[ic+2])
        # loop on all particles in neighbour boxes
        for p in range(iPstart, iPend):
            # index of particle in neighbour box
            l = int(partInCell[p])
            if j != l:
                dx = particles['x'][l] - x
                dy = particles['y'][l] - y
                dz = particles['z'][l] - z
                r = norm(dx, dy, dz)
                if r < 0.001 * rKernel:
                    # impose a minimum distance between particles
                    r = 0.001 * rKernel
                if r < rKernel:
                    hr = rKernel - r
                    dwdr = dfacKernel * hr * hr
                    massl = particles['m'][l]
                    gradhX = gradhX + massl * dwdr * dx / r
                    gradhY = gradhY + massl * dwdr * dy / r
                    gradhZ = gradhZ + massl * dwdr * dz / r

    return gradhX, gradhY,  gradhZ


def getNormal(x, y, Nx, Ny, Nz, csz):
    nx = geoTrans.projectOnRasterRoot(x, y, Nx, csz=csz)
    ny = geoTrans.projectOnRasterRoot(x, y, Ny, csz=csz)
    nz = geoTrans.projectOnRasterRoot(x, y, Nz, csz=csz)
    nx, ny, nz = normalize(nx, ny, nz)
    return nx, ny, nz


def getNormalArray(x, y, Nx, Ny, Nz, csz):
    nrow, ncol = np.shape(Nx)
    nx, _ = geoTrans.projectOnRasterVectRoot(x, y, Nx, csz=csz)
    ny, _ = geoTrans.projectOnRasterVectRoot(x, y, Ny, csz=csz)
    nz, _ = geoTrans.projectOnRasterVectRoot(x, y, Nz, csz=csz)
    nx, ny, nz = normalize(nx, ny, nz)
    return nx, ny, nz


def getNormalVect(z, csz):
    n, m = np.shape(z)
    # first and last row, first and last column are inacurate
    # normal calculation with 4 triangles
    Nx = np.ones((n, m))
    # (Zl - Zr) * csz
    Nx[1:n-1, 1:m-1] = (z[1:n-1, 0:m-2] - z[1:n-1, 2:m]) / csz
    Ny = np.ones((n, m))
    # (Zd - Zu) * csz
    Ny[1:n-1, 1:m-1] = (z[0:n-2, 1:m-1] - z[2:n, 1:m-1]) / csz
    Nz = 2 * np.ones((n, m))

    # # normal calculation with 6 triangles
    # Nx = np.ones((n, m))
    # # (2*(Zl - Zr) - Zur + Zdl + Zu - Zd) * csz
    # Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m])
    #                     - z[2:n, 2:m] + z[0:n-2, 0:m-2]
    #                     + z[2:n, 1:m-1] - z[0:n-2, 1:m-1]) / csz
    # Ny = np.ones((n, m))
    # # (2*(Zd - Zu) + Zur + Zdl - Zu - Zl) * csz
    # Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1])
    #                     + z[2:n, 2:m] + z[0:n-2, 0:m-2]
    #                     - z[2:n, 1:m-1] - z[1:n-1, 0:m-2]) / csz
    # Nz = 6 * np.ones((n, m))

    Nx, Ny, Nz = normalize(Nx, Ny, Nz)

    return Nx, Ny, Nz


def norm(x, y, z):
    norme = np.sqrt(x*x + y*y + z*z)
    return norme


def normalize(x, y, z):
    norme = norm(x, y, z)
    xn = x / norme
    xn = np.where(np.isnan(xn), 0, xn)
    yn = y / norme
    yn = np.where(np.isnan(yn), 0, yn)
    zn = z / norme
    zn = np.where(np.isnan(zn), 0, zn)

    return xn, yn, zn


def SamosATfric(v, p, h):
    Rs = rho * v * v / (p + 0.001)
    div = h / R
    if(div < 1.0):
        div = 1.0
    div = math.log(div) / kappa + B
    tau = p * mu * (1.0 + Rs0 / (Rs0 + Rs)) + rho * v * v / (div * div)
    return tau
