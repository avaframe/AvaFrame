import logging
import time
import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
from avaframe.out3Plot.plotUtils import *
# import avaframe.in2Trans.shpConversion as shpConv
# import avaframe.in3Utils.ascUtils as IOf
# from avaframe.DFAkernel.setParam import *

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def initializeSimulation(cfg, relRaster, dem):
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    massPerPart = cfg.getfloat('massPerPart')
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
    particles['s'] = np.zeros(np.shape(Xpart))
    # adding z component
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')

    particles['m'] = Mpart
    particles['h'] = Hpart
    particles['InCell'] = InCell
    particles['ux'] = np.zeros(np.shape(Xpart))
    particles['uy'] = np.zeros(np.shape(Xpart))
    particles['uz'] = np.zeros(np.shape(Xpart))
    particles['stoppCriteria'] = False
    particles['kineticEne'] = np.sum(0.5 * Mpart * norm(particles['ux'], particles['uy'], particles['uz']))
    particles['potentialEne'] = np.sum(gravAcc * Mpart * particles['z'])

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

    # get particles location (neighbours for sph)
    particles = getNeighbours(particles, dem)
    # update fields (compute grid values)
    t = 0
    particles['t'] = t
    fields = updateFields(cfg, particles, dem, fields)
    # get normal vector of the grid mesh
    Nx, Ny, Nz = getNormalVect(dem['rasterData'], dem['header'].cellsize)
    dem['Nx'] = Nx
    dem['Ny'] = Ny
    dem['Nz'] = Nz

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

    return dem, particles, fields, Cres, Ment


def DFAIterate(cfg, particles, fields, dem, Ment, Cres, Tcpu):
    Tend = cfg.getfloat('Tend')
    dtSave = cfg.getfloat('dtSave')
    dt = cfg.getfloat('dt')

    Particles = [particles.copy()]
    Fields = [fields.copy()]
    nSave = 1
    niter = 0
    iterate = True
    t = particles['t']
    # Start time step computation
    while t < Tend and iterate:
        t = t + dt
        niter = niter + 1
        log.debug('Computing time step t = %f s', t)
        particles['t'] = t

        particles, fields, Tcpu = computeTimeStep(cfg, particles, fields, dem, Ment, Cres, Tcpu)
        iterate = not(particles['stoppCriteria'])
        if t >= nSave * dtSave:
            log.info('Saving results for time step t = %f s', t)
            Particles.append(particles.copy())
            Fields.append(fields.copy())
            nSave = nSave + 1

    Tcpu['niter'] = niter

    return Particles, Fields, Tcpu


def computeTimeStep(cfg, particles, fields, dem, Ment, Cres, Tcpu):
    # get forces
    startTime = time.time()
    force = computeForce(cfg, particles, dem, Ment, Cres)
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce
    # get forces sph
    # forceSPH = tools.computeForceSPH(particles, dem)
    # update velocity and particle position
    startTime = time.time()
    particles = updatePosition(cfg, particles, dem, force)
    tcpuPos = time.time() - startTime
    Tcpu['Pos'] = Tcpu['Pos'] + tcpuPos
    # get particles location (neighbours for sph)
    startTime = time.time()
    particles = getNeighbours(particles, dem)
    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh
    # update fields (compute grid values)
    startTime = time.time()
    fields = updateFields(cfg, particles, dem, fields)
    tcpuField = time.time() - startTime
    Tcpu['Field'] = Tcpu['Field'] + tcpuField

    return particles, fields, Tcpu


def polygon2Raster(demHeader, Line):
    # adim and center dem and polygon
    ncols = demHeader.ncols
    nrows = demHeader.nrows
    xllc = demHeader.xllcenter
    yllc = demHeader.yllcenter
    csz = demHeader.cellsize
    xCoord = (Line['x'] - xllc) / csz
    xCoord0 = np.delete(xCoord, -1)
    yCoord = (Line['y'] - yllc) / csz
    yCoord0 = np.delete(yCoord, -1)
    # get the raster corresponding to the polygon
    mask = geoTrans.poly2maskSimple(xCoord0, yCoord0, ncols, nrows)

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


def computeForce(cfg, particles, dem, Ment, Cres):
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    dt = cfg.getfloat('dt')
    mu = cfg.getfloat('mu')
    entEroEnergy = cfg.getfloat('entEroEnergy')
    rhoEnt = cfg.getfloat('rhoEnt')
    entShearResistance = cfg.getfloat('entShearResistance')
    entDefResistance = cfg.getfloat('entDefResistance')
    hRes = cfg.getfloat('hRes')
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
    TcpuSPH = 0
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

        # Calculating bottom shear and normal stress
        if(effAccNorm > 0.0):
            # if fluid detatched
            log.info('fluid detatched for particle %s', j)
            tau = 0.0
        else:
            # bottom normal stress sigmaB
            sigmaB = - effAccNorm * rho * h
            # SamosAT friction type (bottom shear stress)
            tau = SamosATfric(cfg, uMag, sigmaB, h)
            # coulomb friction type (bottom shear stress)
            # tau = mu * sigmaB

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
        startTime = time.time()
        gradhX, gradhY,  gradhZ = calcGradHSPH(particles, j, ncols)
        tcpuSPH = time.time() - startTime
        TcpuSPH = TcpuSPH + tcpuSPH
        forceY[j] = forceY[j] - gradhY * mass * (-gravAcc) / rho
        forceX[j] = forceX[j] - gradhX * mass * (-gravAcc) / rho
        forceZ[j] = forceZ[j] - gradhZ * mass * (-gravAcc) / rho

    # log.info(('cpu time SPH = %s s' % (TcpuSPH)))
    force = {}
    force['dM'] = dM
    force['forceX'] = forceX
    force['forceY'] = forceY
    force['forceZ'] = forceZ

    return force


def updatePosition(cfg, particles, dem, force):
    dt = cfg.getfloat('dt')
    gravAcc = cfg.getfloat('gravAcc')
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
    kinEne = particles['kineticEne']
    potEne = particles['potentialEne']
    totEne = kinEne + potEne
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
    particles['s'] = particles['s'] + np.sqrt((xNew-x)*(xNew-x) + (yNew-y)*(yNew-y))
    # make sure particle is on the mesh (recompute the z component)
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')

    nx, ny, nz = getNormalArray(xNew, yNew, Nx, Ny, Nz, csz)
    particles['m'] = massNew
    # normal component of the velocity
    uN = uxNew*nx + uyNew*ny + uzNew*nz
    # print(nx, ny, nz)
    # print(norm(ux, uy, uz), uN)
    # remove normal component of the velocity
    particles['ux'] = uxNew - uN * nx
    particles['uy'] = uyNew - uN * ny
    particles['uz'] = uzNew - uN * nz
    # uN = particles['ux']*nx + particles['uy']*ny + particles['uz']*nz
    # print(norm(particles['ux'], particles['uy'], particles['uz']), uN)
    kinEneNew = np.sum(0.5 * massNew * norm(particles['ux'], particles['uy'], particles['uz']))
    potEneNew = np.sum(gravAcc * massNew * particles['z'])
    totEneNew = kinEneNew + potEneNew
    # log.info('total energy variation: %f' % ((totEneNew - totEne) / totEneNew))
    # log.info('kinetic energy variation: %f' % ((kinEneNew - kinEne) / kinEneNew))
    # if (abs(totEneNew - totEne) / totEneNew < 0.01) and (abs(kinEneNew - kinEne) / kinEneNew < 0.01):
    #     log.info('Reached stopping creteria : total energy varied fom less than 1 %')
    #     particles['stoppCriteria'] = True
    particles['kineticEne'] = kinEneNew
    particles['potentialEne'] = potEneNew
    return particles


def updateFields(cfg, particles, dem, fields):
    rho = cfg.getfloat('rho')
    header = dem['header']
    csz = dem['header'].cellsize
    S = csz * csz
    ncols = header.ncols
    nrows = header.nrows
    Npart = particles['Npart']
    m = particles['m']
    x = particles['x']
    y = particles['y']
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
    # startTime = time.time()
    # Mass = geoTrans.pointsToRaster(x, y, m, Mass, csz=csz, interp='bilinear')
    # MomX = geoTrans.pointsToRaster(x, y, m * ux, MomX, csz=csz, interp='bilinear')
    # MomY = geoTrans.pointsToRaster(x, y, m * uy, MomY, csz=csz, interp='bilinear')
    # MomZ = geoTrans.pointsToRaster(x, y, m * uz, MomZ, csz=csz, interp='bilinear')
    # endTime = time.time()
    # log.info(('time = %s s' % (endTime - startTime)))

    # startTime = time.time()
    iC = particles['InCell'][:, 2]
    Mass = Mass.flatten()
    np.add.at(Mass, iC, m)
    MomX = MomX.flatten()
    np.add.at(MomX, iC, m * ux)
    MomY = MomY.flatten()
    np.add.at(MomY, iC, m * uy)
    MomZ = MomZ.flatten()
    np.add.at(MomZ, iC, m * uz)
    Mass = np.reshape(Mass, (nrows, ncols))
    MomX = np.reshape(MomX, (nrows, ncols))
    MomY = np.reshape(MomY, (nrows, ncols))
    MomZ = np.reshape(MomZ, (nrows, ncols))

    # same with a loop
    # # loop on particles
    # for j in range(Npart):
    #     indx, indy, ic = particles['InCell'][j]
    #
    #     Mass[indy, indx] = Mass[indy, indx] + m[j]
    #     MomX[indy, indx] = MomX[indy, indx] + m[j] * ux[j]
    #     MomY[indy, indx] = MomY[indy, indx] + m[j] * uy[j]
    #     MomZ[indy, indx] = MomZ[indy, indx] + m[j] * uz[j]

    # endTime = time.time()
    # log.info(('time = %s s' % (endTime - startTime)))

    # print(np.sum(Mass))
    # print(np.sum(Mass2))
    # plotPosition(particles, dem, Mass)
    # plotPosition(particles, dem, Mass2)

    VX = np.where(Mass > 0, MomX/Mass, MomX)
    VY = np.where(Mass > 0, MomY/Mass, MomY)
    VZ = np.where(Mass > 0, MomZ/Mass, MomZ)
    V = norm(VX, VY, VZ)
    FD = Mass / (S * rho)
    P = V * V * rho
    PV = np.where(V > PV, V, PV)
    PP = np.where(P > PP, P, PP)
    PFD = np.where(FD > PFD, FD, PFD)

    fields['V'] = V
    fields['P'] = P
    fields['FD'] = FD
    fields['PV'] = PV
    fields['PP'] = PP
    fields['PFD'] = PFD

    return fields


def plotPosition(particles, dem, data, Cmap, fig, ax):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    X = PointsX[0, :]
    Y = PointsY[:, 0]
    Z = dem['rasterData']
    x = particles['x'] + xllc
    y = particles['y'] + yllc
    xx = np.arange(ncols) * csz + xllc
    yy = np.arange(nrows) * csz + yllc
    ax.clear()
    ax.set_title('Particles on dem at t=%.2f s' % particles['t'])
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        Cmap, 0.0, np.nanmax(data), continuous=contCmap)
    cmap.set_under(color='w')
    ref0, im = NonUnifIm(ax, xx, yy, data, 'x [m]', 'y [m]',
                         extent=[x.min(), x.max(), y.min(), y.max()],
                         cmap=cmap, norm=norm)
    ax.plot(x, y, 'or', linestyle='None')
    Cp1 = ax.contour(X, Y, Z, levels=10, colors='k')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax)
    plt.pause(0.1)
    # plt.close(fig)
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
    # SPH kernel
    # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
    rKernel = 5
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))
    dfacKernel = -3.0 * facKernel

    indx, indy, _ = particles['InCell'][j]
    indPartInCell = particles['indPartInCell']
    partInCell = particles['partInCell']
    x = particles['x'][j]
    y = particles['y'][j]
    z = particles['z'][j]
    gradhX = 0
    gradhY = 0
    gradhZ = 0
    # gradhX1 = 0
    # gradhY1 = 0
    # gradhZ1 = 0
    index = np.empty((0), dtype=int)
    # find all the neighbour boxes
    for n in range(-1, 2):
        ic = (indx - 1) + ncols * (indy + n)
        iPstart = int(indPartInCell[ic-1])
        iPend = int(indPartInCell[ic+2])
        ind = np.arange(iPstart, iPend, 1)
        index = np.append(index, ind)
        # # loop on all particles in neighbour boxes
        # for p in range(iPstart, iPend):
        #     # index of particle in neighbour box
        #     l = int(partInCell[p])
        #     if j != l:
        #         dx = particles['x'][l] - x
        #         dy = particles['y'][l] - y
        #         dz = particles['z'][l] - z
        #         r = norm(dx, dy, dz)
        #         if r < 0.001 * rKernel:
        #             # impose a minimum distance between particles
        #             r = 0.001 * rKernel
        #         if r < rKernel:
        #             hr = rKernel - r
        #             dwdr = dfacKernel * hr * hr
        #             massl = particles['m'][l]
        #             gradhX1 = gradhX1 + massl * dwdr * dx / r
        #             gradhY1 = gradhY1 + massl * dwdr * dy / r
        #             gradhZ1 = gradhZ1 + massl * dwdr * dz / r
    # no loop option
    L = partInCell[index].astype('int')
    dx = particles['x'][L] - x
    dy = particles['y'][L] - y
    dz = particles['z'][L] - z
    r = norm(dx, dy, dz)
    r = np.where(r < 0.001 * rKernel, 0.001 * rKernel, r)
    hr = rKernel - r
    dwdr = dfacKernel * hr * hr
    massl = particles['m'][L]
    GX = massl * dwdr * dx / r
    GX = np.where(r < rKernel, GX, 0)
    GY = massl * dwdr * dy / r
    GY = np.where(r < rKernel, GY, 0)
    GZ = massl * dwdr * dz / r
    GZ = np.where(r < rKernel, GZ, 0)
    gradhX = gradhX + np.sum(GX)
    gradhY = gradhY + np.sum(GY)
    gradhZ = gradhZ + np.sum(GZ)
    # print(gradhX, gradhX1)
    # print(gradhY, gradhY1)
    # print(gradhZ, gradhZ1)

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


def SamosATfric(cfg, v, p, h):
    rho = cfg.getfloat('rho')
    Rs0 = cfg.getfloat('Rs0')
    mu = cfg.getfloat('mu')
    kappa = cfg.getfloat('kappa')
    B = cfg.getfloat('B')
    R = cfg.getfloat('R')
    Rs = rho * v * v / (p + 0.001)
    div = h / R
    if(div < 1.0):
        div = 1.0
    div = math.log(div) / kappa + B
    tau = p * mu * (1.0 + Rs0 / (Rs0 + Rs)) + rho * v * v / (div * div)
    return tau
