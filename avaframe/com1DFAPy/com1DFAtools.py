"""

    This file is part of Avaframe.

"""

import logging
import time
import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
from avaframe.out3Plot.plotUtils import *
import avaframe.com1DFAPy.timeDiscretizations as tD
# import avaframe.in2Trans.shpConversion as shpConv
# import avaframe.in2Trans.ascUtils as IOf
# from avaframe.DFAkernel.setParam import *

# create local logger
log = logging.getLogger(__name__)
debugPlot = False

# set feature leapfrog time stepping
featLF = True


def initializeMesh(dem):
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    # get normal vector of the grid mesh
    Nx, Ny, Nz = getNormalMesh(dem['rasterData'], csz)
    dem['Nx'] = Nx
    dem['Ny'] = Ny
    dem['Nz'] = Nz

    # get real Area
    Area = getAreaMesh(Nx, Ny, Nz, csz)
    dem['Area'] = Area

    return dem


def initializeSimulation(cfg, relRaster, dem):
    # get simulation parameters
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    massPerPart = cfg.getfloat('massPerPart')
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    A = dem['Area']
    # initialize arrays
    partPerCell = np.zeros(np.shape(relRaster), dtype=np.int64)
    FD = np.zeros((nrows, ncols))
    Npart = 0
    Xpart = np.empty(0)
    Ypart = np.empty(0)
    Mpart = np.empty(0)
    Hpart = np.empty(0)
    InCell = np.empty((0, 3), int)
    # find all non empty cells (meaning release area)
    indY, indX = np.nonzero(relRaster)
    # loop on non empty cells
    for indx, indy in zip(indX, indY):
        # number of particles for this cell
        h = relRaster[indy, indx]
        Vol = A[indy, indx] * h
        mass = Vol * rho
        nPart = np.ceil(mass / massPerPart).astype('int')
        Npart = Npart + nPart
        partPerCell[indy, indx] = nPart
        mPart = mass / nPart
        # TODO make this an independent function
        #######################
        # start ###############
        # place particles randomly in the cell
        xpart = csz * (np.random.rand(nPart) - 0.5 + indx)
        ypart = csz * (np.random.rand(nPart) - 0.5 + indy)
        # if one particle in center of cell
        # xpart = csz * indx
        # ypart = csz * indy
        # end #################
        #######################
        # initialize field Flow depth
        FD[indY, indX] = h
        # initialize particles position mass height...
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
    particles['hNearestNearest'] = Hpart
    particles['hNearestBilinear'] = Hpart
    particles['hBilinearNearest'] = Hpart
    particles['hBilinearBilinear'] = Hpart
    particles['InCell'] = InCell
    particles['ux'] = np.zeros(np.shape(Xpart))
    particles['uy'] = np.zeros(np.shape(Xpart))
    particles['uz'] = np.zeros(np.shape(Xpart))
    particles['stoppCriteria'] = False
    particles['kineticEne'] = np.sum(0.5 * Mpart * norm2(particles['ux'], particles['uy'], particles['uz']))
    particles['potentialEne'] = np.sum(gravAcc * Mpart * particles['z'])

    # initialize entrainment and resistance
    Ment = intializeMassEnt(dem)
    Cres = intializeResistance(dem)

    PV = np.zeros((nrows, ncols))
    PP = np.zeros((nrows, ncols))
    fields = {}
    fields['pv'] = PV
    fields['ppr'] = PP
    fields['pfd'] = FD
    fields['V'] = PV
    fields['P'] = PP
    fields['FD'] = FD

    # get particles location (neighbours for sph)
    # particles = getNeighbours(particles, dem)
    particles = getNeighboursVect(particles, dem)
    # initialize time
    t = 0
    particles['t'] = t

    log.info('Initializted simulation. MTot = %f kg, %s particles in %s cells' % (particles['mTot'], particles['Npart'], np.size(indY)))

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


def intializeMassEnt(dem):
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    Ment = np.zeros((nrows, ncols))
    return Ment


def intializeResistance(dem):
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    Cres = np.zeros((nrows, ncols))
    return Cres


def DFAIterate(cfg, particles, fields, dem, Ment, Cres, Tcpu):
    """ Perform computations and save results for desired interval """

    # Load configuration settings
    Tend = cfg.getfloat('Tend')
    dtSave = cfg.getfloat('dtSave')

    # Initialise Lists to save fields
    Particles = [copy.deepcopy(particles)]
    Fields = [copy.deepcopy(fields)]
    # save Z, S, U, T at each time step for developping purpouses
    Z = np.empty((0, 0))
    S = np.empty((0, 0))
    U = np.empty((0, 0))
    T = np.empty((0, 0))

    # Initialize time and counters
    nSave = 1
    nIter = 0
    iterate = True
    particles['iterate'] = iterate
    t = particles['t']
    # Start time step computation
    while t < Tend and iterate:

        #++++++++++++++++if you want to use cfl time step+++++++++++++++++++
        # CALL TIME STEP:
        # to play around with the courant number
        if featLF:
            cfg['cMax'] = '0.5'
            dtStable = tD.getcflTimeStep(particles, dem, cfg)
        # dt overwrites dt in .ini file, so comment this block if you dont want to use cfl
        #++++++++++++++++++++++++++++++++++++++++++++++
        # get time step
        dt = cfg.getfloat('dt')
        t = t + dt
        nIter = nIter + 1
        log.debug('Computing time step t = %f s', t)
        T = np.append(T, t)
        particles['t'] = t
        Tcpu['nSave'] = nSave

        # Perform computations
        if featLF:
            particles, fields, Tcpu, dt = computeLeapFrogTimeStep(cfg, particles, fields, dt, dem, Ment, Cres, Tcpu)
        else:
            particles, fields, Tcpu = computeTimeStep(cfg, particles, fields, dt, dem, Ment, Cres, Tcpu)
        # Save desired parameters and export to Lists for saving interval
        U = np.append(U, norm(particles['ux'][0], particles['uy'][0], particles['uz'][0]))
        Z = np.append(Z, particles['z'][0])
        S = np.append(S, particles['s'][0])
        iterate = particles['iterate']
        if t >= nSave * dtSave:
            log.info('Saving results for time step t = %f s', t)
            log.info('MTot = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
            Particles.append(copy.deepcopy(particles))
            Fields.append(copy.deepcopy(fields))
            nSave = nSave + 1

    Tcpu['nIter'] = nIter

    return T, U, Z, S, Particles, Fields, Tcpu


def computeTimeStep(cfg, particles, fields, dt, dem, Ment, Cres, Tcpu):
    log.info('Use standard time stepping')
    # get forces
    # loop version of the compute force
    startTime = time.time()
    # forceLoop = computeForce(cfg, particles, dem, Ment, Cres)
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce
    # vectorized version of the compute force
    startTime = time.time()
    force = computeForceVect(cfg, particles, dem, Ment, Cres, dt)
    tcpuForceVect = time.time() - startTime
    Tcpu['ForceVect'] = Tcpu['ForceVect'] + tcpuForceVect
    # compare output from compute force loop and vectorized
    #########################################################
    # print(np.max(np.abs(force['forceX']-forceVect['forceX'])))
    # print(np.max(np.abs(force['forceY']-forceVect['forceY'])))
    # print(np.max(np.abs(force['forceZ']-forceVect['forceZ'])))
    # print(np.allclose(force['forceX'], forceVect['forceX'], atol=0.00001))
    ##########################################################

    # compute lateral force (SPH component of the calculation)
    startTime = time.time()
    particles, force = computeForceSPH(cfg, particles, force, dem)
    tcpuForceSPH = time.time() - startTime
    Tcpu['ForceSPH'] = Tcpu['ForceSPH'] + tcpuForceSPH

    # plot depth computed with different interpolation methods
    nSave = Tcpu['nSave']
    dtSave = cfg.getfloat('dtSave')
    if debugPlot and particles['t'] >= nSave * dtSave:
        hNN = copy.deepcopy(particles['hNearestNearest'])
        hNB = copy.deepcopy(particles['hNearestBilinear'])
        hSPH = copy.deepcopy(particles['hSPH'])
        hBN = copy.deepcopy(particles['hBilinearNearest'])
        hBB = copy.deepcopy(particles['hBilinearBilinear'])
        indexSort = np.argsort(hNN)
        indexSortBB = np.argsort(hBB)

        # print(np.mean(hSPH))
        # fig, ax = plt.subplots(figsize=(2*figW, figH))
        # ax.set_title('flow depth per particle')
        # # ax.plot(hSPH[indexSort], 'b', linewidth=0.5, label='h from SPH')
        # ax.plot(hBB[indexSort], 'y', linewidth=0.5, label='h from bilinear bilinear')
        # # ax.plot(hBN[indexSort], 'g', linewidth=0.5, label='h from bilinear nearest')
        # ax.plot(hNB[indexSort], 'r', linewidth=0.5, label='h from nearest bilinear')
        # ax.plot(hNN[indexSort], 'k', linewidth=0.5, label='h from nearest nearest')
        # plt.legend()
        # plt.show()#block=False)
        # plt.pause(1)
        # plt.close()
        ind = np.where(((particles['y']>905) & (particles['y']<1005)))
        # fig2 = plt.figure()
        # ax2 = fig2.add_subplot(111, projection='3d')
        # ax2.scatter(particles['x'][ind], particles['y'][ind], hNN[ind], 'k')
        # ax2.scatter(particles['x'][ind], particles['y'][ind], hBB[ind], 'b')
        # ax2.set_xlim3d(0, 1000)
        # ax2.set_ylim3d(950,1050)
        # ax2.set_zlim3d(0,1000)

        fig1, ax1 = plt.subplots(figsize=(2*figW, figH))
        ax1.plot(particles['x'][ind], hNN[ind], color='k', marker='.', linestyle = 'None')
        ax1.plot(particles['x'][ind], hBB[ind], color='b', marker='.', linestyle = 'None')
        plt.show()

    # update velocity and particle position
    startTime = time.time()
    particles = updatePosition(cfg, particles, dem, force)
    tcpuPos = time.time() - startTime
    Tcpu['Pos'] = Tcpu['Pos'] + tcpuPos

    # get particles location (neighbours for sph)
    startTime = time.time()
    # particles = getNeighbours(particles, dem)
    particles = getNeighboursVect(particles, dem)
    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

    # update fields (compute grid values)
    startTime = time.time()
    particles, fields = updateFields(cfg, particles, force, dem, fields)
    tcpuField = time.time() - startTime
    Tcpu['Field'] = Tcpu['Field'] + tcpuField

    return particles, fields, Tcpu


def computeLeapFrogTimeStep(cfg, particles, fields, dt, dem, Ment, Cres, Tcpu):
    """ perform all computations that belong to one time step """

    log.info('Use LeapFrog time stepping')
    # start timing
    startTime = time.time()
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce

    # dtK5 is half time step
    dtK5 = 0.5 * dt
    log.info('dt used now is %f' % dt)

    # load required DEM and mesh info
    csz = dem['header'].cellsize
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']

    # particle properties
    mass = particles['m']
    xK = particles['x']
    yK = particles['y']
    zK = particles['z']
    uxK = particles['ux']
    uyK = particles['uy']
    uzK = particles['uz']

    #+++++++++++++Time integration using leapfrog 'Drift-Kick-Drif' scheme+++++
    # first predict position at time t_(k+0.5)
    # 'DRIFT'
    xK5 = xK + dt * 0.5 * uxK
    yK5 = yK + dt * 0.5 * uyK
    zK5 = zK + dt * 0.5 * uzK
    # update position from particles
    particles['x'] = xK5
    particles['y'] = yK5
    # For now z-position is taken from DEM - no detachment enforces...
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')
    # TODO: do we need to update also h from particles?? I think yes! also mass, ent, res
    ##particles['h'] = ?

    # 'KICK'
    # compute velocity at t_(k+0.5)
    # first compute force at t_(k+0.5)
    force = computeForceVect(cfg, particles, dem, Ment, Cres, dtK5)
    particles, force = computeForceSPH(cfg, particles, force, dem)
    mass = particles['m']
    uxNew = uxK + (force['forceX'] + force['forceSPHX']) * dt / mass
    uyNew = uyK + (force['forceY'] + force['forceSPHY']) * dt  / mass
    uzNew = uzK + (force['forceZ'] + force['forceSPHZ']) * dt / mass

    # 'DRIF'
    # now update position at t_(k+ 1)
    xNew = xK5 + dtK5 * uxNew
    yNew = yK5 + dtK5 * uyNew
    zNew = zK5 + dtK5 * uzNew

    # Compute magnitude of velocity
    # print anormally big values
    uNew = norm(uxNew, uyNew, uzNew)
    if np.max(uNew) > 100:
        ind = np.argmax(uNew)
        A = mass[ind] / (h[ind] * rho)
        print('Normal component of velocity: particle index : ', ind)
        print((forceX / mass)[ind])
        print((forceY / mass)[ind])
        print((forceZ / mass)[ind])
        print(uNew[ind])
        print(A)


    #++++++++++++++UPDATE Particle Properties
    # update mass required if entrainment
    massNew = mass + force['dM']
    particles['mTot'] = np.sum(massNew)
    particles['x'] = xNew
    particles['y'] = yNew
    particles['s'] = particles['s'] + np.sqrt((xNew-xK)*(xNew-xK) + (yNew-yK)*(yNew-yK))
    # make sure particle is on the mesh (recompute the z component)
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')

    nx, ny, nz = getNormalArray(xNew, yNew, Nx, Ny, Nz, csz)
    particles['m'] = massNew
    # normal component of the velocity
    uN = uxNew*nx + uyNew*ny + uzNew*nz
    # remove normal component of the velocity
    particles['ux'] = uxNew - uN * nx
    particles['uy'] = uyNew - uN * ny
    particles['uz'] = uzNew - uN * nz

    #################################################################
    # this is dangerous!!!!!!!!!!!!!!
    ###############################################################
    # remove particles that are not located on the mesh any more
    particles = removeOutPart(cfg, particles, dem)

    #++++++++++++++GET particles location (neighbours for sph)
    startTime = time.time()
    particles = getNeighboursVect(particles, dem)
    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

    #++++++++++++++UPDATE FIELDS (compute grid values)
    startTime = time.time()
    particles, fields = updateFields(cfg, particles, force, dem, fields)
    tcpuField = time.time() - startTime
    Tcpu['Field'] = Tcpu['Field'] + tcpuField

    return particles, fields, Tcpu, dt


def prepareArea(releaseLine, dem):
    NameRel = releaseLine['Name']
    StartRel = releaseLine['Start']
    LengthRel = releaseLine['Length']
    relRaster = np.zeros(np.shape(dem['rasterData']))

    for i in range(len(NameRel)):
        name = NameRel[i]
        start = StartRel[i]
        end = start + LengthRel[i]
        avapath = {}
        avapath['x'] = releaseLine['x'][int(start):int(end)]
        avapath['y'] = releaseLine['y'][int(start):int(end)]
        avapath['Name'] = name
        relRaster = polygon2Raster(dem['header'], avapath, relRaster)
    return relRaster


def polygon2Raster(demHeader, Line, Mask):
    # adim and center dem and polygon
    ncols = demHeader.ncols
    nrows = demHeader.nrows
    xllc = demHeader.xllcenter
    yllc = demHeader.yllcenter
    csz = demHeader.cellsize
    xCoord0 = (Line['x'] - xllc) / csz
    yCoord0 = (Line['y'] - yllc) / csz
    if (xCoord0[0] == xCoord0[-1]) and (yCoord0[0] == yCoord0[-1]):
        xCoord = np.delete(xCoord0, -1)
        yCoord = np.delete(yCoord0, -1)
    else:
        xCoord = copy.deepcopy(xCoord0)
        yCoord = copy.deepcopy(yCoord0)
        xCoord0 = np.append(xCoord0, xCoord0[0])
        yCoord0 = np.append(yCoord0, yCoord0[0])

    # get the raster corresponding to the polygon
    mask = geoTrans.poly2maskSimple(xCoord, yCoord, ncols, nrows)
    Mask = Mask + mask
    Mask = np.where(Mask > 0, 1, 0)

    if debugPlot:
        fig, ax = plt.subplots(figsize=(figW, figH))
        ax.set_title('Release area')
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        im = plt.imshow(Mask, cmap, origin='lower')
        ax.plot(xCoord0, yCoord0, 'r', label='release polyline')
        plt.legend()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        plt.show()

    return Mask


def computeForce(cfg, particles, dem, Ment, Cres):
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    dt = cfg.getfloat('dt')
    mu = cfg.getfloat('mu')
    Npart = particles['Npart']
    csz = dem['header'].cellsize
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    # initialize
    Fnormal = np.zeros(Npart)
    forceX = np.zeros(Npart)
    forceY = np.zeros(Npart)
    forceZ = np.zeros(Npart)
    forceSPHX = np.zeros(Npart)
    forceSPHY = np.zeros(Npart)
    forceSPHZ = np.zeros(Npart)
    dM = np.zeros(Npart)
    force = {}
    # loop on particles
    for j in range(Npart):
        mass = particles['m'][j]
        x = particles['x'][j]
        y = particles['y'][j]
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
        ment = Ment[indCellY][indCellX]
        if ment > 0:
            dm, fEntX, fEntY, fEntZ = computeEntMassAndForce(cfg, ment, A, uMag, ux, uy, uz, uxDir, uyDir, uzDir, tau)
            dM[j] = dm
            forceX[j] = forceX[j] + fEntX
            forceY[j] = forceY[j] + fEntY
            forceZ[j] = forceZ[j] + fEntZ

        # adding resistance force du to obstacles
        cres = Cres[indCellY][indCellX]
        if cres > 0:
            fEntX, fEntY, fEntZ = computeResForce(cfg, h, A, rho, cres, uMag, ux, uy, uz)
            forceX[j] = forceX[j] + fEntX
            forceY[j] = forceY[j] + fEntY
            forceZ[j] = forceZ[j] + fEntZ

    # save results
    force['dM'] = dM
    force['forceX'] = forceX
    force['forceY'] = forceY
    force['forceZ'] = forceZ
    force['forceSPHX'] = forceSPHX
    force['forceSPHY'] = forceSPHY
    force['forceSPHZ'] = forceSPHZ
    return force


def computeEntMassAndForce(cfg, ment, A, uMag, ux, uy, uz, uxDir, uyDir, uzDir, tau):
    dt = cfg.getfloat('dt')
    entEroEnergy = cfg.getfloat('entEroEnergy')
    rhoEnt = cfg.getfloat('rhoEnt')
    entShearResistance = cfg.getfloat('entShearResistance')
    entDefResistance = cfg.getfloat('entDefResistance')
    # compute entrained mass
    dm = 0
    fEntX, fEntY, fEntZ = 0, 0, 0
    if ment > 0:
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
            rhoHent = ment
            dm = rhoHent * ABotSwiped
            Aent = rhoHent / rhoEnt
        # adding mass balance contribution
        fEntX = fEntX - dm / dt * ux
        fEntY = fEntY - dm / dt * uy
        fEntZ = fEntZ - dm / dt * uz

        # adding force du to entrained mass
        Fent = width * (entShearResistance + dm / Aent * entDefResistance)
        fEntX = fEntX + Fent * uxDir
        fEntY = fEntY + Fent * uyDir
        fEntZ = fEntZ + Fent * uzDir

    return dm, fEntX, fEntY, fEntZ


def computeResForce(cfg, h, A, rho, cres, uMag, ux, uy, uz):
    hRes = cfg.getfloat('hRes')
    if(h < hRes):
        hResEff = h
    cres = - rho * A * hResEff * cres * uMag
    fResX, fResY, fResZ = cres*ux, cres*uy, cres*uz
    return fResX, fResY, fResZ


def computeForceVect(cfg, particles, dem, Ment, Cres, dt):
    """ Compute forces """

    # Load required parameters
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    mu = cfg.getfloat('mu')
    Npart = particles['Npart']
    csz = dem['header'].cellsize
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']

    # initialize
    Fnormal = np.zeros(Npart)
    forceX = np.zeros(Npart)
    forceY = np.zeros(Npart)
    forceZ = np.zeros(Npart)
    dM = np.zeros(Npart)
    force = {}

    # loop on particles
    mass = particles['m']
    x = particles['x']
    y = particles['y']
    h = particles['h']
    ux = particles['ux']
    uy = particles['uy']
    uz = particles['uz']
    # deduce area
    A = mass / (h * rho)
    # get velocity magnitude and direction
    uMag = norm(ux, uy, uz)
    uxDir, uyDir, uzDir = normalize(ux, uy, uz)
    # get normal at the particle location
    nx, ny, nz = getNormalArray(x, y, Nx, Ny, Nz, csz)
    # get normal at the particle estimated end location
    xEnd = x + dt * ux
    yEnd = y + dt * uy
    nxEnd, nyEnd, nzEnd = getNormalArray(xEnd, yEnd, Nx, Ny, Nz, csz)
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
    Fnormal = np.where(effAccNorm < 0.0, mass * effAccNorm, 0)

    # body forces (tangential component of acceleration of gravity)
    gravAccTangX =          - gravAccNorm * nxAvg
    gravAccTangY =          - gravAccNorm * nyAvg
    gravAccTangZ = -gravAcc - gravAccNorm * nzAvg
    # adding gravity force contribution
    forceX = forceX + gravAccTangX * mass
    forceY = forceY + gravAccTangY * mass
    forceZ = forceZ + gravAccTangZ * mass

    # Calculating bottom shear and normal stress
    # bottom normal stress sigmaB
    sigmaB = - effAccNorm * rho * h
    # SamosAT friction type (bottom shear stress)
    tau = SamosATfric(cfg, uMag, sigmaB, h)
    # coulomb friction type (bottom shear stress)
    # tau = mu * sigmaB
    tau = np.where(effAccNorm > 0.0, 0, tau)

    # adding bottom shear resistance contribution
    forceBotTang = - A * tau
    # print(np.min(np.abs(forceBotTang)))
    # print(np.max(np.abs(forceBotTang)))
    forceX = forceX + forceBotTang * uxDir
    forceY = forceY + forceBotTang * uyDir
    forceZ = forceZ + forceBotTang * uzDir

    force['dM'] = dM
    force['forceX'] = forceX
    force['forceY'] = forceY
    force['forceZ'] = forceZ

    return force


def computeForceSPH(cfg, particles, force, dem):
    """ Compute SPH """

    # Load required parameters
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    Npart = particles['Npart']
    nrows = dem['header'].nrows
    ncols = dem['header'].ncols
    csz = dem['header'].cellsize

    # initialize fields for force
    forceSPHX = np.zeros(Npart)
    forceSPHY = np.zeros(Npart)
    forceSPHZ = np.zeros(Npart)

    H = np.zeros(np.shape(particles['h']))
    # loop on particles
    # TcpuSPH = 0
    # Tcpuadd = 0
    for j in range(Npart):
        mass = particles['m'][j]
        # adding lateral force (SPH component)
        # startTime = time.time()
        # gradhX, gradhY,  gradhZ, _ = calcGradHSPH(particles, j, ncols, nrows, csz)
        h, gradhX, gradhY,  gradhZ, _ = calcGradHSPHVect(particles, j, ncols, nrows, csz)
        H[j] = h / rho
        # tcpuSPH = time.time() - startTime
        # TcpuSPH = TcpuSPH + tcpuSPH
        # startTime = time.time()
        forceSPHX[j] = forceSPHX[j] - gradhX * mass * (-gravAcc) / rho
        forceSPHY[j] = forceSPHY[j] - gradhY * mass * (-gravAcc) / rho
        forceSPHZ[j] = forceSPHZ[j] - gradhZ * mass * (-gravAcc) / rho
        # tcpuadd = time.time() - startTime
        # Tcpuadd = Tcpuadd + tcpuadd

    # log.info(('cpu time SPH = %s s' % (TcpuSPH / Npart)))
    # log.info(('cpu time SPH add = %s s' % (Tcpuadd / Npart)))

    # print(np.min(np.abs(forceSPHX)))
    # print(np.max(np.abs(forceSPHX)))
    force['forceSPHX'] = forceSPHX
    force['forceSPHY'] = forceSPHY
    force['forceSPHZ'] = forceSPHZ
    particles['hSPH'] = H

    return particles, force


def updatePosition(cfg, particles, dem, force):
    dt = cfg.getfloat('dt')
    log.info('dt used now is %f' % dt)
    gravAcc = cfg.getfloat('gravAcc')
    csz = dem['header'].cellsize
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    dM = force['dM']
    forceX = force['forceX']
    forceY = force['forceY']
    forceZ = force['forceZ']
    forceSPHX = force['forceSPHX']
    forceSPHY = force['forceSPHY']
    forceSPHZ = force['forceSPHZ']
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
    uxNew = ux + (forceX + forceSPHX) * dt / mass
    uyNew = uy + (forceY + forceSPHY) * dt / mass
    uzNew = uz + (forceZ + forceSPHZ) * dt / mass
    # print anormally big values
    uNew = norm(uxNew, uyNew, uzNew)
    if np.max(uNew) > 100:
        ind = np.argmax(uNew)
        A = mass[ind] / (h[ind] * rho)
        print('particle index : ', ind)
        print((forceX / mass)[ind])
        print((forceY / mass)[ind])
        print((forceZ / mass)[ind])
        print(uNew[ind])
        print(A)
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

    #################################################################
    # this is dangerous!!!!!!!!!!!!!!
    ###############################################################
    # remove particles that are not located on the mesh any more
    particles = removeOutPart(cfg, particles, dem)
    # uN = particles['ux']*nx + particles['uy']*ny + particles['uz']*nz
    # print(norm(particles['ux'], particles['uy'], particles['uz']), uN)
    # kinEneNew = np.sum(0.5 * massNew * norm2(particles['ux'], particles['uy'], particles['uz']))
    # potEneNew = np.sum(gravAcc * massNew * particles['z'])
    # totEneNew = kinEneNew + potEneNew
    # log.info('total energy variation: %f' % ((totEneNew - totEne) / totEneNew))
    # log.info('kinetic energy variation: %f' % ((kinEneNew - kinEne) / kinEneNew))
    # if (abs(totEneNew - totEne) / totEneNew < 0.01) and (abs(kinEneNew - kinEne) / kinEneNew < 0.01):
    #     log.info('Reached stopping creteria : total energy varied fom less than 1 %')
    #     particles['stoppCriteria'] = True
    # particles['kineticEne'] = kinEneNew
    # particles['potentialEne'] = potEneNew
    return particles


def updateFields(cfg, particles, force, dem, fields):
    rho = cfg.getfloat('rho')
    depMin = cfg.getfloat('depMin')
    header = dem['header']
    Z = dem['rasterData']
    csz = dem['header'].cellsize
    A = dem['Area']
    ncols = header.ncols
    nrows = header.nrows
    Npart = particles['Npart']
    m = particles['m']
    dM = force['dM']
    x = particles['x']
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    uz = particles['uz']
    PV = fields['pv']
    PP = fields['ppr']
    PFD = fields['pfd']

    #########################################
    # Update fields using a SPH approach
    MassSPH = np.zeros((nrows, ncols))
    hSPH = np.ones((nrows, ncols))
    VXSPH = np.zeros((nrows, ncols))
    VYSPH = np.zeros((nrows, ncols))
    VZSPH = np.zeros((nrows, ncols))
    MassSPH = pointsToRasterSPH(particles, rho, Z, m, MassSPH, csz=csz)
    hSPH = pointsToRasterSPH(particles, rho, Z, m, hSPH, csz=csz)
    VXSPH = pointsToRasterSPH(particles, rho, Z, ux, VXSPH, csz=csz)
    VYSPH = pointsToRasterSPH(particles, rho, Z, uy, VYSPH, csz=csz)
    VZSPH = pointsToRasterSPH(particles, rho, Z, uz, VZSPH, csz=csz)
    VSPH = norm(VXSPH, VYSPH, VZSPH)
    FDSPH = hSPH
    # FDSPH = MassSPH / (A * rho)
    PSPH = VSPH * VSPH * rho
    # PV = np.where(VSPH > PV, VSPH, PV)
    # PP = np.where(PSPH > PP, PSPH, PP)
    # PFD = np.where(FDSPH > PFD, FDSPH, PFD)

    #########################################
    # Update fields using a bilinear interpolation
    MassBilinear = np.zeros((nrows, ncols))
    MomBilinearX = np.zeros((nrows, ncols))
    MomBilinearY = np.zeros((nrows, ncols))
    MomBilinearZ = np.zeros((nrows, ncols))
    #
    # # startTime = time.time()
    MassBilinear = geoTrans.pointsToRaster(x, y, m, MassBilinear, csz=csz, interp='bilinear')
    MomBilinearX = geoTrans.pointsToRaster(x, y, m * ux, MomBilinearX, csz=csz, interp='bilinear')
    MomBilinearY = geoTrans.pointsToRaster(x, y, m * uy, MomBilinearY, csz=csz, interp='bilinear')
    MomBilinearZ = geoTrans.pointsToRaster(x, y, m * uz, MomBilinearZ, csz=csz, interp='bilinear')

    # VXBilinear = np.where(MassBilinear > 0, MomBilinearX/MassBilinear, MomBilinearX)
    # VYBilinear = np.where(MassBilinear > 0, MomBilinearY/MassBilinear, MomBilinearY)
    # VZBilinear = np.where(MassBilinear > 0, MomBilinearZ/MassBilinear, MomBilinearZ)
    # VBilinear = norm(VXBilinear, VYBilinear, VZBilinear)
    FDBilinear = MassBilinear / (A * rho)
    # PBilinear = VBilinear * VBilinear * rho
    # PV = np.where(VBilinear > PV, VBilinear, PV)
    # PP = np.where(PBilinear > PP, PBilinear, PP)
    # PFD = np.where(FDBilinear > PFD, FDBilinear, PFD)
    # # endTime = time.time()
    # # log.info(('time = %s s' % (endTime - startTime)))

    #########################################
    # Update fields using a nearest interpolation
    MassNearest = np.zeros((nrows, ncols))
    MomNearestX = np.zeros((nrows, ncols))
    MomNearestY = np.zeros((nrows, ncols))
    MomNearestZ = np.zeros((nrows, ncols))

    # startTime = time.time()
    iC = particles['InCell'][:, 2]
    MassNearest = MassNearest.flatten()
    np.add.at(MassNearest, iC, m)
    MomNearestX = MomNearestX.flatten()
    np.add.at(MomNearestX, iC, m * ux)
    MomNearestY = MomNearestY.flatten()
    np.add.at(MomNearestY, iC, m * uy)
    MomNearestZ = MomNearestZ.flatten()
    np.add.at(MomNearestZ, iC, m * uz)
    MassNearest = np.reshape(MassNearest, (nrows, ncols))
    MomNearestX = np.reshape(MomNearestX, (nrows, ncols))
    MomNearestY = np.reshape(MomNearestY, (nrows, ncols))
    MomNearestZ = np.reshape(MomNearestZ, (nrows, ncols))
    VXNearest = np.where(MassNearest > 0, MomNearestX/MassNearest, MomNearestX)
    VYNearest = np.where(MassNearest > 0, MomNearestY/MassNearest, MomNearestY)
    VZNearest = np.where(MassNearest > 0, MomNearestZ/MassNearest, MomNearestZ)
    VNearest = norm(VXNearest, VYNearest, VZNearest)
    FDNearest = MassNearest / (A * rho)
    PNearest = VNearest * VNearest * rho
    PV = np.where(VNearest > PV, VNearest, PV)
    PP = np.where(PNearest > PP, PNearest, PP)
    PFD = np.where(FDNearest > PFD, FDNearest, PFD)

    # endTime = time.time()
    # log.info(('time = %s s' % (endTime - startTime)))

    #######################################
    print(np.sum(MassNearest))
    print(np.sum(MassBilinear))
    print(np.sum(MassSPH))

    ###################################
    fields['V'] = VNearest
    fields['P'] = PNearest
    fields['FD'] = FDNearest
    fields['pv'] = PV
    fields['ppr'] = PP
    fields['pfd'] = PFD

    hNN, _ = geoTrans.projectOnRasterVectRoot(x, y, FDNearest, csz=csz, interp='nearest')
    particles['hNearestNearest'] = hNN  # np.where(h < depMin, depMin, h)
    hNB, _ = geoTrans.projectOnRasterVectRoot(x, y, FDNearest, csz=csz, interp='bilinear')
    particles['hNearestBilinear'] = hNB  # np.where(h < depMin, depMin, h)
    hBN, _ = geoTrans.projectOnRasterVectRoot(x, y, FDBilinear, csz=csz, interp='nearest')
    particles['hBilinearNearest'] = hBN  # np.where(h2 < depMin, depMin, h2)
    hBB, _ = geoTrans.projectOnRasterVectRoot(x, y, FDBilinear, csz=csz, interp='bilinear')
    particles['hBilinearBilinear'] = hBB  # np.where(h2 < depMin, depMin, h2)

    # choose the interpolation method
    particles['h'] = hNN

    # remove particles that have a too small height
    # particles = removeSmallPart(hmin, particles, dem)

    return particles, fields


def plotPosition(particles, dem, data, Cmap, unit, fig, ax, plotPart=False):
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
    try:
        # Get the images on an axis
        cb = ax.images[-1].colorbar
        if cb:
            cb.remove()
    except IndexError:
        pass

    ax.clear()
    ax.set_title('Particles on dem at t=%.2f s' % particles['t'])
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        Cmap, 0.0, np.nanmax(data), continuous=contCmap)
    cmap.set_under(color='w')
    ref0, im = NonUnifIm(ax, xx, yy, data, 'x [m]', 'y [m]',
                         extent=[x.min(), x.max(), y.min(), y.max()],
                         cmap=cmap, norm=norm)
    if plotPart:
        ax.plot(x, y, 'ok', linestyle='None', markersize=1)
    Cp1 = ax.contour(X, Y, Z, levels=10, colors='k')
    addColorBar(im, ax, ticks, unit)
    plt.pause(0.1)
    # plt.close(fig)
    # ax.set_ylim([510, 530])
    # ax.set_xlim([260, 300])
    return fig, ax


def removeOutPart(cfg, particles, dem):
    header = dem['header']
    nrows = header.nrows
    ncols = header.ncols
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize

    x = particles['x']
    y = particles['y']

    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    Lx = (x - xllc) / csz
    Ly = (y - yllc) / csz
    mask = np.ones(len(x), dtype=bool)
    indOut = np.where(Lx <= -0.5)
    mask[indOut] = False
    indOut = np.where(Ly <= -0.5)
    mask[indOut] = False
    indOut = np.where(Lx >= ncols-0.5)
    mask[indOut] = False
    indOut = np.where(Ly >= nrows-0.5)
    mask[indOut] = False

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they exited the domain' % (nRemove))

    return particles


def removeSmallPart(hmin, particles, dem):

    h = particles['h']

    indOut = np.where(h < hmin)
    mask = np.ones(len(h), dtype=bool)
    mask[indOut] = False

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they were too thin' % (nRemove))
        particles = getNeighboursVect(particles, dem)

    return particles


def removePart(particles, mask, nRemove):

    particles['Npart'] = particles['Npart'] - nRemove
    particles['x'] = particles['x'][mask]
    particles['y'] = particles['y'][mask]
    particles['z'] = particles['z'][mask]
    particles['s'] = particles['s'][mask]
    particles['ux'] = particles['ux'][mask]
    particles['uy'] = particles['uy'][mask]
    particles['uz'] = particles['uz'][mask]
    particles['m'] = particles['m'][mask]
    particles['h'] = particles['h'][mask]
    particles['InCell'] = particles['InCell'][mask, :]
    particles['partInCell'] = particles['partInCell'][mask]

    return particles


def getNeighbours(particles, dem):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    Npart = particles['Npart']
    x = particles['x']
    y = particles['y']
    check = np.zeros((nrows, ncols))
    indPartInCell = np.zeros(ncols*nrows + 1).astype(int)
    partInCell = np.zeros(Npart).astype(int)
    # Count number of particles in each cell
    indx = ((x + csz/2) / csz).astype(int)
    indy = ((y + csz/2) / csz).astype(int)
    check[indy, indx] = check[indy, indx] + 1
    ic = indx + ncols * indy
    ##################################
    # start ########################
    # TODO: test speed between add.at and bincount
    # indPartInCell = np.bincount(ic + 1, minlength=len(indPartInCell))
    # or
    np.add.at(indPartInCell, ic + 1, 1)
    # end ##############################
    ##################################
    indPartInCell = np.cumsum(indPartInCell)

    # make the list of which particles are in which cell
    InCell = np.empty((0, 3), int).astype(int)
    indPartInCell2 = copy.deepcopy(indPartInCell)
    for j in range(Npart):
        indx = int((x[j] + csz/2) / csz)
        indy = int((y[j] + csz/2) / csz)
        ic = indx + ncols * indy
        partInCell[int(indPartInCell2[ic])] = j
        indPartInCell2[ic] = indPartInCell2[ic] + 1
        InCell = np.append(InCell, np.tile(np.array([indx, indy, ic]), (1, 1)), axis=0)

    particles['InCell'] = InCell
    particles['indPartInCell'] = indPartInCell
    particles['partInCell'] = partInCell

    return particles


def getNeighboursVect(particles, dem):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    Npart = particles['Npart']
    x = particles['x']
    y = particles['y']
    check = np.zeros((nrows, ncols))
    indPartInCell = np.zeros(ncols*nrows + 1).astype(int)
    partInCell = np.zeros(Npart).astype(int)
    # Count number of particles in each cell
    indx = ((x + csz/2) / csz).astype(int)
    indy = ((y + csz/2) / csz).astype(int)
    check[indy, indx] = check[indy, indx] + 1
    ic = indx + ncols * indy
    ##################################
    # start ########################
    # TODO: test speed between add.at and bincount
    # indPartInCell = np.bincount(ic + 1, minlength=len(indPartInCell))
    # or
    np.add.at(indPartInCell, ic + 1, 1)
    # end ##############################
    ##################################
    indPartInCell = np.cumsum(indPartInCell)

    # no loop
    partInCell = np.argsort(ic, kind='mergesort')
    InCell = np.vstack((indx, indy))
    InCell = np.vstack((InCell, ic))
    InCell = InCell.T

    particles['InCell'] = InCell
    particles['indPartInCell'] = indPartInCell
    particles['partInCell'] = partInCell

    return particles


def calcGradHSPH(particles, j, ncols, nrows, csz):
    # SPH kernel
    # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
    rKernel = csz
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))
    dfacKernel = -3.0 * facKernel

    indx, indy, _ = particles['InCell'][j]
    indPartInCell = particles['indPartInCell']
    partInCell = particles['partInCell']
    x = particles['x'][j]
    y = particles['y'][j]
    z = particles['z'][j]
    # With loop
    gradhX = 0
    gradhY = 0
    gradhZ = 0
    # startTime = time.time()
    L = np.empty((0), dtype=int)
    # check if we are on the bottom ot top row!!!
    lInd = -1
    rInd = 2
    if indy == 0:
        lInd = 0
    if indy == nrows - 1:
        rInd = 1
    for n in range(lInd, rInd):
        ic = (indx - 1) + ncols * (indy + n)
        # make sure not to take particles from the other edge
        iPstart = indPartInCell[max(ic, ncols * (indy + n))]
        iPend = indPartInCell[min(ic+3, ncols * (indy + n + 1))]
        # loop on all particles in neighbour boxes
        for p in range(iPstart, iPend):
            # index of particle in neighbour box
            l = int(partInCell[p])
            if j != l:
                L = np.append(L, l)
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
    # print((time.time() - startTime)/1)

    return gradhX, gradhY,  gradhZ, L


def calcGradHSPHVect(particles, j, ncols, nrows, csz):
    # SPH kernel
    # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
    rKernel = csz
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
    h = 0
    # no loop option

    # startTime = time.time()
    # find all the neighbour boxes
    # index of the left box
    # check if we are on the bottom ot top row!!!
    lInd = -1
    rInd = 2
    if indy == 0:
        lInd = 0
    if indy == nrows - 1:
        rInd = 1
    ic = (indx - 1) + ncols * (indy + np.arange(lInd, rInd, 1))
    # go throught all index from left middle and right box and get the
    # particles in those boxes
    # make sure not to take particles from the other edge
    iPstart = indPartInCell[np.maximum(ic, ncols * (indy + np.arange(lInd, rInd, 1)))]
    iPend = indPartInCell[np.minimum(ic + 3, ncols * (indy + np.arange(lInd, rInd, 1) + 1))]

    # gather all the particles
    index = np.concatenate([np.arange(x, y) for x, y in zip(iPstart, iPend)])

    # compute SPH gradient
    L = partInCell[index]
    # make sure to remove the j particle
    indJ = np.where(L == j)
    L = np.delete(L, indJ)
    dx = particles['x'][L] - x
    dy = particles['y'][L] - y
    dz = 0 * (particles['z'][L] - z)
    r = norm(dx, dy, dz)
    # impose a minimum distance between particles
    r = np.where(r < 0.001 * rKernel, 0.001 * rKernel, r)
    hr = rKernel - r
    w = facKernel * hr * hr * hr
    dwdr = dfacKernel * hr * hr
    massl = particles['m'][L]
    # GX = massl * dwdr * dx / r
    # GX = np.where(r < rKernel, GX, 0)
    # GY = massl * dwdr * dy / r
    # GY = np.where(r < rKernel, GY, 0)
    # GZ = massl * dwdr * dz / r
    # GZ = np.where(r < rKernel, GZ, 0)
    # ------------------------------
    indOut = np.where(r >= rKernel)
    massl[indOut] = 0
    mdwdrr = massl * dwdr / r
    H = massl * w
    GX = mdwdrr * dx
    GY = mdwdrr * dy
    GZ = mdwdrr * dz
    # -----------------------------
    gradhX = gradhX + np.sum(GX)
    gradhY = gradhY + np.sum(GY)
    gradhZ = gradhZ + np.sum(GZ)
    h = h + np.sum(H)
    # leInd = len(index)
    # print((time.time() - startTime)/leInd)

    return h, gradhX, gradhY,  gradhZ, L


def pointsToRasterSPH(particles, rho, Z, f, F, csz=1, xllc=0, yllc=0):
    """
    Vectorized version of projectOnRaster
    Projects the points Points on Raster using a bilinear or nearest
    interpolation and returns the z coord
    Input :
    Points: (x, y) coord of the pointsi
    Output:
    PointsZ: z coord of the points
             ioob number of out of bounds indexes
    """
    # SPH kernel
    # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
    rKernel = csz
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))

    x = particles['x']
    y = particles['y']
    z = particles['z']
    m = particles['m']

    nrow, ncol = np.shape(F)

    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    Lx = (x - xllc) / csz
    Ly = (y - yllc) / csz

    # find coordinates of the 4 nearest cornes on the raster
    Lx0 = np.floor(Lx).astype('int')
    Ly0 = np.floor(Ly).astype('int')
    Lx1 = Lx0 + 1
    Ly1 = Ly0 + 1

    dx = (Lx - Lx0) * csz
    dy = (Ly - Ly0) * csz

    dz = z - Z[Lx1, Ly0]
    R10 = norm(csz-dx, dy, dz)
    dz = z - Z[Lx0, Ly1]
    R01 = norm(dx, csz-dy, dz)
    dz = z - Z[Lx1, Ly1]
    R11 = norm(csz-dx, csz-dy, dz)

    F = F.flatten()
    Z = Z.flatten()

    # add sph contribution to bottom left grid cell
    ic00 = Lx0 + ncol * Ly0
    dz = z - Z[ic00]
    R00 = norm(dx, dy, dz)
    indOut = np.where(R00 >= rKernel)
    R00[indOut] = 0
    hr00 = rKernel - R00
    W00 = facKernel * hr00 * hr00 * hr00
    f00 = f * m * W00 / rho
    np.add.at(F, ic00, f00)

    # add sph contribution to bottom right grid cell
    ic10 = Lx1 + ncol * Ly0
    dz = z - Z[ic10]
    R10 = norm(dx, dy, dz)
    indOut = np.where(R10 >= rKernel)
    R10[indOut] = 0
    hr10 = rKernel - R10
    W10 = facKernel * hr10 * hr10 * hr10
    f10 = f * m * W10 / rho
    np.add.at(F, ic10, f10)

    # add sph contribution to top left grid cell
    ic01 = Lx0 + ncol * Ly1
    dz = z - Z[ic01]
    R01 = norm(dx, dy, dz)
    indOut = np.where(R01 >= rKernel)
    R01[indOut] = 0
    hr01 = rKernel - R01
    W01 = facKernel * hr01 * hr01 * hr01
    f01 = f * m * W01 / rho
    np.add.at(F, ic01, f01)

    # add sph contribution to top right grid cell
    ic11 = Lx1 + ncol * Ly1
    dz = z - Z[ic11]
    R11 = norm(dx, dy, dz)
    indOut = np.where(R11 >= rKernel)
    R11[indOut] = 0
    hr11 = rKernel - R10
    W11 = facKernel * hr11 * hr11 * hr11
    f11 = f * m * W11 / rho
    np.add.at(F, ic11, f11)

    F = np.reshape(F, (nrow, ncol))

    return F


def getNormal(x, y, Nx, Ny, Nz, csz):
    # by default bilinear interpolation of the Nx, Ny, Nz of the grid
    nx = geoTrans.projectOnRasterRoot(x, y, Nx, csz=csz)
    ny = geoTrans.projectOnRasterRoot(x, y, Ny, csz=csz)
    nz = geoTrans.projectOnRasterRoot(x, y, Nz, csz=csz)
    nx, ny, nz = normalize(nx, ny, nz)
    return nx, ny, nz


def getNormalArray(x, y, Nx, Ny, Nz, csz):
    nrow, ncol = np.shape(Nx)
    # by default bilinear interpolation of the Nx, Ny, Nz of the grid
    nx, _ = geoTrans.projectOnRasterVectRoot(x, y, Nx, csz=csz)
    ny, _ = geoTrans.projectOnRasterVectRoot(x, y, Ny, csz=csz)
    nz, _ = geoTrans.projectOnRasterVectRoot(x, y, Nz, csz=csz)
    nx, ny, nz = normalize(nx, ny, nz)
    return nx, ny, nz


def getNormalMesh(z, csz, num=4):
    n, m = np.shape(z)
    Nx = np.ones((n, m))
    Ny = np.ones((n, m))
    Nz = np.ones((n, m))
    # first and last row, first and last column are inacurate
    if num == 4:
        # normal calculation with 4 triangles
        # (Zl - Zr) * csz
        Nx[1:n-1, 1:m-1] = (z[1:n-1, 0:m-2] - z[1:n-1, 2:m]) / csz
        # (Zd - Zu) * csz
        Ny[1:n-1, 1:m-1] = (z[0:n-2, 1:m-1] - z[2:n, 1:m-1]) / csz
        Nz = 2 * Nz

    if num == 6:
        # normal calculation with 6 triangles
        # (2*(Zl - Zr) - Zur + Zdl + Zu - Zd) * csz
        Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m])
                            - z[2:n, 2:m] + z[0:n-2, 0:m-2]
                            + z[2:n, 1:m-1] - z[0:n-2, 1:m-1]) / csz
        # (2*(Zd - Zu) - Zur + Zdl - Zl + Zr) * csz
        Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1])
                            - z[2:n, 2:m] + z[0:n-2, 0:m-2]
                            - z[1:n-1, 0:m-2] + z[1:n-1, 2:m]) / csz
        Nz = 6 * Nz
    if num == 8:
        # normal calculation with 8 triangles
        # (2*(Zl - Zr) + Zul - Zur + Zdl - Zdr) * csz
        Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m])
                            + z[2:n, 0:m-2] - z[2:n, 2:m]
                            + z[0:n-2, 0:m-2] - z[0:n-2, 2:m]) / csz
        # (2*(Zd - Zu) - Zul - Zur + Zdl + Zdr) * csz
        Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1])
                            - z[2:n, 0:m-2] - z[2:n, 2:m]
                            + z[0:n-2, 0:m-2] + z[0:n-2, 2:m]) / csz
        Nz = 8 * Nz

    Nx, Ny, Nz = normalize(Nx, Ny, Nz)

    return Nx, Ny, Nz


def getAreaMesh(Nx, Ny, Nz, csz):
    A = 1/(Nz*Nz) * np.sqrt((Nz*Nz + Nx*Nx) * (Nz*Nz + Ny*Ny)) * csz*csz
    return A


def norm(x, y, z):
    norme = np.sqrt(x*x + y*y + z*z)
    return norme


def norm2(x, y, z):
    norme2 = (x*x + y*y + z*z)
    return norme2


def normalize(x, y, z):
    # TODO : avoid error message when input vector is zero and make sure
    # to return zero
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
    div = np.where(div < 1.0, 1, div)
    # if(div < 1.0):
    #     div = 1.0
    div = np.log(div) / kappa + B
    tau = p * mu * (1.0 + Rs0 / (Rs0 + Rs)) + rho * v * v / (div * div)
    return tau
