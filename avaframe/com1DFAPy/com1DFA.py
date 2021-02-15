"""
    Main functions for python DFA kernel
"""

import logging
import time
import os
import numpy as np
import glob
import math
import copy
import pickle
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.makePalette as makePalette
import avaframe.com1DFAPy.timeDiscretizations as tD
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.DFAfunctionsCython as DFAfunC
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils

#######################################
# Set flags here
#######################################
# create local logger
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')
# set feature flag for initial particle distribution
# particles are homegeneosly distributed with a little random variation
flagSemiRand = False
# particles are randomly distributed
flagRand = True
# set feature flag for flow deth calculation
# use SPH to get the particles flow depth
flagFDSPH = False
# set feature leapfrog time stepping
featLF = False
featCFL = False
featCFLConstrain = False


def initializeMesh(dem, num=4):
    """ Create rectangular mesh

    Reads the DEM information, computes the normal vector field and
    boundries to the DEM

    Parameters
    ----------
    dem : dict
        dictionary with dem information
    num : int
        chose between 4, 6 or 8 (using then 4, 6 or 8 triangles) or
        1 to use the simple cross product method

    Returns
    -------
    dem : dict
        dictionary completed with normal field and boundaries
    """
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    # get normal vector of the grid mesh
    Nx, Ny, Nz = DFAtls.getNormalMesh(dem['rasterData'], csz, num=num)
    dem['Nx'] = np.where(np.isnan(Nx), 0, Nx)
    dem['Ny'] = np.where(np.isnan(Ny), 0, Ny)
    # build no data mask (used to find out of dem particles)
    bad = np.where(Nz > 1, True, False)
    dem['Nz'] = np.where(Nz > 1, 0, Nz)
    dem['Bad'] = bad
    if debugPlot:
        x = np.arange(ncols) * csz
        y = np.arange(nrows) * csz
        fig = plt.figure(figsize=(pU.figW, pU.figH))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.plot(x, dem['Nx'][int(nrows/2+1), :], '--k', label='Nx(x)')
        ax1.plot(x, dem['Ny'][int(nrows/2+1), :], '--b', label='Ny(x)')
        ax1.plot(x, dem['Nz'][int(nrows/2+1), :], '--r', label='Nz(x)')
        ax2.plot(y, dem['Nx'][:, int(ncols/2+1)], ':k', label='Nx(y)')
        ax2.plot(y, dem['Ny'][:, int(ncols/2+1)], ':b', label='Ny(y)')
        ax2.plot(y, dem['Nz'][:, int(ncols/2+1)], ':r', label='Nz(y)')
        ax1.legend()
        ax2.legend()
        plt.show()
        IOf.writeResultToAsc(dem['header'], dem['Nx'], 'Nx.asc')
        IOf.writeResultToAsc(dem['header'], dem['Ny'], 'Ny.asc')
        IOf.writeResultToAsc(dem['header'], dem['Nz'], 'Nz.asc')

    # get real Area
    Area = DFAtls.getAreaMesh(Nx, Ny, Nz, csz)
    dem['Area'] = Area
    log.info('Largest cell area: %f m²' % (np.amax(Area)))
    log.debug('Projected Area :', ncols * nrows * csz * csz)
    log.debug('Total Area :', np.sum(Area))

    return dem


def initializeSimulation(cfg, relRaster, dem):
    """ Initialize DFA simulation

    Create particles and fields dictionary according to config parameters
    release raster and dem

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    relRaster: 2D numpy array
        release depth raster
    dem : dict
        dictionary with dem information

    Returns
    -------
    particles : dict
        particles dictionary at initial time step
    fields : dict
        fields dictionary at initial time step
    Cres : 2D numpy array
        resistance raster
    Ment : 2D numpy array
        entrained mass raster
    """
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
    Mraster = np.sum(A*relRaster*rho)
    # initialize arrays
    partPerCell = np.zeros(np.shape(relRaster), dtype=np.int64)
    FD = np.zeros((nrows, ncols))
    Npart = 0
    NPPC = np.empty(0)
    Apart = np.empty(0)
    Xpart = np.empty(0)
    Ypart = np.empty(0)
    Mpart = np.empty(0)
    Hpart = np.empty(0)
    InCell = np.empty((0), int)
    IndX = np.empty((0), int)
    IndY = np.empty((0), int)
    # find all non empty cells (meaning release area)
    indRelY, indRelX = np.nonzero(relRaster)
    # loop on non empty cells
    for indRelx, indRely in zip(indRelX, indRelY):
        # compute number of particles for this cell
        h = relRaster[indRely, indRelx]
        Vol = A[indRely, indRelx] * h
        mass = Vol * rho
        xpart, ypart, mPart, nPart = placeParticles(mass, indRelx, indRely, csz, massPerPart)
        Npart = Npart + nPart
        partPerCell[indRely, indRelx] = nPart
        # initialize field Flow depth
        FD[indRely, indRelx] = h
        # initialize particles position, mass, height...
        NPPC = np.append(NPPC, nPart*np.ones(nPart))
        Apart = np.append(Apart, A[indRely, indRelx]*np.ones(nPart)/nPart)
        Xpart = np.append(Xpart, xpart)
        Ypart = np.append(Ypart, ypart)
        Mpart = np.append(Mpart, mPart * np.ones(nPart))
        Hpart = np.append(Hpart, h * np.ones(nPart))
        ic = indRelx + ncols * indRely
        IndX = np.append(IndX, np.ones(nPart)*indRelx)
        IndY = np.append(IndY, np.ones(nPart)*indRely)
        InCell = np.append(InCell, np.ones(nPart)*ic)

    Hpart, _ = geoTrans.projectOnRasterVectRoot(Xpart, Ypart, relRaster, csz=csz, interp='bilinear')
    Mpart = rho * Hpart * Apart
    # create dictionnary to store particles properties
    particles = {}
    particles['Npart'] = Npart
    particles['NPPC'] = NPPC
    particles['x'] = Xpart
    particles['y'] = Ypart
    particles['s'] = np.zeros(np.shape(Xpart))
    # adding z component
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')
    # readjust mass
    mTot = np.sum(Mpart)
    particles['m'] = Mpart*Mraster/mTot
    particles['mTot'] = np.sum(particles['m'])
    particles['h'] = Hpart
    particles['hNearestNearest'] = Hpart
    particles['hNearestBilinear'] = Hpart
    particles['hBilinearNearest'] = Hpart
    particles['hBilinearBilinear'] = Hpart
    particles['hSPH'] = Hpart
    particles['GHX'] = np.zeros(np.shape(Xpart))
    particles['GHY'] = np.zeros(np.shape(Xpart))
    particles['GHZ'] = np.zeros(np.shape(Xpart))
    particles['InCell'] = InCell
    particles['indX'] = IndX
    particles['indY'] = IndY
    particles['ux'] = np.zeros(np.shape(Xpart))
    particles['uy'] = np.zeros(np.shape(Xpart))
    particles['uz'] = np.zeros(np.shape(Xpart))
    particles['stoppCriteria'] = False
    kineticEne = np.sum(0.5 * Mpart * DFAtls.norm2(particles['ux'], particles['uy'], particles['uz']))
    particles['kineticEne'] = kineticEne
    particles['potentialEne'] = np.sum(gravAcc * Mpart * particles['z'])
    particles['peakKinEne'] = kineticEne

    # initialize entrainment and resistance
    Ment = intializeMassEnt(dem)
    Cres = intializeResistance(dem)
    PFV = np.zeros((nrows, ncols))
    PP = np.zeros((nrows, ncols))
    fields = {}
    fields['pfv'] = PFV
    fields['ppr'] = PP
    fields['pfd'] = FD
    fields['FV'] = PFV
    fields['P'] = PP
    fields['FD'] = FD
    fields['Vx'] = PFV
    fields['Vy'] = PFV
    fields['Vz'] = PFV

    # get particles location (neighbours for sph)
    # particles = getNeighbours(particles, dem)

    particles = DFAfunC.getNeighboursC(particles, dem)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)

    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    indX = (particles['indX']).astype('int')
    indY = (particles['indY']).astype('int')
    H, C, W = DFAfunC.computeFDC(cfg, particles, header, Nx, Ny, Nz, indX, indY)
    H = np.asarray(H)
    # particles['h'] = H
    # H, W = SPHC.computeFDC(cfg, particles, header, Nx, Ny, Nz, indX, indY)
    # particles['h'] = hh
    # H = np.asarray(H)
    W = np.asarray(W)
    particles['hSPH'] = H/W
    # initialize time
    t = 0
    particles['t'] = t

    log.info('Expeced mass. Mexpected = %f kg.' % (Mraster))
    log.info('Initializted simulation. MTot = %f kg, %s particles in %s cells. %f kg per particle' %
             (particles['mTot'], particles['Npart'], np.size(indRelY), particles['mTot']/particles['Npart']))

    if debugPlot:
        x = np.arange(ncols) * csz
        y = np.arange(nrows) * csz
        # X, Y = np.meshgrid(x, y)
        # fig = plt.figure()
        # ax = fig.gca(projection='3d')
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        ref0, im = pU.NonUnifIm(ax, x, y, A, 'x [m]', 'y [m]',
                                extent=[x.min(), x.max(), y.min(), y.max()],
                                cmap=cmap, norm=None)

        ax.plot(Xpart, Ypart, 'or', linestyle='None')
        pU.addColorBar(im, ax, None, 'm²')

        # ax.quiver(x, 2500*np.ones(nrows), dem['rasterData'][500, :], Nx[500, :], Ny[500, :], Nz[500, :], length=50, colors='b')
        # ax.quiver(2500*np.ones(ncols), y, dem['rasterData'][:, 500], Nx[:, 500], Ny[:, 500], Nz[:, 500], length=50, colors='r')
        # ax.set_xlim([2000, 3000])
        # ax.set_ylim([2000, 3000])
        # ax.set_zlim([0, 500])
        # # Label each axis
        # ax.set_xlabel('x')
        # ax.set_ylabel('y')
        # ax.set_zlabel('z')

        plt.show()

    return particles, fields, Cres, Ment


def placeParticles(mass, indx, indy, csz, massPerPart):
    """ Create particles in given cell

    Compute number of particles to create in a given cell.
    Place particles in cell according to the chosen pattern (random semirandom
    or ordered)

    Parameters
    ----------
    mass: float
        mass of snow in cell
    indx: int
        column index of the cell
    indy: int
        row index of the cell
    csz : float
        cellsize
    massPerPart : float
        maximum mass per particle

    Returns
    -------
    xpart : 1D numpy array
        x position of particles
    ypart : 1D numpy array
        y position of particles
    mPart : 1D numpy array
        mass of particles
    nPart : int
        number of particles created
    """
    if flagRand:
        nPart = (np.round(mass / massPerPart) + round(np.random.rand())).astype('int')
        # enshure that there is at last one particle
        nPart = np.maximum(nPart, 1)
    else:
        n = (np.floor(np.sqrt(mass / massPerPart)) + 1).astype('int')
        nPart = n*n
        d = csz/n
        pos = np.linspace(0, csz-d, n) + d/2
        x, y = np.meshgrid(pos, pos)
        x = x.flatten()
        y = y.flatten()

    mPart = mass / nPart
    # TODO make this an independent function
    #######################
    # start ###############
    if flagSemiRand:
        # place particles equaly distributed with a small variation
        xpart = csz * (- 0.5 + indx) + x + (np.random.rand(nPart) - 0.5) * d
        ypart = csz * (- 0.5 + indy) + y + (np.random.rand(nPart) - 0.5) * d
    elif flagRand:
        # place particles randomly in the cell
        xpart = csz * (np.random.rand(nPart) - 0.5 + indx)
        ypart = csz * (np.random.rand(nPart) - 0.5 + indy)
    else:
        # place particles equaly distributed
        xpart = csz * (- 0.5 + indx) + x
        ypart = csz * (- 0.5 + indy) + y
    return xpart, ypart, mPart, nPart


def intializeMassEnt(dem):
    """ Intialize mass for entrainment

    Parameters
    ----------
    dem: dict
        dem dictionary

    Returns
    -------
    Ment : 2D numpy array
        raster of available mass for entrainment
    """
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    Ment = np.zeros((nrows, ncols))
    return Ment


def intializeResistance(dem):
    """ Intialize resistance matrix

    Parameters
    ----------
    dem: dict
        dem dictionary

    Returns
    -------
    Cres : 2D numpy array
        raster of resistance coefficients
    """
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    Cres = np.zeros((nrows, ncols))
    return Cres


def DFAIterate(cfg, particles, fields, dem, Ment, Cres):
    """ Perform time loop for DFA simulation

     Save results at desired intervals

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at initial time step
    fields : dict
        fields dictionary at initial time step
    dem : dict
        dictionary with dem information
    Ment : 2D numpy array
        entrained mass raster
    Cres : 2D numpy array
        resistance raster
    Tcpu : dict
        computation time dictionary

    Returns
    -------
    Particles : list
        list of particles dictionary
    Fields : list
        list of fields dictionary (for each time step saved)
    Tcpu : dict
        computation time dictionary
    """

    # Initialise cpu timing
    Tcpu = {}
    Tcpu['Force'] = 0.
    Tcpu['ForceVect'] = 0.
    Tcpu['ForceSPH'] = 0.
    Tcpu['Pos'] = 0.
    Tcpu['Neigh'] = 0.
    Tcpu['Field'] = 0.

    # Load configuration settings
    Tend = cfg.getfloat('Tend')
    dtSave = cfg.getfloat('dtSave')
    sphOption = cfg.getint('sphOption')
    log.info('using sphOption %s:' % sphOption)

    # Initialise Lists to save fields
    Particles = [copy.deepcopy(particles)]
    Fields = [copy.deepcopy(fields)]
    # save Z, S, U, T at each time step for developping purpouses
    Z = np.empty((0, 0))
    S = np.empty((0, 0))
    U = np.empty((0, 0))
    T = np.empty((0, 0))
    Tsave = [0]

    if featLF:
        log.info('Use LeapFrog time stepping')
    else:
        log.info('Use standard time stepping')
    # Initialize time and counters
    nSave = 1
    Tcpu['nSave'] = nSave
    nIter = 0
    nIter0 = 0
    iterate = True
    particles['iterate'] = iterate
    t = particles['t']
    # Start time step computation
    while t < Tend and iterate:
        # ++++++++++++++++if you want to use cfl time step+++++++++++++++++++
        # CALL TIME STEP:
        # to play around with the courant number
        if featCFL:
            dtStable = tD.getcflTimeStep(particles, dem, cfg)
        elif featCFLConstrain:
            dtStable = tD.getcfldTwithConstraints(particles, dem, cfg)

        # dt overwrites dt in .ini file, so comment this block if you dont want to use cfl
        # ++++++++++++++++++++++++++++++++++++++++++++++
        # get time step
        dt = cfg.getfloat('dt')
        t = t + dt
        nIter = nIter + 1
        nIter0 = nIter0 + 1
        log.debug('Computing time step t = %f s', t)

        # Perform computations
        if featLF:
            particles, fields, Tcpu, dt = computeLeapFrogTimeStep(
                cfg, particles, fields, dt, dem, Ment, Cres, Tcpu)
        else:
            particles, fields, Tcpu = computeEulerTimeStep(
                cfg, particles, fields, dt, dem, Ment, Cres, Tcpu)


        T = np.append(T, t)
        Tcpu['nSave'] = nSave
        particles['t'] = t
        # Save desired parameters and export to Lists for saving interval
        U = np.append(U, DFAtls.norm(particles['ux'][0], particles['uy'][0], particles['uz'][0]))
        Z = np.append(Z, particles['z'][0])
        S = np.append(S, particles['s'][0])
        iterate = particles['iterate']
        if t >= nSave * dtSave:
            Tsave.append(t)
            log.info('Saving results for time step t = %f s', t)
            log.info('MTot = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
            log.info(('cpu time Force = %s s' % (Tcpu['Force'] / nIter)))
            log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / nIter)))
            log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / nIter)))
            log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / nIter)))
            log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / nIter)))
            log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / nIter)))
            Particles.append(copy.deepcopy(particles))
            Fields.append(copy.deepcopy(fields))
            nSave = nSave + 1

    Tcpu['nIter'] = nIter
    log.info('Ending computation at time t = %f s', t)
    log.info('Saving results for time step t = %f s', t)
    log.info('MTot = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
    log.info(('cpu time Force = %s s' % (Tcpu['Force'] / nIter)))
    log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / nIter)))
    log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / nIter)))
    log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / nIter)))
    log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / nIter)))
    log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / nIter)))
    Tsave.append(t)
    Particles.append(copy.deepcopy(particles))
    Fields.append(copy.deepcopy(fields))

    return Tsave, T, U, Z, S, Particles, Fields, Tcpu


def computeEulerTimeStep(cfg, particles, fields, dt, dem, Ment, Cres, Tcpu):
    """ compute next time step using an euler forward scheme


    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at t
    fields : dict
        fields dictionary at t
    dt : float
        time step
    dem : dict
        dictionary with dem information
    Ment : 2D numpy array
        entrained mass raster
    Cres : 2D numpy array
        resistance raster
    Tcpu : dict
        computation time dictionary

    Returns
    -------
    particles : dict
        particles dictionary at t + dt
    fields : dict
        fields dictionary at t + dt
    Tcpu : dict
        computation time dictionary
    """
    sphOption = cfg.getint('sphOption')
    # get forces
    startTime = time.time()
    # loop version of the compute force
    particles, force = DFAfunC.computeForceC(cfg, particles, fields, dem, Ment, Cres, dt)
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce


    # compute lateral force (SPH component of the calculation)
    startTime = time.time()
    particles, force = DFAfunC.computeForceSPHC(cfg, particles, force, dem, SPHOption=sphOption, gradient=0)
    tcpuForceSPH = time.time() - startTime
    Tcpu['ForceSPH'] = Tcpu['ForceSPH'] + tcpuForceSPH

    nSave = Tcpu['nSave']
    dtSave = cfg.getfloat('dtSave')
    hmin = cfg.getfloat('hmin')

    # update velocity and particle position
    startTime = time.time()
    # particles = updatePosition(cfg, particles, dem, force)
    particles = DFAfunC.updatePositionC(cfg, particles, dem, force)
    tcpuPos = time.time() - startTime
    Tcpu['Pos'] = Tcpu['Pos'] + tcpuPos

    # get particles location (neighbours for sph)
    startTime = time.time()
    particles = DFAfunC.getNeighboursC(particles, dem)

    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

    # update fields (compute grid values)
    startTime = time.time()
    # particles, fields = updateFields(cfg, particles, force, dem, fields)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)
    tcpuField = time.time() - startTime
    Tcpu['Field'] = Tcpu['Field'] + tcpuField

    # get SPH flow depth

    if flagFDSPH:
        # particles = SPH.computeFlowDepth(cfg, particles, dem)
        header = dem['header']
        Nx = dem['Nx']
        Ny = dem['Ny']
        Nz = dem['Nz']
        indX = (particles['indX']).astype('int')
        indY = (particles['indY']).astype('int')
        H, C, W = DFAfunC.computeFDC(cfg, particles, header, Nx, Ny, Nz, indX, indY)
        H = np.asarray(H)
        # particles['h'] = H
        # H, W = SPHC.computeFDC(cfg, particles, header, Nx, Ny, Nz, indX, indY)
        # particles['h'] = hh
        # H = np.asarray(H)
        W = np.asarray(W)
        H = np.where(W > 0, H/W, H)
        particles['hSPH'] = np.where(H <= hmin, hmin, H)
        particles['h'] = particles['hSPH']

        # print(np.min(particles['h']))

    return particles, fields, Tcpu


def computeLeapFrogTimeStep(cfg, particles, fields, dt, dem, Ment, Cres, Tcpu):
    """ compute next time step using a Leap Frog scheme


    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at t
    fields : dict
        fields dictionary at t
    dt : float
        time step
    dem : dict
        dictionary with dem information
    Ment : 2D numpy array
        entrained mass raster
    Cres : 2D numpy array
        resistance raster
    Tcpu : dict
        computation time dictionary

    Returns
    -------
    particles : dict
        particles dictionary at t + dt
    fields : dict
        fields dictionary at t + dt
    Tcpu : dict
        computation time dictionary
    dt : float
        time step
    """

    # start timing
    startTime = time.time()
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce

    # dtK5 is half time step
    dtK5 = 0.5 * dt
    # cfg['dt'] = str(dtK5)
    log.debug('dt used now is %f' % dt)

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

    # +++++++++++++Time integration using leapfrog 'Drift-Kick-Drif' scheme+++++
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
    # particles['h'] = ?

    # 'KICK'
    # compute velocity at t_(k+0.5)
    # first compute force at t_(k+0.5)
    startTime = time.time()
    # TODO check  effect of artificial viscosity - update of velocity works here too
    particles, force = DFAfunC.computeForceC(cfg, particles, fields, dem, Ment, Cres, dtK5)
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce
    # force = computeForceVect(cfg, particles, dem, Ment, Cres, dtK5)
    startTime = time.time()
    particles, force = DFAfunC.computeForceSPHC(cfg, particles, force, dem)
    tcpuForceSPH = time.time() - startTime
    Tcpu['ForceSPH'] = Tcpu['ForceSPH'] + tcpuForceSPH
    # particles, force = computeForceSPH(cfg, particles, force, dem)
    mass = particles['m']
    uxNew = uxK + (force['forceX'] + force['forceSPHX']) * dt / mass
    uyNew = uyK + (force['forceY'] + force['forceSPHY']) * dt / mass
    uzNew = uzK + (force['forceZ'] + force['forceSPHZ']) * dt / mass

    # 'DRIF'
    # now update position at t_(k+ 1)
    xNew = xK5 + dtK5 * uxNew
    yNew = yK5 + dtK5 * uyNew
    zNew = zK5 + dtK5 * uzNew

    # ++++++++++++++UPDATE Particle Properties
    # update mass required if entrainment
    massNew = mass + force['dM']
    particles['mTot'] = np.sum(massNew)
    particles['x'] = xNew
    particles['y'] = yNew
    particles['s'] = particles['s'] + np.sqrt((xNew-xK)*(xNew-xK) + (yNew-yK)*(yNew-yK))
    # make sure particle is on the mesh (recompute the z component)
    particles, _ = geoTrans.projectOnRasterVect(dem, particles, interp='bilinear')

    nx, ny, nz = DFAtls.getNormalArray(xNew, yNew, Nx, Ny, Nz, csz)
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

    # ++++++++++++++GET particles location (neighbours for sph)
    startTime = time.time()
    particles = DFAfunC.getNeighboursC(particles, dem)
    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

    # ++++++++++++++UPDATE FIELDS (compute grid values)
    # update fields (compute grid values)
    startTime = time.time()
    # particles, fields = updateFields(cfg, particles, force, dem, fields)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)
    tcpuField = time.time() - startTime
    Tcpu['Field'] = Tcpu['Field'] + tcpuField

    return particles, fields, Tcpu, dt


def prepareArea(releaseLine, dem):
    """ convert shape file polygon to raster

    Parameters
    ----------
    releaseLine: dict
        line dictionary
    dem : dict
        dictionary with dem information
    Returns
    -------

    Raster : 2D numpy array
        raster
    """
    NameRel = releaseLine['Name']
    StartRel = releaseLine['Start']
    LengthRel = releaseLine['Length']
    Raster = np.zeros(np.shape(dem['rasterData']))

    for i in range(len(NameRel)):
        name = NameRel[i]
        start = StartRel[i]
        end = start + LengthRel[i]
        avapath = {}
        avapath['x'] = releaseLine['x'][int(start):int(end)]
        avapath['y'] = releaseLine['y'][int(start):int(end)]
        avapath['Name'] = name
        Raster = polygon2Raster(dem['header'], avapath, Raster)
    return Raster


def polygon2Raster(demHeader, Line, Mask):
    """ convert line to raster

    Parameters
    ----------
    demHeader: dict
        dem header dictionary
    Line : dict
        line dictionary
    Mask : 2D numpy array
        raster to update
    Returns
    -------

    Mask : 2D numpy array
        updated raster
    """
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
        x = np.arange(ncols) * csz
        y = np.arange(nrows) * csz
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        ax.set_title('Release area')
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        ref0, im = pU.NonUnifIm(ax, x, y, Mask, 'x [m]', 'y [m]',
                                extent=[x.min(), x.max(), y.min(), y.max()],
                                cmap=cmap, norm=None)
        ax.plot(xCoord0 * csz, yCoord0 * csz, 'r', label='release polyline')
        plt.legend()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        plt.show()

    return Mask


def plotPosition(fig, ax, particles, dem, data, Cmap, unit, plotPart=False, last=False):
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
    ax.set_title('t=%.2f s' % particles['t'])
    cmap, _, lev, norm, ticks = makePalette.makeColorMap(
        Cmap, 0.0, np.nanmax(data), continuous=True)
    cmap.set_under(color='w')
    ref0, im = pU.NonUnifIm(ax, xx, yy, data, 'x [m]', 'y [m]',
                         extent=[x.min(), x.max(), y.min(), y.max()],
                         cmap=cmap, norm=norm)

    Cp1 = ax.contour(X, Y, Z, levels=10, colors='k')
    pU.addColorBar(im, ax, ticks, unit)
    if plotPart:
        # ax.plot(x, y, '.b', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 1], y[NPPC == 1], '.c', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 4], y[NPPC == 4], '.b', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 9], y[NPPC == 9], '.r', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 16], y[NPPC == 16], '.m', linestyle='None', markersize=1)
        # load variation colormap
        variable = particles['h']
        cmap, _, _, norm, ticks = makePalette.makeColorMap(
            pU.cmapDepth, np.amin(variable), np.amax(variable), continuous=True)
        # set range and steps of colormap
        cc = variable
        sc = ax.scatter(x, y, c=cc, cmap=cmap, marker='.')

        if last:
            pU.addColorBar(sc, ax, ticks, 'm', 'Flow Depth')

    plt.pause(0.1)
    return fig, ax


def plotContours(fig, ax, particles, dem, data, Cmap, unit, last=False):
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
    try:
        # Get the images on an axis
        cb = ax.images[-1].colorbar
        if cb:
            cb.remove()
    except IndexError:
        pass

    ax.clear()
    ax.set_title('t=%.2f s' % particles['t'])
    cmap, _, lev, norm, ticks = makePalette.makeColorMap(
        Cmap, 0.0, np.nanmax(data), continuous=True)
    cmap.set_under(color='w')

    CS = ax.contour(X, Y, data, levels=8, origin='lower', cmap=cmap,
                    linewidths=2)
    lev = CS.levels

    if last:
        # pU.addColorBar(im, ax, ticks, unit, 'Flow Depth')
        CB = fig.colorbar(CS)
        ax.clabel(CS, inline=1, fontsize=8)

    plt.pause(0.1)
    return fig, ax, cmap, lev


def removeOutPart(cfg, particles, dem):
    """ find and remove out of raster particles

    Parameters
    ----------
    cfg : configparser
        DFA parameters
    particles : dict
        particles dictionary
    dem : dict
        dem dictionary

    Returns
    -------
    particles : dict
        particles dictionary
    """
    dt = cfg.getfloat('dt')
    header = dem['header']
    nrows = header.nrows
    ncols = header.ncols
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    Bad = dem['Bad']

    x = particles['x']
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    indX = particles['indX']
    indY = particles['indY']
    x = x + ux*dt
    y = y + uy*dt

    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    Lx = (x - xllc) / csz
    Ly = (y - yllc) / csz
    mask = np.ones(len(x), dtype=bool)
    indOut = np.where(Lx <= 1.5)
    mask[indOut] = False
    indOut = np.where(Ly <= 1.5)
    mask[indOut] = False
    indOut = np.where(Lx >= ncols-1.5)
    mask[indOut] = False
    indOut = np.where(Ly >= nrows-1.5)
    mask[indOut] = False

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they exited the domain' % (nRemove))

    x = particles['x']
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    mask = np.ones(len(x), dtype=bool)
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    indX = particles['indX']
    indY = particles['indY']
    indOut = np.where(Bad[indY, indX], False, True)
    mask = np.logical_and(mask, indOut)
    indOut = np.where(Bad[indY+np.sign(uy).astype('int'), indX], False, True)
    mask = np.logical_and(mask, indOut)
    indOut = np.where(Bad[indY, indX+np.sign(ux).astype('int')], False, True)
    mask = np.logical_and(mask, indOut)
    indOut = np.where(Bad[indY+np.sign(uy).astype('int'), indX+np.sign(ux).astype('int')], False, True)
    mask = np.logical_and(mask, indOut)

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they exited the domain' % (nRemove))

    return particles


def removeSmallPart(hmin, particles, dem):
    """ find and remove too small particles

    Parameters
    ----------
    hmin : float
        minimum depth
    particles : dict
        particles dictionary
    dem : dict
        dem dictionary

    Returns
    -------
    particles : dict
        particles dictionary
    """
    h = particles['h']

    indOut = np.where(h < hmin)
    mask = np.ones(len(h), dtype=bool)
    mask[indOut] = False

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they were too thin' % (nRemove))
        particles = DFAfunC.getNeighboursC(particles, dem)

    return particles


def removePart(particles, mask, nRemove):
    """ remove given particles

    Parameters
    ----------
    particles : dict
        particles dictionary
    mask : 1D numpy array
        particles to keep
    nRemove : int
        number of particles removed

    Returns
    -------
    particles : dict
        particles dictionary
    """
    particles['Npart'] = particles['Npart'] - nRemove
    particles['NPPC'] = particles['NPPC'][mask]
    particles['x'] = particles['x'][mask]
    particles['y'] = particles['y'][mask]
    particles['z'] = particles['z'][mask]
    particles['s'] = particles['s'][mask]
    particles['ux'] = particles['ux'][mask]
    particles['uy'] = particles['uy'][mask]
    particles['uz'] = particles['uz'][mask]
    particles['m'] = particles['m'][mask]
    particles['h'] = particles['h'][mask]
    particles['InCell'] = particles['InCell'][mask]
    particles['indX'] = particles['indX'][mask]
    particles['indY'] = particles['indY'][mask]
    particles['partInCell'] = particles['partInCell'][mask]

    return particles


def savePartToPickle(dictList, outDir):
    """ Save each dictionary from a list to a pickle in outDir; works also for one dictionary instead of list

        Parameters
        ---------
        dictList: list or dict
            list of dictionaries or single dictionary
        outDir: str
            path to output directory
    """
    if len(dictList) > 1:
        for dict in dictList:
            pickle.dump(dict, open(os.path.join(outDir, "particles%f.p" % dict['t']), "wb"))
    else:
        pickle.dump(dictList, open(os.path.join(outDir, "particles%f.p" % dictList['t']), "wb"))


def readPartFromPickle(inDir, flagAvaDir=False):
    """ Read pickles within a directory and return List of dicionaries read from pickle

        Parameters
        -----------
        inDir: str
            path to input directory
        flagAvaDir: bool
            if True inDir corresponds to an avalanche directory and pickles are
            read from avaDir/Outputs/com1DFAPy/particles 
    """

    if flagAvaDir:
        inDir = os.path.join(inDir, 'Outputs', 'com1DFAPy', 'particles')

    # search for all pickles within directory
    PartDicts = glob.glob(os.path.join(inDir, '*.p'))

    # initialise list of particle dictionaries
    Particles = []
    TimeStepInfo = []
    for particles in PartDicts:
        particles = pickle.load(open(particles, "rb"))
        Particles.append(particles)
        TimeStepInfo.append(particles['t'])

    return Particles, TimeStepInfo


def exportFields(cfgGen, Tsave, Fields, relFile, demOri, outDir):
    """ export result fields to Outputs directory according to result parameters and time step
        that can be specified in the configuration file

        Parameters
        -----------
        cfgGen: dict
            configurations
        Tsave: list
            list of time step that corresponds to each dict in Fields
        Fields: list
            list of Fields for each dtSave
        relFile: str
            path to release area shapefile
        outDir: str
            outputs Directory


        Returns
        --------
        exported peak fields are saved in Outputs/com1DFAPy/peakFiles

    """

    resTypesString = cfgGen['resType']
    resTypes = resTypesString.split('_')
    tSteps = fU.getTimeIndex(cfgGen, Fields)
    for tStep in tSteps:
        finalFields = Fields[tStep]
        for resType in resTypes:
            resField = finalFields[resType]
            if resType == 'ppr':
                resField = resField * 0.001
            relName = os.path.splitext(os.path.basename(relFile))[0]
            dataName = relName + '_' + 'null' + '_' + 'dfa' + '_' + '0.155' + '_' + resType + '_'  + 't%.2f' % (Tsave[tStep]) +'.asc'
            # create directory
            outDirPeak = os.path.join(outDir, 'peakFiles')
            fU.makeADir(outDirPeak)
            outFile = os.path.join(outDirPeak, dataName)
            IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
            if tStep == -1:
                log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f - FINAL time step ' % (resType,Tsave[tStep]))
            else:
                log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f ' % (resType,Tsave[tStep]))


def analysisPlots(Particles, Fields, cfg, demOri, dem, outDir):
    """ create analysis plots during simulation run """

    cfgGen = cfg['GENERAL']
    partRef = Particles[0]
    Z0 = partRef['z'][0]
    rho = cfgGen.getfloat('rho')
    gravAcc = cfgGen.getfloat('gravAcc')
    mu = cfgGen.getfloat('mu')
    repeat = True
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
            fig, ax = plotPosition(
                fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True)
            fig.savefig(os.path.join(outDir, 'particlest%f.%s' % (part['t'], pU.outputFormat)))

        fig, ax = plotPosition(
                fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True, last=True)
        fig.savefig(os.path.join(outDir, 'particlesFinal.%s' % (pU.outputFormat)))
        value = input("[y] to repeat:\n")
        if value != 'y':
            repeat = False

    fieldEnd = Fields[-1]
    partEnd = Particles[-1]
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig3, ax3 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig1, ax1 = plotPosition(
        fig1, ax1, partEnd, demOri, fieldEnd['FD'], pU.cmapPres, 'm', plotPart=False)
    fig2, ax2 = plotPosition(
        fig2, ax2, partEnd, demOri, fieldEnd['FV'], pU.cmapPres, 'm/s', plotPart=False)
    fig3, ax3 = plotPosition(
        fig3, ax3, partEnd, demOri, fieldEnd['P']/1000, pU.cmapPres, 'kPa', plotPart=False)
    plt.show()
