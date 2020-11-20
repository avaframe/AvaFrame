"""
    function regarding time discretization and time stepping for com1DFA

    This file is part of Avaframe.
"""

# Load modules
import copy
import logging
import math
import numpy as np

# Local imports
import avaframe.com1DFAPy.DFAtools as DFAtls


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


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
                r = DFAtls.norm(dx, dy, dz)
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


def calcGradHSPHVect(particles, j, ncols, nrows, csz, nx, ny, nz):
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
    dz = (particles['z'][L] - z)
    dn = nx*dx + ny*dy + nz*dz
    dx = dx - dn*nx
    dy = dy - dn*ny
    dz = dz - dn*nz
    r = DFAtls.norm(dx, dy, dz)
    # impose a minimum distance between particles
    dx = np.where(r < 0.0001 * rKernel, 0.0001 * rKernel * dx, dx)
    dy = np.where(r < 0.0001 * rKernel, 0.0001 * rKernel * dy, dy)
    dz = np.where(r < 0.0001 * rKernel, 0.0001 * rKernel * dz, dz)
    r = np.where(r < 0.0001 * rKernel, 0.0001 * rKernel, r)
    hr = rKernel - r
    dwdr = dfacKernel * hr * hr
    massl = particles['m'][L]
    # ------------------------------
    indOut = np.where(r >= rKernel)
    massl[indOut] = 0
    mdwdrr = massl * dwdr / r
    GX = mdwdrr * dx
    GY = mdwdrr * dy
    GZ = mdwdrr * dz
    # -----------------------------
    gradhX = gradhX + np.sum(GX)
    gradhY = gradhY + np.sum(GY)
    gradhZ = gradhZ + np.sum(GZ)
    # leInd = len(index)
    # print((time.time() - startTime))

    return gradhX, gradhY,  gradhZ, L


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
    R10 = DFAtls.norm(csz-dx, dy, dz)
    dz = z - Z[Lx0, Ly1]
    R01 = DFAtls.norm(dx, csz-dy, dz)
    dz = z - Z[Lx1, Ly1]
    R11 = DFAtls.norm(csz-dx, csz-dy, dz)

    F = F.flatten()
    Z = Z.flatten()

    # add sph contribution to bottom left grid cell
    ic00 = Lx0 + ncol * Ly0
    dz = z - Z[ic00]
    R00 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R00 >= rKernel)
    R00[indOut] = 0
    hr00 = rKernel - R00
    W00 = facKernel * hr00 * hr00 * hr00
    f00 = f * m * W00 / rho
    np.add.at(F, ic00, f00)

    # add sph contribution to bottom right grid cell
    ic10 = Lx1 + ncol * Ly0
    dz = z - Z[ic10]
    R10 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R10 >= rKernel)
    R10[indOut] = 0
    hr10 = rKernel - R10
    W10 = facKernel * hr10 * hr10 * hr10
    f10 = f * m * W10 / rho
    np.add.at(F, ic10, f10)

    # add sph contribution to top left grid cell
    ic01 = Lx0 + ncol * Ly1
    dz = z - Z[ic01]
    R01 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R01 >= rKernel)
    R01[indOut] = 0
    hr01 = rKernel - R01
    W01 = facKernel * hr01 * hr01 * hr01
    f01 = f * m * W01 / rho
    np.add.at(F, ic01, f01)

    # add sph contribution to top right grid cell
    ic11 = Lx1 + ncol * Ly1
    dz = z - Z[ic11]
    R11 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R11 >= rKernel)
    R11[indOut] = 0
    hr11 = rKernel - R10
    W11 = facKernel * hr11 * hr11 * hr11
    f11 = f * m * W11 / rho
    np.add.at(F, ic11, f11)

    F = np.reshape(F, (nrow, ncol))

    return F


def computeFD(cfg, particles, dem):
    rho = cfg.getfloat('rho')
    Npart = particles['Npart']
    nrows = dem['header'].nrows
    ncols = dem['header'].ncols
    csz = dem['header'].cellsize
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    # initialize
    H = np.zeros(np.shape(particles['h']))
    # loop on particles
    # TcpuSPH = 0
    # Tcpuadd = 0
    for j in range(Npart):
        # adding lateral force (SPH component)
        # startTime = time.time()
        # gradhX, gradhY,  gradhZ, _ = calcGradHSPH(particles, j, ncols, nrows, csz)

        # get normal at the particle location
        x = particles['x'][j]
        y = particles['y'][j]
        nx, ny, nz = DFAtls.getNormal(x, y, Nx, Ny, Nz, csz)

        h, _ = calcHSPHVect(particles, j, ncols, nrows, csz, nx, ny, nz)
        H[j] = h / rho

    particles['hSPH'] = H

    return particles


def calcHSPHVect(particles, j, ncols, nrows, csz, nx, ny, nz):
    # SPH kernel
    # use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
    rKernel = csz
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))

    indx, indy, _ = particles['InCell'][j]
    indPartInCell = particles['indPartInCell']
    partInCell = particles['partInCell']
    x = particles['x'][j]
    y = particles['y'][j]
    z = particles['z'][j]
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
    dz = (particles['z'][L] - z)
    dn = nx*dx + ny*dy + nz*dz
    dx = dx - dn*nx
    dy = dy - dn*ny
    dz = dz - dn*nz
    r = DFAtls.norm(dx, dy, dz)
    # impose a minimum distance between particles
    r = np.where(r < 0.0001 * rKernel, 0.0001 * rKernel, r)
    hr = rKernel - r
    w = facKernel * hr * hr * hr
    massl = particles['m'][L]
    indOut = np.where(r >= rKernel)
    massl[indOut] = 0
    H = massl * w
    # -----------------------------
    h = h + np.sum(H)
    # leInd = len(index)
    # print((time.time() - startTime)/leInd)

    return h, L
