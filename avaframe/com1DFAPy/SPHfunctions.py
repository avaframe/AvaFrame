"""
    function related to SPH calculations in com1DFA
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

# how to compute SPH gradient:
# 1) like in SAMOS AT just project on the plane
# 2) project on the plane, compute the gradient in the local coord sys related
# to the local plane and flow direction (orthogonal coord sys)
# enables to choose earth pressure coefficients
# 3) project on the plane, compute the gradient in the local coord sys related
# to the local plane (tau1, tau2, n) non orthogonal coord sys
SPHoption = 2


def getNeighbours(particles, dem):
    """ Locate particles in cell for SPH computation (for loop implementation)

    Ĺocate each particle in a grid cell and build the indPartInCell and
    partInCell arrays. See issue #200 and documentation for details

    Parameters
    ----------
    particles : dict
    dem : dict

    Returns
    -------
    particles : dict
      updated particles dictionary with indPartInCell and partInCell arrays
    """
    # get grid information
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    # get particle location
    Npart = particles['Npart']
    x = particles['x']
    y = particles['y']

    # initialize outputs
    indPartInCell = np.zeros(ncols*nrows + 1).astype(int)
    partInCell = np.zeros(Npart).astype(int)
    # Count number of particles in each cell
    indx = ((x + csz/2) / csz).astype(int)
    indy = ((y + csz/2) / csz).astype(int)
    # get index of cell containing the particle
    ic = indx + ncols * indy
    # count particles in each cell
    ##################################
    # start ########################
    # TODO: test speed between add.at and bincount
    # indPartInCell = np.bincount(ic + 1, minlength=len(indPartInCell))
    # or
    np.add.at(indPartInCell, ic + 1, 1)
    # end ##############################
    ##################################
    # make the cumulative sum out of it
    indPartInCell = np.cumsum(indPartInCell)

    # make the list of which particles are in which cell
    indX = np.empty((0), int).astype(int)
    indY = np.empty((0), int).astype(int)
    InCell = np.empty((0), int).astype(int)
    indPartInCell2 = copy.deepcopy(indPartInCell)
    for j in range(Npart):
        indx = int((x[j] + csz/2) / csz)
        indy = int((y[j] + csz/2) / csz)
        ic = indx + ncols * indy
        partInCell[int(indPartInCell2[ic])] = j
        indPartInCell2[ic] = indPartInCell2[ic] + 1
        indX = np.append(indX, indx)
        indY = np.append(indY, indy)
        InCell = np.append(InCell, ic)
        # InCell = np.append(InCell, np.tile(np.array([indx, indy, ic]), (1, 1)), axis=0)

    particles['indX'] = indX
    particles['indY'] = indY
    particles['InCell'] = InCell
    particles['indPartInCell'] = indPartInCell
    particles['partInCell'] = partInCell

    return particles


def getNeighboursVect(particles, dem):
    """ Locate particles in cell for SPH computation (no loop implementation)

    Ĺocate each particle in a grid cell and build the indPartInCell and
    partInCell arrays. See issue #200 and documentation for details.
    Without using loops.

    Parameters
    ----------
    particles : dict
    dem : dict

    Returns
    -------
    particles : dict
      updated particles dictionary with indPartInCell and partInCell arrays
    """
    # get grid information
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    # get particle location
    Npart = particles['Npart']
    x = particles['x']
    y = particles['y']

    # initialize outputs
    indPartInCell = np.zeros(ncols*nrows + 1).astype(int)
    partInCell = np.zeros(Npart).astype(int)
    # Count number of particles in each cell
    indx = ((x + csz/2) / csz).astype(int)
    indy = ((y + csz/2) / csz).astype(int)
    # get index of cell containing the particle
    ic = indx + ncols * indy
    # count particles in each cell
    ##################################
    # start ########################
    # TODO: test speed between add.at and bincount
    # indPartInCell = np.bincount(ic + 1, minlength=len(indPartInCell))
    # or
    np.add.at(indPartInCell, ic + 1, 1)
    # end ##############################
    ##################################
    # make the cumulative sum out of it
    indPartInCell = np.cumsum(indPartInCell)

    # make the list of which particles are in which cell
    partInCell = np.argsort(ic, kind='mergesort')

    particles['indX'] = indx
    particles['indY'] = indy
    particles['InCell'] = ic
    particles['indPartInCell'] = indPartInCell
    particles['partInCell'] = partInCell

    return particles


def calcGradHSPH(particles, j, ncols, nrows, csz, minRKern):
    """ Compute gradient of Flow Depth using SPH (for loop implementation)

    Parameters
    ----------
    particles : dict
    j: int
        index of particle under consideration
    ncols: int
        number of columns of the DEM
    nrows: int
        number of rows of the DEM
    csz  : float
        cellsize of the DEM

    Returns
    -------
    gradhX: float
        x coordinate of the gradient of the flow depth at particle j location
    gradhY: float
        x coordinate of the gradient of the flow depth at particle j location
    gradhZ: float
        x coordinate of the gradient of the flow depth at particle j location
    L: 1D numpy array
        index of particles within the kernel function radius
    """
    # SPH kernel
    # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
    rKernel = csz
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))
    dfacKernel = -3.0 * facKernel

    indx = particles['indX'][j]
    indy = particles['indY'][j]
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
                if r < minRKern * rKernel:
                    # impose a minimum distance between particles
                    r = minRKern * rKernel
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
    """ Compute gradient of Flow Depth using SPH (no loop implementation)

    Compute the gradient of the flow depth at the location of patricle j
    using SPH method.

    Parameters
    ----------
    particles : dict
    j: int
        index of particle under consideration
    ncols: int
        number of columns of the DEM
    nrows: int
        number of rows of the DEM
    csz  : float
        cellsize of the DEM

    Returns
    -------
    gradhX: float
        x coordinate of the gradient of the flow depth at particle j location
    gradhY: float
        x coordinate of the gradient of the flow depth at particle j location
    gradhZ: float
        x coordinate of the gradient of the flow depth at particle j location
    L: 1D numpy array
        index of particles within the kernel function radius
    """
    # Define SPH kernel
    # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
    rKernel = csz
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))
    dfacKernel = -3.0 * facKernel

    # get particle j information
    indx = particles['indX'][j]
    indy = particles['indY'][j]
    indPartInCell = particles['indPartInCell']
    partInCell = particles['partInCell']
    x = particles['x'][j]
    y = particles['y'][j]
    z = particles['z'][j]

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
    icMin = ncols * (indy + np.arange(lInd, rInd, 1))
    icMax = ncols * (indy + np.arange(lInd, rInd, 1) + 1)
    iPstart = indPartInCell[np.maximum(ic, icMin)]
    iPend = indPartInCell[np.minimum(ic + 3, icMax)]

    # gather all the particles indexes
    index = np.concatenate([np.arange(x, y) for x, y in zip(iPstart, iPend)])
    # gather all particles
    L = partInCell[index]
    # make sure to remove the j particle
    indJ = np.where(L == j)
    L = np.delete(L, indJ)

    # compute SPH gradient
    # get xj - xl
    dx = particles['x'][L] - x
    dy = particles['y'][L] - y
    dz = (particles['z'][L] - z)

    if SPHoption == 1 or SPHoption == 3:
        # remove the normal part (make sure that r = xj - xl lies in the plane
        # defined by the normal at xj)
        dn = nx*dx + ny*dy + nz*dz
        dx = dx - dn*nx
        dy = dy - dn*ny
        dz = dz - dn*nz

        # get norm of r = xj - xl
        r = DFAtls.norm(dx, dy, dz)
        # impose a minimum distance between particles
        dx = np.where(r < minRKern * rKernel, minRKern * rKernel * dx, dx)
        dy = np.where(r < minRKern * rKernel, minRKern * rKernel * dy, dy)
        dz = np.where(r < minRKern * rKernel, minRKern * rKernel * dz, dz)
        r = np.where(r < minRKern * rKernel, minRKern * rKernel, r)

    elif SPHoption == 2:
        # Option 2
        # projecting onto the tengent plane and taking the change of coordinates into account
        # the coord sysem used is the orthogonal coord system related to the flow
        # (e1, e2, n=e3)
        # this way it is possible to include the earth pressure coeficients

        # get velocity magnitude and direction
        ux = particles['ux'][j]
        uy = particles['uy'][j]
        uz = particles['uz'][j]
        uMag = DFAtls.norm(ux, uy, uz)
        if uMag < 0.1:
            ax = 1
            ay = 0
            uxDir = ax
            uyDir = ay
            uzDir = -(ax*nx + ay*ny) / nz
            uxDir, uyDir, uzDir = DFAtls.normalize(uxDir, uyDir, uzDir)
        else:
            # TODO check if direction is non zero, if it is define another u1 direction
            uxDir, uyDir, uzDir = DFAtls.normalize(ux, uy, uz)

        uxOrtho, uyOrtho, uzOrtho = DFAtls.croosProd(nx, ny, nz, uxDir, uyDir, uzDir)
        uxOrtho, uyOrtho, uzOrtho = DFAtls.normalize(uxOrtho, uyOrtho, uzOrtho)

        v1 = np.array([[uxDir, uyDir, uzDir]])
        v2 = np.array([[uxOrtho, uyOrtho, uzOrtho]])
        v3 = np.array([[nx, ny, nz]])
        # build the transformation matrix from (e1, e2, e3) to (ex, ey, ez)
        M = np.concatenate((v1.T, v2.T), axis=1)
        M = np.concatenate((M, v3.T), axis=1)
        # compute the transformation matrix from (ex, ey, ez) to (e1, e2, e3)
        MM1 = M.T  # because M is orthogonal, it inverse is its transpose !!! np.linalg.inv(M).T

        # now take into accout the fact that we are on the surface so the r3 or x3
        # component is not independent from the 2 other ones!!
        MM1[0, 0] = MM1[0, 0] - M[2, 0]*M[0, 2]/M[2, 2]
        MM1[0, 1] = MM1[0, 1] - M[2, 0]*M[1, 2]/M[2, 2]
        MM1[1, 0] = MM1[1, 0] - M[2, 1]*M[0, 2]/M[2, 2]
        MM1[1, 1] = MM1[1, 1] - M[2, 1]*M[1, 2]/M[2, 2]

        # buil the matrix that transforms the gradient in (r1, r2, r3) to the one in (x1, x2, x3)
        MMGrad = MM1.T

        # get coordinates in local coord system
        r1 = DFAtls.scalProd(dx, dy, dz, uxDir, uyDir, uzDir)
        r2 = DFAtls.scalProd(dx, dy, dz, uxOrtho, uyOrtho, uzOrtho)
        # impse r3=0 even if the particle is not exactly on the tengent plane
        # get norm of r = xj - xl
        r = DFAtls.norm(r1, r2, 0)
        # impose a minimum distance between particles
        r1 = np.where(r < minRKern * rKernel, minRKern * rKernel * r1, r1)
        r2 = np.where(r < minRKern * rKernel, minRKern * rKernel * r2, r2)
        r = np.where(r < minRKern * rKernel, minRKern * rKernel, r)

    # compute derivative of kernel function
    hr = rKernel - r
    dwdr = dfacKernel * hr * hr
    massl = particles['m'][L]
    # keep only the particles in the support of the kernel function
    indOut = np.where(r >= rKernel)
    massl[indOut] = 0
    mdwdrr = massl * dwdr / r

    # ----------------------------------------
    if SPHoption == 1:
        # Option 1
        # only projecting dr on the tangent plane
        GX = mdwdrr * dx
        GY = mdwdrr * dy
        GZ = mdwdrr * dz
        GX = np.sum(GX)
        GY = np.sum(GY)
        GZ = np.sum(GZ)

    elif SPHoption == 2:
        # Option 2
        # projecting onto the tengent plane and taking the change of coordinates into account
        # the coord sysem used is the orthogonal coord system related to the flow
        # (e1, e2, n=e3)
        # this way it is possible to include the earth pressure coeficients
        K1, K2 = 1, 1
        G1 = mdwdrr * K1*r1
        G2 = mdwdrr * K2*r2
        G3 = 0

        g1 = nx/(nz)
        g2 = ny/(nz)

        GX1 = MMGrad[0, 0]*G1 + MMGrad[0, 1]*G2
        GY1 = MMGrad[1, 0]*G1 + MMGrad[1, 1]*G2
        GZ1 = (- g1*GX1 - g2*GY1)

        GX = np.sum(GX1)
        GY = np.sum(GY1)
        GZ = np.sum(GZ1)

    elif SPHoption == 3:
        # Option 3
        # No proof yet....
        # projecting onto the tengent plane and taking the change of coordinates into account
        # the coord sysem used is the non orthogonal coord system related to the surface
        # (Tau1, Tau2, n)
        g1 = nx/(nz)
        g2 = ny/(nz)
        g12 = g1*g2
        # g33 = (1 + g1*g1 + g2*g2)
        g11 = 1 + g1*g1
        g22 = 1 + g2*g2

        GX = mdwdrr * (g11*dx + g12*dy)
        GY = mdwdrr * (g22*dy + g12*dx)
        # the gradient hat to be tangent to the plane...
        GZ = (- g1*GX - g2*GY)

        GX = np.sum(GX)
        GY = np.sum(GY)
        GZ = np.sum(GZ)

    # -----------------------------
    gradhX = GX
    gradhY = GY
    gradhZ = GZ

    # leInd = len(index)
    # print((time.time() - startTime))

    return gradhX, gradhY,  gradhZ, L


def pointsToRasterSPH(particles, dem, rho, f):
    """ SPH interpolation from particles (unstructured grid) to grid

    Interpolates the scalar field f associated to the particles in particles
    (equivalent of unstructured grid) on a structured grid (associated to
    the dem). Values of f on the structured grid are called F

    Parameters
    ----------
    particles : dict
    dem : dict
    rho: float
        fluid density
    f: 1D numpy array
        values of scalar field at particle location

    Returns
    -------
    F: 2d numpy array
    values of scalar field at grid point location
    """

    # get mesh information
    Z = dem['rasterData']
    header = dem['header']
    csz = header.cellsize
    xllc = header.xllcenter
    yllc = header.yllcenter
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    # get particle information (location and mass)
    x = particles['x']
    y = particles['y']
    z = particles['z']
    m = particles['m']

    # Define SPH kernel
    # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
    rKernel = csz
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))

    nrows, ncols = np.shape(Nx)
    F = np.zeros((nrows, ncols))

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
    Nx = Nx.flatten()
    Ny = Ny.flatten()
    Nz = Nz.flatten()

    for lx, ly in zip([Lx0, Lx0, Lx1, Lx1], [Ly0, Ly1, Ly0, Ly1]):
        dx = (lx - Lx) * csz
        dy = (ly - Ly) * csz
        ic00 = lx + ncols * ly
        dz = z - Z[ic00]
        # remove the normal part (make sure that r = xj - xl lies in the plane
        # defined by the normal at xj)
        dn = Nx*dx + Ny*dy + Nz*dz
        dx = dx - dn*Nx
        dy = dy - dn*Ny
        dz = dz - dn*Nz
        # add sph contribution to grid cell (lx, ly)
        R00 = DFAtls.norm(dx, dy, dz)
        indOut = np.where(R00 >= rKernel)
        R00[indOut] = 0
        hr00 = rKernel - R00
        W00 = facKernel * hr00 * hr00 * hr00
        f00 = f * m * W00 / rho
        np.add.at(F, ic00, f00)

    # add sph contribution to bottom left grid cell
    ic00 = Lx0 + ncols * Ly0
    dz = z - Z[ic00]
    R00 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R00 >= rKernel)
    R00[indOut] = 0
    hr00 = rKernel - R00
    W00 = facKernel * hr00 * hr00 * hr00
    f00 = f * m * W00 / rho
    np.add.at(F, ic00, f00)

    # add sph contribution to bottom right grid cell
    ic10 = Lx1 + ncols * Ly0
    dz = z - Z[ic10]
    R10 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R10 >= rKernel)
    R10[indOut] = 0
    hr10 = rKernel - R10
    W10 = facKernel * hr10 * hr10 * hr10
    f10 = f * m * W10 / rho
    np.add.at(F, ic10, f10)

    # add sph contribution to top left grid cell
    ic01 = Lx0 + ncols * Ly1
    dz = z - Z[ic01]
    R01 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R01 >= rKernel)
    R01[indOut] = 0
    hr01 = rKernel - R01
    W01 = facKernel * hr01 * hr01 * hr01
    f01 = f * m * W01 / rho
    np.add.at(F, ic01, f01)

    # add sph contribution to top right grid cell
    ic11 = Lx1 + ncols * Ly1
    dz = z - Z[ic11]
    R11 = DFAtls.norm(dx, dy, dz)
    indOut = np.where(R11 >= rKernel)
    R11[indOut] = 0
    hr11 = rKernel - R10
    W11 = facKernel * hr11 * hr11 * hr11
    f11 = f * m * W11 / rho
    np.add.at(F, ic11, f11)

    F = np.reshape(F, (nrows, ncols))

    return F


def computeFlowDepth(cfg, particles, dem):
    """ Compute Flow Depth using SPH (no loop implementation)

    Compute the flow depth at the location of patricle j
    using SPH method.

    Parameters
    ----------
    cfg: configParser
    paricles: dict
    dem: dict

    Returns
    -------
    paricles: dict
        particles dictionnary updated with the SPH flow depth (hSPH)
    """
    rho = cfg.getfloat('rho')
    minRKern = cfg.getfloat('minRKern')
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
        # minRKern = cfg.getfloat('minRKern')
        # gradhX, gradhY,  gradhZ, _ = calcGradHSPH(particles, j, ncols, nrows, csz, minRKern)

        # get normal at the particle location
        x = particles['x'][j]
        y = particles['y'][j]
        nx, ny, nz = DFAtls.getNormal(x, y, Nx, Ny, Nz, csz)

        h, _ = calcHSPHVect(particles, j, ncols, nrows, csz, nx, ny, nz, minRKern)
        H[j] = h / rho

    particles['hSPH'] = H

    return particles


def calcHSPHVect(particles, j, ncols, nrows, csz, nx, ny, nz, minRKern):
    """ Compute Flow Depth using SPH (no loop implementation)

    Compute the flow depth at the location of patricle j
    using SPH method.

    Parameters
    ----------
    particles : dict
    j: int
        index of particle under consideration
    ncols: int
        number of columns of the DEM
    nrows: int
        number of rows of the DEM
    csz  : float
        cellsize of the DEM
    nx  : float
        x component of the surface normal at particle j position
    ny  : float
        y component of the surface normal at particle j position
    nz  : float
        z component of the surface normal at particle j position

    Returns
    -------
    h: float
        flow depth at particle j location
    L: 1D numpy array
        index of particles within the kernel function radius
    """
    # SPH kernel
    # use "spiky" kernel: w = (rKernel - r)**3 * 10/(pi*rKernel**5)
    rKernel = csz
    facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))

    indx = particles['indX'][j]
    indy = particles['indY'][j]
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
    r = np.where(r < minRKern * rKernel, minRKern * rKernel, r)
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
