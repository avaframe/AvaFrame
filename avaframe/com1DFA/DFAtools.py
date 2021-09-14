"""
    Basic tools for getting grid normals, area and working with vectors.
"""

# Load modules
import logging
import numpy as np
import numbers

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getNormalArray(x, y, Nx, Ny, Nz, csz):
    """ Interpolate vector field from grid to unstructures points

        Originaly created to get the normal vector at location (x,y) given the
        normal vector field on the grid. Grid has its origin in (0,0).
        Can be used to interpolate any vector field.
        Interpolation using a bilinear interpolation

        Parameters
        ----------
            x: numpy array
                location in the x location of desiered interpolation
            y: numpy array
                location in the y location of desiered interpolation
            Nx: 2D numpy array
                x component of the vector field at the grid nodes
            Ny: 2D numpy array
                y component of the vector field at the grid nodes
            Nz: 2D numpy array
                z component of the vector field at the grid nodes
            csz: float
                cellsize of the grid

        Returns
        -------
            nx: numpy array
                x component of the interpolated vector field at position (x, y)
            ny: numpy array
                y component of the interpolated vector field at position (x, y)
            nz: numpy array
                z component of the interpolated vector field at position (x, y)
    """
    nrow, ncol = np.shape(Nx)
    # by default bilinear interpolation of the Nx, Ny, Nz of the grid
    nx, _ = geoTrans.projectOnGrid(x, y, Nx, csz=csz)
    ny, _ = geoTrans.projectOnGrid(x, y, Ny, csz=csz)
    nz, _ = geoTrans.projectOnGrid(x, y, Nz, csz=csz)
    return nx, ny, nz


def getNormalMesh(dem, num):
    """ Compute normal to surface at grid points

        Get the normal vectors to the surface defined by a DEM.
        Either by adding the normal vectors of the adjacent triangles for each
        points (using 4, 6 or 8 adjacent triangles). Or use the next point in
        x direction and the next in y direction to define two vectors and then
        compute the cross product to get the normal vector

        Parameters
        ----------
            dem: dict
                header :
                    dem header (cellsize, ncols, nrows)
                rasterData : 2D numpy array
                    elevation at grid points
            num: int
                chose between 4, 6 or 8 (using then 4, 6 or 8 triangles) or
                1 to use the simple cross product method (with the diagonals)

        Returns
        -------
            Nx: 2D numpy array
                x component of the normal vector field on grid points
            Ny: 2D numpy array
                y component of the normal vector field on grid points
            Nz: 2D numpy array
                z component of the normal vector field on grid points
    """
    # read dem header
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    csz = header['cellsize']
    # read rasterData
    z = dem['rasterData']
    n, m = np.shape(z)
    Nx = np.ones((n, m))
    Ny = np.ones((n, m))
    Nz = np.ones((n, m))
    # first and last row, first and last column are inacurate
    if num == 4:
        # filling the inside of the matrix
        # normal calculation with 4 triangles
        # (Zl - Zr) / csz
        Nx[1:n-1, 1:m-1] = (z[1:n-1, 0:m-2] - z[1:n-1, 2:m]) / csz
        # (Zd - Zu) * csz
        Ny[1:n-1, 1:m-1] = (z[0:n-2, 1:m-1] - z[2:n, 1:m-1]) / csz
        Nz = 2 * Nz
        # filling the first col of the matrix
        # -2*(Zr - Zp) / csz
        Nx[1:n-1, 0] = - 2*(z[1:n-1, 1] - z[1:n-1, 0]) / csz
        # (Zd - Zu) / csz
        Ny[1:n-1, 0] = (z[0:n-2, 0] - z[2:n, 0]) / csz
        # filling the last col of the matrix
        # 2*(Zl - Zp) / csz
        Nx[1:n-1, m-1] = 2*(z[1:n-1, m-2] - z[1:n-1, m-1]) / csz
        # (Zd - Zu) / csz
        Ny[1:n-1, m-1] = (z[0:n-2, m-1] - z[2:n, m-1]) / csz
        # filling the first row of the matrix
        # (Zl - Zr) / csz
        Nx[0, 1:m-1] = (z[0, 0:m-2] - z[0, 2:m]) / csz
        # -2*(Zu - Zp) / csz
        Ny[0, 1:m-1] = - 2*(z[1, 1:m-1] - z[0, 1:m-1]) / csz
        # filling the last row of the matrix
        # (Zl - Zr) / csz
        Nx[n-1, 1:m-1] = (z[n-1, 0:m-2] - z[n-1, 2:m]) / csz
        # 2*(Zd - Zp) / csz
        Ny[n-1, 1:m-1] = 2*(z[n-2, 1:m-1] - z[n-1, 1:m-1]) / csz
        # filling the corners of the matrix
        Nx[0, 0] = -(z[0, 1] - z[0, 0]) / csz
        Ny[0, 0] = -(z[1, 0] - z[0, 0]) / csz
        Nz[0, 0] = 1
        Nx[n-1, 0] = -(z[n-1, 1] - z[n-1, 0]) / csz
        Ny[n-1, 0] = (z[n-2, 0] - z[n-1, 0]) / csz
        Nz[n-1, 0] = 1
        Nx[0, m-1] = (z[0, m-2] - z[0, m-1]) / csz
        Ny[0, m-1] = -(z[1, m-1] - z[0, m-1]) / csz
        Nz[0, m-1] = 1
        Nx[n-1, m-1] = (z[n-1, m-2] - z[n-1, m-1]) / csz
        Ny[n-1, m-1] = (z[n-2, m-1] - z[n-1, m-1]) / csz
        Nz[n-1, m-1] = 1

    if num == 6:
        # filling the inside of the matrix
        # normal calculation with 6 triangles
        # (2*(Zl - Zr) - Zur + Zdl + Zu - Zd) / csz
        Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m])
                            - z[2:n, 2:m] + z[0:n-2, 0:m-2]
                            + z[2:n, 1:m-1] - z[0:n-2, 1:m-1]) / csz
        # (2*(Zd - Zu) - Zur + Zdl - Zl + Zr) / csz
        Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1])
                            - z[2:n, 2:m] + z[0:n-2, 0:m-2]
                            - z[1:n-1, 0:m-2] + z[1:n-1, 2:m]) / csz
        Nz = 6 * Nz
        # filling the first col of the matrix
        # (- 2*(Zr - Zp) + Zu - Zur ) / csz
        Nx[1:n-1, 0] = (- 2*(z[1:n-1, 1] - z[1:n-1, 0]) + z[2:n, 0] - z[2:n, 1]) / csz
        # (Zd - Zu + Zr - Zur) / csz
        Ny[1:n-1, 0] = (z[0:n-2, 0] - z[2:n, 0] + z[1:n-1, 1] - z[2:n, 1]) / csz
        Nz[1:n-1, 0] = 3
        # filling the last col of the matrix
        # (2*(Zl - Zp) + Zdl - Zd) / csz
        Nx[1:n-1, m-1] = (2*(z[1:n-1, m-2] - z[1:n-1, m-1]) + z[0:n-2, m-2] - z[0:n-2, m-1]) / csz
        # (Zd - Zu + Zdl - Zl) / csz
        Ny[1:n-1, m-1] = (z[0:n-2, m-1] - z[2:n, m-1] + z[0:n-2, m-2] - z[1:n-1, m-2]) / csz
        Nz[1:n-1, m-1] = 3
        # filling the first row of the matrix
        # (Zl - Zr + Zu - Zur) / csz
        Nx[0, 1:m-1] = (z[0, 0:m-2] - z[0, 2:m] + z[1, 1:m-1] - z[1, 2:m]) / csz
        # (-2*(Zu - Zp) + Zr - Zur) / csz
        Ny[0, 1:m-1] = (- 2*(z[1, 1:m-1] - z[0, 1:m-1]) + z[0, 2:m] - z[1, 2:m]) / csz
        Nz[0, 1:m-1] = 3
        # filling the last row of the matrix
        # (Zl - Zr + Zdl - Zd) / csz
        Nx[n-1, 1:m-1] = (z[n-1, 0:m-2] - z[n-1, 2:m] + z[n-2, 0:m-2] - z[n-2, 1:m-1]) / csz
        # (2*(Zd - Zp) + Zdl - Zl) / csz
        Ny[n-1, 1:m-1] = (2*(z[n-2, 1:m-1] - z[n-1, 1:m-1]) + z[n-2, 0:m-2] - z[n-1, 0:m-2]) / csz
        Nz[n-1, 1:m-1] = 3
        # filling the corners of the matrix
        Nx[0, 0] = (z[1, 0] - z[1, 1] - (z[0, 1] - z[0, 0])) / csz
        Ny[0, 0] = (z[0, 1] - z[1, 1] - (z[1, 0] - z[0, 0])) / csz
        Nz[0, 0] = 2
        Nx[n-1, 0] = -(z[n-1, 1] - z[n-1, 0]) / csz
        Ny[n-1, 0] = (z[n-2, 0] - z[n-1, 0]) / csz
        Nz[n-1, 0] = 1
        Nx[0, m-1] = (z[0, m-2] - z[0, m-1]) / csz
        Ny[0, m-1] = -(z[1, m-1] - z[0, m-1]) / csz
        Nz[0, m-1] = 1
        Nx[n-1, m-1] = (z[n-1, m-2] - z[n-1, m-1] + z[n-2, m-2] - z[n-2, m-1]) / csz
        Ny[n-1, m-1] = (z[n-2, m-1] - z[n-1, m-1] + z[n-2, m-2] - z[n-1, m-2]) / csz
        Nz[n-1, m-1] = 2

    if num == 8:
        # filling the inside of the matrix
        # normal calculation with 8 triangles
        # (2*(Zl - Zr) + Zul - Zur + Zdl - Zdr) / csz
        Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m]) + z[2:n, 0:m-2] - z[2:n, 2:m]
                            + z[0:n-2, 0:m-2] - z[0:n-2, 2:m]) / csz
        # (2*(Zd - Zu) - Zul - Zur + Zdl + Zdr) / csz
        Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1]) - z[2:n, 0:m-2] - z[2:n, 2:m]
                            + z[0:n-2, 0:m-2] + z[0:n-2, 2:m]) / csz
        Nz = 8 * Nz
        # filling the first col of the matrix
        # (- 2*(Zr - Zp) + Zu - Zur + Zd - Zdr) / csz
        Nx[1:n-1, 0] = (- 2*(z[1:n-1, 1] - z[1:n-1, 0]) + z[2:n, 0] - z[2:n, 1] + z[0:n-2, 0] - z[0:n-2, 1]) / csz
        # (Zd - Zu + Zdr - Zur) / csz
        Ny[1:n-1, 0] = (z[0:n-2, 0] - z[2:n, 0] + z[0:n-2, 1] - z[2:n, 1]) / csz
        Nz[1:n-1, 0] = 4
        # filling the last col of the matrix
        # (2*(Zl - Zp) + Zdl - Zd + Zul - Zu) / csz
        Nx[1:n-1, m-1] = (2*(z[1:n-1, m-2] - z[1:n-1, m-1]) + z[0:n-2, m-2]
                          - z[0:n-2, m-1] + z[2:n, m-2] - z[2:n, m-1]) / csz
        # (Zd - Zu + Zdl - Zul) / csz
        Ny[1:n-1, m-1] = (z[0:n-2, m-1] - z[2:n, m-1] + z[0:n-2, m-2] - z[2:n, m-2]) / csz
        Nz[1:n-1, m-1] = 4
        # filling the first row of the matrix
        # (Zl - Zr + Zul - Zur) / csz
        Nx[0, 1:m-1] = (z[0, 0:m-2] - z[0, 2:m] + z[1, 0:m-2] - z[1, 2:m]) / csz
        # (-2*(Zu - Zp) + Zr - Zur + Zl - Zul) / csz
        Ny[0, 1:m-1] = (- 2*(z[1, 1:m-1] - z[0, 1:m-1]) + z[0, 2:m] - z[1, 2:m] + z[0, 0:m-2] - z[1, 0:m-2]) / csz
        Nz[0, 1:m-1] = 4
        # filling the last row of the matrix
        # (Zl - Zr + Zdl - Zdr) / csz
        Nx[n-1, 1:m-1] = (z[n-1, 0:m-2] - z[n-1, 2:m] + z[n-2, 0:m-2] - z[n-2, 2:m]) / csz
        # (2*(Zd - Zp) + Zdl - Zl + Zdr - Zr) / csz
        Ny[n-1, 1:m-1] = (2*(z[n-2, 1:m-1] - z[n-1, 1:m-1]) + z[n-2, 0:m-2] - z[n-1, 0:m-2] + z[n-2, 2:m] - z[n-1, 2:m]) / csz
        Nz[n-1, 1:m-1] = 4
        # filling the corners of the matrix
        Nx[0, 0] = (z[1, 0] - z[1, 1] - (z[0, 1] - z[0, 0])) / csz
        Ny[0, 0] = (z[0, 1] - z[1, 1] - (z[1, 0] - z[0, 0])) / csz
        Nz[0, 0] = 2
        Nx[n-1, 0] = (-(z[n-1, 1] - z[n-1, 0]) + z[n-2, 0] - z[n-2, 1]) / csz
        Ny[n-1, 0] = (z[n-2, 1] - z[n-1, 1] + z[n-2, 0] - z[n-1, 0]) / csz
        Nz[n-1, 0] = 2
        Nx[0, m-1] = (z[1, m-2] - z[1, m-1] + z[0, m-2] - z[0, m-1]) / csz
        Ny[0, m-1] = (-(z[1, m-1] - z[0, m-1]) + z[0, m-2] - z[1, m-2]) / csz
        Nz[0, m-1] = 2
        Nx[n-1, m-1] = (z[n-1, m-2] - z[n-1, m-1]
                        + z[n-2, m-2] - z[n-2, m-1]) / csz
        Ny[n-1, m-1] = (z[n-2, m-1] - z[n-1, m-1]
                        + z[n-2, m-2] - z[n-1, m-2]) / csz
        Nz[n-1, m-1] = 2

    if num == 1:
        # using the simple cross product
        z1 = np.append(z, z[:, -2].reshape(n, 1), axis=1)
        n1, m1 = np.shape(z1)
        z2 = np.append(z1, z1[-2, :].reshape(1, m1), axis=0)
        n2, m2 = np.shape(z2)

        Nx = - ((z2[0:n2-1, 1:m2] - z2[1:n2, 0:m2-1]) + (z2[1:n2, 1:m2] - z2[0:n2-1, 0:m2-1])) * csz
        Ny = - ((z2[1:n2, 1:m2] - z2[0:n2-1, 0:m2-1]) - (z2[0:n2-1, 1:m2] - z2[1:n2, 0:m2-1])) * csz
        Nz = 2 * Nz * csz * csz

        # Nx = - (z2[0:n2-1, 1:m2] - z2[0:n2-1, 0:m2-1]) / csz
        # Ny = - (z2[1:n2, 0:m2-1] - z2[0:n2-1, 0:m2-1]) / csz
        Ny[n-1, 0:m-1] = -Ny[n-1, 0:m-1]
        Nx[0:n-1, m-1] = -Nx[0:n-1, m-1]
        Ny[n-1, m-1] = -Ny[n-1, m-1]
        Nx[n-1, m-1] = -Nx[n-1, m-1]
        # TODO, Try to replicate samosAT notmal computation
        # if method num=1 is used, the normals are computed at com1DFA (original) cell center
        # this corresponds to our cell vertex
        # Create com1DFA (original) vertex grid
        x = np.linspace(-csz/2., (ncols-1)*csz - csz/2., ncols)
        y = np.linspace(-csz/2., (nrows-1)*csz - csz/2., nrows)
        X, Y = np.meshgrid(x, y)
        # interpolate the normal from com1DFA (original) center to his vertex
        # this means from our vertex to our centers
        Nx, Ny, NzCenter = getNormalArray(X, Y, Nx, Ny, Nz, csz)
        # this is for tracking mesh cell with actual data
        NzCenter = np.where(np.isnan(Nx), Nz, NzCenter)
        Nz = NzCenter

    return 0.5*Nx, 0.5*Ny, 0.5*Nz


def getAreaMesh(Nx, Ny, Nz, csz, num):
    """ Get area of grid cells.

        Parameters
        ----------
            Nx: 2D numpy array
                x component of the normal vector field on grid points
            Ny: 2D numpy array
                y component of the normal vector field on grid points
            Nz: 2D numpy array
                z component of the normal vector field on grid points
            csz: float
                cellsize of the grid

        Returns
        -------
            A: 2D numpy array
                Area of grid cells
    """
    # see documentation and issue 202
    if num == 1:
        A = norm(Nx, Ny, Nz)
    else:
        _, _, NzNormed = normalize(Nx, Ny, Nz)
        A = csz * csz / NzNormed
    return A


def removePart(particles, mask, nRemove, reasonString=''):
    """ remove given particles

    Parameters
    ----------
    particles : dict
        particles dictionary
    mask : 1D numpy array
        particles to keep
    nRemove : int
        number of particles removed
    reasonString: str
        reason why removing particles - for log message

    Returns
    -------
    particles : dict
        particles dictionary
    """

    log.info('removed %s particles %s' % (nRemove, reasonString))
    nPart = particles['Npart']
    for key in particles:
        if key == 'Npart':
            particles['Npart'] = particles['Npart'] - nRemove
        elif type(particles[key]).__module__ == np.__name__:
            if np.size(particles[key]) == nPart:
                particles[key] = particles[key][mask]

    particles['mTot'] = np.sum(particles['m'])
    return particles


def splitPart(particles):
    """Split big particles

    Split particles bigger than 1.5 times the massPerPart

    Parameters
    ----------
    particles : dict
        particles dictionary

    Returns
    -------
    particles : dict
        particles dictionary

    """
    massPerPart = particles['massPerPart']
    m = particles['m']
    nSplit = np.round(m/massPerPart)
    Ind = np.where(nSplit > 1)[0]
    if np.size(Ind) > 0:
        # loop on particles to split
        for ind in Ind:
            # get old values
            nPart = particles['Npart']
            nID = particles['nID']
            # compute new mass and particles to add
            mNew = m[ind] / nSplit[ind]
            nAdd = (nSplit[ind]-1).astype('int')
            # update total number of particles and number of IDs used so far
            particles['Npart'] = particles['Npart'] + nAdd
            particles['nID'] = nID + nAdd
            log.debug('Spliting particle %s in %s' % (ind, nAdd+1))
            for key in particles:
                # update splitted particle mass
                particles['m'][ind] = mNew
                # add new particles at the end of the arrays
                if type(particles[key]).__module__ == np.__name__:
                    # create unique ID for the new particles
                    if key == 'ID':
                        particles['ID'] = np.append(particles['ID'], np.arange(nID, nID + nAdd, 1))
                    # set the parent properties to new particles due to splitting
                    elif np.size(particles[key]) == nPart:
                        particles[key] = np.append(particles[key], particles[key][ind]*np.ones((nAdd)))

    particles['mTot'] = np.sum(particles['m'])
    return particles


def mergeParticleDict(particles1, particles2):
    """Merge two particles dictionary

    Parameters
    ----------
    particles1 : dict
        first particles dictionary
    particles2 : dict
        second particles dictionary

    Returns
    -------
    particles : dict
        merged particles dictionary

    """
    particles = {}
    nPart1 = particles1['Npart']
    # loop on the keys from particles1 dicionary
    for key in particles1:
        # deal with specific cases
        # Npart: just sum them up
        if key == 'Npart':
            particles['Npart'] = particles1['Npart'] + particles2['Npart']
        # massPerPart, should stay unchanged. If ever they are different take
        # the minimum
        # ToDo: are we sure we want the minimum?
        elif key == 'massPerPart':
            particles['massPerPart'] = min(particles1['massPerPart'], particles2['massPerPart'])
        # now if the value is a numpy array and this key is also in particles2
        elif (key in particles2) and (type(particles1[key]).__module__ == np.__name__):
            # deal with the specific cases:
            # in the case of ID or 'parentID' we assume that bot particles1 and
            # particles2 were initialized with an ID and parentID starting at 0
            # here whene we merge the 2 arrays we make sure to shift the value
            # of particles2 so that the ID stays a unique identifier and
            # that the parentID is consistent with this shift.
            if (key == 'ID') or (key == 'parentID'):
                particles[key] = np.append(particles1[key], particles2[key] + particles1['nID'])
            # general case where the key value is an array with as many elements
            # as particles
            elif np.size(particles1[key]) == nPart1:
                particles[key] = np.append(particles1[key], particles2[key])
            # if the array is of size one, (potential energy, mTot...) we just
            # sum the 2 values
            else:
                particles[key] = particles1[key] + particles2[key]
        # the key is in both dictionaries, it is not an array but it is a
        # number (int, double, float) then we sum the 2 values
        elif (key in particles2) and (isinstance(particles1[key], numbers.Number)):
            particles[key] = particles1[key] + particles2[key]
        # finaly, if the key is only in particles1 then we give this value to
        # the new particles
        else:
            particles[key] = particles1[key]
    return particles


def findParticles2Track(particles, center, radius):
    '''Find particles within a circle arround a given point

    Parameters
    ----------
    particles : dict
        particles dictionary (with the 'parentID' array)
    center : dict
        point dictionary:
            x : x coordinate
            y : y coordinate
            z : z coordinate
    radius : float
        radius of the circle around point

    Returns
    -------
    particles2Track : numpy array
        array with Parent ID of particles to track
    track: boolean
        False if no particles are tracked
    '''
    track = True
    x = particles['x']
    y = particles['y']
    z = particles['z']
    xc = center['x']
    yc = center['y']
    zc = center['z']
    r = norm(x-xc, y-yc, z-zc)
    index = np.where(r <= radius)
    particles2Track = particles['parentID'][index]
    log.info('Tracking %d particles' % len(index[0]))
    if len(index[0]) < 1:
        log.warning('Found No particles to track ')
        track = False

    return particles2Track, track


def getTrackedParticles(particlesList, particles2Track):
    '''Track particles along time given the parentID of the particles to track

    Parameters
    ----------
    particlesList : list
        list of particles dictionaries (with the 'parentID' array)
    particles2Track : numpy array
        array with the parentID of the particles to track

    Returns
    -------
    particlesList : list
        list of particles dictionaries updated with the 'trackedParticles'
        array (in the array, the ones correspond to the particles that
        are tracked)
    nPartTracked : int
        total number of tracked particles
    '''
    nPartTracked = np.size(particles2Track)
    # add trackedParticles array to the particles dictionary for every saved time step
    for particles in particlesList:
        # find index of particles to track
        index = [ind for ind, parent in enumerate(particles['parentID']) if parent in particles2Track]
        nPartTracked = max(nPartTracked, len(index))
        trackedParticles = np.zeros(particles['Npart'])
        trackedParticles[index] = 1
        particles['trackedParticles'] = trackedParticles
    return particlesList, nPartTracked


def getTrackedParticlesProperties(particlesList, nPartTracked, properties):
    '''Get the desired properties for the tracked particles

    Parameters
    ----------
    particlesList : list
        list of particles dictionaries (with the 'parentID' array)
    nPartTracked : int
        total number of tracked particles
    properties : list
        list of strings

    Returns
    -------
    trackedPartProp : dict
        dictionary with 2D numpy arrays corresponding to the time series of the
        properties for the tracked particles (for example if
        properties = ['x', 'y'], the dictionary will have the keys 'time',
        'x' and 'y'. trackedPartProp['x'] will be a 2D numpy array, each line
        corresponds to the 'x' time series of a tracked particle)
    '''
    # buid time series for desired properties of tracked particles
    nTimeSteps = len(particlesList)
    trackedPartProp = {}
    trackedPartProp['time'] = np.zeros(nTimeSteps)
    newProperties = []
    # initialize
    for key in properties:
        if key in particlesList[0]:
            trackedPartProp[key] = np.zeros((nTimeSteps, nPartTracked))
            newProperties.append(key)
        else:
            log.warning('%s is not a particle property' % key)

    # extract wanted properties and build the time series
    trackedPartID = []
    for particles, nTime in zip(particlesList, range(nTimeSteps)):
        trackedParticles = particles['trackedParticles']
        trackedPartProp['time'][nTime] = particles['t']
        index = np.where(trackedParticles == 1)
        for ind, id in zip(index[0], particles['ID'][index]):
            if id not in trackedPartID:
                trackedPartID.append(id)
            indCol = trackedPartID.index(id)
            for key in newProperties:
                trackedPartProp[key][nTime, indCol] = particles[key][ind]

    return trackedPartProp

##############################################################################
# ###################### Vectorial functions #################################
##############################################################################


def norm(x, y, z):
    """ Compute the Euclidean norm of the vector (x, y, z).

    (x, y, z) can be numpy arrays.

    Parameters
    ----------
        x: numpy array
            x component of the vector
        y: numpy array
            y component of the vector
        z: numpy array
            z component of the vector

    Returns
    -------
        norme: numpy array
            norm of the vector
    """
    norme = np.sqrt(x*x + y*y + z*z)
    return norme


def norm2(x, y, z):
    """ Compute the square of the Euclidean norm of the vector (x, y, z).

    (x, y, z) can be numpy arrays.

    Parameters
    ----------
        x: numpy array
            x component of the vector
        y: numpy array
            y component of the vector
        z: numpy array
            z component of the vector

    Returns
    -------
        norme2: numpy array
            square of the norm of the vector
    """
    norme2 = (x*x + y*y + z*z)
    return norme2


def normalize(x, y, z):
    """ Normalize vector (x, y, z) for the Euclidean norm.

    (x, y, z) can be np arrays.

    Parameters
    ----------
        x: numpy array
            x component of the vector
        y: numpy array
            y component of the vector
        z: numpy array
            z component of the vector

    Returns
    -------
        x: numpy array
            x component of the normalized vector
        y: numpy array
            y component of the normalized vector
        z: numpy array
            z component of the normalized vector
    """
    norme = norm(x, y, z)
    ind = np.where(norme > 0)
    x[ind] = x[ind] / norme[ind]
    y[ind] = y[ind] / norme[ind]
    z[ind] = z[ind] / norme[ind]

    return x, y, z


def crossProd(ux, uy, uz, vx, vy, vz):
    """ Compute cross product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
    """
    wx = uy*vz - uz*vy
    wy = uz*vx - ux*vz
    wz = ux*vy - uy*vx

    return wx, wy, wz


def scalProd(ux, uy, uz, vx, vy, vz):
    """ Compute scalar product of vector u = (ux, uy, uz) and v = (vx, vy, vz).
    """
    scal = ux*vx + uy*vy + uz*vz

    return scal
