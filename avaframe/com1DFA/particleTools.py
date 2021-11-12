"""
    Tools for handling particles, splitting, merging and tracking.
"""

# Load modules
import logging
import numpy as np
import numbers
import math

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


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
    if reasonString != '':
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


def addParticles(particles, nAdd, ind, mNew, xNew, yNew, zNew):
    """ add particles

    Parameters
    ----------
    particles : dict
        particles dictionary
    nAdd : int
        number of particles added (one particles is modified, nAdd are added)
    ind : int
        index of particle modified
    mNew: float
        new mass of the particles
    xNew: numpy array
        new x position of the particles
    yNew: numpy array
        new y position of the particles
    zNew: numpy array
        new z position of the particles

    Returns
    -------
    particles : dict
        particles dictionary with modified particle and new ones
    """
    # get old values
    nPart = particles['Npart']
    nID = particles['nID']
    # update total number of particles and number of IDs used so far
    particles['Npart'] = particles['Npart'] + nAdd
    particles['nID'] = nID + nAdd
    # log.info('Spliting particle %s in %s' % (ind, nSplit))
    for key in particles:
        # update splitted particle mass
        # first update the old particle
        particles['m'][ind] = mNew
        particles['x'][ind] = xNew[0]
        particles['y'][ind] = yNew[0]
        particles['z'][ind] = zNew[0]
        # add new particles at the end of the arrays
        if type(particles[key]).__module__ == np.__name__:
            # create unique ID for the new particles
            if key == 'ID':
                particles['ID'] = np.append(particles['ID'], np.arange(nID, nID + nAdd, 1))
            elif key == 'x':
                particles[key] = np.append(particles[key], xNew[1:])
            elif key == 'y':
                particles[key] = np.append(particles[key], yNew[1:])
            elif key == 'z':
                particles[key] = np.append(particles[key], zNew[1:])
            # set the parent properties to new particles due to splitting
            elif np.size(particles[key]) == nPart:
                particles[key] = np.append(particles[key], particles[key][ind]*np.ones((nAdd)))
            # ToDo: maybe also update the h smartly
    return particles


def splitPart(particles, cfg):
    """Split big particles

    Split particles bigger than thresholdMassSplit x massPerPart
    place the new particle in the flow direction at distance epsilon x rPart
    (this means splitting happens only if particles grow -> entrainment)

    Parameters
    ----------
    particles : dict
        particles dictionary
    cfg : configParser
        GENERAL configuration for com1DFA

    Returns
    -------
    particles : dict
        particles dictionary

    """
    rho = cfg.getfloat('rho')
    thresholdMassSplit = cfg.getfloat('thresholdMassSplit')
    distSplitPart = cfg.getfloat('distSplitPart')
    massPerPart = particles['massPerPart']
    mPart = particles['m']
    # decide which particles to split
    nSplit = np.ceil(mPart/(massPerPart*thresholdMassSplit))
    Ind = np.where(nSplit > 1)[0]
    # loop on particles to split
    for ind in Ind:
        # compute new mass (split particle in 2)
        mNew = mPart[ind] / 2  # nSplit[ind]
        nAdd = 1  # (nSplit[ind]-1).astype('int')
        xNew, yNew, zNew = getSplitPartPositionSimple(particles, rho, distSplitPart, ind)
        # add new particles
        particles = addParticles(particles, nAdd, ind, mNew, xNew, yNew, zNew)

    particles['mTot'] = np.sum(particles['m'])
    return particles


def testSplitPart(particles, cfg, dem):
    """Split big particles

    Split particles to keep enough particles within the kernel radius.
    place the new particle in the flow direction at distance epsilon x rPart

    Parameters
    ----------
    particles : dict
        particles dictionary
    cfg : configParser
        GENERAL configuration for com1DFA
    dem : dict
        dem dictionary

    Returns
    -------
    particles : dict
        particles dictionary

    """
    # get cfg info
    rho = cfg.getfloat('rho')
    sphKernelRadius = cfg.getfloat('sphKernelRadius')
    cMinNPPK = cfg.getfloat('cMinNPPK')
    cMinMass = cfg.getfloat('cMinMass')
    nSplit = cfg.getint('nSplit')
    # get dem info
    csz = dem['header']['cellsize']
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']
    # get the threshold area over which we split the particle
    massPerPart = particles['massPerPart']
    nPPK = particles['nPPK']
    aMax = math.pi * sphKernelRadius**2 / (cMinNPPK * nPPK)
    mMin = massPerPart * cMinMass
    # get particle area
    mPart = particles['m']
    hPart = particles['h']
    aPart = mPart/(rho*hPart)
    # find particles to split
    tooBig = np.where((aPart > aMax) & (mPart/nSplit > mMin))[0]
    # count new particles
    nNewPart = 0
    # loop on particles to split
    for ind in tooBig:
        # compute new mass and particles to add
        mNew = mPart[ind] / nSplit
        nAdd = nSplit-1
        nNewPart = nNewPart + nAdd
        # get position of new particles
        xNew, yNew, zNew = getSplitPartPosition(cfg, particles, aPart, Nx, Ny, Nz, csz, nSplit, ind)
        # add new particles
        particles = addParticles(particles, nAdd, ind, mNew, xNew, yNew, zNew)
    log.debug('Added %s because of splitting' % (nNewPart))

    particles['mTot'] = np.sum(particles['m'])
    return particles


def testMergePart(particles, cfg, dem):
    """merge small particles

    merge particles to avoid too many particles within the kernel radius.
    place the new merge particle between the two old ones. The new position and velocity are the
    mass averaged ones

    Parameters
    ----------
    particles : dict
        particles dictionary
    cfg : configParser
        GENERAL configuration for com1DFA
    dem : dict
        dem dictionary

    Returns
    -------
    particles : dict
        particles dictionary

    """
    # get cfg info
    rho = cfg.getfloat('rho')
    sphKernelRadius = cfg.getfloat('sphKernelRadius')
    cMaxNPPK = cfg.getfloat('cMaxNPPK')
    # get the threshold area under which we merge the particle
    nPPK = particles['nPPK']
    aMin = math.pi * sphKernelRadius**2 / (cMaxNPPK * nPPK)
    # get particle area
    mPart = particles['m']
    hPart = particles['h']
    xPart = particles['x']
    yPart = particles['y']
    zPart = particles['z']
    aPart = mPart/(rho*hPart)
    # find particles to merge
    tooSmall = np.where(aPart < aMin)[0]
    keepParticle = np.ones((particles['Npart']), dtype=bool)
    nRemoved = 0
    # loop on particles to merge
    for ind in tooSmall:
        if keepParticle[ind]:
            # find nearest particle
            rMerge, neighbourInd = getClosestNeighbour(xPart, yPart, zPart, ind, sphKernelRadius, keepParticle)
            # only merge a particle if it is closer thant the kernel radius
            if rMerge < sphKernelRadius:
                # remove neighbourInd from tooSmall if possible
                keepParticle[neighbourInd] = False
                nRemoved = nRemoved + 1
                # compute new mass and particles to add
                mNew = mPart[ind] + mPart[neighbourInd]
                # compute mass averaged values
                for key in ['x', 'y', 'z', 'ux', 'uy', 'uz']:
                    particles[key][ind] = (mPart[ind]*particles[key][ind] +
                                           mPart[neighbourInd]*particles[key][neighbourInd]) / mNew
                particles['m'][ind] = mNew
                # ToDo: mabe also update h

    particles = removePart(particles, keepParticle, nRemoved, reasonString='')  # 'because of colocation')
    return particles


def getClosestNeighbour(xPartArray, yPartArray, zPartArray, ind, sphKernelRadius, keepParticle):
    """ find closest neighbour

    Parameters
    ----------
    xPartArray: numpy array
        x position of the particles
    yPartArray: numpy array
        y position of the particles
    zPartArray: numpy array
        z position of the particles
    ind : int
        index of particle modified
    sphKernelRadius: float
        kernel radius
    keepParticle: numpy array
        boolean array telling if particles are kept or merged

    Returns
    -------
    rMerge : float
        distance to the closest neighbour
    neighbourInd : int
        index of closest neighbour
    """
    r = DFAtls.norm(xPartArray-xPartArray[ind], yPartArray-yPartArray[ind], zPartArray-zPartArray[ind])
    # make sure you don't find the particle itself
    r[ind] = 2*sphKernelRadius
    # make sure you don't find a particle that has already been merged
    r[~keepParticle] = 2*sphKernelRadius
    # find nearest particle
    neighbourInd = np.argmin(r)
    rMerge = r[neighbourInd]
    return rMerge, neighbourInd


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


def getSplitPartPosition(cfg, particles, aPart, Nx, Ny, Nz, csz, nSplit, ind):
    """Compute the new particle potion due to splitting

    Parameters
    ----------
    cfg : configParser
        GENERAL configuration for com1DFA
    particles : dict
        particles dictionary
    aPart : numpy array
        particle area array
    Nx : numpy 2D array
        x component of the normal vector on the grid
    Ny : numpy 2D array
        y component of the normal vector on the grid
    Nz : numpy 2D array
        z component of the normal vector on the grid
    csz : float
        grid cell size
    nSplit : int
        in how many particles do we split?
    ind : int
        index of the particle to split

    Returns
    -------
    xNew : numpy array
        x components of the splitted particles
    yNew : numpy array
        y components of the splitted particles
    zNew : numpy array
        z components of the splitted particles
    """
    rng = np.random.default_rng(int(cfg['seed']))
    x = particles['x']
    y = particles['y']
    z = particles['z']
    ux = particles['ux']
    uy = particles['uy']
    uz = particles['uz']
    rNew = np.sqrt(aPart[ind] / (math.pi * nSplit))
    alpha = 2*math.pi*(np.arange(nSplit)/nSplit + rng.random(1))
    cos = rNew*np.cos(alpha)
    sin = rNew*np.sin(alpha)
    # nx, ny, nz = getNormalArray(np.array([particles['x'][ind]]), np.array([particles['y'][ind]]), Nx, Ny, Nz, csz)
    nx, ny, nz = DFAtls.getNormalArray(np.array([x[ind]]), np.array([y[ind]]), Nx, Ny, Nz, csz)
    e1x, e1y, e1z, e2x, e2y, e2z = getTangenVectors(nx, ny, nz, np.array([ux[ind]]), np.array([uy[ind]]), np.array([uz[ind]]))
    xNew = x[ind] + cos * e1x + sin * e2x
    yNew = y[ind] + cos * e1y + sin * e2y
    zNew = z[ind] + cos * e1z + sin * e2z
    # toDo: do we need to reproject the particles on the dem?
    return xNew, yNew, zNew


def getSplitPartPositionSimple(particles, rho, distSplitPart, ind):
    """Compute the new particle potion due to splitting

    Parameters
    ----------
    particles : dict
        particles dictionary
    rho : float
        density
    distSplitPart : float
        distance coefficient
    ind : int
        index of the particle to split

    Returns
    -------
    xNew : numpy array
        x components of the splitted particles
    yNew : numpy array
        y components of the splitted particles
    zNew : numpy array
        z components of the splitted particles
    """
    mPart = particles['m'][ind]
    hPart = particles['h'][ind]
    xPart = particles['x'][ind]
    yPart = particles['y'][ind]
    zPart = particles['z'][ind]
    uxPart = particles['ux'][ind]
    uyPart = particles['uy'][ind]
    uzPart = particles['uz'][ind]
    # get the area of the particle as well as the distance expected between the old and new particle
    # note that if we did not update the particles FD, we use here the h from the previous time step
    aPart = mPart/(rho*hPart)
    rNew = distSplitPart * np.sqrt(aPart/math.pi)
    cos = rNew*np.array([-1, 1])
    # compute velocity mag to get the direction of the flow (e_1)
    uMag = DFAtls.norm(uxPart, uyPart, uzPart)
    xNew = xPart + cos * uxPart/uMag
    yNew = yPart + cos * uyPart/uMag
    zNew = zPart + cos * uzPart/uMag
    # toDo: do we need to reproject the particles on the dem?
    return xNew, yNew, zNew


def getTangenVectors(nx, ny, nz, ux, uy, uz):
    """Compute the tangent vector to the surface

    If possible, e1 is in the velocity direction, if not possible,
    use the tangent vector in x direction for e1 (not that any other u vector could be provided,
    it does not need to be the velocity vector, it only needs to be in the tangent plane)
    Parameters
    ----------
    nx : float
        x component of the normal vector
    ny : float
        y component of the normal vector
    nz : float
        z component of the normal vector
    ux : float
        x component of the velocity vector
    uy : float
        y component of the velocity vector
    uz : float
        z component of the velocity vector

    Returns
    -------
    e1x : float
        x component of the first tangent vector
    e1y : float
        y component of the first tangent vector
    e1z : float
        z component of the first tangent vector
    e2x : float
        x component of the second tangent vector
    e2y : float
        y component of the second tangent vector
    e2z : float
        z component of the second tangent vector
    """
    # compute the velocity magnitude
    velMag = DFAtls.norm(ux, uy, uz)
    if velMag > 0:
        e1x = ux / velMag
        e1y = uy / velMag
        e1z = uz / velMag
    else:
        # if vector u is zero use the tangent vector in x direction for e1
        e1x = np.array([1])
        e1y = np.array([0])
        e1z = -nx/nz
        e1x, e1y, e1z = DFAtls.normalize(e1x, e1y, e1z)
    # compute the othe tengent vector
    e2x, e2y, e2z = DFAtls.crossProd(nx, ny, nz, e1x, e1y, e1z)
    e2x, e2y, e2z = DFAtls.normalize(e2x, e2y, e2z)
    return e1x, e1y, e1z, e2x, e2y, e2z


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
    r = DFAtls.norm(x-xc, y-yc, z-zc)
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
