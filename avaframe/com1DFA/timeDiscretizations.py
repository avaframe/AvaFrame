"""
    Functions regarding time discretization and time stepping for com1DFA
"""

# Load modules
import logging
import numpy as np

# Local imports
import avaframe.com1DFA.DFAtools as DFAtls


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getcflTimeStep(particles, dem, cfg):
    """ Compute cfl time step  """

    # determine max velocity of particles
    vmagnitude = DFAtls.norm(particles['ux'], particles['uy'], particles['uz'])
    vMax = np.amax(vmagnitude)

    # get cell size
    cszDEM = dem['header']['cellsize']
    cszNeighbourGrid = dem['headerNeighbourGrid']['cellsize']
    # use the smallest of those two values
    csz = min(cszDEM, cszNeighbourGrid)

    # courant number
    cMax = float(cfg['cMax'])

    # compute stable time step
    # if velocity is zero - divided by zero error so to avoid:
    if vMax <= (cMax * csz)/float(cfg['maxdT']):
        dtStable = float(cfg['maxdT'])
    else:
        dtStable = (cMax * csz) / vMax
        if cfg.getboolean('constrainCFL'):
            if dtStable < float(cfg['mindT']):
                dtStable = float(cfg['mindT'])

    log.debug('dtStable is with cMAX=%.1f is: %.4f with vMax:%.2f' % (cMax, dtStable, vMax))

    # return stable time step
    return dtStable


def getSphKernelRadiusTimeStep(dem, cfg):
    """ Compute the time step  given the sph kernel radius and the cMax coefficient
    This is based on the article from Ben Moussa et Vila
    DOI:10.1137/S0036142996307119
    Parameters
    -----------
    dem: dict
        dem dictionary (with info about sph kernel radius and mesh size)
    cfg: configparser
        the cfg cith cMax
    Returns
    --------
    dtStable: float
        corresponding time step
    """
    # get cell size
    cszDEM = dem['header']['cellsize']
    cszNeighbourGrid = dem['headerNeighbourGrid']['cellsize']
    # use the minimum of those two values
    csz = min(cszDEM, cszNeighbourGrid)

    # courant number
    cMax = float(cfg['cMax'])

    dtStable = cMax * csz
    log.debug('dtStable is with cMAX=%.1f is: %.4f' % (cMax, dtStable))

    # return stable time step
    return dtStable
