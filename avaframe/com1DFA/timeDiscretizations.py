"""
    Functions regarding time discretization and time stepping for com1DFA
"""

# Load modules
import logging

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


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
