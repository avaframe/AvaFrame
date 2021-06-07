"""
    function regarding time discretization and time stepping for com1DFA
    This file is part of Avaframe.
"""

# Load modules
import logging
import numpy as np

# Local imports
import avaframe.com1DFAPy.DFAtools as DFAtls


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getcflTimeStep(particles, dem, cfg):
    """ Compute cfl time step  """

    # determine max velocity of particles
    vmagnitude = DFAtls.norm(particles['ux'], particles['uy'], particles['uz'])
    vmax = np.amax(vmagnitude)

    # get cell size
    cszDEM = dem['header'].cellsize
    cszNeighbourGrid = dem['headerNeighbourGrid'].cellsize
    # use the smallest of those two values
    csz = min(cszDEM, cszNeighbourGrid)

    # courant number
    cMax = float(cfg['cMax'])

    # compute stable time step
    # if velocity is zero - divided by zero error so to avoid:
    if vmax <= (cMax * csz)/float(cfg['maxdT']):
        dtStable = float(cfg['maxdT'])
    else:
        dtStable = (cMax * csz) / vmax
        if cfg.getboolean('constrainCFL'):
            if dtStable < float(cfg['mindT']):
                dtStable = float(cfg['mindT'])

    log.debug('dtStable is with cMAX=%.1f is: %.4f with vmax:%.2f' % (cMax, dtStable, vmax))

    # return stable time step
    return dtStable
