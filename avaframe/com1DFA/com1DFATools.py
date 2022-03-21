"""
    Tools specific to the com1DFA computational kernel
"""

# Load modules
import logging
import math

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getPartInitMethod(cfg, csz, relThForPart):
    """ Get particle initialization parameters

    Get the massPerPart and nPPK corresponding to the desired initialization method

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    csz: float
        cell size
    relThForPart: float
        relTh value

    Returns
    -------
    massPerPart : float
        mass per particle desired
    nPPK : int
        number of particles per kernel radius desired
    """
    rho = cfg.getfloat('rho')
    massPerParticleDeterminationMethod = cfg['massPerParticleDeterminationMethod']
    nPPK = 0
    # derive mass per particle to define number of particles per cell:
    if massPerParticleDeterminationMethod == 'MPPDIR':
        massPerPart = cfg.getfloat('massPerPart')
        log.debug('Number of particles defined by: mass per particle %s' % cfg['massPerPart'])
    elif massPerParticleDeterminationMethod == 'MPPDH':
        deltaTh = cfg.getfloat('deltaTh')
        ds = min(csz, cfg.getfloat('sphKernelRadius'))
        massPerPart = rho * ds * ds * deltaTh
        log.debug('Number of particles defined by: release thickness per particle: %s' % cfg['deltaTh'])
        log.debug('mass per particle is %.2f' % massPerPart)
    elif massPerParticleDeterminationMethod == 'MPPKR':
        sphKernelRadius = cfg.getfloat('sphKernelRadius')
        cszMin = min(csz, sphKernelRadius)
        nPPK0 = cfg.getfloat('nPPK0')
        sphKR0 = cfg.getfloat('sphKR0')
        aPPK = cfg.getfloat('aPPK')
        nPPK = round(nPPK0 * (cszMin/sphKR0)**aPPK)
        massPerPart = rho * math.pi * cszMin * cszMin * relThForPart / nPPK

    return massPerPart, nPPK
