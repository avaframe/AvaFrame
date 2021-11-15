"""
    Tools specific to the com1DFA computational kernel
"""

# Load modules
import logging
import math

# Local imports


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getPartInitMethod(cfg, csz):
    rho = cfg.getfloat('rho')
    massPerParticleDeterminationMethod = cfg['massPerParticleDeterminationMethod']
    nPPK =  0
    # derive mass per particle to define number of particles per cell:
    if massPerParticleDeterminationMethod == 'MPPDIR':
        log.debug('Number of particles defined by: mass per particle %s' % cfg['massPerPart'])
    elif massPerParticleDeterminationMethod == 'MPPDH':
        deltaTh = cfg.getfloat('deltaTh')
        ds = min(csz, cfg.getfloat('sphKernelRadius'))
        massPerPart = rho * ds * ds * deltaTh
        log.debug('Number of particles defined by: release thickness per particle: %s' % cfg['deltaTh'])
        log.debug('mass per particle is %.2f' % massPerPart)
    elif massPerParticleDeterminationMethod == 'MPPKR':
        nPPK = cfg.getfloat('nPPK')
        relTh = cfg.getfloat('relTh')
        ds = min(csz, cfg.getfloat('sphKernelRadius'))
        massPerPart = rho * math.pi * ds * ds * relTh / nPPK

    return massPerPart, nPPK
