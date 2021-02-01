"""
    function regarding time discretization and time stepping for com1DFA
    This file is part of Avaframe.
"""

# Load modules
import os
import glob
import logging
import numpy as np

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFAPy.DFAtools as DFAtls


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getcflTimeStep(particles, dem, cfgDict):
    """ Compute cfl time step  """

    # determine max velocity of particles
    vmagnitude = DFAtls.norm(particles['ux'], particles['uy'], particles['uz'])
    vmax = np.amax(vmagnitude)

    # get cell size
    csz = dem['header'].cellsize
    # use the smoothing length of the kernel instead, for now
    # rKernel = csz

    # courant number
    cMax = float(cfgDict['cMax'])

    # compute stable time step
    # if velocity is zero - divided by zero error so to avoid:
    if vmax == 0.0:
        dtStable = float(cfgDict['mindT'])
    else:
        dtStable = (cMax * csz) / vmax

    # 'overwrite' dt that is read from cfg ini file
    cfgDict['dt'] = str(dtStable)
    log.info('dtStable is with cMAX=%.1f is: %.4f with vmax:%.2f' % (cMax, dtStable, vmax))

    # return stable time step
    return dtStable


def getcfldTwithConstraints(particles, dem, cfgDict):
    """ Compute cfl time step  """

    # determine max velocity of particles
    vmagnitude = DFAtls.norm(particles['ux'], particles['uy'], particles['uz'])
    vmax = np.amax(vmagnitude)

    # get cell size
    csz = dem['header'].cellsize
    # use the smoothing length of the kernel instead, for now
    # rKernel = csz

    # courant number
    cMax = float(cfgDict['cMax'])

    # compute stable time step
    # if velocity is zero - divided by zero error so to avoid:
    if vmax == 0.0:
        dtStable = float(cfgDict['mindT'])
    else:
        dtStable = (cMax * csz) / vmax
        if dtStable < float(cfgDict['mindT']):
            dtStable = float(cfgDict['mindT'])
        elif dtStable > float(cfgDict['maxdT']):
            dtStable = float(cfgDict['maxdT'])

    # 'overwrite' dt that is read from cfg ini file
    cfgDict['dt'] = str(dtStable)
    log.debug('dtStable is with cMAX=%.1f is: %.4f with vmax:%.2f' % (cMax, dtStable, vmax))

    # return stable time step
    return dtStable
