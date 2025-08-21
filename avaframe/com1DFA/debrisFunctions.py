"""
Functions that are used specifically for modeling debris flows (within DebrisFrame).
"""

# Load modules
import logging
import numpy as np
import copy

import avaframe.com1DFA.DFAfunctionsCython as DFAfunC
import avaframe.com1DFA.particleTools as particleTools
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.com1DFA.com1DFA as com1DFA

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def updateParticlesHydrograph(cfg, inputSimLines, particles, fields, dem, zPartArray0, t):
    """
    Update particles with "new" particles initialised by a hydrograph.

    Parameters
    ---------
    cfg: dict
        configuration settings
    inputSimLines : dict
        dictionary with input data dictionaries (releaseLine, hydrographLine,...)
    particles : dict
        particles dictionary at t that are in the flow already
    fields: dict
        fields dictionary at t
    dem: dict
        dictionary with info on DEM data
    zPartArray0: dict
        dictionary containing z - value of particles at timestep 0
    t: float
        timestep of iteration

    Returns
    ---------
    particles: dict
        particles dictionary at t including the hydrograph particles
    fields: dict
        updated fields dictionary at t including the hydrograph particles
    zPartArray0: dict
        dictionary containing z - value of particles at timestep 0
    """
    hydrValues = inputSimLines["hydrographLine"]["values"]
    if round(t, 1) in hydrValues["timeStep"]:
        i = np.where(hydrValues["timeStep"] == round(t, 1))
        log.info(
            "add hydrograph at timestep: %f with thickness %s and velocity %s"
            % (t, hydrValues["thickness"][i], hydrValues["velocity"][i])
        )
        # similar workflow to secondary release!
        particles, zPartArray0 = addHydrographParticles(
            cfg,
            particles,
            inputSimLines,
            hydrValues["thickness"][i],
            hydrValues["velocity"][i],
            dem,
            zPartArray0,
        )
        particles = DFAfunC.getNeighborsC(particles, dem)
        # update fields (compute grid values)
        if fields["computeTA"]:
            particles = DFAfunC.computeTrajectoryAngleC(particles, zPartArray0)
        particles, fields = DFAfunC.updateFieldsC(cfg["GENERAL"], particles, dem, fields)

    return particles, fields, zPartArray0


def addHydrographParticles(cfg, particles, inputSimLines, thickness, velocityMag, dem, zPartArray0):
    """
    add new particles initialized by a hydrograph to particles that are in the flow already

    Parameters
    ---------
    cfg: dict
        configuration settings
    particles : dict
        particles dictionary at t that are in the flow already
    inputSimLines : dict
        dictionary with input data dictionaries (releaseLine, hydrographLine,...)
    thickness: float
        thickness of incoming hydrograph
    velocityMag: float
        velocity of incoming hydrograph
    dem: dict
        dictionary with info on DEM data
    zPartArray0: dict
        dictionary containing z - value of particles at timestep 0

    Returns
    ---------
    particles: dict
        particles dictionary at t including the hydrograph particles
    zPartArray0: dict
        dictionary containing z - value of particles at timestep 0
    """
    hydrLine = inputSimLines["hydrographLine"]
    hydrLine["header"] = dem["originalHeader"].copy()
    hydrLine = geoTrans.prepareArea(
        hydrLine,
        dem,
        np.sqrt(2),
        thList=[thickness],
        combine=True,
        checkOverlap=False,
    )
    particlesHydrograph = com1DFA.initializeParticles(
        cfg["GENERAL"],
        hydrLine,
        dem,
    )
    particlesHydrograph = DFAfunC.updateInitialVelocity(
        cfg["GENERAL"], particlesHydrograph, dem, velocityMag
    )
    particles = particleTools.mergeParticleDict(particles, particlesHydrograph)
    # save initial z position for travel angle computation
    zPartArray0 = np.append(zPartArray0, copy.deepcopy(particlesHydrograph["z"]))
    return particles, zPartArray0


def checkHydrograph(cfgGen, hydrographValues, hydrCsv):
    """
    check if hydrograph satisfied some requirements
    Parameters
    -----------
    hydrCsv: str
        directory to csv table containing hydrograph values
    cfgGen: configparser
        configuration settings, part "GENERAL"
    hydrographValues: dict
        contains hydrograph values: timestep, thickness, velocity
    """
    # check if timesteps are unique
    timeStepUnique = np.unique(hydrographValues["timeStep"])
    if timeStepUnique.ndim == 0:
        if timeStepUnique != hydrographValues["timeStep"]:
            message = "The provided hydrograph time steps in %s are not unique" % (hydrCsv)
    elif len(timeStepUnique) != len(hydrographValues["timeStep"]):
        message = "The provided hydrograph timesteps in %s are not unique" % (hydrCsv)
        log.error(message)
        raise ValueError(message)

    if cfgGen.getboolean("hydrograph") and cfgGen.getboolean("noRelArea"):
        if 0 not in timeStepUnique:
            message = (
                "If no release area is released, a thickness needs to be provided for  time step 0 s in %s"
                % (hydrCsv)
            )
            log.error(message)
            raise ValueError(message)
    for th in hydrographValues["thickness"]:
        if th <= 0:
            message = "For every release time step a thickness > 0 needs to be provided in %s" % (hydrCsv)
            log.error(message)
            raise ValueError(message)
