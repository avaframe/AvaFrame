"""
    functions for initialising particle distribution
"""

import logging
import time
import pathlib
import numpy as np
import copy
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

# Local imports
import avaframe.com1DFA.DFAfunctionsCython as DFAfunC
import avaframe.out3Plot.outDebugPlots as debPlot
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.com1DFA.particleTools as particleTools

# create local logger
log = logging.getLogger(__name__)

def getIniPosition(cfg, particles, dem, fields, inputSimLines, relThField):
    """ Redistribute particles so that SPH force reduces with fixed particles as boundaries

        Parameters
        ------------
        cfg: configparser object
            configuration settings
        particles: dict
            dictionary with particle properties
        dem: dict
            dictionary with dem header and data
        fields: dict
            dictionary with fields of result types
        inputSimLines: dict
            info on input files
        relThField: numpy array
            release thickness field

        Returns
        --------
        particles: dict
            updated particles dict with new positions
        fields: fields
            updated fields dictionary
        """

    iterate = True
    particles['iterate'] = True
    particlesList = [particles.copy()]
    countIterations = 0
    # load relRaster from buffered release line
    relRaster = inputSimLines['releaseLineBuffer']['rasterData']

    # redistribute particles to reduce SPH force
    while iterate and countIterations < cfg['GENERAL'].getint('maxIterations'):

        # compute artificial viscosity effect on velocity
        particles, force = DFAfunC.computeIniMovement(cfg['GENERAL'], particles, dem, cfg['GENERAL'].getfloat('dtIni'),
            fields)

        # compute SPH force
        particles, force = DFAfunC.computeForceSPHC(cfg['GENERAL'], particles, force, dem,
            cfg['GENERAL'].getint('sphOptionIni'), gradient=0)

        # update position as a result of SPH force and artifical viscosity
        particles = DFAfunC.updatePositionC(cfg['GENERAL'], particles, dem, force, cfg['GENERAL'].getfloat('dtIni'),
            typeStop=1)

        particlesList.append(particles.copy())

        # compute neighbours and update fields
        particles = DFAfunC.getNeighborsC(particles, dem)
        particles, fields = DFAfunC.updateFieldsC(cfg['GENERAL'], particles, dem, fields)

        iterate = particles['iterate']

        # count iterations
        countIterations = countIterations + 1

    # reset iterate for performing ava simulations
    particles['iterate'] = True

    # remove particles that lie outside of release polygons
    if len(relThField) != 0:
        relRaster = relThField
    particles = resetMassPerParticle(cfg, particles, dem, relRaster)
    particles = com1DFA.checkParticlesInRelease(particles, inputSimLines['releaseLine'],
            cfg['GENERAL'].getfloat('thresholdPointInPoly'))

    # adjust mass of particles in order to match good final mass
    if len(relThField) != 0:
        particles['m'] = particles['mIni']
        particles['mTot'] = np.sum(particles['m'])
    else:
        newMass = np.sum(particles['mIni'])
        oldMass = np.sum(particles['m'])
        particles['m'] = particles['m'] * (newMass / oldMass)
        particles['mTot'] = np.sum(particles['m'])
        log.info('oldMass: %.2f and newMass: %.2f - mass factor: %.2f - total mass: %.2f' %
                (oldMass, newMass, newMass/oldMass, particles['mTot']))

    # reset particles IDs
    particles['ID'] = np.arange(particles['nPart'])
    particles['nID'] = particles['nPart']
    particles['parentID'] = np.arange(particles['nPart'])
    # reset particle properties
    nPart = particles['nPart']
    particles['ux'] = np.zeros(nPart)
    particles['uy'] = np.zeros(nPart)
    particles['uz'] = np.zeros(nPart)
    particles['s'] = np.zeros(nPart)
    particles['l'] = np.zeros(nPart)
    particles['stoppCriteria'] = False
    particles['kineticEne'] = 0.0
    particles['peakKinEne'] = 0.0
    particles['peakMassFlowing'] = 0.0
    particles['potentialEne'] = np.sum(cfg['GENERAL'].getfloat('gravAcc') * particles['m'] * particles['z'])
    # TODO: note particle flow depth is not updated- this is done in updateFieldsC in the next step as the flow
    # depth is currently computed from the mass and an interpolation on the grid

    # for final configuration get neighbors and update fields
    particles = DFAfunC.getNeighborsC(particles, dem)
    particles, fields = DFAfunC.updateFieldsC(cfg['GENERAL'], particles, dem, fields)

    fields['pfv'] = fields['FD']
    fields['ppr'] = fields['P']
    fields['pfd'] = fields['FV']

    # save particles to file for visualisation
    avaDir = pathlib.Path(cfg['GENERAL']['avalancheDir'])
    outDir = avaDir / 'Outputs' / 'com1DFA' / 'particlesIni'
    particleTools.savePartToCsv(cfg['VISUALISATION']['particleProperties'], particlesList, outDir)
    log.info('Initialising particles finalized, total mass: %.2f, number of particles: %d' % (np.sum(particles['m']), particles['nPart']))

    return particles, fields


def resetMassPerParticle(cfg, particles, dem, relRaster):
    """ recompute mass of particles according to their location with respect to relRaster

        Parameters
        ------------
        cfg: configparser object
            configuration settings
        particles: dict
            dictionary with particles properties
        dem: dict
            dictionary with info on dem
        relRaster: np.array
            raster of release thickness values

        Returns
        ---------
        particles: dict
            updated particles dictionary with new mass
    """

    ncols = dem['header']['ncols']
    nrows = dem['header']['nrows']
    indPartInCell = particles['indPartInCell']
    partInCell = particles['partInCell']
    rho = cfg['GENERAL'].getfloat('rho')

    indRelY, indRelX = np.nonzero(relRaster)
    particles['mIni'] = np.zeros(particles['nPart'])
    for indRelx, indRely in zip(indRelX, indRelY):
        # compute number of particles for this cell
        hCell = relRaster[indRely, indRelx]
        volCell = dem['areaRaster'][indRely, indRelx] * hCell
        massCell = volCell * rho
        ic = indRelx + ncols * indRely
        iStart = indPartInCell[ic]
        iStop = indPartInCell[ic+1]
        indParts = partInCell[iStart:iStop]
        mNew = massCell / len(indParts)
        particles['mIni'][indParts] = mNew

    return particles


def createReleaseBuffer(cfg, inputSimLines):
    """ add a buffer around release polygons to get boundary particles

        Parameters
        -----------
        cfg: configparser object
            configuration settings
        inputSimLines: dict
            dictionary with input data info

        Returns
        --------
        inputSimLines: dict
            updated inputSimLines with releaseLineBuffer
    """

    # get start indices and lengths of release polygons
    lengthRels = inputSimLines['releaseLine']['Length']
    startLines = inputSimLines['releaseLine']['Start']
    count = 0
    xBuffered = np.empty(0)
    yBuffered = np.empty(0)
    lengthArray = []
    startArray = [int(startLines[0])]
    for m in range(len(lengthRels)):

        # get coordinates of release line for each feature
        indStart = int(startLines[count])
        indStop = int(startLines[count] + lengthRels[count])
        xRel = inputSimLines['releaseLine']['x'][indStart:indStop]
        yRel = inputSimLines['releaseLine']['y'][indStart:indStop]
        # create list of point tuples
        points = []
        for m in range(len(xRel)):
            pointsT = xRel[m], yRel[m]
            points.append(pointsT)

        # create polygon
        pol1 = Polygon(points)
        # compute bufferZone
        bufferZone = cfg['GENERAL'].getfloat('sphKernelRadius') * cfg['GENERAL'].getfloat('bufferZoneFactor')
        polBuffered = pol1.buffer(bufferZone)
        xnew, ynew = polBuffered.exterior.coords.xy
        xnew = np.asarray(xnew)
        ynew = np.asarray(ynew)
        xBuffered = np.append(xBuffered, xnew)
        yBuffered = np.append(yBuffered, ynew)
        lengthArray.append(len(xnew))
        startArray.append(int(startArray[-1]+len(xnew)))
        count = count + 1

    # add buffered release lines to dictionary
    releaseLineBuffer = inputSimLines['releaseLine'].copy()
    releaseLineBuffer['x'] = xBuffered
    releaseLineBuffer['y'] = yBuffered
    releaseLineBuffer['Length'] = np.asarray(lengthArray)
    releaseLineBuffer['Start'] = np.asarray(startArray)
    releaseLineBuffer['z'] = np.zeros(len(xBuffered))
    inputSimLines['releaseLineBuffer'] = releaseLineBuffer

    # plt.plot(inputSimLines['releaseLine']['x'], inputSimLines['releaseLine']['y'], 'g')
    # plt.plot(xBuffered, yBuffered, 'b')
    # plt.show()

    return inputSimLines
