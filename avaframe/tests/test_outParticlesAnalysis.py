''' Tests for module outParticlesAnalysis  '''
import numpy as np
import numpy.ma as ma
import pandas as pd
import pathlib
import configparser
import matplotlib.pyplot as plt

# Local imports
import avaframe.out3Plot.outParticlesAnalysis as oA
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


def test_velocityEnvelope():
    """ test computing envelope of particlesTimeArrays """

    # setup required input
    particlesTimeArrays = {'t': [1., 2., 3., 4., 5.]}
    velMag = np.zeros((5, 10))
    velMag[1:,1:] = 1.
    velMag[1:,1] = 5.
    velMag[2, 3] = 7
    trajectoryLengthXYZ = np.zeros((5, 10))
    trajectoryLengthXYZ[1,:] = 2.
    trajectoryLengthXYZ[2,:] = 3.
    trajectoryLengthXYZ[3,:] = 5.
    trajectoryLengthXYZ[4,:] = 7.

    uAcc = np.zeros((5, 10))
    uAcc[1:,1:] = 4.
    particlesTimeArrays['velocityMag'] = velMag
    particlesTimeArrays['trajectoryLengthXYZ'] = trajectoryLengthXYZ
    particlesTimeArrays['uAcc'] = uAcc

    # call function to be tested
    dictVelEnvelope = oA.velocityEnvelope(particlesTimeArrays)

    # TODO check what if a particle has been merged or removed what happens to array
    b = 19./10.
    assert np.array_equal(np.asarray(dictVelEnvelope['Time']), np.asarray([1., 2., 3., 4., 5.]))
    assert np.array_equal(dictVelEnvelope['Velocity'], velMag)
    assert np.array_equal(dictVelEnvelope['Acc'], uAcc)
    assert np.array_equal(dictVelEnvelope['Mean'], np.asarray([0, 1.3, b, 1.3, 1.3]))
    assert np.array_equal(dictVelEnvelope['Max'], np.asarray([0, 5., 7., 5., 5.]))
    assert np.array_equal(dictVelEnvelope['Min'], np.asarray([0., 0., 0., 0., 0.]))
    assert np.array_equal(dictVelEnvelope['Median'], np.asarray([0., 1., 1., 1., 1.]))
    assert np.array_equal(dictVelEnvelope['SxyzMax'], np.asarray([0., 2., 3., 5., 7.]))
    assert np.array_equal(dictVelEnvelope['SxyzMax'], np.asarray([0., 2., 3., 5., 7.]))


def test_velocityEnvelopeThalweg():
    """ test computing envelope  along the thalweg of particlesTimeArrays """

    # setup required input
    particlesTimeArrays = {'t': [1., 2., 3., 4., 5.]}
    velMag = np.zeros((3, 5))
    velMag[1:,1:] = 1.
    velMag[1:,1] = 5.
    velMag[2, 3] = 7
    trajectoryLengthXYZ = np.zeros((3, 5))
    trajectoryLengthXYZ[1,:] = np.asarray([0., 2., 3., 4., 5.])
    trajectoryLengthXYZ[2,:] = np.asarray([0., 3., 4., 6., 7.])

    sAimec = np.zeros((3, 5))
    sAimec[1,:] = np.asarray([0., 1., 2., 3., 4.])
    sAimec[2,:] = np.asarray([0., 4., 5., 7., 8.])


    elevation = np.zeros((3, 5))
    elevation[1,:] = np.asarray([0., 10., 20., 30., 40.])
    elevation[2,:] = np.asarray([0., 20., 30., 40., 50.])


    uAcc = np.zeros((3, 5))
    uAcc[1:,1:] = 4.
    particlesTimeArrays['velocityMag'] = velMag
    particlesTimeArrays['trajectoryLengthXYZ'] = trajectoryLengthXYZ
    particlesTimeArrays['uAcc'] = uAcc
    particlesTimeArrays['sBetaPoint'] = 4.
    particlesTimeArrays['sAimec'] = sAimec
    particlesTimeArrays['z'] = elevation

    # call function to be tested
    dictVelAltThalweg = oA.velocityEnvelopeThalweg(particlesTimeArrays)

    velMagMax = np.asarray([0.0, 5.0, 1., 1., 5., 1., 7., 1.])
    velMagMin = np.asarray([0.0, 5.0, 1., 1., 1., 1., 7., 1.])
    velMagMean = np.asarray([0.0, 5.0, 1., 1., 3., 1., 7., 1.])
    sxyzMax = np.asarray([0., 2., 3., 4., 5., 4., 6., 7.])
    sxyzMin = np.asarray([0., 2., 3., 4., 3., 4., 6., 7.])

    assert np.array_equal(dictVelAltThalweg['maxVelocity'], velMagMax)
    assert np.array_equal(dictVelAltThalweg['minVelocity'], velMagMin)
    assert np.array_equal(dictVelAltThalweg['meanVelocity'], velMagMean)
    assert np.array_equal(dictVelAltThalweg['medianVelocity'], velMagMean)
    assert np.array_equal(dictVelAltThalweg['maxSxyz'], sxyzMax)
    assert np.array_equal(dictVelAltThalweg['minSxyz'], sxyzMin)
