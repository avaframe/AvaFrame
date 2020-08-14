#!/usr/bin/env python
# coding: utf-8
""" Main file for module com2AB - Alpha Beta
"""

import os
import glob
import pickle
import logging
import numpy as np
import matplotlib.pyplot as plt

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def setEqParameters(smallAva, customParam):
    """Set alpha beta equation parameters to
    - standard (default)
    - small avalanche
    - custom
    TODO: test
    """

    eqParameters = {}

    if smallAva is True:
        log.debug('Using small Avalanche Setup')
        eqParameters['k1'] = 0.933
        eqParameters['k2'] = 0.0
        eqParameters['k3'] = 0.0088
        eqParameters['k4'] = -5.02
        eqParameters['SD'] = 2.36

        ParameterSet = "Small avalanches"

    elif customParam:
        log.debug('Using custom Avalanche Setup')
        eqParameters['k1'] = customParam['k1']
        eqParameters['k2'] = customParam['k2']
        eqParameters['k3'] = customParam['k3']
        eqParameters['k4'] = customParam['k4']
        eqParameters['SD'] = customParam['SD']

        ParameterSet = "Custom"

    else:
        log.debug('Using standard Avalanche Setup')
        eqParameters['k1'] = 1.05
        eqParameters['k2'] = -3130.0
        eqParameters['k3'] = 0.0
        eqParameters['k4'] = -2.38
        eqParameters['SD'] = 1.25

        ParameterSet = "Standard"
    eqParameters['ParameterSet'] = ParameterSet
    return eqParameters


def com2ABMain(dem, Avapath, splitPoint, saveOutPath, cfgsetup):
    """ Loops on the given Avapath and runs com2AB to compute AlpahBeta model
    Inputs : dem header and rater (as np array),
            Avapath and Split points .shp file,
            optional output save path,
            avalanche type,
            reamplind lenght for the Avapath
    Outputs : writes raw results to saveOutPath
    """
    smallAva = cfgsetup.getboolean('smallAva')
    customParam = cfgsetup.getboolean('customParam')
    # customParam = str(customParam or None)
    distance = float(cfgsetup['distance'])

    NameAva = Avapath['Name']
    StartAva = Avapath['Start']
    LengthAva = Avapath['Length']

    for i in range(len(NameAva)):
        name = NameAva[i]
        start = StartAva[i]
        end = start + LengthAva[i]
        avapath = {}
        avapath['x'] = Avapath['x'][int(start):int(end)]
        avapath['y'] = Avapath['y'][int(start):int(end)]
        avapath['Name'] = name
        com2AB(dem, avapath, splitPoint, saveOutPath,
               smallAva, customParam, distance)


def com2AB(dem, avapath, splitPoint, OutPath,
           smallAva, customParam, distance):
    """ Computes the AlphaBeta model given an input raster (of the dem),
    an avalanche path and split points
    Inputs : dem header and rater (as np array),
            single avapath as np array,
            Split points as np array,
            output save path,
            avalanche type,
            resamplind lenght for the Avapath
    Outputs : writes raw results to OutPath
    """
    name = avapath['Name']
    abVersion = '4.1'
    log.info('Running Alpha Beta %s on: %s ', abVersion, name)
    eqParams = setEqParameters(smallAva, customParam)

    # TODO: make rest work with dict

    # read inputs, ressample ava path
    # make pofile and project split point on path
    AvaProfile, projSplitPoint, splitPoint = geoTrans.prepareLine(
        dem, avapath, splitPoint, distance)

    # Sanity check if first element of AvaProfile[3,:]
    # (i.e z component) is highest:
    # if not, flip all arrays
    projSplitPoint, AvaProfile = geoTrans.checkProfile(AvaProfile, projSplitPoint)

    AvaProfile['indSplit'] = projSplitPoint['indSplit']  # index of split point

    eqOut = calcAB(AvaProfile, eqParams)
    savename = name + '_com2AB_eqparam.pickle'
    save_file = os.path.join(OutPath, savename)
    with open(save_file, 'wb') as handle:
        pickle.dump(eqParams, handle, protocol=pickle.HIGHEST_PROTOCOL)
    savename = name + '_com2AB_eqout.pickle'
    save_file = os.path.join(OutPath, savename)
    with open(save_file, 'wb') as handle:
        pickle.dump(eqOut, handle, protocol=pickle.HIGHEST_PROTOCOL)


def readABinputs(cfgAva):

    cfgPath = {}

    profileLayer = glob.glob(cfgAva + '/Inputs/LINES/*AB*.shp')
    cfgPath['profileLayer'] = ''.join(profileLayer)

    demSource = glob.glob(cfgAva + '/Inputs/*.asc')
    try:
        assert len(demSource) == 1, 'There should be exactly one topography .asc file in ' + \
            cfgAva + '/Inputs/'
    except AssertionError:
        raise

    cfgPath['demSource'] = ''.join(demSource)

    splitPointSource = glob.glob(cfgAva + '/Inputs/POINTS/*.shp')
    cfgPath['splitPointSource'] = ''.join(splitPointSource)

    saveOutPath = os.path.join(cfgAva, 'Outputs/com2AB/')
    if not os.path.exists(saveOutPath):
        # log.info('Creating output folder %s', saveOutPath)
        os.makedirs(saveOutPath)
    cfgPath['saveOutPath'] = saveOutPath

    defaultName = str(cfgAva).split('/')[-1]
    cfgPath['defaultName'] = defaultName

    return cfgPath


def find_10Point(tmp, delta_ind):
    """ find the beta point: first point under 10°
     (make sure that the delta_ind next indexes are also under 10°)
     otherwise keep looking
     """
    i = 0
    while True:
        ind = tmp[0][i]
        condition = True
        for j in range(delta_ind):
            condition = condition and (tmp[0][i+j+1] == ind+j+1)
            if not condition:
                i = i + j + 1
                break
        if condition:
            ids_10Point = ind - 1
            break
    return ids_10Point


def calcAB(AvaProfile, eqParameters):
    """
    Calculate Alpha Beta for data in eqInput according to chosen eqParameters
    """
    log.debug("Calculating alpha beta")
    k1 = eqParameters['k1']
    k2 = eqParameters['k2']
    k3 = eqParameters['k3']
    k4 = eqParameters['k4']
    SD = eqParameters['SD']

    s = AvaProfile['s']
    z = AvaProfile['z']
    distance = s[1] - s[0]
    delta_ind = max(int(np.floor(30/distance)), 1)
    indSplit = AvaProfile['indSplit']
    ds = np.abs(s - np.roll(s, 1))
    dz = np.abs(z - np.roll(z, 1))
    ds[0] = 0.0
    dz[0] = 0.0
    angle = np.rad2deg(np.arctan2(dz, ds))
    CuSplit = s[indSplit]
    # TODO SPLIT POINT READING
    # get all values where Angle < 10 but >0
    # get index of first occurance and go one back to get previous value
    # (i.e. last value above 10 deg)
    # tmp = x[(angle < 10.0) & (angle > 0.0) & (x > 450)]

    tmp = np.where((angle < 10.0) & (angle > 0.0) & (s > CuSplit))

    # find the beta point: first point under 10°
    # (make sure that the 30 next meters are also under 10°)
    ids_10Point = find_10Point(tmp, delta_ind)
    if debugPlot:
        plt.figure(figsize=(10, 6))
        plt.plot(s, angle)
        plt.plot(s[ids_10Point], angle[ids_10Point], 'or')
        plt.axhline(y=10, color='0.8',
                    linewidth=1, linestyle='-.', label='10^\circ line')
        plt.show()

    # Do a quadtratic fit and get the polynom for 2nd derivative later
    zQuad = np.polyfit(s, z, 2)
    poly = np.poly1d(zQuad)
    # Get H0: max - min for parabola
    H0 = max(poly(s)) - min(poly(s))
    # get beta
    dz_beta = z[0] - z[ids_10Point]
    beta = np.rad2deg(np.arctan2(dz_beta, s[ids_10Point]))
    # get Alpha
    alpha = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4

    # get Alpha standard deviations
    SDs = [SD, -1*SD, -2*SD]
    alphaSD = k1 * beta + k2 * poly.deriv(2)[0] + k3 * H0 + k4 + SDs

    AvaProfile['CuSplit'] = CuSplit
    AvaProfile['ids_10Point'] = ids_10Point
    AvaProfile['poly'] = poly
    AvaProfile['beta'] = beta
    AvaProfile['alpha'] = alpha
    AvaProfile['SDs'] = SDs
    AvaProfile['alphaSD'] = alphaSD
    return AvaProfile
