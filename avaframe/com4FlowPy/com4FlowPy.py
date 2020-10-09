#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:23:00 2018

@author: Michael Neuhauser
"""
# import standard libraries
import os
import glob
import sys
import psutil
import numpy as np
from datetime import datetime
from multiprocessing import cpu_count
import multiprocessing as mp
import logging

# Flow-Py Libraries
import avaframe.com4FlowPy.raster_io as io
import avaframe.com4FlowPy.flow_core_gui as fc

# create local logger
log = logging.getLogger(__name__)


def readFlowPyinputs(cfgAva):
    cfgPath = {}
    # read release pixels
    releasePath = glob.glob(os.path.join(cfgAva, 'Inputs', 'FlowPy', '*rel*.tif'))
    try:
        message = 'There should be exactly one *rel*.tif file containing the release area in ' + cfgAva + '/Inputs/flowPy/'
        assert len(releasePath) == 1, message
    except AssertionError:
        raise
    cfgPath['releasePath'] = ''.join(releasePath)
    # read infra pixel
    infraPath = glob.glob(cfgAva + '/Inputs/FlowPy/*infra*.tif')
    if len(infraPath) == 0:
        infraPath = None
    cfgPath['infraPath'] = ''.join(infraPath)

    # read DEM
    demSource = glob.glob(cfgAva + '/Inputs/*.asc')
    try:
        assert len(demSource) == 1, 'There should be exactly one topography .asc file in ' + \
            cfgAva + '/Inputs/'
    except AssertionError:
        raise
    cfgPath['demSource'] = ''.join(demSource)

    # make output path
    saveOutPath = os.path.join(cfgAva, 'Outputs/com4FlowPy/')
    if not os.path.exists(saveOutPath):
        # log.info('Creating output folder %s', saveOutPath)
        os.makedirs(saveOutPath)
    cfgPath['saveOutPath'] = saveOutPath

    return cfgPath


def com4FlowPyMain(cfgPath, cfgSetup):

    alpha = float(cfgSetup['alpha'])
    exp = float(cfgSetup['exp'])
    process = cfgSetup['process']
    saveOutPath = cfgPath['saveOutPath']
    dem_path = cfgPath['demSource']
    release_path = cfgPath['releasePath']
    infra_path = cfgPath['infraPath']

    start = datetime.now().replace(microsecond=0)
    calc_bool = False
    # Create result directory
    time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
    # Start of Calculation
    log.info('Start Calculation')
    log.info('Alpha Angle: {}'.format(alpha))
    log.info('Exponent: {}'.format(exp))
    # Read in raster files
    try:
        dem, header = io.read_raster(dem_path)
        log.info('DEM File: {}'.format(dem_path))
    except FileNotFoundError:
        log.error("DEM: Wrong filepath or filename")
        return

    try:
        release, release_header = io.read_raster(release_path)
        log.info('Release File: {}'.format(release_path))
    except FileNotFoundError:
        log.error("Wrong filepath or filename")
        return

    # Check if Layers have same size!!!
    if header['ncols'] == release_header['ncols'] and header['nrows'] == release_header['nrows']:
        log.info("DEM and Release Layer ok!")
    else:
        log.error("Error: Release Layer doesn't match DEM!")
        return

    try:
        infra, infra_header = io.read_raster(infra_path)
        if header['ncols'] == infra_header['ncols'] and header['nrows'] == infra_header['nrows']:
            log.info("Infra Layer ok!")
            calc_bool = True
            log.info('Infrastructure File: {}'.format(infra_path))
        else:
            log.error("Error: Infra Layer doesn't match DEM!")
            return
    except:
        infra = np.zeros_like(dem)

    log.info('Files read in')

    log.info('Process: {}'.format(process))
    z_delta = np.zeros_like(dem)
    susc = np.zeros_like(dem)
    cell_counts = np.zeros_like(dem)
    z_delta_sum = np.zeros_like(dem)
    backcalc = np.zeros_like(dem)
    fp_ta = np.zeros_like(dem)
    sl_ta = np.zeros_like(dem)

    avaiable_memory = psutil.virtual_memory()[1]
    needed_memory = sys.getsizeof(dem)

    max_number_procces = int(avaiable_memory / (needed_memory * 10))

    log.info(
        "There are {} Bytes of Memory avaiable and {} Bytes needed per process. Max. Nr. of Processes = {}".format(
            avaiable_memory, needed_memory*10, max_number_procces))

    # Calculation
    log.info('Multiprocessing starts, used cores: {}'.format(cpu_count()))

    if calc_bool:
        release_list = fc.split_release(release, release_header, min(mp.cpu_count() * 2, max_number_procces))

        log.info("{} Processes started.".format(len(release_list)))
        pool = mp.Pool(len(release_list))
        results = pool.map(fc.calculation,
                           [[dem, header, infra, process, release_pixel, alpha, exp]
                            for release_pixel in release_list])
        pool.close()
        pool.join()
    else:
        release_list = fc.split_release(release, release_header, min(mp.cpu_count() * 4, max_number_procces))

        log.info("{} Processes started.".format(len(release_list)))
        pool = mp.Pool(mp.cpu_count())
        # results = pool.map(gc.calculation, iterable)
        results = pool.map(fc.calculation_effect,
                           [[dem, header, process, release_pixel, alpha, exp] for
                            release_pixel in release_list])
        pool.close()
        pool.join()

    z_delta_list = []
    susc_list = []
    cc_list = []
    z_delta_sum_list = []
    backcalc_list = []
    fp_ta_list = []
    sl_ta_list = []
    for i in range(len(results)):
        res = results[i]
        res = list(res)
        z_delta_list.append(res[0])
        susc_list.append(res[1])
        cc_list.append(res[2])
        z_delta_sum_list.append(res[3])
        backcalc_list.append(res[4])
        fp_ta_list.append(res[5])
        sl_ta_list.append(res[6])

    log.info('Calculation finished, getting results.')
    for i in range(len(z_delta_list)):
        z_delta = np.maximum(z_delta, z_delta_list[i])
        susc = np.maximum(susc, susc_list[i])
        cell_counts += cc_list[i]
        z_delta_sum += z_delta_sum_list[i]
        backcalc = np.maximum(backcalc, backcalc_list[i])
        fp_ta = np.maximum(fp_ta, fp_ta_list[i])
        sl_ta = np.maximum(sl_ta, sl_ta_list[i])

    if process == 'Avalanche':
        proc = 'ava'
    if process == 'Rockfall':
        proc = 'rf'
    if process == 'Soil Slides':
        proc = 'ds'
    # time_string = datetime.now().strftime("%Y%m%d_%H%M%S")
    log.info('Writing Output Files')
    output_format = '.tif'
    io.output_raster(dem_path,
                     saveOutPath + "susceptibility_{}{}".format(proc, output_format),
                     susc)
    io.output_raster(dem_path,
                     saveOutPath + "z_delta_{}{}".format(proc, output_format),
                     z_delta)
    io.output_raster(dem_path,
                     saveOutPath + "FP_travel_angle_{}{}".format(proc, output_format),
                     fp_ta)
    io.output_raster(dem_path,
                     saveOutPath + "SL_travel_angle_{}{}".format(proc, output_format),
                     sl_ta)
    if not calc_bool:  # if no infra
        io.output_raster(dem_path,
                         saveOutPath + "cell_counts_{}{}".format(proc, output_format),
                         cell_counts)
        io.output_raster(dem_path,
                         saveOutPath + "z_delta_sum_{}{}".format(proc, output_format),
                         z_delta_sum)
    if calc_bool:  # if infra
        io.output_raster(dem_path,
                         saveOutPath + "backcalculation_{}{}".format(proc, output_format),
                         backcalc)

    log.info("Calculation finished")
    end = datetime.now().replace(microsecond=0)
    log.info('Calculation needed: ' + str(end - start) + ' seconds')
