#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import logging
import os
import platform
if os.name == 'nt':
    from multiprocessing.pool import Pool as Pool    
elif platform.system() == 'Darwin':
    from multiprocessing.pool import ThreadPool as Pool
else:
    from multiprocessing import Pool

from avaframe.com4FlowPy.flowClass import Cell


def get_start_idx(dem, release):
    """Sort Release Pixels by altitude and return the result as lists for the
    Rows and Columns, starting with the highest altitude

    Input parameters:
        dem         Digital Elevation Model to gain information about altitude
        release     The release layer, release pixels need int value > 0

    Output parameters:
        row_list    Row index of release pixels sorted by altitude
        col_list    Column index of release pixels sorted by altitude
        """
    row_list, col_list = np.where(release > 0)  # Gives back the indices of the release areas
    if len(row_list) > 0:
        altitude_list = []
        for i in range(len(row_list)):
            altitude_list.append(dem[row_list[i], col_list[i]])
        altitude_list, row_list, col_list = list(zip(*sorted(zip(altitude_list, row_list, col_list), reverse=True)))
        # Sort this lists by altitude
    return row_list, col_list

def split_release(release, pieces):
    """Split the release layer in several tiles, the number is depending on the
    available CPU Cores, so every Core gets one tile. The area is determined by
    the number of release pixels in it, so that every tile has the same amount
    of release pixels in it. Splitting in x(Columns) direction.
    The release tiles have still the size of the original layer, so no split
    for the DEM is needed.

    Input parameters:
        release         the release layer with release pixels as int > 0
        header_release  the header of the release layer to identify the
                        noDataValue

    Output parameters:
        release_list    A list with the tiles(arrays) in it [array0, array1, ..]
        """

    release[release < 0] = 0
    release[release > 1] = 1
    sumRelease = np.sum(release) # Count number of release pixels
    sum_per_split = sumRelease/pieces  # Divide the number by avaiable Cores
    release_list = []
    #2022-09-06 - AH: Note - maybe we need to think of sth. smarter here!!
    #i.e. not slicing in columns??    
    breakpoint_x = 0
    for i in range(release.shape[1]):
        if (len(release_list) == (pieces - 1)) or (np.sum(release[:, i:]) <= sum_per_split):
            c = np.zeros_like(release)
            c[:, breakpoint_x:] = release[:, breakpoint_x:]
            release_list.append(c)
            break
        if np.sum(release[:, breakpoint_x:i]) < sum_per_split:
            continue
        else:
            c = np.zeros_like(release)
            c[:, breakpoint_x:i] = release[:, breakpoint_x:i]
            release_list.append(c)
            print("Release Split from {} to {}".format(breakpoint_x, i))
            breakpoint_x = i
    
    return release_list

def back_calculation(back_cell):
    """Here the back calculation from a run out pixel that hits a infrastructure
    to the release pixel is performed.

    Input parameters:
        hit_cell_list        All cells that hit a Infrastructure

    Output parameters:
        Back_list   List of pixels that are on the way to the start cell
                    Maybe change it to array like DEM?
    """
    #start = time.time()
    #if len(hit_cell_list) > 1:
        #hit_cell_list.sort(key=lambda cell: cell.altitude, reverse=False)
        #print("{} Elements sorted!".format(len(hit_cell_list)))
    back_list = []
    for parent in back_cell.parent:
        if parent not in back_list:
            back_list.append(parent)
    for cell in back_list:
        for parent in cell.parent:
            # Check if parent already in list
            if parent not in back_list:
                back_list.append(parent)
    #end = time.time()
    #print('\n Backcalculation needed: ' + str(end - start) + ' seconds')
    return back_list

def run(optTuple):
    
    log = logging.getLogger(__name__)
    tempDir = optTuple[8]
    infraBool = optTuple[9]
    nCPU = optTuple[10]

    dem = np.load(tempDir / ("dem_%s_%s.npy" % (optTuple[0], optTuple[1])))
    release = np.load(tempDir / ("init_%s_%s.npy" % (optTuple[0], optTuple[1])))
    if infraBool:
        infra = np.load(tempDir / ("infra_%s_%s.npy" % (optTuple[0], optTuple[1])))
    else:
        infra = np.zeros_like(dem)

    alpha = float(optTuple[2])
    exp = float(optTuple[3])
    cellsize = float(optTuple[4])
    nodata = float(optTuple[5])
    flux_threshold = float(optTuple[6])
    max_z_delta = float(optTuple[7])

    log.info("Multiprocessing starts, used cores: %i" % (nCPU))
        
    release_list = split_release(release, nCPU)
    
    with Pool(processes=nCPU) as pool:
        results = pool.map(calculation,[[dem, infra, release_sub, alpha, exp, flux_threshold, max_z_delta, nodata, cellsize, infraBool]
                            for release_sub in release_list])
        pool.close()
        pool.join()

    z_delta_array = np.zeros_like(dem)
    flux_array = np.zeros_like(dem)
    count_array = np.zeros_like(dem)
    z_delta_sum = np.zeros_like(dem)
    backcalc = np.zeros_like(dem)
    fp_travelangle_array = np.zeros_like(dem)
    sl_travelangle_array = np.zeros_like(dem)

    z_delta_list = []
    flux_list = []
    cc_list = []
    z_delta_sum_list = []
    backcalc_list = []
    fp_ta_list = []
    sl_ta_list = []
    
    for i in range(len(results)):
        res = results[i]
        res = list(res)
        z_delta_list.append(res[0])
        flux_list.append(res[1])
        cc_list.append(res[2])
        z_delta_sum_list.append(res[3])
        backcalc_list.append(res[4])
        fp_ta_list.append(res[5])
        sl_ta_list.append(res[6])

    logging.info('Calculation finished, getting results.')
    for i in range(len(z_delta_list)):
        z_delta_array = np.maximum(z_delta_array, z_delta_list[i])
        flux_array = np.maximum(flux_array, flux_list[i])
        count_array += cc_list[i]
        z_delta_sum += z_delta_sum_list[i]
        backcalc = np.maximum(backcalc, backcalc_list[i])
        fp_travelangle_array = np.maximum(fp_travelangle_array, fp_ta_list[i])
        sl_travelangle_array = np.maximum(sl_travelangle_array, sl_ta_list[i])        
        
    # Save Calculated tiles
    np.save(tempDir / ("res_z_delta_%s_%s" % (optTuple[0], optTuple[1])), z_delta_array)
    np.save(tempDir / ("res_z_delta_sum_%s_%s" % (optTuple[0], optTuple[1])), z_delta_sum)
    np.save(tempDir / ("res_flux_%s_%s" % (optTuple[0], optTuple[1])), flux_array)
    np.save(tempDir / ("res_count_%s_%s" % (optTuple[0], optTuple[1])), count_array)
    np.save(tempDir / ("res_fp_%s_%s" % (optTuple[0], optTuple[1])), fp_travelangle_array)
    np.save(tempDir / ("res_sl_%s_%s" % (optTuple[0], optTuple[1])), sl_travelangle_array)
    if infraBool:
        np.save(tempDir / ("res_backcalc_%s_%s" % (optTuple[0], optTuple[1])), backcalc)

def calculation(args):
    """This is the core function where all the data handling and calculation is
    done.

    Input parameters:
        dem         The digital elevation model
        header      The header of the elevation model
        infra       The infra layer
        release     The list of release arrays
        alpha
        exp
        flux_threshold
        max_z_delta

    Output parameters:
        z_delta     Array like DEM with the max. kinetic Energy Height for every
                    pixel
        flux_array  Array with max. concentration factor saved
        count_array Array with the number of hits for every pixel
        elh_sum     Array with the sum of Energy Line Height
        back_calc   Array with back calculation, still to do!!!
        """
    
    dem = args[0]
    infra = args[1]
    release = args[2]
    alpha = args[3]
    exp = args[4]
    flux_threshold = args[5]
    max_z_delta=args[6]
    nodata = args[7]
    cellsize = args[8]
    infraBool = args[9]

    z_delta_array = np.zeros_like(dem, dtype=np.float32)
    z_delta_sum = np.zeros_like(dem, dtype=np.float32)
    flux_array = np.zeros_like(dem, dtype=np.float32)
    count_array = np.zeros_like(dem, dtype=np.int32)
    
    fp_travelangle_array = np.zeros_like(dem, dtype=np.float32)  # fp = Flow Path
    sl_travelangle_array = np.zeros_like(dem, dtype=np.float32) * 90  # sl = Straight Line
    
    backcalc = np.zeros_like(dem, dtype=np.int32)
    
    if infraBool:        
        back_list = []

    # Core
    #start = datetime.now().replace(microsecond=0)
    row_list, col_list = get_start_idx(dem, release)

    startcell_idx = 0
    while startcell_idx < len(row_list):

        sys.stdout.write('\r' "Calculating Startcell: " + str(startcell_idx + 1) + " of " + str(len(row_list)) + " = " + str(
            round((startcell_idx + 1) / len(row_list) * 100, 2)) + "%" '\r')
        sys.stdout.flush()

        cell_list = []
        row_idx = row_list[startcell_idx]
        col_idx = col_list[startcell_idx]
        dem_ng = dem[row_idx - 1:row_idx + 2, col_idx - 1:col_idx + 2]  # neighbourhood DEM
        if (nodata in dem_ng) or np.size(dem_ng) < 9:
            startcell_idx += 1
            continue

        startcell = Cell(row_idx, col_idx, dem_ng, cellsize, 1, 0, None,
                         alpha, exp, flux_threshold, max_z_delta, startcell=True)
        # If this is a startcell just give a Bool to startcell otherwise the object startcell

        cell_list.append(startcell)

        for idx, cell in enumerate(cell_list):
            row, col, flux, z_delta = cell.calc_distribution()

            if len(flux) > 0:
                # mass, row, col  = list(zip(*sorted(zip( mass, row, col), reverse=False)))

                z_delta, flux, row, col = list(zip(*sorted(zip(z_delta, flux, row, col), reverse=False)))
                # Sort this lists by elh, to start with the highest cell

            for i in range(idx, len(cell_list)):  # Check if Cell already exists
                k = 0
                while k < len(row):
                    if row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex:
                        cell_list[i].add_os(flux[k])
                        cell_list[i].add_parent(cell)
                        if z_delta[k] > cell_list[i].z_delta:
                            cell_list[i].z_delta = z_delta[k]
                        row = np.delete(row, k)
                        col = np.delete(col, k)
                        flux = np.delete(flux, k)
                        z_delta = np.delete(z_delta, k)
                    else:
                        k += 1

            for k in range(len(row)):
                dem_ng = dem[row[k] - 1:row[k] + 2, col[k] - 1:col[k] + 2]  # neighbourhood DEM
                if (nodata in dem_ng) or np.size(dem_ng) < 9:
                    continue
                cell_list.append(
                    Cell(row[k], col[k], dem_ng, cellsize, flux[k], z_delta[k], cell, alpha, exp, flux_threshold, max_z_delta, startcell))

            z_delta_array[cell.rowindex, cell.colindex] = max(z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
            flux_array[cell.rowindex, cell.colindex] = max(flux_array[cell.rowindex, cell.colindex], cell.flux)
            count_array[cell.rowindex, cell.colindex] += int(1)
            z_delta_sum[cell.rowindex, cell.colindex] += cell.z_delta
            fp_travelangle_array[cell.rowindex, cell.colindex] = max(fp_travelangle_array[cell.rowindex, cell.colindex], cell.max_gamma)
            sl_travelangle_array[cell.rowindex, cell.colindex] = max(sl_travelangle_array[cell.rowindex, cell.colindex], cell.sl_gamma)

            #Backcalculation
            if infraBool:
                if infra[cell.rowindex, cell.colindex] > 0:
                    #backlist = []
                    back_list = back_calculation(cell)
    
                    for back_cell in back_list:
                        backcalc[back_cell.rowindex, back_cell.colindex] = max(backcalc[back_cell.rowindex, back_cell.colindex],
                                                                               infra[cell.rowindex, cell.colindex])
        if infraBool:
            release[z_delta_array > 0] = 0
            # Check if i hit a release Cell, if so set it to zero and get again the indexes of release cells
            row_list, col_list = get_start_idx(dem, release)
        startcell_idx += 1
    #end = datetime.now().replace(microsecond=0)
    return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array


