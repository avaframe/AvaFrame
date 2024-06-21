#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import logging
import os
import platform
import gc
import psutil
import time

if os.name == "nt":
    from multiprocessing.pool import Pool as Pool
elif platform.system() == "Darwin":
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
    """Split the release layer in several tiles. The area is determined by
    the number of release pixels in it, so that every tile has the same amount
    of release pixels in it.

    In this version the split is performed along a flattened 2D-array to ensure
    a more even splitting of release pixels than just along the x-Axis ...

    NOTE: TO DO: Ideally a 'greedy' algorithm would let idle CPU cores 'snatch' any
    un-processed release cell until all releaseCells are handled -- this would
    ensure that the total workload is distributed evenly along all CPUs (which
    becomes an important factor for bigger model areas) !!!

    The release tiles have still the size of the original layer, so no split
    for the DEM is needed.

    Parameters
    -----------
    release: np.array - assumes a binary 0|1 array with  release pixels designated by '1'
    pieces:  int - number of chunck in which the release layer should be split

    Returns
    -----------
    release_list:    A list with the tiles(arrays) in it [array0, array1, ..]
    """

    # Flatten the array and compute the cumulative sum
    flat_release = release.flatten()
    cumulative_sum = np.cumsum(flat_release)

    total_sum = cumulative_sum[-1]
    sum_per_split = total_sum / pieces

    release_list = []
    start_index = 0

    for i in range(1, pieces):
        # Find the split point in the flattened array
        split_index = np.searchsorted(cumulative_sum, sum_per_split * i)

        # Create a new array for this split
        split_flat = np.zeros_like(flat_release)
        split_flat[start_index:split_index] = flat_release[start_index:split_index]

        # Reshape the flat array back to 2D and add to the list
        release_list.append(split_flat.reshape(release.shape))

        start_index = split_index

    # Handle the last piece
    split_flat = np.zeros_like(flat_release)
    split_flat[start_index:] = flat_release[start_index:]
    release_list.append(split_flat.reshape(release.shape))

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

    back_list = []
    for parent in back_cell.lOfParents:
        if parent not in back_list:
            back_list.append(parent)
    for cell in back_list:
        for parent in cell.lOfParents:
            # Check if parent already in list
            if parent not in back_list:
                back_list.append(parent)

    return back_list


def run(optTuple):
    """This is a wrapper around calculation() for performing model runs for a single tile of the model domain
    using multiprocessing across multiple CPUs

    Parameters:
    ----------
    optTuple: tuple with all necessary model information
        - optTuple[0] - int - i index of the processed tile (for loading correct data for tiles)
        - optTuple[1] - int - j index of the processed tile (for loading correct data for tiles)

        - optTuple[2] - {} - dict containing modelParameters
        - optTuple[3] - {} - dict containing modelPaths
        - optTuple[4] - {} - dict containing rasterAttributes
        - optTuple[5] - {} - dict containing forestParameters
        - optTuple[6] - {} - dict containing MPOptions

    Returns:
    ----------
    Nothing --> saves results for processed tile to temp folder (modelPaths["tempFolder"])
    """

    log = logging.getLogger(__name__)

    # Flow-Py parameters
    alpha = float(optTuple[2]["alpha"])
    exp = float(optTuple[2]["exp"])
    flux_threshold = float(optTuple[2]["flux_threshold"])
    max_z_delta = float(optTuple[2]["max_z"])
    infraBool = optTuple[2]["infraBool"]
    forestBool = optTuple[2]["forestBool"]

    # Temp-Dir (all input files are located here and results are written back in here)
    tempDir = optTuple[3]["tempDir"]

    # raster-layer Attributes
    cellsize = float(optTuple[4]["cellsize"])
    nodata = float(optTuple[4]["nodata"])

    MPOptions = optTuple[6]  # CPU, Multiprocessing options ...

    dem = np.load(tempDir / ("dem_%s_%s.npy" % (optTuple[0], optTuple[1])))
    release = np.load(tempDir / ("init_%s_%s.npy" % (optTuple[0], optTuple[1])))
    if infraBool:
        infra = np.load(tempDir / ("infra_%s_%s.npy" % (optTuple[0], optTuple[1])))
    else:
        infra = np.zeros_like(dem)

    # if forestBool == 'True'
    # --> load forestFile
    # --> read parametersOfForestExtension
    # NOTE-TODO: this is a quick work-around to simply include forest information - should probably be handled more
    # elegantly/explicitly AND error handling should probably be included
    if forestBool:
        forestArray = np.load(tempDir / ("forest_%s_%s.npy" % (optTuple[0], optTuple[1])))
        forestParams = optTuple[5]

    # convert release areas to binary (0: no release areas, 1: release areas)
    # every positive value >0 is interpreted as release area
    release[release < 0] = 0
    release[release > 0] = 1

    nRel = np.sum(release)
    log.info("Number of release cells: %i" % nRel)

    nProcesses, nChunks = calculateMultiProcessingOptions(
        nRel,
        MPOptions["nCPU"],
        procPerCPU=MPOptions["procPerCPU"],
        maxChunks=MPOptions["maxChunks"],
        chunkSize=MPOptions["chunkSize"],
    )

    release_list = split_release(release, nChunks)
    log.info("Multiprocessing starts, used Cores/Processes/Chunks: %i/%i/%i" % (MPOptions["nCPU"], nProcesses, nChunks))

    if forestBool:
        with Pool(processes=nProcesses) as pool:
            results = pool.map(
                calculation,
                [
                    [
                        dem, infra, release_sub,
                        alpha, exp, flux_threshold, max_z_delta,
                        nodata, cellsize,
                        infraBool, forestBool,
                        forestArray, forestParams,
                    ]
                    for release_sub in release_list
                ],
            )
            pool.close()
            pool.join()
    else:
        with Pool(processes=nProcesses) as pool:
            results = pool.map(
                calculation,
                [
                    [
                        dem, infra, release_sub,
                        alpha, exp, flux_threshold, max_z_delta,
                        nodata, cellsize,
                        infraBool, forestBool,
                    ]
                    for release_sub in release_list
                ],
            )
            pool.close()
            pool.join()

    # TODO - provide option in com4FlowPyCfg.ini file for which output layers to write
    # e.g.: default  [zDelta, cellCounts, fpTravelAngle, (backCalc if Infra)]
    #       optional [flux, slTravelAngle]
    # TODO - NOTE:
    # also move this part into a separate function
    # initializing arrays for storing the results from the multiprocessing step
    z_delta_array = np.zeros_like(dem, dtype=np.float32)
    flux_array = np.zeros_like(dem, dtype=np.float32)
    count_array = np.zeros_like(dem, dtype=np.int32)
    z_delta_sum = np.zeros_like(dem, dtype=np.float32)
    backcalc = np.zeros_like(dem, dtype=np.int32)
    fp_travelangle_array = np.zeros_like(dem, dtype=np.float32)
    sl_travelangle_array = np.zeros_like(dem, dtype=np.float32)
    travel_length_array = np.zeros_like(dem, dtype=np.float32)

    z_delta_list = []
    flux_list = []
    cc_list = []
    z_delta_sum_list = []
    backcalc_list = []
    fp_ta_list = []
    sl_ta_list = []
    travel_length_list = []

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
        travel_length_list.append(res[7])

    logging.info("Calculation finished, getting results.")
    for i in range(len(z_delta_list)):
        z_delta_array = np.maximum(z_delta_array, z_delta_list[i])
        flux_array = np.maximum(flux_array, flux_list[i])
        count_array += cc_list[i]
        z_delta_sum += z_delta_sum_list[i]
        backcalc = np.maximum(backcalc, backcalc_list[i])
        fp_travelangle_array = np.maximum(fp_travelangle_array, fp_ta_list[i])
        sl_travelangle_array = np.maximum(sl_travelangle_array, sl_ta_list[i])
        travel_length_array = np.maximum(travel_length_array, travel_length_list[i])

    # Save Calculated tiles
    np.save(tempDir / ("res_z_delta_%s_%s" % (optTuple[0], optTuple[1])), z_delta_array)
    np.save(tempDir / ("res_z_delta_sum_%s_%s" % (optTuple[0], optTuple[1])), z_delta_sum)
    np.save(tempDir / ("res_flux_%s_%s" % (optTuple[0], optTuple[1])), flux_array)
    np.save(tempDir / ("res_count_%s_%s" % (optTuple[0], optTuple[1])), count_array)
    np.save(tempDir / ("res_fp_%s_%s" % (optTuple[0], optTuple[1])), fp_travelangle_array)
    np.save(tempDir / ("res_sl_%s_%s" % (optTuple[0], optTuple[1])), sl_travelangle_array)
    np.save(tempDir / ("res_travel_length_%s_%s" % (optTuple[0], optTuple[1])), travel_length_array)
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

    # check if there's enough RAM available (default value set to 5%)
    # if not, wait for 30 secs and check again
    # should prevent the occurence of broken pipe errors or similar issues related
    # to RAM overflow
    handleMemoryAvailability()

    dem = args[0]
    infra = args[1]
    release = args[2]
    alpha = args[3]
    exp = args[4]
    flux_threshold = args[5]
    max_z_delta = args[6]
    nodata = args[7]
    cellsize = args[8]
    infraBool = args[9]
    forestBool = args[10]

    if forestBool:
        forestArray = args[11]
        forestParams = args[12]

    z_delta_array = np.zeros_like(dem, dtype=np.float32)
    z_delta_sum = np.zeros_like(dem, dtype=np.float32)
    flux_array = np.zeros_like(dem, dtype=np.float32)
    count_array = np.zeros_like(dem, dtype=np.int32)

    fp_travelangle_array = np.zeros_like(dem, dtype=np.float32)  # fp = Flow Path
    sl_travelangle_array = np.zeros_like(dem, dtype=np.float32) * 90  # sl = Straight Line

    travel_length_array = np.zeros_like(dem, dtype=np.float32)

    # NOTE-TODO maybe also include a switch for INFRA (like Forest) and not implicitly always use an empty infra array ?
    backcalc = np.zeros_like(dem, dtype=np.int32)

    if infraBool:
        back_list = []

    # Core
    # start = datetime.now().replace(microsecond=0)
    # NOTE-TODO: row_list, col_list are tuples - rethink variable naming
    row_list, col_list = get_start_idx(dem, release)

    startcell_idx = 0
    while startcell_idx < len(row_list):

        cell_list = []
        row_idx = row_list[startcell_idx]
        col_idx = col_list[startcell_idx]
        dem_ng = dem[row_idx - 1: row_idx + 2, col_idx - 1: col_idx + 2]  # neighbourhood DEM

        if (nodata in dem_ng) or np.size(dem_ng) < 9:
            startcell_idx += 1
            continue

        if forestBool:
            startcell = Cell(
                row_idx, col_idx,
                dem_ng, cellsize,
                1, 0, None,
                alpha, exp, flux_threshold, max_z_delta,
                startcell=True,
                FSI=forestArray[row_idx, col_idx],
                forestParams=forestParams,
            )
            # If this is a startcell just give a Bool to startcell otherwise the object startcell
        else:
            startcell = Cell(
                row_idx, col_idx,
                dem_ng, cellsize,
                1, 0, None,
                alpha, exp, flux_threshold, max_z_delta,
                startcell=True
            )
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
                dem_ng = dem[row[k] - 1: row[k] + 2, col[k] - 1: col[k] + 2]  # neighbourhood DEM

                # This bit handles edge cases and noData-values in the DEM!! this is an important piece of code, since
                # no-data handling is expected (by some users/applications) to behave like here:
                # i.e. if nodata in the 3x3 neighbourhood --> no calculation
                if (nodata in dem_ng) or np.size(dem_ng) < 9:
                    continue

                if forestBool:
                    cell_list.append(
                        Cell(row[k], col[k],
                             dem_ng, cellsize,
                             flux[k], z_delta[k],
                             cell,
                             alpha, exp, flux_threshold, max_z_delta,
                             startcell,
                             FSI=forestArray[row[k], col[k]],
                             forestParams=forestParams,
                             )
                    )
                else:
                    cell_list.append(
                        Cell(row[k], col[k],
                             dem_ng, cellsize,
                             flux[k], z_delta[k],
                             cell,
                             alpha, exp, flux_threshold, max_z_delta,
                             startcell,
                             )
                    )

            z_delta_array[cell.rowindex, cell.colindex] = max(z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
            flux_array[cell.rowindex, cell.colindex] = max(flux_array[cell.rowindex, cell.colindex], cell.flux)
            count_array[cell.rowindex, cell.colindex] += int(1)
            z_delta_sum[cell.rowindex, cell.colindex] += cell.z_delta
            fp_travelangle_array[cell.rowindex, cell.colindex] = max(fp_travelangle_array[cell.rowindex, cell.colindex],
                                                                     cell.max_gamma)
            sl_travelangle_array[cell.rowindex, cell.colindex] = max(sl_travelangle_array[cell.rowindex, cell.colindex],
                                                                     cell.sl_gamma)
            travel_length_array[cell.rowindex, cell.colindex] = max(travel_length_array[cell.rowindex, cell.colindex], 
                                                                    cell.min_distance)

            # Backcalculation
            if infraBool:
                # NOTE-TODO:
                # just store 'affected' infrastructure cells (row,index-colindex) here and
                # do backcalculation after the path calculation is finished
                if infra[cell.rowindex, cell.colindex] > 0:
                    # backlist = []
                    backList = back_calculation(cell)

                    for bCell in backList:
                        backcalc[bCell.rowindex, bCell.colindex] = max(backcalc[bCell.rowindex, bCell.colindex],
                                                                       infra[cell.rowindex, cell.colindex])

        if infraBool:
            release[z_delta_array > 0] = 0
            # Check if i hit a release Cell, if so set it to zero and get again the indexes of release cells
            row_list, col_list = get_start_idx(dem, release)

        del cell_list

        startcell_idx += 1
    # end = datetime.now().replace(microsecond=0)
    gc.collect()
    return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array


def enoughMemoryAvailable(limit=0.05):
    """simple function to monitor memory(RAM) availability during parallel processing
    of calculation() inside run(). utilizing psutil

    Parameters
    -----------
    limit: float (between 0 and 1) - default at 0.05 (i.e. 5%)

    Returns
    -----------
    'True' if more than the defined memory-limit is still available
    'False' if less than the defined memory-limit is available
    """

    log = logging.getLogger(__name__)
    availableMemory = psutil.virtual_memory().available / psutil.virtual_memory().total

    if availableMemory >= limit:
        # log.info('RAM availability o.k. -- %.2f %% of %.2f GiB'%
        # (availableMemory*100,psutil.virtual_memory().total/(1024.**3)))
        return True
    else:
        log.info(
            "RAM availability at limit -- %.2f %% of %.2f GiB - maybe recheck multiProcessing/Tiling settings"
            % (availableMemory * 100, psutil.virtual_memory().total / (1024.0**3))
        )
        return False


def calculateMultiProcessingOptions(nRel, nCPU, procPerCPU=1, maxChunks=500, chunkSize=50):
    """compute required options for multiprocessing of calulation() function inside run() and accompanied splitting of
    release cells into chunks in split_release().append
    The general idea is to make good use of available CPU resources to speed up calculations while not getting into
    trouble with RAM issues ...

    NOTE: this is still a quick'n'dirty hack, it might make sense to have a more sophisticated approach for optimization
          of CPU and RAM resource usage during multiprocessing depending on e.g.:
              - size of the numpy arrays that are processed (depending on tileSize and rasterResolution)
              - density of release areas in the tile
              - total available RAM and CPUs on the machine
              - (other com4FlowPy Parameterization)

    Parameters
    -----------
    nRel: int - number of release Pixels inside the tile (i.e. all cells/pixels with values >=1 in 'release')
    nCPU: int - number of available CPUs (as defined in the .ini files)
    procPerCPU: int - number of processes to be spawned per CPU (default = 1) - might be set higher for increased
                      performance

    maxChunks: int - hard limit to the maximum number of chunks that is used --> a larger number of chunks will very
               probably increase performance in terms of maximising CPU workload (especially with large numbers of
               nCPU) but also cause higher RAM consumption (in the current multiprocessing implementation)
    chunkSize: int - default number of release pixels per chunk in cases where the chunk-size is not constrained by
               nCPU*procPerCPU or maxChunks
    Returns
    -----------
    nChunks: int - the number of chunks into which the release layer/array is split for multiprocessing
    nProcesses: int - the number of processes used in Pool.map() inside run()
    """
    nProcesses = int(nCPU * procPerCPU)
    # check if release is empty - if so, there's no reason to split
    if nRel == 0:
        nChunks = 1
    # if the number of release cells is smaller/equal than/to the number of processes
    # each single release cell is assigned to a different process
    elif nRel <= nProcesses:
        nChunks = nRel
    # if there are more release cells than number of Processes available (this is the main case!)
    # then either divide release cells equally to the available processes - however limit the size of single chunks
    # to chunkSize if possible ...
    else:
        _nChunks = max(nProcesses, int(nRel / chunkSize))
        nChunks = min(_nChunks, maxChunks)

    return int(nProcesses), int(nChunks)


def handleMemoryAvailability(recheckInterval=30):
    """function is called at the start of each subProcess for parallel processing to check if enough memory is available
    and handle the situation if not

    NOTE: currently only time.sleep() is called to delay the subprocess for a defined time and then re-check
          memory availability.
          other possible options:
              - log message and abort model run?
    NOTE: The implementation with time.sleep() can cause an "infinite loop" in time-sleep if for some
          reason memory is not freed after a sensible amount of time.altzone
          Memory consumption is dependend on tile-sizes and number of Chunks/tile ...
    Parameters
    -----------
    recheckInterval: int - delay time for the process after which memory availability is re-checked
    """
    while not enoughMemoryAvailable():
        time.sleep(recheckInterval)
