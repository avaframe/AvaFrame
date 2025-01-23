#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Calculation functions (raster level)
"""

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

    Parameters
    -----------
    dem: numpy array
        Digital Elevation Model to gain information about altitude
    release: numpy array
        The release layer, release pixels need int value > 0

    Returns
    -----------
    row_list: list
        Row indices of release pixels sorted by altitude
    col_list: list
        Column indices of release pixels sorted by altitude
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
    release: np.array
        a binary 0|1 array with release pixels designated by '1'
    pieces:  int
        number of chunck in which the release layer should be split

    Returns
    -----------
    release_list:   list
        contains the tiles(arrays) [array0, array1, ..]
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

    Parameters
    -----------
    back_cell:  class-object cell
        cell that hit a Infrastructure

    Returns
    -----------
    back_list: list
        List of cells that are on the way to the start cell TODO:Maybe change it to array like DEM?
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
    using multiprocessing across multiple CPUs and for saving results for processed tile to temp folder
    (modelPaths["tempFolder"])

    Parameters
    -----------
    optTuple: tuple
        with all necessary model information:

        - optTuple[0] (int) - i index of the processed tile (for loading correct data for tiles)
        - optTuple[1] (int) - j index of the processed tile (for loading correct data for tiles)
        - optTuple[2] (dict) - containing modelParameters
        - optTuple[3] (dict) - containing modelPaths
        - optTuple[4] (dict) - containing rasterAttributes
        - optTuple[5] (dict) - containing forestParameters
        - optTuple[6] (dict) - containing MPOptions

    """

    log = logging.getLogger(__name__)

    # Flow-Py parameters
    modelParameters = optTuple[2]
    infraBool = optTuple[2]["infraBool"]
    forestBool = optTuple[2]["forestBool"]
    forestInteraction = optTuple[2]["forestInteraction"]
    varUmaxBool = optTuple[2]["varUmaxBool"]
    varAlphaBool = optTuple[2]["varAlphaBool"]
    varExponentBool = optTuple[2]["varExponentBool"]

    # Temp-Dir (all input files are located here and results are written back in here)
    tempDir = optTuple[3]["tempDir"]

    # raster-layer Attributes
    rasterAttributes = optTuple[4]

    MPOptions = optTuple[6]  # CPU, Multiprocessing options ...

    dem = np.load(tempDir / ("dem_%s_%s.npy" % (optTuple[0], optTuple[1])))
    release = np.load(tempDir / ("init_%s_%s.npy" % (optTuple[0], optTuple[1])))
    if infraBool:
        infra = np.load(tempDir / ("infra_%s_%s.npy" % (optTuple[0], optTuple[1])))
    else:
        infra = None

    # if forestBool == 'True'
    # --> load forestFile
    # --> read parametersOfForestExtension
    # NOTE-TODO: this is a quick work-around to simply include forest information - should probably be handled more
    # elegantly/explicitly AND error handling should probably be included
    if forestBool:
        forestArray = np.load(tempDir / ("forest_%s_%s.npy" % (optTuple[0], optTuple[1])))
        forestParams = optTuple[5]
        forestParams["forestInteraction"] = forestInteraction
    else:
        forestParams = None
        forestArray = None
    modelParameters["forestParams"] = forestParams

    if varUmaxBool:
        varUmaxArray = np.load(tempDir / ("varUmax_%s_%s.npy" % (optTuple[0], optTuple[1])))
        if optTuple[2]["varUmaxType"].lower() == 'umax':
            varUmaxArray[varUmaxArray > 0] = varUmaxArray[varUmaxArray > 0] ** 2 / 2 / 9.81
        elif optTuple[2]["varUmaxType"].lower() != 'zdeltalim':
            log.error("PLease provide the type of the uMax Limit: 'uMax' (in m/s) or zDeltaMax (in m)!")
    else:
        varUmaxArray = None

    if varAlphaBool:
        varAlphaArray = np.load(tempDir / ("varAlpha_%s_%s.npy" % (optTuple[0], optTuple[1])))
    else:
        varAlphaArray = None

    if varExponentBool:
        varExponentArray = np.load(tempDir / ("varExponent_%s_%s.npy" % (optTuple[0], optTuple[1])))
    else:
        varExponentArray = None

    modelData = {
        "dem": dem,
        "infra": infra,
        "forestArray": forestArray,
        'varUmaxArray': varUmaxArray,
        'varAlphaArray': varAlphaArray,
        'varExponentArray': varExponentArray,
    }
    del infra, forestArray, varUmaxArray, varAlphaArray, varExponentArray
    outputFileList = set(optTuple[3]['outputFileList'])
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

    with Pool(processes=nProcesses) as pool:
        results = pool.map(
            calculation,
            [
                [   # TODO: write in dicts:
                    modelParameters,
                    modelData,
                    rasterAttributes,
                    outputFileList,
                    release_sub,
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

    computeSaveTileResults(outputFileList, results, dem, tempDir, optTuple[0], optTuple[1])

def calculation(args):
    """This is the core function where all the data handling and calculation is
    done.

    Parameters
    -----------
    args: list
        contains the following model input data:

        - args[0] (dict) - contains all parameters for the model
        - args[1] (dict) - contains necessary input data
        - args[2] (rasterAttributes) - contains raster attributes od the input data
        - args[3] (list) - contains the names for the outputs
        - args[4] (np.array) - cells with value 1 are PRAs

    Returns
    -----------
    outputArray: dict
        contains the numpy arrays for the output
    """
    log = logging.getLogger(__name__)
    # check if there's enough RAM available (default value set to 5%)
    # if not, wait for 30 secs and check again
    # should prevent the occurence of broken pipe errors or similar issues related
    # to RAM overflow
    handleMemoryAvailability()
    modelParameters, modelData, rasterAttributes, outputFileList, release_sub = args
    dem = modelData["dem"]
    infra = modelData["infra"]
    release = release_sub
    alpha = float(modelParameters["alpha"])
    exp = float(modelParameters["exp"])
    flux_threshold = float(modelParameters["flux_threshold"])
    max_z_delta = float(modelParameters["max_z"])
    nodata = rasterAttributes["nodata"]
    cellsize = rasterAttributes["cellsize"]
    infraBool = modelParameters["infraBool"]
    forestBool = modelParameters["forestBool"]
    varUmaxBool = modelParameters["varUmaxBool"]
    varUmaxArray = modelData["varUmaxArray"]
    varAlphaBool = modelParameters["varAlphaBool"]
    varAlphaArray = modelData["varAlphaArray"]
    varExponentBool = modelParameters["varExponentBool"]
    varExponentArray = modelData["varExponentArray"]
    fluxDistOldVersionBool = modelParameters["fluxDistOldVersionBool"]

    if forestBool:
        forestArray = modelData["forestArray"]
        forestParams = modelParameters["forestParams"]
        forestInteraction = forestParams["forestInteraction"]
    else:
        forestInteraction = False
        forestArray = None
        forestParams = None

    if 'zDelta' in outputFileList:
        zDeltaArray = np.zeros_like(dem, dtype=np.float32)
    else:
        zDeltaArray = None
    if 'zDeltaSum' in outputFileList:
        zDeltaSumArray = np.zeros_like(dem, dtype=np.float32)
        zDeltaPathList = []
    else:
        zDeltaSumArray = None
    if 'routFluxSum' in outputFileList:
        routFluxSumArray = np.zeros_like(dem, dtype=np.float32)
    else:
        routFluxSumArray = None
    if 'depFluxSum' in outputFileList:
        depFluxSumArray = np.zeros_like(dem, dtype=np.float32)
    else:
        depFluxSumArray = None
    if 'flux' in outputFileList:
        fluxArray = np.zeros_like(dem, dtype=np.float32)
    else:
        fluxArray = None
    if 'cellCounts' in outputFileList:
        cellCountsArray = np.zeros_like(dem, dtype=np.int32)
    else:
        cellCountsArray = None
    if 'fpTravelAngle' in outputFileList:
        fpTravelAngleArray = np.zeros_like(dem, dtype=np.float32)  # fp = Flow Path
    else:
        fpTravelAngleArray = None
    if 'slTravelAngle' in outputFileList:
        slTravelAngleArray = np.zeros_like(dem, dtype=np.float32)  # sl = Straight Line
    else:
        slTravelAngleArray = None
    if 'travelLength' in outputFileList:
        travelLengthArray = np.zeros_like(dem, dtype=np.float32)
    else:
        travelLengthArray = None
    if 'forestInteraction' in outputFileList:
        if forestInteraction:
            forestIntArray = np.ones_like(dem, dtype=np.float32) * -9999
        else:
            log.error("forestInteraction must set to True, when 'forestInteraction' is in outputs (.ini - File)")
            raise ValueError("forestInteraction must set to True, when 'forestInteraction' is in outputs (.ini - File)")
    else:
        forestIntArray = None

    # NOTE-TODO maybe also include a switch for INFRA (like Forest) and not implicitly always use an empty infra array ?
    if 'backcalc' in outputFileList:
        if infraBool:
            backcalc = np.zeros_like(dem, dtype=np.int32)
        else:
            log.error("infra must set to True, when 'backcalc' is in outputs (.ini - File)")
            raise ValueError("infra must set to True, when 'backcalc' is in outputs (.ini - File)")
    else:
        backcalc = None
    if infraBool:
        back_list = []

    # Core
    # start = datetime.now().replace(microsecond=0)
    # NOTE-TODO: row_list, col_list are tuples - rethink variable naming
    row_list, col_list = get_start_idx(dem, release)

    startcell_idx = 0
    while startcell_idx < len(row_list):

        processedCells = {}  # dictionary of cells that have been processed already
        zDeltaPathArray = np.zeros_like(dem, dtype=np.float32)
        cell_list = []
        row_idx = row_list[startcell_idx]
        col_idx = col_list[startcell_idx]
        dem_ng = dem[row_idx - 1: row_idx + 2, col_idx - 1: col_idx + 2]  # neighbourhood DEM
        if varUmaxBool and varUmaxArray is not None:
            if varUmaxArray[row_idx, col_idx] > 0 and varUmaxArray[row_idx, col_idx] <= 8848:
                max_z_delta = varUmaxArray[row_idx, col_idx]
        if varAlphaBool and varAlphaArray is not None:
            if varAlphaArray[row_idx, col_idx] > 0 and varAlphaArray[row_idx, col_idx] <= 90:
                alpha = varAlphaArray[row_idx, col_idx]
        if varExponentBool and varExponentArray is not None:
            if varExponentArray[row_idx, col_idx] > 0:
                exp = varExponentArray[row_idx, col_idx]

        if (nodata in dem_ng) or np.size(dem_ng) < 9:
            startcell_idx += 1
            continue

        startcell = Cell(
            row_idx, col_idx,
            dem_ng, cellsize,
            1, 0, None,
            alpha, exp, flux_threshold, max_z_delta,
            startcell=True, fluxDistOldVersionBool=fluxDistOldVersionBool,
            FSI=forestArray[row_idx, col_idx] if isinstance(forestArray, np.ndarray) else None,
            forestParams=forestParams,
        )

        # dictionary of all the cells that have been processed and the number of times the cell has been visited
        processedCells[(startcell.rowindex, startcell.colindex)] = 1
        # list of flowClass.Cell() Objects that is contains the "path" for each release-cell
        cell_list.append(startcell)

        for idx, cell in enumerate(cell_list):

            row, col, flux, z_delta = cell.calc_distribution()

            if len(flux) > 0:
                # mass, row, col  = list(zip(*sorted(zip( mass, row, col), reverse=False)))
                z_delta, flux, row, col = list(zip(*sorted(zip(z_delta, flux, row, col), reverse=False)))
                # Sort this lists by elh, to start with the highest cell

            # check if cell already exists
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

                # if the current child cell is already in processedCells
                # just add +1 to the visit-counter, else add it to the
                # processedCells dictionary with visit-count = 1
                if (row[k], col[k]) in processedCells:
                    processedCells[(row[k], col[k])] += 1
                else:
                    processedCells[(row[k], col[k])] = 1

                cell_list.append(Cell(
                            row[k], col[k],
                            dem_ng, cellsize,
                            flux[k], z_delta[k],
                            cell,
                            alpha, exp, flux_threshold, max_z_delta,
                            startcell, fluxDistOldVersionBool=fluxDistOldVersionBool,
                            FSI=forestArray[row[k], col[k]] if isinstance(forestArray, np.ndarray) else None,
                            forestParams=forestParams,
                                     ))
            if 'zDelta' in outputFileList:
                zDeltaArray[cell.rowindex, cell.colindex] = max(zDeltaArray[cell.rowindex, cell.colindex], cell.z_delta)
            if 'flux' in outputFileList:
                fluxArray[cell.rowindex, cell.colindex] = max(fluxArray[cell.rowindex, cell.colindex], cell.flux)
            if 'routFluxSum' in outputFileList:
                routFluxSumArray[cell.rowindex, cell.colindex] += cell.flux
            if 'depFluxSum' in outputFileList:
                depFluxSumArray[cell.rowindex, cell.colindex] += cell.fluxDep
            if 'zDeltaSum' in outputFileList:
                zDeltaPathArray[cell.rowindex, cell.colindex] = max(zDeltaPathArray[cell.rowindex, cell.colindex], cell.z_delta)
            if 'fpTravelAngle' in outputFileList:
                if (not cell.is_start):
                    cell.calc_fp_travelangle()
                fpTravelAngleArray[cell.rowindex, cell.colindex] = max(fpTravelAngleArray[cell.rowindex, cell.colindex],
                                                                   cell.max_gamma)
            if 'slTravelAngle' in outputFileList:
                if (not cell.is_start):
                    cell.calc_sl_travelangle()
                slTravelAngleArray[cell.rowindex, cell.colindex] = max(slTravelAngleArray[cell.rowindex, cell.colindex],
                                                                   cell.sl_gamma)
            if 'travelLength' in outputFileList:
                if (not cell.is_start):
                    cell.calc_fp_travelangle()
                travelLengthArray[cell.rowindex, cell.colindex] = max(travelLengthArray[cell.rowindex, cell.colindex],
                                                                  cell.min_distance)
            if 'cellCounts' in outputFileList:
                if processedCells[(cell.rowindex, cell.colindex)] == 1:
                    cellCountsArray[cell.rowindex, cell.colindex] += int(1)
            if 'forestInteraction' in outputFileList:
                if forestIntArray[cell.rowindex, cell.colindex] >= 0 and cell.forestIntCount >= 0:
                    forestIntArray[cell.rowindex, cell.colindex] = min(forestIntArray[cell.rowindex, cell.colindex],
                                                                       cell.forestIntCount)
                else:
                    forestIntArray[cell.rowindex, cell.colindex] = max(forestIntArray[cell.rowindex, cell.colindex],
                                                                       cell.forestIntCount)

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
            release[zDeltaArray > 0] = 0
            # Check if i hit a release Cell, if so set it to zero and get again the indexes of release cells
            row_list, col_list = get_start_idx(dem, release)
        if 'zDeltaSum' in outputFileList:
            zDeltaPathList.append(zDeltaPathArray)
            del zDeltaPathArray
        del cell_list, processedCells

        startcell_idx += 1
    # end = datetime.now().replace(microsecond=0)
    if 'zDeltaSum' in outputFileList:
        for zDeltaPathArray in zDeltaPathList:
            zDeltaSumArray += zDeltaPathArray
    gc.collect()
    outputArray = {}
    _arrayNames = ['zDelta', 'flux', 'cellCounts', 'zDeltaSum', 'backcalc', 'fpTravelAngle', 'slTravelAngle', \
                  'travelLength', 'routFluxSum', 'depFluxSum', 'forestInteraction']
    _output = zDeltaArray, fluxArray, cellCountsArray, zDeltaSumArray, backcalc, fpTravelAngleArray, slTravelAngleArray, travelLengthArray, routFluxSumArray, depFluxSumArray, forestIntArray
    for name, array in zip(_arrayNames, _output):
        outputArray[name] = array
    del _output, _arrayNames
    return outputArray


def enoughMemoryAvailable(limit=0.05):
    """simple function to monitor memory(RAM) availability during parallel processing
    of calculation() inside run(). utilizing psutil

    Parameters
    -----------
    limit: float
        available RAM memory limit (between 0 and 1) - default at 0.05 (i.e. 5%)

    Returns
    -----------
    bool
        'True' if more than the defined memory-limit is still available;
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
    nRel: int
        number of release Pixels inside the tile (i.e. all cells/pixels with values >=1 in 'release')
    nCPU: int
        number of available CPUs (as defined in the .ini files)
    procPerCPU: int
        number of processes to be spawned per CPU (default = 1) - might be set higher for increased performance
    maxChunks: int
        hard limit to the maximum number of chunks that is used --> a larger number of chunks will very
        probably increase performance in terms of maximising CPU workload (especially with large numbers of
        nCPU) but also cause higher RAM consumption (in the current multiprocessing implementation)
    chunkSize: int
        default number of release pixels per chunk in cases where the chunk-size is not constrained by
        nCPU*procPerCPU or maxChunks

    Returns
    -----------
    nChunks: int
        the number of chunks into which the release layer/array is split for multiprocessing
    nProcesses: int
        the number of processes used in Pool.map() inside run()
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
    recheckInterval: int
        delay time (in seconds) for the process after which memory availability is re-checked
    """
    while not enoughMemoryAvailable():
        time.sleep(recheckInterval)


def computeSaveTileResults(outputFileList, results, dem, tempDir, indI, indJ):
    """function computes the output rasters for the processed tile and saves the rasters
    in the temporary fodler

    Parameters
    -----------
    outputFileList: list
        contains the names of the outputs
    results: dict
        contains output rasters computed by every CPU core
    dem: numpy array
        digital elevation model of the processed tile
    tempDir: str
        path to temporary directory
    indI: int
        i index of the processed tile (for loading correct data for tiles)
    indJ: int
        j index of the processed tile (for loading correct data for tiles)

    """

    if 'zDelta' in outputFileList:
        zDeltaArray = np.zeros_like(dem, dtype=np.float32)
    if 'flux' in outputFileList:
        fluxArray = np.zeros_like(dem, dtype=np.float32)
    if 'cellCounts'  in outputFileList:
        cellCountsArray = np.zeros_like(dem, dtype=np.int32)
    if 'zDeltaSum'  in outputFileList:
        zDeltaSumArray = np.zeros_like(dem, dtype=np.float32)
    if 'routFluxSum'  in outputFileList:
        routFluxSumArray = np.zeros_like(dem, dtype=np.float32)
    if 'depFluxSum'  in outputFileList:
        depFluxSumArray = np.zeros_like(dem, dtype=np.float32)
    if 'backcalc'  in outputFileList:
        backcalc = np.zeros_like(dem, dtype=np.int32)
    if 'fpTravelAngle'  in outputFileList:
        fpTravelAngleArray = np.zeros_like(dem, dtype=np.float32)
    if 'slTravelAngle'  in outputFileList:
        slTravelAngleArray = np.zeros_like(dem, dtype=np.float32)
    if 'travelLength'  in outputFileList:
        travelLengthArray = np.zeros_like(dem, dtype=np.float32)
    if 'forestInteraction'  in outputFileList:
        forestIntArray = np.ones_like(dem, dtype=np.int32) * -9999

    for i in range(len(results)):
        res = results[i]
        if 'zDelta' in outputFileList:
            zDeltaArray = np.maximum(zDeltaArray, res['zDelta'])
        if 'flux' in outputFileList:
            fluxArray = np.maximum(fluxArray, res['flux'])
        if 'cellCounts'  in outputFileList:
            cellCountsArray += res['cellCounts']
        if 'zDeltaSum'  in outputFileList:
            zDeltaSumArray += res['zDeltaSum']
        if 'backcalc'  in outputFileList:
            backcalc = np.maximum(backcalc, res['backcalc'])
        if 'fpTravelAngle'  in outputFileList:
            fpTravelAngleArray = np.maximum(fpTravelAngleArray, res['fpTravelAngle'])
        if 'slTravelAngle'  in outputFileList:
            slTravelAngleArray = np.maximum(slTravelAngleArray, res['slTravelAngle'])
        if 'travelLength'  in outputFileList:
            travelLengthArray = np.maximum(travelLengthArray, res['travelLength'])
        if 'routFluxSum'  in outputFileList:
            routFluxSumArray += res['routFluxSum']
        if 'depFluxSum'  in outputFileList:
            depFluxSumArray += res['depFluxSum']
        if 'forestInteraction'  in outputFileList:
            forestIntArray = np.where((forestIntArray >= 0) & (res['forestInteraction'] >= 0),
                                    np.minimum(forestIntArray, res['forestInteraction']),
                                    np.maximum(forestIntArray, res['forestInteraction']))

    logging.info("Calculation finished, getting results.")
    """
    for i in range(len(zDeltaList)):
        zDeltaArray = np.maximum(zDeltaArray, zDeltaList[i])
        fluxArray = np.maximum(fluxArray, fluxList[i])
        cellCountsArray += ccList[i]
        zDeltaSumArray += zDeltaSumList[i]
        routFluxSumArray += routFluxSumList[i]
        depFluxSumArray += depFluxSumList[i]
        backcalc = np.maximum(backcalc, backcalcList[i])
        fpTravelAngleArray = np.maximum(fpTravelAngleArray, fpTravelAngleList[i])
        slTravelAngleArray = np.maximum(slTravelAngleArray, slTravelAngleList[i])
        travelLengthArray = np.maximum(travelLengthArray, travelLengthList[i])
        if forestInteraction:
            forestIntArray = np.where((forestIntArray >= 0) & (forestIntList[i] >= 0),
                                    np.minimum(forestIntArray, forestIntList[i]),
                                    np.maximum(forestIntArray, forestIntList[i]))
    """

    # Save Calculated tiles
    if 'zDelta' in outputFileList:
        np.save(tempDir / ("res_z_delta_%s_%s" % (indI, indJ)), zDeltaArray)
    if 'zDeltaSum'  in outputFileList:
        np.save(tempDir / ("res_z_delta_sum_%s_%s" % (indI, indJ)), zDeltaSumArray)
    if 'routFluxSum'  in outputFileList:
        np.save(tempDir / ("res_rout_flux_sum_%s_%s" % (indI, indJ)), routFluxSumArray)
    if 'depFluxSum'  in outputFileList:
        np.save(tempDir / ("res_dep_flux_sum_%s_%s" % (indI, indJ)), depFluxSumArray)
    if 'flux' in outputFileList:
        np.save(tempDir / ("res_flux_%s_%s" % (indI, indJ)), fluxArray)
    if 'cellCounts'  in outputFileList:
        np.save(tempDir / ("res_count_%s_%s" % (indI, indJ)), cellCountsArray)
    if 'fpTravelAngle'  in outputFileList:
        np.save(tempDir / ("res_fp_%s_%s" % (indI, indJ)), fpTravelAngleArray)
    if 'slTravelAngle'  in outputFileList:
        np.save(tempDir / ("res_sl_%s_%s" % (indI, indJ)), slTravelAngleArray)
    if 'travelLength'  in outputFileList:
        np.save(tempDir / ("res_travel_length_%s_%s" % (indI, indJ)), travelLengthArray)
    if 'backcalc'  in outputFileList:
        np.save(tempDir / ("res_backcalc_%s_%s" % (indI, indJ)), backcalc)
    if 'forestInteraction'  in outputFileList:
        np.save(tempDir / ("res_forestInt_%s_%s" % (indI, indJ)), forestIntArray)