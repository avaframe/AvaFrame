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
#Paula
from avaframe.com4FlowPy.flowPath import Path
import matplotlib.pyplot as plt
#end paula

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

def path_calc_analysis(path_list):
    # PAULA
    path_travel_lengths = np.empty(0)
    path_altitude = np.empty(0)
    path_z_delta_sum = np.empty(0)
    path_z_delta_area_1 = np.empty(0)
    path_z_delta_area_mean = np.empty(0)
    path_z_delta_max = np.empty(0)
    path_area = np.empty(0)

    for path in path_list: # calculate for every path
        path_travel_lengths = np.append(path_travel_lengths, max(path.s_coE)) # travel length of the whole path 
        path_altitude = np.append(path_altitude, max(path.altitude_coE)-min(path.altitude_coE)) #drop height
        path_z_delta_sum = np.append(path_z_delta_sum, sum(path.z_delta_coE)) #sum of z_delta
        path_z_delta_max = np.append(path_z_delta_max, max(path.z_delta_coE)) #max of z_delta
        path_area = np.append(path_area, path.path_area)

        #idea for calculating area between zdelta and terrain
        distance = [0]
        z_delta_mean = np.empty(0)

        for i,s in enumerate(path.s_coE):
            if i < (len(path.s_coE)-1):
                dist = path.s_coE[i+1]-s
                distance.append(dist)

                z_mean = np.mean([path.z_delta_coE[i+1],path.z_delta_coE[i]])
                z_delta_mean = np.append(z_delta_mean, z_mean)

        # Distanz * z_delta an jeweiligem Ort
        path_z_delta_area_1 = np.append(path_z_delta_area_1, np.sum(np.array(distance) * np.array(path.z_delta_coE)))
        # Distanz * mitterlwert von z_delta zwischen zwei punkten
        path_z_delta_area_mean = np.append(path_z_delta_area_mean, np.sum(np.array(distance[1:]) * np.array(z_delta_mean)))
                    
        #print(f'Area between z_delta and terrain: calculating with local z_delta: {path_z_delta_area_1},
        #calculated with mean: {path_z_delta_area_mean}')

        ''' EXPENSIVE!
        fig = path.plot_pathanaylsis() #this line takes lot of time!!
        fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/plot_pathlist_col{path_test.start_col},row{path_test.start_row}.png')
        plt.close(fig)
        '''
    return path_travel_lengths, path_altitude, path_z_delta_sum, path_z_delta_area_mean, path_area, path_z_delta_max

def path_plot_analysis(path_analysis_list):
    # get path variables from path_analysis_list
    path_travel_lengths = []
    path_altitude = []
    path_z_delta_sum = []
    path_z_delta_area_mean = []    
    path_area = []
    path_z_delta_max = []


    for var_processes in path_analysis_list:
        path_travel_lengths.extend(var_processes[0])
        path_altitude.extend(var_processes[1])
        path_z_delta_sum.extend(var_processes[2])
        path_z_delta_area_mean.extend(var_processes[3])
        path_area.extend(var_processes[4])
        path_z_delta_max.extend(var_processes[5])

    # Histograms
    fig,ax = plt.subplots()
    ax.hist(path_travel_lengths)
    plt.xlabel('max. travel length of coE path [m]')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/hist_travel_length.png')
    # HARDCODED!!!
    plt.close(fig)

    fig,ax = plt.subplots()
    ax.hist(path_altitude)
    plt.xlabel('drop height of coE path [m]')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/hist_altitude.png')
    plt.close(fig)

    fig,ax = plt.subplots()
    ax.hist(path_z_delta_sum)
    plt.xlabel('sum of Z^{\delta} along coE path')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/hist_z_delta_sum.png')
    plt.close(fig)

    fig,ax = plt.subplots()
    ax.hist(path_z_delta_max)
    plt.xlabel('maximum of Z_delta along coE path')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/hist_z_delta_max.png')
    plt.close(fig)

    fig,ax = plt.subplots()
    ax.hist(path_z_delta_area_mean)
    plt.xlabel('area between Z^{\delta} and topography')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/hist_z_delta_area.png')
    plt.close(fig)

    fig,ax = plt.subplots()
    ax.hist(path_area)
    plt.xlabel('path area [km²]')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/path_area.png')
    plt.close(fig)

    # BOxplots
    fig,ax = plt.subplots()
    ax.boxplot(path_travel_lengths)
    plt.ylabel('max. travel length of coE path [m]')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/boxpl_travel_length.png')
    # HARDCODED!!!
    plt.close(fig)   

    #Scatterplot
    fig,ax = plt.subplots(1,2)
    ax[0].scatter(path_altitude, path_z_delta_area_mean)
    ax[0].set(xlabel = 'Drop height [m]')
    ax[0].set(ylabel = 'area between Z^{\delta} and topography')

    ax[1].scatter(path_z_delta_max, path_z_delta_area_mean)
    ax[1].set(xlabel = 'max z_delta')
    ax[1].set(ylabel = 'area between Z^{\delta} and topography')
    fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/scatter_area.png')
    # HARDCODED!!!
    plt.close(fig)   
    


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
    print('pool closed')
    z_delta_array = np.zeros_like(dem)
    flux_array = np.zeros_like(dem)
    count_array = np.zeros_like(dem)
    z_delta_sum = np.zeros_like(dem)
    backcalc = np.zeros_like(dem)
    fp_travelangle_array = np.zeros_like(dem)
    sl_travelangle_array = np.zeros_like(dem)
    #Chris
    travel_length_array = np.zeros_like(dem)
    #ende chris
    #Paula
    flow_energy_array = np.zeros_like(dem)
    #ende paula

    z_delta_list = []
    flux_list = []
    cc_list = []
    z_delta_sum_list = []
    backcalc_list = []
    fp_ta_list = []
    sl_ta_list = []
    #Chris
    travel_length_list = []
    #chris ende
    #Paula
    flow_energy_list = []
    path_analysis_list = []
    #ende paula

    print('start for loop for results')
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
        #Chris
        travel_length_list.append(res[7])
        #chris ende
        #Paula
        flow_energy_list.append(res[8])
        path_analysis_list.append(res[9])
    print('end for loop for results')
    
    print('start path plot')
    path_plot_analysis(path_analysis_list)
    print('end path plot')
        #ende paula

    logging.info('Calculation finished, getting results.')
    for i in range(len(z_delta_list)):
        z_delta_array = np.maximum(z_delta_array, z_delta_list[i])
        flux_array = np.maximum(flux_array, flux_list[i])
        count_array += cc_list[i]
        z_delta_sum += z_delta_sum_list[i]
        backcalc = np.maximum(backcalc, backcalc_list[i])
        fp_travelangle_array = np.maximum(fp_travelangle_array, fp_ta_list[i])
        sl_travelangle_array = np.maximum(sl_travelangle_array, sl_ta_list[i])   
        #Chris
        travel_length_array = np.maximum(travel_length_array, travel_length_list[i])
        #ende chris    
        #Paula
        flow_energy_array = np.maximum(flow_energy_array, flow_energy_list[i])
        #ende paula
        
    # Save Calculated tiles (Zwischenspeicher)
    np.save(tempDir / ("res_z_delta_%s_%s" % (optTuple[0], optTuple[1])), z_delta_array)
    np.save(tempDir / ("res_z_delta_sum_%s_%s" % (optTuple[0], optTuple[1])), z_delta_sum)
    np.save(tempDir / ("res_flux_%s_%s" % (optTuple[0], optTuple[1])), flux_array)
    np.save(tempDir / ("res_count_%s_%s" % (optTuple[0], optTuple[1])), count_array)
    np.save(tempDir / ("res_fp_%s_%s" % (optTuple[0], optTuple[1])), fp_travelangle_array)
    np.save(tempDir / ("res_sl_%s_%s" % (optTuple[0], optTuple[1])), sl_travelangle_array)
    #Chris
    np.save(tempDir / ("res_travel_length_%s_%s" % (optTuple[0], optTuple[1])), travel_length_array)
    #chris ende
    #Paula
    np.save(tempDir / ("res_flow_energy_%s_%s" % (optTuple[0], optTuple[1])), flow_energy_array)
    #np.save(tempDir / ("res_path_list_%s_%s" % (optTuple[0], optTuple[1])), path_list_list)
    #ende paula
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

    #Chris
    travel_length_array = np.zeros_like(dem, dtype=np.float32)
    #ende chris
    #Paula
    flow_energy_array = np.zeros_like(dem, dtype=np.float32)
    path_list = []
    #ende paula
    

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

        # Michi generation
        #cell_list.append(startcell)
        cell_list = [startcell] # list of parents for current iteration
        gen_list = [cell_list]  # list of all cells (which are calculated), oreganised in generations
        child_list = []         # list of childs of the current iteration

        #for idx, cell in enumerate(cell_list):
        for gen, cell_list in enumerate(gen_list):
            flux_sum = 0
            for cell in cell_list:
        #ende michi
                row, col, flux, z_delta = cell.calc_distribution()

                # Michi generation
                #if len(flux) > 0:
                    #z_delta, flux, row, col = list(zip(*sorted(zip(z_delta, flux, row, col), reverse=False)))
                    # Sort this lists by elh, to start with the highest cell
                if len(row) > 1:  # if there are more than 1 element in list, sort it by z_delta, lowest -> highest
                    z_delta, flux, row, col = list(zip(*sorted(zip(z_delta, flux, row, col), reverse=False)))  # reverse = True == descending
                    row = list(row)
                    col = list(col)
                    flux = list(flux)
                    z_delta = list(z_delta)

                #for i in range(idx, len(cell_list)):  # Check if Cell already exists
                for i in range(len(cell_list)):
                #ende michi
                    k = 0
                    while k < len(row):
                        if row[k] == cell_list[i].rowindex and col[k] == cell_list[i].colindex:
                            cell_list[i].add_os(flux[k])
                            cell_list[i].add_parent(cell)
                            if z_delta[k] > cell_list[i].z_delta:
                                cell_list[i].z_delta = z_delta[k]

                            #row = np.delete(row, k)
                            #col = np.delete(col, k)
                            #flux = np.delete(flux, k)
                            #z_delta = np.delete(z_delta, k)

                            # MICHI generation
                            row.pop(k)
                            col.pop(k)
                            flux.pop(k)
                            z_delta.pop(k)
                            #ende michi
                        else:
                            k += 1

                # MICHI generation
                for i in range(len(child_list)):  # Check if Cell already exists in child_list
                    k = 0
                    while k < len(row):
                        if row[k] == child_list[i].rowindex and col[k] == child_list[i].colindex:
                            child_list[i].add_os(flux[k])
                            child_list[i].add_parent(cell)
                            if z_delta[k] > child_list[i].z_delta:
                                child_list[i].z_delta = z_delta[k]

                            row.pop(k)
                            col.pop(k)
                            flux.pop(k)
                            z_delta.pop(k)
                        else:
                            k += 1
                #ende michi

                for k in range(len(row)):
                    dem_ng = dem[row[k] - 1:row[k] + 2, col[k] - 1:col[k] + 2]  # neighbourhood DEM
                    if (nodata in dem_ng) or np.size(dem_ng) < 9:
                        continue
                    # Michi generation
                    #cell_list.append(
                    child_list.append(
                        Cell(row[k], col[k], dem_ng, cellsize, flux[k], z_delta[k], cell, alpha, exp, flux_threshold, max_z_delta, startcell))

            # prepare lists for next iteration
            if len(child_list) > 0:
                cell_list = child_list               
                gen_list.append(cell_list)
                child_list = []
            
        #PAULA
        #list with all paths (every startcell has one path)
        path_list.append(Path(dem, row_list[startcell_idx], col_list[startcell_idx], gen_list))
        path_list[-1].calc_all_analysis()
        #ende paula

            #Michi generation
        for gen, cell_list in enumerate(gen_list):
            for cell in cell_list:
        #ende michi

                z_delta_array[cell.rowindex, cell.colindex] = max(z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
                flux_array[cell.rowindex, cell.colindex] = max(flux_array[cell.rowindex, cell.colindex], cell.flux)
                count_array[cell.rowindex, cell.colindex] += int(1)
                z_delta_sum[cell.rowindex, cell.colindex] += cell.z_delta
                fp_travelangle_array[cell.rowindex, cell.colindex] = max(fp_travelangle_array[cell.rowindex, cell.colindex], cell.max_gamma)
                sl_travelangle_array[cell.rowindex, cell.colindex] = max(sl_travelangle_array[cell.rowindex, cell.colindex], cell.sl_gamma)
                #Chris
                travel_length_array[cell.rowindex, cell.colindex] = max(travel_length_array[cell.rowindex, cell.colindex], cell.min_distance)
                #ende chris
                #PAula
                flow_energy_array[cell.rowindex, cell.colindex] = max(flow_energy_array[cell.rowindex, cell.colindex], cell.flow_energy)
                #ende paula

            
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
    #return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array
    
    #Chris/Paula
    print('start calc path analysis')
    res_path_data = path_calc_analysis(path_list)
    print('end calc path analysis')

    return z_delta_array, flux_array, count_array, z_delta_sum, backcalc, fp_travelangle_array, sl_travelangle_array, travel_length_array, flow_energy_array, res_path_data
    #ende 

