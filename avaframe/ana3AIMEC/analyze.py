# -*- coding: utf-8 -*-

# ==============================================================================
# packages
# ==============================================================================
import os
import logging
import numpy as np
import scipy as sp
import copy
import operator
from multiprocessing import Pool
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


# create local logger
log = logging.getLogger(__name__)


# ==============================================================================
# analyzeDocu
# ==============================================================================
def analyzeDocu(p_lim, fnames, rasterInd, pressureData,
                depthData, dokuData, dhmData, with_depth, with_dhm, with_doku):
    """
    Vergleich Simulation - Dokumentation

    Andreas Kofler, 2013

    this function is used to compare each simulation with a given doku-polyline
    (or if isn't existing a documatation: simulation#1 (ref-sim) --> documentation)

    input: p_lim(doku),names of sim-files,deskewed_rasterInd,pressureData,depthData,dokuData,with_depth=0/1
    ouput: structure{teilfächen(4),mean-values for pressure and depth(4)} +
    runout-length of doku
    """
#    initialisieren
#   fnames:       full name of data file including path and extension
#   polyfnames:   file with polylinie
    n_topo = len(fnames)
    avalData = np.array([[None for m in range(n_topo)] for n in range(8)])

    if with_doku:
        log.info('Compare pressure data with doku data')
        s_coordinate_doku = rasterInd['s_coord']
        l_coordinate_doku = rasterInd['l_coord']
        new_mask = dokuData
    else:
        log.warning('No doku found. Use ref-sim: %s instead' % (fnames[0].split('/')[-1]))
        # if no docu-polyline --> ref.sim. ≃ doku
        s_coordinate_doku = rasterInd['s_coord']
        l_coordinate_doku = rasterInd['l_coord']
        new_mask = copy.deepcopy(pressureData[0])

        doku_cross = np.array((np.nanmax(new_mask, 1),
                               np.nanmean(new_mask, 1)))
#        doku(s)
#        search in doku values
        lindex = np.nonzero(doku_cross[0] > p_lim)[0]
        if lindex.any():
            cupper_doku = min(lindex)
            clower_doku = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % p_lim)
            cupper_doku = 0
            clower_doku = 0
        runout_doku = s_coordinate_doku[clower_doku]
#    analyze shape of the avalanche front (last frontshape m)
#    from runout point frontshapelength back and measure medium??? width
#    find the point of frontal length before runout
        fs_i, fs_v = min(enumerate(abs(s_coordinate_doku - runout_doku + 200)),
                         key=operator.itemgetter(1))
#        fs_i, fs_v = min(enumerate(abs(s_coordinate_doku - runout_doku + frontshape_length)),
#                         key=operator.itemgetter(1))
        new_mask[0:fs_i] = 0
        new_mask[np.where(np.nan_to_num(new_mask) < p_lim)] = 0
        new_mask[np.where(np.nan_to_num(new_mask) >= p_lim)] = 1

#   comparison rasterdata with mask
    log.info('Sim number\tTP\tFN\tFP\tTN')

#     rasterinfo
    n_start, m_start = np.nonzero(np.nan_to_num(new_mask))
    n_start = min(n_start)

    n_total = len(rasterInd['s_coord'])
    m_total = len(rasterInd['l_coord'])
    cellarea = rasterInd['absRasterData']

    for i in range(n_topo):
        rasterdata = pressureData[i]
        # include for depth
        if with_depth:
            rasterdata_depth = depthData[i]

        """
        area
        # true positive: reality(mask)=1, model(rasterdata)=1
        # false negative: reality(mask)=1, model(rasterdata)=0
        # false positive: reality(mask)=0, model(rasterdata)=1
        # true negative: reality(mask)=0, model(rasterdata)=0
        """
#    for each pressure-file p_lim is introduced (1/3/.. kPa), where the avalanche has stopped
        new_rasterdata = copy.deepcopy(rasterdata)
        new_rasterdata[np.where(np.nan_to_num(new_rasterdata) < p_lim)] = 0
        new_rasterdata[np.where(np.nan_to_num(new_rasterdata) >= p_lim)] = 1

        tpInd = np.where((new_mask[n_start:n_total+1] == True) &
                         (new_rasterdata[n_start:n_total+1] == True))
        fpInd = np.where((new_mask[n_start:n_total+1] == False) &
                         (new_rasterdata[n_start:n_total+1] == True))
        fnInd = np.where((new_mask[n_start:n_total+1] == True) &
                         (new_rasterdata[n_start:n_total+1] == False))
        tnInd = np.where((new_mask[n_start:n_total+1] == False) &
                         (new_rasterdata[n_start:n_total+1] == False))

#     Teilrasterpunkte
        tpCount = len(tpInd[0])
        fpCount = len(fpInd[0])
        fnCount = len(fnInd[0])
        tnCount = len(tnInd[0])

#     subareas
        tp = sum(cellarea[tpInd[0] + n_start, tpInd[1]])
        fp = sum(cellarea[fpInd[0] + n_start, fpInd[1]])
        fn = sum(cellarea[fnInd[0] + n_start, fnInd[1]])
        tn = sum(cellarea[tnInd[0] + n_start, tnInd[1]])

# for mean-pressure and mean-depth over simualtion(s)
        if tpCount == 0:
            tp_pressure_mean = 0
            tp_depth_mean = 0
        else:
            tp_pressure_mean = sum(rasterdata[tpInd[0] + n_start, tpInd[1]]) / tpCount
            if with_depth:
                tp_depth_mean = sum(rasterdata_depth[tpInd[0] + n_start, tpInd[1]]) / tpCount
            else:
                tp_depth_mean = 0
        if fpCount == 0:
            fp_pressure_mean = 0
            fp_depth_mean = 0
        else:
            fp_pressure_mean = sum(rasterdata[fpInd[0] + n_start, fpInd[1]]) / fpCount
            if with_depth:
                fp_depth_mean = sum(rasterdata_depth[fpInd[0] + n_start, fpInd[1]]) / fpCount
            else:
                fp_depth_mean = 0

        avalData[0][i] = tp
        avalData[1][i] = fn
        avalData[2][i] = fp
        avalData[3][i] = tn
        avalData[4][i] = tp_depth_mean
        avalData[5][i] = fp_depth_mean
        avalData[6][i] = tp_pressure_mean
        avalData[7][i] = fp_pressure_mean

        area_sum = tp + fn + fp + tn
        log.info('%s\t %f\t %f\t %f\t %f' % (i+1, tp/area_sum, fn/area_sum,
                                             fp/area_sum, tn/area_sum))
        # runout
        # doku(s)
        doku_cross = np.array((np.nanmax(new_mask, 1),
                               np.nanmean(new_mask, 1)))
        # search in doku values
#        lindex = np.nonzero(doku_cross[0] == p_lim)[0]
        lindex = np.nonzero(doku_cross[0] == 1)[0]
        if lindex.any():
            cupper_doku = min(lindex)
            clower_doku = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % p_lim)
            cupper_doku = 0
            clower_doku = 0

    runout_doku = s_coordinate_doku[clower_doku]

    if with_dhm:
        # if dhm delta h analysis
        # Achtung Fehler in SamosAT: Druckraster und DHM-Raster stimmen nicht exakt überein!
        # Eventuell shift in assignData berücksichtigen
        new_dhm = dhmData[0]
        # find first cells that have flow - (randomly) choose simulation
        # P(s)
        p_cross = np.array((np.nanmax(rasterdata, 1),
                            np.nanmean(rasterdata, 1)))

        # search in max values
        lindex = np.nonzero(p_cross[0] > p_lim)[0]
        if lindex.any():
            cupper = min(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % p_lim)
            cupper = 0
        elevRel = new_dhm[cupper, int(m_total/2)]
        deltah = new_dhm[cupper, int(m_total/2)] - new_dhm[clower_doku, int(m_total/2)]
    else:
        deltah = 0

    return avalData, runout_doku, deltah, elevRel

# ==============================================================================
# analyzeImpact
# ==============================================================================


def analyzeImpact(fnames, pressureData, depthData, damagesData, with_depth):
    """
    Vergleich Simulation - Dokumentation von Schäden

    Andreas Kofler, 2014

    this function is used to compare each simulation with a given doku-polyline

    input: p_lim(doku),names of sim-files,deskewed_rasterInd,pressureData,depthData,dokuData,with_depth=0/1
    ouput: mean-pressure over damaged docu cells
    """

#    initialisieren
#   fnames:       full name of data file including path and extension
#   polyfnames:   file with polylinie
    n_topo = len(fnames)
    damages_mean = np.array([0 for m in range(n_topo)])
    damages_max = np.array([0 for m in range(n_topo)])

    #   comparison rasterdata with mask
    for i in range(n_topo):
        rasterdata = pressureData[i]
        # include for depth
        if with_depth:
            rasterdata_depth = depthData[i]

        tpInd = np.nonzero(damagesData)
        tpCount = len(tpInd[0])

# for mean-pressure and max-pressure over simualtion(s)
        if tpCount == 0:
            tp_damages_mean = 0
            tp_damages_max = 0
        else:
            tp_damages_mean = sum(rasterdata[tpInd[0], tpInd[1]]) / tpCount
            tp_damages_max = max(rasterdata[tpInd[0], tpInd[1]])
#            if with_depth:
#                tp_depth_mean = sum(rasterdata_depth[tpInd[0] + n_start, tpInd[1]]) / tpCount
#            else: tp_depth_mean = 0

        damages_mean[i] = tp_damages_mean
        damages_max[i] = tp_damages_max

    return damages_mean, damages_max

# ==============================================================================
# analyzeDataWithDepth
# ==============================================================================


def analyzeDataWithDepth(rasterInd, p_lim, fname, data, data_depth, visu, outpath, out):
    """
    ANALYZEData

    JT Fischer BFW 2010 - 2012
    """

    if data_depth == None:
        use_depth = None
    else:
        use_depth = 1

#    initialize Arrays
    n_topo = len(fname)

    runout = np.zeros((n_topo))
    runout_mean = np.zeros((n_topo))
    ampp = np.zeros((n_topo))
    mmpp = np.zeros((n_topo))
    amd = np.zeros((n_topo))
    mmd = np.zeros((n_topo))
    frontal_shape = np.zeros((n_topo))

    s_coordinate = rasterInd['s_coord']
    l_coordinate = rasterInd['l_coord']

    p_cross_all = np.zeros((n_topo, len(s_coordinate)))
    log.info('Sim number \t rRunout \t rampp \t ramd \t FS')
#    For each data set
    for i in range(n_topo):
        rasterdata = data[i]
#    include for depth
        if use_depth:
            rasterdata_depth = data_depth[i]
#    Maximum and averaged Values of pressure and depth
        # P(s)
        p_cross = np.array((np.nanmax(rasterdata, 1),
                            np.nanmean(rasterdata, 1)))
#        # P(l)
#        p_long = np.array((np.amax(rasterdata, 0),
#                           np.mean(rasterdata, 0)))
        # D(s)
        if use_depth:
            d_cross = np.array((np.nanmax(rasterdata_depth, 1),
                                np.nanmean(rasterdata_depth, 1)))

#    Determine runout according to maximum and averaged values
        # search in max values
        lindex = np.nonzero(p_cross[0] > p_lim)[0]
        if lindex.any():
            cupper = min(lindex)
            clower = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % p_lim)
            cupper = 0
            clower = 0
        #        search in mean values
        lindex = np.nonzero(p_cross[1] > p_lim)[0]
        if lindex.any():
            cupper_m = min(lindex)
            clower_m = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % p_lim)
            cupper_m = 0
            clower_m = 0
    #    Mean max dpp of Cross-Section
        ampp[i] = np.mean(p_cross[0][cupper:clower+1])
        mmpp[i] = max(p_cross[0][cupper:clower+1])
        if use_depth:
            amd[i] = np.mean(d_cross[0][cupper:clower+1])
            mmd[i] = max(d_cross[0][cupper:clower+1])
    #    Runout
        runout[i] = s_coordinate[clower]
        runout_mean[i] = s_coordinate[clower_m]

# analyze shape of the avalanche front (last frontshape m)
# from runout point frontshapelength back and measure medium??? width
# find the point of frontal length before runout
#        fs_i, fs_v = min(enumerate(abs(s_coordinate - runout[i] + frontshape_length)),
#                         key=operator.itemgetter(1))
# determine the width of the frontal zone
#        p_frontlong = np.array((np.amax(rasterdata[fs_i:clower+1, :], 0),
#                                   np.mean(rasterdata[fs_i:clower+1, :], 0)))
#       # P(X) % maximum value and averaged value
#        frontalIndex = np.where(p_frontlong[0] > p_lim)[0]
#       # p_long only in frontal area
#
#        if frontalIndex.any():
#            frleft = min(frontalIndex)
#            frright = max(frontalIndex)
#        else:
#            print('[DATA]: No values > p_lim found for frontal width determination. p_lim = %10.4f, too high?' % p_lim
#            frleft = 0
##            frright = len(p_long[0])-1
#            frright = 0
#
#        front_width = abs(l_coordinate[frright]-l_coordinate[frleft])
#        frontal_shape[i] = front_width/frontshape_length

        log.info('%s\t%10.4f\t%10.4f\t%10.4f' % (i+1, runout[i], ampp[i], amd[i]))

    # visu
        figure_width = 3*5
        figure_height = 3*4

        if i == 0:
            fig = plt.figure(figsize=(figure_width, figure_height), dpi=150)
            ax1 = plt.subplot(121)
            ref1 = ax1.plot([l_coordinate[0], l_coordinate[len(l_coordinate)-1]],
                            [s_coordinate[clower], s_coordinate[clower]], 'g-')
            ref2 = ax1.plot([l_coordinate[0], l_coordinate[len(l_coordinate)-1]],
                            [s_coordinate[clower_m], s_coordinate[clower_m]], 'r-')
#            ref3 = ax1.plot([l_coordinate[frleft], l_coordinate[frleft]],
#                            [s_coordinate[fs_i], s_coordinate[clower]],'k-')
#            ref4 = ax1.plot([l_coordinate[frright], l_coordinate[frright]],
#                         [s_coordinate[fs_i], s_coordinate[clower]],'k-')
            isosurf = copy.deepcopy(rasterdata)
            xx, yy = np.meshgrid(l_coordinate, s_coordinate)
            masked_array = np.ma.masked_where(isosurf == 0, isosurf)
            cmap = copy.copy(matplotlib.cm.jet)
            cmap.set_bad('w', 1.)
            ref0 = ax1.imshow(masked_array, vmin=isosurf.min(),
                              vmax=isosurf.max(), origin='lower', cmap=cmap,
                              extent=[xx.min(), xx.max(), yy.min(), yy.max()],
                              aspect='auto', label='pressure data')
            ax1.set_xlim([xx.min(), xx.max()])
            ax1.set_ylim([yy.min(), yy.max()])
            ax1.set_xlabel('l [m]')
            ax1.set_ylabel('s [m]')
            cbh = plt.colorbar(ref0, use_gridspec=True)
            cbh.set_label('peak pressure [kPa]')
#            ax1.legend((ref1[0], ref2[0], ref3[0]),
#                       ('runout max', 'runout mean', 'frontal width'), loc=0)
            ax1.legend((ref1[0], ref2[0]), ('runout max', 'runout mean'), loc=0)

        p_cross_all[i] = p_cross[0]

    ax2 = plt.subplot(122)
#    p_mean = p_cross_all.mean(axis=0)
    p_median = np.median(p_cross_all, axis=0)
    p_percentile = sp.percentile(p_cross_all, [2.5, 50, 97.5], axis=0)
    ax2.fill_betweenx(s_coordinate, p_percentile[2], p_percentile[0],
                      facecolor=[.8, .8, .8], alpha=0.5)
    ref1 = mpatches.Patch(alpha=0.5, color=[.8, .8, .8])
    ax2.plot(p_median, s_coordinate, color='r')
    ref2 = mlines.Line2D([], [], color='r', linewidth=2)
#    ax2.plot(p_mean, s_coordinate, color='b')
#    ref3 = mlines.Line2D([], [], color='b', linewidth=2)
    ax2.set_ylabel('s [m]')
    ax2.set_ylim([yy.min(), yy.max()])
    ax2.set_xlim(auto=True)
    ax2.set_xlabel('Pmax(s) [kPa]')
    ax2.legend((ref1, ref2),
               ('quantiles', 'median'), loc=0)

    fig.tight_layout()
    if out:
        pro_name = fname[0].split('/')[-3]
        outname_fin = ''.join([outpath, '/pics/', pro_name, '_dptr',
                               str(int(p_lim)), '_simulationsl', '.pdf'])
        if not os.path.exists(os.path.dirname(outname_fin)):
            os.makedirs(os.path.dirname(outname_fin))
        fig.savefig(outname_fin, transparent=True)

    if visu:
        plt.show()
    else:
        plt.ioff()

    return runout, runout_mean, ampp, mmpp, amd, mmd


# ==============================================================================
#  read_write
# ==============================================================================
def read_write(fname_ent):
    #    load data
    #    time, total mass, entrained mass
    mass_time = np.loadtxt(fname_ent, skiprows=2)
    maxind, maxval = max(enumerate(mass_time[:, 1]),
                         key=operator.itemgetter(1))
    timeResults = [mass_time[0, 0], mass_time[maxind, 0], mass_time[-1, 0]]
    totMassResults = [mass_time[0, 1], mass_time[maxind, 1], mass_time[-1, 1]]
    entMassResults = [mass_time[0, 2], mass_time[maxind, 2], mass_time[-1, 2]]
    relMass = totMassResults[0]
    entMass = entMassResults[2]
#   growth results
#    growth index = maxmasse / anfang
#    growth index =  endmasse / anfangs
#    growthIndex = totMassResults[1]/totMassResults[0]
    growthIndex = totMassResults[2]/totMassResults[0]
#    (MASS_max - MASS_start) / T(MASS_max) - T(0)
#    (MASS_end - MASS_start) / T(MASS_end) - T(0)
#    growthGrad  = (totMassResults[1]-totMassResults[0])/(timeResults[1]-timeResults[0])
    growthGrad = (totMassResults[2] - totMassResults[0]) / (timeResults[2] - timeResults[0])
#    print('[ENT] %s\t %s\t %s\t %f\t %f\t' % (fname_ent.split('/')[-1],
#                                              relMass, entMass, growthIndex, growthGrad)
    return relMass, entMass, growthIndex, growthGrad

# ==============================================================================
# analyzeEntrainmentdata
# ==============================================================================


def analyzeEntrainmentdata(fnames):
    """
    Jan-Thomas Fischer BFW 2012

    header: lauf, q, es, ed, eb, rhorel, rhoSnowCover

    folgende hierachie ist zu beachten:
        ist ein entrainment dhm gesetzt wird das daraus berechnete q_DHM
        verwendet, nicht das hier gesetzte!
        wenn e_b gesetzt ist sind e_s und e_d irrelevant!
        rho_rel = 150 kg/m³ rhosc = 50 kg/m² standard...
        der zusammenhang von e_d und e_s ist uber q gegeben.. e_s = q* e_d
        wenn e_b klein ist ist q_er = q
    """

    grIndex = np.ones(len(fnames))
    grGrad = np.ones(len(fnames))
    releaseMass = np.ones(len(fnames))
    entrainedMass = np.ones(len(fnames))

    log.info('Sim number\t GI\t Ggrad')

    for i in range(len(fnames)):
        releaseMass[i], entrainedMass[i], grIndex[i], grGrad[i] = read_write(fnames[i])

    if (releaseMass == releaseMass[0]).any():
        releaseMass = releaseMass[0]
    else:
        log.warning('Release masses differs between simulations!')

    return releaseMass, entrainedMass, grIndex, grGrad
