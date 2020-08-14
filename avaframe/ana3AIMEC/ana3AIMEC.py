import sys
import os
import logging
import glob
import math
import numpy as np
import scipy as sp
import pandas as pd
import copy
import functools
import operator
from multiprocessing import Pool
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import cm

# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in3Utils.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)

# -----------------------------------------------------------
# Aimec read inputs tools
# -----------------------------------------------------------

def readAIMECinputs(avalancheDir):
    """
    Reads the requiered files for AIMEC postpocessing
    given an avalanche directory
    """
    cfgPath = {}
    pathPressure = avalancheDir + '/Outputs/com1DFA/dfa_pressure'
    pathFlowHeight = avalancheDir + '/Outputs/com1DFA/dfa_depth'
    pathMassBalance = avalancheDir + '/Outputs/com1DFA/dfa_mass_balance'

    if not os.path.exists(pathMassBalance):
        os.makedirs(pathMassBalance)

    profileLayer = glob.glob(avalancheDir + '/Inputs/LINES/*aimec*.shp')
    cfgPath['profileLayer'] = ''.join(profileLayer)

    demSource = glob.glob(avalancheDir + '/Inputs/*.asc')
    try:
        assert len(demSource) == 1, 'There should be only one and only one DEM .asc file in ' + \
            avalancheDir + '/Inputs/'
    except AssertionError as e:
        raise
    cfgPath['demSource'] = ''.join(demSource)
    pressurefileList = [pathPressure +
                        '/' +
                        str(name) for name in
                        sorted(os.listdir(pathPressure)) if os.path.isfile(os.path.join(pathPressure, name))]
    cfgPath['pressurefileList'] = pressurefileList

    depthfileList = [str(pathFlowHeight) +
                     '/' +
                     str(name) for name in
                     sorted(os.listdir(pathFlowHeight)) if os.path.isfile(os.path.join(pathFlowHeight, name))]
    cfgPath['depthfileList'] = depthfileList


    massfileList = [str(pathMassBalance) +
                     '/' +
                     str(name) for name in
                     sorted(os.listdir(pathMassBalance)) if os.path.isfile(os.path.join(pathMassBalance, name))]
    cfgPath['massfileList'] = massfileList

    pathResult = avalancheDir + '/Outputs/AimecResults'
    cfgPath['pathResult'] = pathResult

    defaultName = str(avalancheDir).split('/')[-1]
    cfgPath['defaultName'] = defaultName

    set_name = pressurefileList[0].split('/')[-4]
    cfgPath['set_name'] = set_name
    project_name = str(profileLayer).split('/')[-4]
    cfgPath['project_name'] = project_name
    path_name = str(profileLayer[0]).split('/')[-1]
    cfgPath['path_name'] = path_name

    return cfgPath


# -----------------------------------------------------------
# Aimec processing tools
# -----------------------------------------------------------

def processDataInd(cfgPath, domainWidth, cfgFlags):
    """
    process data ind
    this function is used to process the rasterdata such that it can be
    analysed with the methods for a regular grid
    data given in a regulare grid is projected on a nonuniform grid given by
    a polyline

    JT Fischer, Uwe Schlifkowitz BFW 2010-2012
    AK BFW 2014

    input: names of rasterfiles, poly names, path width
    ouput: structure{x coordinate along new raster, y coordinate, rasterdata}
    """

    rasterSource = cfgPath['pressurefileList'][0]
    ProfileLayer = cfgPath['profileLayer']
    w = domainWidth
    outpath = cfgPath['pathResult']
    DefaultName = cfgPath['defaultName']

    aval_data = {}

    m = 0
    n = 0
    m_total = 0
    n_total = 0
    m_alt = 0

    log.info('Data-file %s analysed' % rasterSource)
    # read data
    dem = IOf.readRaster(rasterSource)
    header = dem['header']
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize

    rasterdata = dem['rasterData']

    Avapath = shpConv.readLine(ProfileLayer, DefaultName, dem['header'])
    x_path = Avapath['x']
    y_path = Avapath['y']
    z_path = np.zeros(np.shape(Avapath['x']))

    log.info('Creating new raster along polyline: %s' % ProfileLayer)

#     erzeugung der eckpuntke der segmente des unregelmaesigen gitters,
#     Domain Boundaries DB
#     input: mittlerer path
#     output: eckpunkte für punkt entlang der linie
    DB_x_rl, DB_y_rl, DB_x_csz, DB_y_csz = geoTrans.path2domain(x_path, y_path,
                                                       w/2., cellsize)

#     Shift path to raster because raster is not in global grid
    DB_x_rl -= xllcenter
    DB_y_rl -= yllcenter
#    for calculation cellsize
    DB_x_csz -= xllcenter
    DB_y_csz -= yllcenter

#    use bresemham algorithm to determine the new raster
    for i in range(len(DB_x_rl[0])):
        #   for each segment check the number of CROSS-cells
        n = geoTrans.bresenham(DB_x_rl[0, i], DB_y_rl[0, i],
                      DB_x_rl[1, i], DB_y_rl[1, i], cellsize)
        n_total = max(n_total, len(n))
#   number of raster cells of edges parallel to polyline
    m = np.zeros(len(DB_x_rl[0])-1).astype('int')
    for i in range(len(DB_x_rl[0])-1):
        # left edge
        zl = geoTrans.bresenham(DB_x_rl[0, i], DB_y_rl[0, i],
                       DB_x_rl[0, i+1], DB_y_rl[0, i+1], cellsize)
        # right edge
        zr = geoTrans.bresenham(DB_x_rl[1, i], DB_y_rl[1, i],
                       DB_x_rl[1, i+1], DB_y_rl[1, i+1], cellsize)
        m[i] = max(len(zl), len(zr))
#    delete the lines that are double at ech segment connection
    m_total = sum(m) - (len(DB_x_rl[0])-2)

#    Calculation of segments
    new_raster = np.zeros((2, m_total, n_total)) + np.NaN
    # new raster filled with NaN

#    Each dataset needs its own s_coord, l_coord. This is saved to a cell array
#    and returned from this function for use in other parts of the
#    program package

#    l_coord is distance from polyline
    l_coord = np.linspace(-w/2, w/2, n_total)
    s_coord = np.zeros(m_total)
    ds1 = 0
    log.info('Transferring data from old to new raster')
    for i in range(len(DB_x_rl[0])-1):  # for each segment
        # Division of side edges in n segments
        # x-/y-values of lines are linspace(x0,x1,m) etc.
        # DB_x_rl, DB_y_rl are values from polyline2path
        bxl = np.linspace(DB_x_rl[0][i], DB_x_rl[0][i+1], m[i])  # left
        byl = np.linspace(DB_y_rl[0][i], DB_y_rl[0][i+1], m[i])

        bxr = np.linspace(DB_x_rl[1][i], DB_x_rl[1][i+1], m[i])  # right
        byr = np.linspace(DB_y_rl[1][i], DB_y_rl[1][i+1], m[i])

        # => generation of grid points for each segment
        # grid points can be assigned to original raster data, using
        # rasterize function
        new_rastersegment = np.zeros((2, m[i], n_total))

        for j in range(m[i]):
            x = np.linspace(bxl[j], bxr[j], n_total)  # line coordinates x
            y = np.linspace(byl[j], byr[j], n_total)  # line coordinates y
            for k in range(n_total):
                # x,y-Koordinaten of cells on line
                # xy_coord = bresenham(x(k),y(k),x(k),y(k),cellsize);
                xy_coord = [round(x[k]/cellsize) * cellsize,
                            round(y[k]/cellsize) * cellsize]
                # cell coordinates of new raster
                xy_ind = [xy_coord[0]/cellsize + 1, xy_coord[1]/cellsize + 1]
                # translate coordinate of cell to cell index
                # THIS IS THE NEAREST NEIGHBOUR APPROXIMATION
                # Assign pressure data to loc
                new_rastersegment[:, j, k] = [xy_ind[0], xy_ind[1]]

#        For each segment following the first we must delete the first
#        line since it is identical to the last line of the previous
#        segment.

#        s_coord = x-Coordinate along Polyline.
#        % Start of Polylinie at x = 0
        m_neu = m[i]

        if (i == 0):
            # Distance from starting point of current segment
            # from beginning of polyline
            ds0 = 0
        else:
            ds0 += ds1

        ds1 = math.sqrt((x_path[i+1]-x_path[i])**2 +
                        (y_path[i+1]-y_path[i])**2)
        new_raster[:, m_alt:m_alt+m_neu, :] = [new_rastersegment[0],
                                               new_rastersegment[1]]

        s_coord[m_alt:m_alt+m_neu] = np.linspace(ds0, ds0+ds1, m[i])
        m_alt = m_alt+m_neu-1
        if (i == 0):
            s_coordmin = s_coord[0]
        s_coord -= s_coordmin

#    calclation of cellsize (for area)
    new_raster_area = np.zeros((m_total, n_total)) + np.NaN
    sum_mi = 0

#    if m_total % 2 == 0:
#        seg_boundary_lines = len(DB_x_csz)-2
#    else:
    seg_boundary_lines = len(DB_x_csz)-1
    for i in range(seg_boundary_lines):  # for offset segment boundary lines
        # DB_x_rl DB_y_rl for i --> Koordinaten
        x_DB_i = np.linspace(DB_x_csz[i][0], DB_x_csz[i][1], n_total+1)
        y_DB_i = np.linspace(DB_y_csz[i][0], DB_y_csz[i][1], n_total+1)
        # DB_x_rl DB_y_rl for i+1 --> Koordinaten
        x_DB_ii = np.linspace(DB_x_csz[i+1][0], DB_x_csz[i+1][1], n_total+1)
        y_DB_ii = np.linspace(DB_y_csz[i+1][0], DB_y_csz[i+1][1], n_total+1)

        if i % 2 == 0:  # i gerade
            for j in range(n_total):
                x_seg_j = [x_DB_i[j], x_DB_ii[j]]
                y_seg_j = [y_DB_i[j], y_DB_ii[j]]

                x_seg_jj = [x_DB_i[j+1], x_DB_ii[j+1]]
                y_seg_jj = [y_DB_i[j+1], y_DB_ii[j+1]]

                k = 0
                a = sum_mi
                new_raster_area[a][j] = 1./2 * ((y_seg_j[k]-y_seg_jj[k+1]) *
                                                (x_seg_jj[k]-x_seg_j[k+1]) +
                                                (y_seg_j[k+1]-y_seg_jj[k]) *
                                                (x_seg_j[k]-x_seg_jj[k+1]))

            sum_mi = sum_mi+1
        else:  # i ungerade
            m_i = int(np.floor((i/2)))  # m for each segment
            for j in range(n_total):
                x_seg_j = np.linspace(x_DB_i[j], x_DB_ii[j], m[m_i]-1)
                y_seg_j = np.linspace(y_DB_i[j], y_DB_ii[j], m[m_i]-1)

                x_seg_jj = np.linspace(x_DB_i[j+1], x_DB_ii[j+1], m[m_i]-1)
                y_seg_jj = np.linspace(y_DB_i[j+1], y_DB_ii[j+1], m[m_i]-1)

                for k in range(m[m_i]-2):
                    a = sum_mi+k
                    # print(np.shape(new_raster_area))
                    # print(a,j)
                    new_raster_area[a, j] = 1./2*((y_seg_j[k]-y_seg_jj[k+1]) *
                                                  (x_seg_jj[k]-x_seg_j[k+1]) +
                                                  (y_seg_j[k+1]-y_seg_jj[k]) *
                                                  (x_seg_j[k]-x_seg_jj[k+1]))

            sum_mi = a + 1

    log.info('Size of rasterdata- old: %d x %d - new: %d x %d' % (
        np.size(rasterdata, 0), np.size(rasterdata, 1),
        np.size(new_raster, 1), np.size(new_raster, 2)))

    aval_data['header'] = header
    aval_data['s_coord'] = s_coord
    aval_data['l_coord'] = l_coord
    aval_data['rasterData'] = new_raster
    aval_data['absRasterData'] = abs(new_raster_area)

    # visu
    figure_width = 2*5
    figure_height = 2*4
    lw = 1

    fig = plt.figure(figsize=(figure_width, figure_height), dpi=150)

#    for figure: referenz-simulation bei p_lim=1
    new_rasterdata = rasterdata
    masked_array = np.ma.masked_where(new_rasterdata == 0, new_rasterdata)
    cmap = copy.copy(matplotlib.cm.jet)
    cmap.set_bad('w', 1.)

    n, m = np.shape(new_rasterdata)
    xx, yy = np.meshgrid(np.arange(m), np.arange(n))

    ref1 = plt.imshow(masked_array, vmin=new_rasterdata.min(),
                      vmax=new_rasterdata.max(),
                      origin='lower',
                      cmap=cmap,
                      label='pressure data',
                      aspect='auto',
                      extent=[xx.min()*cellsize+xllcenter, xx.max()*cellsize+xllcenter,
                              yy.min()*cellsize+yllcenter, yy.max()*cellsize+yllcenter])
    plt.autoscale(False)
    ref2 = plt.plot(x_path, y_path,
                    'b-', linewidth=lw, label='flow path')
    ref3 = plt.plot(DB_x_rl+xllcenter, DB_y_rl+yllcenter,
                    'g-', linewidth=lw, label='domain')
    ref3 = plt.plot(DB_x_rl.T+xllcenter, DB_y_rl.T+yllcenter,
                    'g-', linewidth=lw, label='domain')
    refs = [ref2[0], ref3[0]]

    labels = ['flow path', 'domain']

    plt.legend(refs, labels, loc=0)
    plt.xlim([xx.min()*cellsize+xllcenter, xx.max()*cellsize+xllcenter])
    plt.ylim([yy.min()*cellsize+yllcenter, yy.max()*cellsize+yllcenter])
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    cbh = plt.colorbar()
    cbh.set_label('peak pressure [kPa]')
    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    if cfgFlags.getboolean('savePlot'):
        pro_name = rasterSource.split('/')[-3]  # CoSiCa-samos-structure
#        pro_name = fnames[0].split('/')[-5] + '_' + fnames[0].split('/')[-2] # DAKUMO_structure
        outname_fin = ''.join([outpath, '/pics/', pro_name,
                               '_simulationxy', '.pdf'])
        if not os.path.exists(os.path.dirname(outname_fin)):
            os.makedirs(os.path.dirname(outname_fin))
        fig.savefig(outname_fin, transparent=True)

    plt.close(fig)

    return aval_data


def transform(fname, rasterIndData):
    """
    transform
    this function is used to process the rasterdata such that it can be
    analysed with the methods for a regular grid
    data given in a regular grid is projected on a nonuniform grid given
    a transformation raster

    JT Fischer, Uwe Schlifkowitz BFW 2010-2012
    AK BFW 2014

    input: names of rasterfiles, transformation raster
    ouput: structure{x coordinate along new raster, y coordinate, rasterdata_xind, rasterdata_yind}

    fnames:       full name of data file including path and extension
    rasterIndData: raster of new shape containing indices from corresponding points of old raster
    """
    name = fname.split('/')

    # xy_oldind = rasterIndData[3].astype('int')
    xy_oldind = rasterIndData['rasterData'].astype('int')

    dem = IOf.readRaster(fname)
    header = dem['header']
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize

    rasterdata = dem['rasterData']

#        out of bounds counter
    i_oob = 0
    i_ib = 0

    new_raster = np.zeros((len(rasterIndData['s_coord']), len(rasterIndData['l_coord'])))

    for x_ind in range(new_raster.shape[0]):
        for y_ind in range(new_raster.shape[1]):
            i_ib += 1
            try:
                new_raster[x_ind, y_ind] = rasterdata[xy_oldind[1]
                                                      [x_ind, y_ind]][xy_oldind[0][x_ind, y_ind]]
            except:
                i_oob += 1
                new_raster[x_ind, y_ind] = np.NaN

    log.info('Data-file: %s - %d raster values transferred - %d out of original raster bounds!' %
                 (name[-1], i_ib-i_oob, i_oob))

    return new_raster


def assignData(fnames, rasterIndData):
    """
    assignData

    this function is used to process the rasterdata such that it can be
    analysed with the methods for a regular grid
    data given in a regular grid is projected on a nonuniform grid given
    a transformation raster

    JT Fischer, Uwe Schlifkowitz BFW 2010-2012
    AK BFW 2014

    input: names of rasterfiles, transformation raster
    ouput: structure{x coordinate along new raster, y coordinate, rasterdata_xind, rasterdata_yind}

    fnames:       full name of data file including path and extension
    rasterIndData: raster of new shape containing indices from corresponding points of old raster
    """
# if started without arguments
#    if (nargin == 0):
#        return None

    maxtopo = len(fnames)
#    aval_data = np.array([[None for m in xrange(4)] for n in xrange(maxtopo)])
    aval_data = np.array(([None] * maxtopo))

    log.info('Transfer data of %d file(s) from old to new raster' % maxtopo)

    pool = Pool()
    aval_data = pool.map(functools.partial(
        transform, rasterIndData=rasterIndData), fnames)
    pool.close()
    pool.join()

    return aval_data


# -----------------------------------------------------------
# Aimec analysis tools
# -----------------------------------------------------------

def analyzeDocu(p_lim, fnames, rasterInd, pressureData,
                depthData, dokuData, dhmData, with_doku):
    """
    Compare Simulation - Documentation

    Andreas Kofler, 2013

    this function is used to compare each simulation with a given doku-polyline
    (or if isn't existing a documatation: simulation#1 (ref-sim) --> documentation)

    input: p_lim(doku),names of sim-files,deskewed_rasterInd,pressureData,depthData,dokuData
    ouput: structure{teilfächen(4),mean-values for pressure and depth(4)} +
    runout-length of doku
    """


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
            tp_depth_mean = sum(rasterdata_depth[tpInd[0] + n_start, tpInd[1]]) / tpCount
        if fpCount == 0:
            fp_pressure_mean = 0
            fp_depth_mean = 0
        else:
            fp_pressure_mean = sum(rasterdata[fpInd[0] + n_start, fpInd[1]]) / fpCount
            fp_depth_mean = sum(rasterdata_depth[fpInd[0] + n_start, fpInd[1]]) / fpCount

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

    return avalData, runout_doku, deltah, elevRel


def analyzeDataWithDepth(rasterInd, p_lim, data, data_depth, cfgPath, cfgFlags):
    """
    ANALYZEData

    JT Fischer BFW 2010 - 2012
    """
    fname = cfgPath['pressurefileList']
    outpath = cfgPath['pathResult']

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

    if cfgFlags.getboolean('savePlot'):
        pro_name = fname[0].split('/')[-3]
        outname_fin = ''.join([outpath, '/pics/', pro_name, '_dptr',
                               str(int(p_lim)), '_simulationsl', '.pdf'])
        if not os.path.exists(os.path.dirname(outname_fin)):
            os.makedirs(os.path.dirname(outname_fin))
        fig.savefig(outname_fin, transparent=True)

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()

    return runout, runout_mean, ampp, mmpp, amd, mmd

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



# -----------------------------------------------------------
# Aimec out tools
# -----------------------------------------------------------

def result_write(data_name, data, outfile, header):
    """
    This function is used as standart output file

    example path: 'log/doublestep_5_5m_calibration_results.txt'
    example header: 'runout, maxvel, rho, mu, tau_0, R_s^0, kappa, R, B, q \n'

    INPUT: data, path, header
    """

#    chekc if folder exists / create
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    output = data
    fid = open(outfile, 'w')
    fid.write(header)
    for i in range(len(output[0])):
        tmp = os.path.basename(data_name[i])
#        name = tmp.split('.')[0] # CoSiCa-Samos
        name = os.path.splitext(tmp)[0]  # DAKUMO
        fid.write('%s' % name)
        for j in range(len(output)):
            try:
                fid.write(',%5.2f' % output[j][i])
            except:
                fid.write(',NaN')
        fid.write('\n')
    fid.close()

    log.info('File written: %s' % outfile)


def colorvar(k, k_end, colorflag, disp=0):
    """
    jt colorvariation editor - JT 2012
    determine how color changes from runnumber = 1 to runnumber = runlength
    input: runnumber,runlength,colorflag,verbose
    output: [R G B]

    possible colorflags:
    'ry': red to yellow
    'bb': blue to light blue
    'pw': pink to white
    'kg' black to green
    'bp' blue to pink
    'pb' pink to blue
    'gy' green to yellow
    'cw' cyan to white
    'kr' black to red
    'gb' green to blue
    'rp' red to pink
    'yw' yellow to white
    'kb' black to blue
    'kc' black to cyan
    'kp' black to pink
    'kw' black to white
    """

    colors = {
        'ry': [1., k/k_end, 0.],  # rot zu gelb
        'bb': [0., k/k_end, 1.],  # blau zu hellbalu
        'pw': [0.8, k/k_end, 0.8],  # pink zu weiss
        'kg': [0., k/k_end, 0.],  # schwarz zu gruen
        'bp': [k/k_end, 0., 1.],  # blau zu pink
        'pb': [1.-k/k_end, 1., 0.],  # blau zu pink
        'gy': [k/k_end, 1., 0.],  # green zu yellow
        'cw': [k/k_end, 1., 1.],  # cyan zu weiss
        'kr': [k/k_end, 1., 1.],  # black to red
        'gb': [0., 1., k/k_end],  # gruen zu blau
        'rp': [1., 0., k/k_end],  # rot zu pink
        'yw': [1., 1., k/k_end],  # yellow to white
        'kb': [0., 0., k/k_end],  # black tp blue
        'kc': [0., k/k_end, k/k_end],  # black zu cyan
        'kp': [k/k_end, 0., k/k_end],  # black zu pink
        'kw': [1.-k/k_end, 1.-k/k_end, 1.-k/k_end]  # black to white
    }

    colornames = {
        'ry': 'red to yellow',
        'bb': 'blue to cyan',
        'pw': 'pink to white',
        'kg': 'black to green',
        'bp': 'blue to pink',
        'pb': 'blue to pink',
        'gy': 'green to yellow',
        'cw': 'cyan to white',
        'kr': 'black to red',
        'gb': 'green to blue',
        'rp': 'rot to pink',
        'yw': 'yellow to white',
        'kb': 'black to blue',
        'kc': 'black to cyan',
        'kp': 'black to pink',
        'kw': 'black to white'
    }

    if colorflag.lower() in colors:
        farbe = colors.get(colorflag.lower())
        if k == 0:
            log.info('Color is: %s' % colornames.get(colorflag.lower()))
    else:
        farbe = [0, 0, 0]
        if k == 0:
            log.info('Color is black')

    return farbe


def result_visu(fnames, rasterSource, ProfileLayer, runout, mean_max_dpp,
                max_max_dpp, doku, GI, dpp_threshold, outpath, DefaultName=None):
    """
    Visualize results in a nice way
    Jan-Thomas Fischer BFW 2010-2012
    AK BFW 2014-2015
    """
    cvar = ['ry', 'bb', 'pw', 'gy']
    colorflag = cvar[0]

    figure_width = 7*2
    figure_height = 4*2
    fs = 20
    mks = 10
    lw = 2
    # includes flag for y axis -
    # 1 = rddp
    # 2 = frontal shape
    # 3 = groth index
    # 4 = runout for intrapraevent
    # 5 = pressure data
    flag = 3
    if (len(fnames) > 100):
        plot_density = 1
    else:
        plot_density = 0

    if flag == 1:
        log.info('Visualizing pressure data')
        tipo = 'rapp'
        data = mean_max_dpp / mean_max_dpp[0]
        yaxis_label = 'rAPP [-]'
        ytick_increment = 0.25
        ymax = 3
    elif flag == 2:
        log.info('Visualizing EGU growth index data')
        tipo = 'GI'
        data = GI
        yaxis_label = 'growth index [GI]'
        ytick_increment = 2
    elif flag == 3:
        log.info('Visualizing pressure data')
        tipo = 'rmpp'
        data = max_max_dpp / max_max_dpp[0]
        yaxis_label = 'rMPP [-]'
        ytick_increment = 0.1
#        ymax = 0.12
#        ymin = 0.08
#        ymax = 0.5
#        ymin = 0.3
        ymax = max(data[1:])+(max(data[1:])-min(data[1:]))*0.1
        ymin = min(data[1:])-(max(data[1:])-min(data[1:]))*0.1
    else:
        log.error('Wrong flag')
        return None

    # read data
    dem = IOf.readRaster(rasterSource)
    header = dem['header']
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize

    rasterdata = dem['rasterData']

    Avapath = shpConv.readLine(ProfileLayer, DefaultName, dem['header'])
    AvaProfile, SplitPoint, splitPoint = geoTrans.prepareLine(
        dem, Avapath, splitPoint=None, distance=10)
    x_path = AvaProfile['x']
    y_path = AvaProfile['y']
    z_path = AvaProfile['z']
    s_path = AvaProfile['s']

    xlim_prof_axis = max(s_path) + 50

    # Final result diagram - z_profile+data
    fig = plt.figure(figsize=(figure_width, figure_height), dpi=300)

    markers = ['+', 'o', 'x', '*', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '.',
               '^', 'v', '>', '<', 'p', 'h', '.']
    mk = 0

#    show flow path
    ax1 = fig.add_subplot(111)
#    plt.xlim([0, xlim_prof_axis])
#    plt.ylim([0, math.ceil(max(data)+0.25)])
#    plt.ylim([0, ymax])
#    plt.yticks(np.arange([0, math.ceil(max(data)+0.25), ytick_increment]))
    ax1.set_ylabel(yaxis_label, color='b', fontsize=2*fs)
    ax1.set_xlabel(''.join(['s [m] - runout with ', str(dpp_threshold),
                            ' kPa threshold']), color='black', fontsize=2*fs)
    if plot_density:  # estimate 2D histogram --> create pcolormesh
        nbins = 100
        H, xedges, yedges = np.histogram2d(runout, data, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.Blues)
#        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.cool)
        cbar = plt.colorbar(data_density, orientation='horizontal')
        cbar.ax.set_ylabel('Counts')
    ax2 = ax1.twinx()
    ax2.set_ylabel('z [m]', color='g', fontsize=2*fs)
    ax2.plot(s_path, z_path, color='green', label='path', linestyle='--', linewidth=2*lw)
    plt.xlim([0, xlim_prof_axis])
    plt.ylim([math.floor(min(z_path)/10)*10, math.ceil(max(z_path)/10)*10])
    if not plot_density:
        for k in range(len(runout)):
            topo_name = fnames[k].split('/')[-1]
            pfarbe = colorvar(float(k), len(runout), colorflag)
            if k == 0:
                ax1.plot(runout[k], data[k], marker='+',
                         markersize=2*mks, color='g', label=topo_name)
    #            plt.yticks(np.arange([0,5000,250]))
                # Make the y-tick labels of first axes match the line color.
                for tl in ax1.get_yticklabels():
                    tl.set_color('b')
            else:
                ax1.plot(runout[k], data[k], label=topo_name, marker=markers[mk],
                         markersize=mks, color=pfarbe, linewidth=lw)
            mk = mk+1
            if mk == len(markers):
                mk = 1
    plt.grid('on')
#    plt.legend()
#    ax1.legend(loc=0)
#    ax2.legend(loc=0)

    pro_name = fnames[0].split('/')[-3]
    outname_fin = ''.join([outpath, '/pics/', pro_name, '_dptr',
                           str(int(dpp_threshold)), '_', tipo, '.pdf'])

    if not os.path.exists(os.path.dirname(outname_fin)):
        os.makedirs(os.path.dirname(outname_fin))
    fig.savefig(outname_fin, transparent=True)

    plt.close(fig)

    # Final result diagram - roc-plots
    rTP = (np.array(doku[0]) / (float(doku[0][0]) + float(doku[1][0]))).astype(float)
    rFP = (np.array(doku[2]) / (float(doku[2][0]) + float(doku[3][0]))).astype(float)


#    rFP = (np.array(doku[2]) / (float(doku[0][0]) + float(doku[1][0]))).astype(float)

    fig = plt.figure(figsize=(figure_width, figure_height), dpi=300)

    mk = 0
    ax1 = fig.add_subplot(111)
    ax1.set_ylabel('True positive rate', fontsize=2*fs)
    ax1.set_xlabel('False positive rate', fontsize=2*fs)
    if plot_density:  # estimate 2D histogram --> create pcolormesh
        nbins = 100
        H, xedges, yedges = np.histogram2d(rFP, rTP, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
#        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.Blues)
        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.cool)
        cbar = plt.colorbar(data_density, orientation='horizontal')
        cbar.ax.set_ylabel('hit rate density')
    if not plot_density:
        for k in range(len(rTP)):
            topo_name = fnames[k].split('/')[-1]
            pfarbe = colorvar(float(k), len(rTP), colorflag)
            ax1.plot(rFP[k], rTP[k], label=topo_name, marker=markers[mk],
                     markersize=mks, color=pfarbe, linewidth=lw)
            mk = mk+1
            if mk == len(markers):
                mk = 0
    plt.xlim([0, max(1, max(rFP))])
    plt.ylim([0, 1])
    plt.grid('on')

    outname_fin = ''.join([outpath, '/pics/', pro_name, '_dptr',
                           str(int(dpp_threshold)), '_ROC.pdf'])

    if not os.path.exists(os.path.dirname(outname_fin)):
        os.makedirs(os.path.dirname(outname_fin))
    fig.savefig(outname_fin, transparent=True)

    plt.close(fig)

    return
