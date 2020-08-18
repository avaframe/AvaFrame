import sys
import os
import time
import logging
import glob
import math
import numpy as np
import scipy as sp
import copy
import operator
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
debugPlotFlag = False

def readAIMECinputs(avalancheDir, dirName = 'com1DFA'):
    """
    Reads the requiered files for AIMEC postpocessing
    given an avalanche directory
    """
    cfgPath = {}
    pathPressure = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName,'dfa_pressure')
    pathFlowHeight = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName,'dfa_depth')
    pathMassBalance = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName,'dfa_mass_balance')

    if not os.path.exists(pathMassBalance):
        os.makedirs(pathMassBalance)

    profileLayer = glob.glob(os.path.join(avalancheDir, 'Inputs', 'LINES', '*aimec*.shp'))
    cfgPath['profileLayer'] = ''.join(profileLayer)

    demSource = glob.glob(os.path.join(avalancheDir, 'Inputs','*.asc'))
    try:
        assert len(demSource) == 1, 'There should be exactly one topography .asc file in ' + \
            avalancheDir + '/Inputs/'
    except AssertionError:
        raise
    cfgPath['demSource'] = ''.join(demSource)

    cfgPath['pressurefileList'] = getFileList(pathPressure)
    cfgPath['depthfileList'] = getFileList(pathFlowHeight)
    cfgPath['massfileList'] = getFileList(pathMassBalance)

    pathResult = os.path.join(avalancheDir, 'Outputs', 'AimecResults')
    cfgPath['pathResult'] = pathResult

    project_name = os.path.basename(avalancheDir)
    cfgPath['project_name'] = project_name
    path_name = os.path.basename(profileLayer[0])
    cfgPath['path_name'] = path_name
    cfgPath['dirName']  = 'com1DFA'

    return cfgPath


def getFileList(path2Folder):
    """ Get list of all files in folder """
    fileList = [path2Folder +
                os.path.sep +
                str(name) for name in
                sorted(os.listdir(path2Folder)) if os.path.isfile(os.path.join(path2Folder, name))]
    return fileList

# -----------------------------------------------------------
# Aimec main
# -----------------------------------------------------------


def mainAIMEC(cfgPath, cfg):
    """
    Main logic for AIMEC postprocessing
    """

    # Extract input parameters
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    domainWidth = float(cfgSetup['domainWidth'])
    pressureLimit = float(cfgSetup['pressureLimit'])

    log.info('Prepare data for post-ptocessing')
    # create new raster + preparing new raster assignment function
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    deskewedRasterInd = processDataInd(cfgPath, domainWidth, cfgFlags)

    # transform pressure_data and depth_data in new raster

    # assign pressure data
    log.info("Assigning pressure data to deskewed raster")
    deskewedRasterPressure = assignData(cfgPath['pressurefileList'], deskewedRasterInd)

    # assign depth data
    log.info("Assigning depth data to deskewed raster")
    deskewedRasterDepth = assignData(cfgPath['depthfileList'], deskewedRasterInd)

    # assign dem data
    log.info("Assigning dem data to deskewed raster")
    deskewedRasterDEM = assignData([cfgPath['demSource']], deskewedRasterInd)
    dem_name = os.path.basename(cfgPath['demSource'])

    # Analyze data
    log.info('Analyzing data')
    # analyze doku
    log.info('Comparing data to reference')
    doku, runout_doku, delta_h, elevRel = analyzeDocu(pressureLimit,
                                                      cfgPath['pressurefileList'],
                                                      deskewedRasterInd,
                                                      deskewedRasterPressure,
                                                      deskewedRasterDepth,
                                                      deskewedRasterPressure[0],
                                                      deskewedRasterDEM,
                                                      with_doku=False)

    # analyze mass / entrainment
    log.info('Analyzing entrainment data')
    # determine growth index from entrainment data
    # [relMass, entMass, gr_index, gr_grad] = analyze.analyzeEntrainmentdata(cfgPath['massfileList'])
    gr_index = 0
    relMass = 0
    entMass = 0
    # analyze pressure_data and depth_data
    # determine runount, AMPP, AMD, FS,
    log.info('Analyzing data in path coordinate system')
    [runout, runout_mean, AMPP, MMPP, AMD, MMD] = analyzeDataWithDepth(deskewedRasterInd,
                                                                       pressureLimit,
                                                                       deskewedRasterPressure,
                                                                       deskewedRasterDepth,
                                                                       cfgPath,
                                                                       cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of results')
    result_visu(cfgPath, runout, AMPP, MMPP, doku, gr_index, pressureLimit)

    # -----------------------------------------------------------
    # write results to file
    # -----------------------------------------------------------
    log.info('Writing results to file')
    doku_name = None
    dam_mean = []
    out_header = ''.join(['project_name: ',  cfgPath['project_name'], '\n',
                          'path: ', cfgPath['path_name'], '\n',
                          'docu: ', str(doku_name), '\n',
                          'dhm: ', str(dem_name), '\n',
                          'domain_width: ', str(domainWidth), '\n',
                          'pressure_limit: ', str(pressureLimit), '\n',
                          'runout_doku: ', str(runout_doku), '\n',
                          'fall_height: ', str(delta_h), '\n',
                          'release_mass: ', str(relMass), '\n',
                          'elevation_release: ', str(elevRel),  '\n',
                          'filenr, runout, AMPP, MMPP, entMass, growth_index, AMD, MMD, TPs, FNs, FPs, TNs, TP_depth, TP_pressure, damages_mean (%i)\n' % len(dam_mean)])
    outname = ''.join([cfgPath['pathResult'], os.path.sep,
                       'Results_pl', str(int(pressureLimit)), '_w', str(int(domainWidth)), '.txt'])

    log.info('write output file: %s' % outname)
    resfile = [runout, AMPP, MMPP, entMass, gr_index, AMD, MMD,
               doku[0], doku[1], doku[2], doku[3], doku[4], doku[6]]
    resfile.extend(dam_mean)
    result_write(cfgPath['pressurefileList'], resfile, outname, out_header)

# -----------------------------------------------------------
# Aimec processing tools
# -----------------------------------------------------------
def processDataInd(cfgPath, domainWidth, cfgFlags):
    """
    process data ind
    this function is used to process the rasterdata such that it can be
    analysed with the methods for a regular grid
    data given in a regular grid is projected on a nonuniform grid given by
    a polyline

    JT Fischer, Uwe Schlifkowitz BFW 2010-2012
    AK BFW 2014

    input: names of rasterfiles, poly names, path width
    ouput: structure{x coordinate along new raster, y coordinate, rasterdata}
    """
    #Read input parameters
    rasterSource = cfgPath['pressurefileList'][0]
    ProfileLayer = cfgPath['profileLayer']
    w = domainWidth
    outpath = cfgPath['pathResult']
    DefaultName = cfgPath['project_name']

    log.info('Data-file %s analysed' % rasterSource)
    # read data
    # read raster data
    sourceData = IOf.readRaster(rasterSource)
    header = sourceData['header']
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize
    rasterdata = sourceData['rasterData']
    # read avaPath
    Avapath = shpConv.readLine(ProfileLayer, DefaultName, sourceData['header'])

    log.info('Creating new raster along polyline: %s' % ProfileLayer)
    # Initialize transformation dictionary
    raster_transfo = {}

    # Get new Domain Boundaries DB
    # input: ava path
    # output: Left and right side points for the domain
    DB = geoTrans.path2domain(Avapath, w, header)



    ### Make transformation matrix
    raster_transfo = makeTransfoMat(raster_transfo, DB, w, cellsize)

    # calculate the real area of the new cells as well as the s_coord
    raster_transfo = getSArea(raster_transfo)

    log.info('Size of rasterdata- old: %d x %d - new: %d x %d' % (
            np.size(rasterdata, 0), np.size(rasterdata, 1),
            np.size(raster_transfo['grid_x'], 0), np.size(raster_transfo['grid_x'], 1)))

    # affect values
    raster_transfo['header'] = header
    # put back scale and origin
    raster_transfo['s_coord'] = raster_transfo['s_coord']*cellsize
    raster_transfo['l_coord'] = raster_transfo['l_coord']*cellsize
    raster_transfo['grid_x'] = raster_transfo['grid_x']*cellsize + header.xllcorner
    raster_transfo['grid_y'] = raster_transfo['grid_y']*cellsize + header.yllcorner
    raster_transfo['rasterArea'] = raster_transfo['rasterArea']*cellsize*cellsize

    aval_data = transform(rasterSource,raster_transfo)
    # visu
    input_data = {}
    input_data['aval_data'] = aval_data
    input_data['sourceData'] = sourceData
    input_data['Avapath'] = Avapath
    input_data['DB'] = DB

    visu_transfo(raster_transfo, input_data, cfgPath, cfgFlags)


    return raster_transfo


def split_section(DB, i):
    """
    Splits the domain DB in the s direction (direction of the path)
    """
    # left edge
    xl0 = DB['DB_x_l'][i]
    xl1 = DB['DB_x_l'][i+1]
    yl0 = DB['DB_y_l'][i]
    yl1 = DB['DB_y_l'][i+1]
    dxl = xl1 - xl0
    dyl = yl1 - yl0
    Vl = np.array((dxl,dyl))
    zl = np.linalg.norm(Vl)

    # right edge
    xr0 = DB['DB_x_r'][i]
    xr1 = DB['DB_x_r'][i+1]
    yr0 = DB['DB_y_r'][i]
    yr1 = DB['DB_y_r'][i+1]
    dxr = xr1 - xr0
    dyr = yr1 - yr0
    Vr = np.array((dxr,dyr))
    zr = np.linalg.norm(Vr)

    # number of segments
    m = int(max(np.ceil(zl), np.ceil(zr))+1)
    # make left segment
    bxl = np.linspace(xl0, xl1, m)
    byl = np.linspace(yl0, yl1, m)
    # make right segment
    bxr = np.linspace(xr0, xr1, m)
    byr = np.linspace(yr0, yr1, m)

    return bxl, byl, bxr, byr, m

def makeTransfoMat(raster_transfo, DB, w, cellsize):
    """ Make transformation matrix.
        Takes a Domain Boundary and finds the (x,y) coordinates of the new raster
        (the one following the path)
    """
    # number of points describing the avaPath
    n_pnt = np.shape(DB['DB_x_r'])[0]
    ## Working with no dimentions (the cellsize scaling will be readded at the end)
    # l_coord is the distance from polyline (cross section)
    # maximum step should be smaller then the cellsize
    n_total = np.ceil(w/cellsize)
    # take the next odd integer. This ensure that the l_coord = o exists
    n_total = int(n_total+1) if ((n_total % 2) == 0) else int(n_total)
    n_2tot = int(np.floor(n_total/2))
    l_coord = np.linspace(-n_2tot, n_2tot, n_total) # this way, 0 is in l_coord

    # initialize new_rasters
    new_grid_raster_x = np.array([]) # x_coord of the points of the new raster
    new_grid_raster_y = np.array([]) # y_coord of the points of the new raster
    # loop on each section of the path
    for i in range(n_pnt-1):
        # split edges in segments
        bxl, byl, bxr, byr, m = split_section(DB, i)
        # bxl, byl, bxr, byr reprensent the s direction (olong path)
        # loop on segments of section
        for j in range(m-1):
            # this is the cross section segment (l direction)
            x = np.linspace(bxl[j], bxr[j], n_total)  # line coordinates x
            y = np.linspace(byl[j], byr[j], n_total)  # line coordinates y
            # save x and y coordinates of the new raster points
            if i==0 and j==0:
                new_grid_raster_x = x.reshape(1,n_total)
                new_grid_raster_y = y.reshape(1,n_total)
            else:
                new_grid_raster_x = np.append(new_grid_raster_x,x.reshape(1,n_total), axis=0)
                new_grid_raster_y = np.append(new_grid_raster_y,y.reshape(1,n_total), axis=0)

    # add last column
    x = np.linspace(bxl[m-1], bxr[m-1], n_total)  # line coordinates x
    y = np.linspace(byl[m-1], byr[m-1], n_total)  # line coordinates y
    new_grid_raster_x = np.append(new_grid_raster_x,x.reshape(1,n_total), axis=0)
    new_grid_raster_y = np.append(new_grid_raster_y,y.reshape(1,n_total), axis=0)

    raster_transfo['l_coord'] = l_coord
    raster_transfo['grid_x'] = new_grid_raster_x
    raster_transfo['grid_y'] = new_grid_raster_y

    return raster_transfo

def getSArea(raster_transfo):
    """
    Find the s_coord corresponding to the transformation and the Area of
    the cells of the new raster
    """
    x_coord = raster_transfo['grid_x']
    y_coord = raster_transfo['grid_y']
    # add ghost lines and columns to the coord matrix
    # in order to perform dx and dy calculation
    n, m =np.shape(x_coord)
    x_coord = np.append(x_coord, x_coord[:,-2].reshape(n, 1), axis=1)
    y_coord = np.append(y_coord, y_coord[:,-2].reshape(n, 1), axis=1)
    n, m =np.shape(x_coord)
    x_coord = np.append(x_coord, x_coord[-2,:].reshape(1, m), axis=0)
    y_coord = np.append(y_coord, y_coord[-2,:].reshape(1, m), axis=0)
    n, m =np.shape(x_coord)
    # calculate dx and dy for each point in the l direction
    dxl = x_coord[0:n-1,1:m]-x_coord[0:n-1,0:m-1]
    dyl = y_coord[0:n-1,1:m]-y_coord[0:n-1,0:m-1]
    # deduce the distance in l direction
    Vl = np.sqrt(dxl*dxl + dyl*dyl)
    # calculate dx and dy for each point in the s direction
    dxs = x_coord[1:n,0:m-1]-x_coord[0:n-1,0:m-1]
    dys = y_coord[1:n,0:m-1]-y_coord[0:n-1,0:m-1]
    # deduce the distance in s direction
    Vs = np.sqrt(dxs*dxs + dys*dys)

    # calculate area of each cell
    new_area_raster = np.abs(Vl*Vs)
    raster_transfo['rasterArea'] = new_area_raster
    # get s_coord
    ds = Vs[:,int(np.floor(m/2))-1]
    s_coord = np.cumsum(ds)-ds[0]
    raster_transfo['s_coord'] = s_coord

    return raster_transfo

def visu_transfo(raster_transfo, input_data, cfgPath, cfgFlags):
    """
    Plot and save the domain transformation figure
    """
    # read paths
    pathResult = cfgPath['pathResult']
    project_name = cfgPath['dirName']
    # read rasterdata
    sourceData = input_data['sourceData']
    header = sourceData['header']
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize
    rasterdata = sourceData['rasterData']
    # read avaPath with scale
    Avapath = input_data['Avapath']
    x_path = Avapath['x']*cellsize+xllcenter
    y_path = Avapath['y']*cellsize+yllcenter
    # read domain boundarries with scale
    DB =input_data['DB']
    DB_x_l = DB['DB_x_l']*cellsize+xllcenter
    DB_x_r = DB['DB_x_r']*cellsize+xllcenter
    DB_y_l = DB['DB_y_l']*cellsize+yllcenter
    DB_y_r = DB['DB_y_r']*cellsize+yllcenter

    figure_width = 2*10
    figure_height = 2*5
    lw = 1

    fig = plt.figure(figsize=(figure_width, figure_height), dpi=150)

#    for figure: referenz-simulation bei p_lim=1
    ax1 = plt.subplot(121)
    new_rasterdata = rasterdata
    masked_array = np.ma.masked_where(new_rasterdata == 0, new_rasterdata)
    cmap = copy.copy(matplotlib.cm.jet)
    cmap.set_bad('w', 1.)

    n, m = np.shape(new_rasterdata)
    xx, yy = np.meshgrid(np.arange(m)*cellsize+xllcenter, np.arange(n)*cellsize+yllcenter)
    ref1 = ax1.imshow(masked_array, vmin=new_rasterdata.min(),
                      vmax=new_rasterdata.max(),
                      origin='lower',
                      cmap=cmap,
                      label='pressure data',
                      aspect='auto',
                      extent=[xx.min(), xx.max(),
                              yy.min(), yy.max()])
    plt.autoscale(False)
    ref2 = plt.plot(x_path, y_path,
                    'b-', linewidth=lw, label='flow path')
    ref3 = plt.plot(DB_x_l, DB_y_l,
                    'g-', linewidth=lw, label='domain')
    ref3 = plt.plot(DB_x_r, DB_y_r,
                    'g-', linewidth=lw, label='domain')
    ref3 = plt.plot([DB_x_l, DB_x_r],[DB_y_l, DB_y_r],
                    'g-', linewidth=lw, label='domain')
    refs = [ref2[0], ref3[0]]

    labels = ['flow path', 'domain']
    ax1.title.set_text('XY Domain')
    ax1.legend(refs, labels, loc=0)
    ax1.set_xlim([xx.min(), xx.max()])
    ax1.set_ylim([yy.min(), yy.max()])
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    cbh = plt.colorbar(ref1, use_gridspec=True)
    cbh.set_label('peak pressure [kPa]')

    ax2 = plt.subplot(122)
    ax2.title.set_text('sl Domain')
    isosurf = copy.deepcopy(input_data['aval_data'])
    xx, yy = np.meshgrid(raster_transfo['l_coord'], raster_transfo['s_coord'] )
    masked_array = np.ma.masked_where(isosurf == 0, isosurf)
    cmap = copy.copy(matplotlib.cm.jet)
    cmap.set_bad('w', 1.)
    # ref0 = plt.pcolormesh(xx,yy,masked_array, vmin=isosurf.min(),
    #                   vmax=isosurf.max(), cmap=cmap,
    #                   #extent=[xx.min(), xx.max(), yy.min(), yy.max()],
    #                   label='pressure data')
    ref0 = ax2.imshow(masked_array, vmin=isosurf.min(),
                      vmax=isosurf.max(), origin='lower', cmap=cmap,
                      extent=[xx.min(), xx.max(), yy.min(), yy.max()],
                      aspect='auto', label='pressure data')
    ax2.set_xlim([xx.min(), xx.max()])
    ax2.set_ylim([yy.min(), yy.max()])
    ax2.set_xlabel('l [m]')
    ax2.set_ylabel('s [m]')
    cbh = plt.colorbar(ref0, use_gridspec=True)
    cbh.set_label('peak pressure [kPa]')

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    if cfgFlags.getboolean('savePlot'):
        outname_fin = ''.join([pathResult, '/pics/', project_name,
                               '_domTransfo', '.pdf'])
        if not os.path.exists(os.path.dirname(outname_fin)):
            os.makedirs(os.path.dirname(outname_fin))
        fig.savefig(outname_fin, transparent=True)

    plt.close(fig)

def transform(fname, raster_transfo):

    """
    Affect value to the points of the new raster (after domain transormation)
    input:
            -fname = name of rasterfile to transform
            -raster_transfo = transformation info
    ouput:
            -new_data = z, pressure or depth... corresponding to fname on the new raster
    """
    name = os.path.basename(fname)
    data = IOf.readRaster(fname)

    # read tranformation info
    new_grid_raster_x = raster_transfo['grid_x']
    new_grid_raster_y = raster_transfo['grid_y']

    n, m = np.shape(new_grid_raster_x)
    xx = new_grid_raster_x
    yy = new_grid_raster_y
    Points = {}
    Points['x'] = xx.flatten()
    Points['y'] = yy.flatten()
    Points = geoTrans.projectOnRaster_Vect(data, Points)
    new_data = Points['z'].reshape(n,m)
    log.info('Data-file: %s - raster values transferred' % (name))

    return new_data

def assignData(fnames, raster_transfo):
    """
    Affect value to the points of the new raster (after domain transormation)
    input:
            -fnames = list of names of rasterfiles to transform
            -raster_transfo = transformation info
    ouput: aval_data = z, pressure or depth... corresponding to fnames on the new rasters
    """

    maxtopo = len(fnames)
    aval_data = np.array(([None] * maxtopo))

    log.info('Transfer data of %d file(s) from old to new raster' % maxtopo)
    for i in range(maxtopo):
        fname = fnames[i]
        aval_data[i] = transform(fname,raster_transfo)

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
    cellarea = rasterInd['rasterArea']

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
            # ref0 = plt.pcolormesh(xx,yy,masked_array, vmin=isosurf.min(),
            #                   vmax=isosurf.max(), cmap=cmap,
            #                   #extent=[xx.min(), xx.max(), yy.min(), yy.max()],
            #                   label='pressure data')
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


def result_visu(cfgPath, runout, mean_max_dpp, max_max_dpp, doku, GI, dpp_threshold):
    """
    Visualize results in a nice way
    Jan-Thomas Fischer BFW 2010-2012
    AK BFW 2014-2015
    """

    fnames = cfgPath['pressurefileList']
    rasterSource = cfgPath['demSource']
    ProfileLayer = cfgPath['profileLayer']
    outpath = cfgPath['pathResult']
    DefaultName = cfgPath['project_name']

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
    AvaProfile, SplitPoint = geoTrans.prepareLine(dem, Avapath, distance=10)
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
