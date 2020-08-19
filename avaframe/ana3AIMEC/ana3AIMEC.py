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


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in3Utils.ascUtils as IOf
import avaframe.out3SimpPlot.outAIMEC as outAimec

# create local logger
log = logging.getLogger(__name__)

# -----------------------------------------------------------
# Aimec read inputs tools
# -----------------------------------------------------------
debugPlotFlag = False


def readAIMECinputs(avalancheDir, dirName='com1DFA'):
    """
    Reads the requiered files for AIMEC postpocessing
    given an avalanche directory
    """
    cfgPath = {}
    pathPressure = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName, 'dfa_pressure')
    pathFlowHeight = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName, 'dfa_depth')
    pathMassBalance = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName, 'dfa_mass_balance')

    if not os.path.exists(pathMassBalance):
        os.makedirs(pathMassBalance)

    profileLayer = glob.glob(os.path.join(avalancheDir, 'Inputs', 'LINES', '*aimec*.shp'))
    cfgPath['profileLayer'] = ''.join(profileLayer)

    demSource = glob.glob(os.path.join(avalancheDir, 'Inputs', '*.asc'))
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
    cfgPath['dirName'] = 'com1DFA'

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
    raster_transfo = processDataInd(cfgPath, domainWidth, cfgFlags)

    # transform pressure_data and depth_data in new raster
    newRasters = {}
    # assign pressure data
    log.info("Assigning pressure data to deskewed raster")
    newRasterPressure = assignData(cfgPath['pressurefileList'], raster_transfo)
    newRasters['newRasterPressure'] = newRasterPressure
    # assign depth data
    log.info("Assigning depth data to deskewed raster")
    newRasterDepth = assignData(cfgPath['depthfileList'], raster_transfo)
    newRasters['newRasterDepth'] = newRasterDepth
    # assign dem data
    log.info("Assigning dem data to deskewed raster")
    newRasterDEM = assignData([cfgPath['demSource']], raster_transfo)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.info('Analyzing data')

    # analyze mass / entrainment
    log.info('Analyzing entrainment data')
    # determine growth index from entrainment data
    # [relMass, entMass, gr_index, gr_grad] = analyzeEntrainmentdata(cfgPath['massfileList'])
    gr_index = 0
    relMass = 0
    entMass = 0
    # analyze pressure_data and depth_data
    # determine runount, AMPP, AMD, FS,
    log.info('Analyzing data in path coordinate system')
    resAnalysis = analyzeData(raster_transfo, pressureLimit, newRasters, cfgPath, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    # log.info('Visualisation of results')
    # outAimec.result_visu(cfgPath, resAnalysis, doku, gr_index, pressureLimit)

    # -----------------------------------------------------------
    # write results to file
    # -----------------------------------------------------------
    log.info('Writing results to file')

    outAimec.result_write(cfgPath, cfgSetup, resAnalysis)

# -----------------------------------------------------------
# Aimec processing tools
# -----------------------------------------------------------


def processDataInd(cfgPath, domainWidth, cfgFlags):
    """
    this function is used to process the rasterdata such that it can be
    analysed with the methods for a regular grid
    data given in a regular grid is projected on a nonuniform grid given by
    a polyline
    This function returns the information about this transformation

    input: cfgPath, domainWidth, cfgFlags
    ouput: raster_transfo =
            -(x,y) coordinates of the points of the new raster
            -(s,l) new coordinate System
            -rasterArea, real area of the cells of the new raster
    """
    # Read input parameters
    rasterSource = cfgPath['pressurefileList'][0] #cfgPath['demSource'] #
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

    # Make transformation matrix
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

    aval_data = transform(rasterSource, raster_transfo)
    # visu
    input_data = {}
    input_data['aval_data'] = aval_data
    input_data['sourceData'] = sourceData
    input_data['Avapath'] = Avapath
    input_data['DB'] = DB

    outAimec.visu_transfo(raster_transfo, input_data, cfgPath, cfgFlags)

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
    Vl = np.array((dxl, dyl))
    zl = np.linalg.norm(Vl)

    # right edge
    xr0 = DB['DB_x_r'][i]
    xr1 = DB['DB_x_r'][i+1]
    yr0 = DB['DB_y_r'][i]
    yr1 = DB['DB_y_r'][i+1]
    dxr = xr1 - xr0
    dyr = yr1 - yr0
    Vr = np.array((dxr, dyr))
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
    # Working with no dimentions (the cellsize scaling will be readded at the end)
    # l_coord is the distance from polyline (cross section)
    # maximum step should be smaller then the cellsize
    n_total = np.ceil(w/cellsize)
    # take the next odd integer. This ensure that the l_coord = o exists
    n_total = int(n_total+1) if ((n_total % 2) == 0) else int(n_total)
    n_2tot = int(np.floor(n_total/2))
    l_coord = np.linspace(-n_2tot, n_2tot, n_total)  # this way, 0 is in l_coord

    # initialize new_rasters
    new_grid_raster_x = np.array([])  # x_coord of the points of the new raster
    new_grid_raster_y = np.array([])  # y_coord of the points of the new raster
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
            if i == 0 and j == 0:
                new_grid_raster_x = x.reshape(1, n_total)
                new_grid_raster_y = y.reshape(1, n_total)
            else:
                new_grid_raster_x = np.append(new_grid_raster_x, x.reshape(1, n_total), axis=0)
                new_grid_raster_y = np.append(new_grid_raster_y, y.reshape(1, n_total), axis=0)

    # add last column
    x = np.linspace(bxl[m-1], bxr[m-1], n_total)  # line coordinates x
    y = np.linspace(byl[m-1], byr[m-1], n_total)  # line coordinates y
    new_grid_raster_x = np.append(new_grid_raster_x, x.reshape(1, n_total), axis=0)
    new_grid_raster_y = np.append(new_grid_raster_y, y.reshape(1, n_total), axis=0)

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
    n, m = np.shape(x_coord)
    x_coord = np.append(x_coord, x_coord[:, -2].reshape(n, 1), axis=1)
    y_coord = np.append(y_coord, y_coord[:, -2].reshape(n, 1), axis=1)
    n, m = np.shape(x_coord)
    x_coord = np.append(x_coord, x_coord[-2, :].reshape(1, m), axis=0)
    y_coord = np.append(y_coord, y_coord[-2, :].reshape(1, m), axis=0)
    n, m = np.shape(x_coord)
    # calculate dx and dy for each point in the l direction
    dxl = x_coord[0:n-1, 1:m]-x_coord[0:n-1, 0:m-1]
    dyl = y_coord[0:n-1, 1:m]-y_coord[0:n-1, 0:m-1]
    # deduce the distance in l direction
    Vl = np.sqrt(dxl*dxl + dyl*dyl)
    # calculate dx and dy for each point in the s direction
    dxs = x_coord[1:n, 0:m-1]-x_coord[0:n-1, 0:m-1]
    dys = y_coord[1:n, 0:m-1]-y_coord[0:n-1, 0:m-1]
    # deduce the distance in s direction
    Vs = np.sqrt(dxs*dxs + dys*dys)

    # calculate area of each cell
    new_area_raster = np.abs(Vl*Vs)
    raster_transfo['rasterArea'] = new_area_raster
    # get s_coord
    ds = Vs[:, int(np.floor(m/2))-1]
    s_coord = np.cumsum(ds)-ds[0]
    raster_transfo['s_coord'] = s_coord

    return raster_transfo


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
    Points, i_ib, i_oob = geoTrans.projectOnRaster_Vect(data, Points, interp='bilinear')
    new_data = Points['z'].reshape(n, m)
    log.info('Data-file: %s - %d raster values transferred - %d out of original raster bounds!' % (name, i_ib-i_oob, i_oob))

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
        aval_data[i] = transform(fname, raster_transfo)

    return aval_data


# -----------------------------------------------------------
# Aimec analysis tools
# -----------------------------------------------------------

def analyzeData(raster_transfo, p_lim, newRasters, cfgPath, cfgFlags):
    """
    ANALYZEData
    """
    fname = cfgPath['pressurefileList']
    fname_mass = cfgPath['massfileList']
    outpath = cfgPath['pathResult']

    dataPressure = newRasters['newRasterPressure']
    dataDepth = newRasters['newRasterDepth']
    dataDEM = newRasters['newRasterDEM']
    s_coord = raster_transfo['s_coord']
    l_coord = raster_transfo['l_coord']
    rasterArea = raster_transfo['rasterArea']

    resAnalysis = {}

    # initialize Arrays
    n_topo = len(fname)
    runout = np.zeros((n_topo))
    runout_mean = np.zeros((n_topo))
    ampp = np.zeros((n_topo))
    mmpp = np.zeros((n_topo))
    amd = np.zeros((n_topo))
    mmd = np.zeros((n_topo))
    elevRel = np.zeros((n_topo))
    deltaH = np.zeros((n_topo))
    grIndex = np.ones((n_topo))
    grGrad = np.ones((n_topo))
    releaseMass = np.ones((n_topo))
    entrainedMass = np.ones((n_topo))

    n = np.shape(l_coord)[0]
    p_cross_all = np.zeros((n_topo, len(s_coord)))
    # log.info('Sim number \t rRunout \t rampp \t ramd \t FS')
    log.info('{: <15} {: <15} {: <15} {: <15}'.format(
        'Sim number ', 'rRunout ', 'rampp ', 'ramd ', 'FS'))
    # For each data set
    for i in range(n_topo):
        rasterdataPres = dataPressure[i]
        rasterdataDepth = dataDepth[i]

        # get mean max for each cross section for pressure
        presCrossMean = np.nansum(rasterdataPres*rasterArea, axis=1)/np.nansum(rasterArea, axis=1)
        presCrossMax = np.nanmax(rasterdataPres, 1)
        # also get the Area corresponding to those cells
        ind_presCrossMax = np.nanargmax(rasterdataPres, 1)
        ind_1 = np.arange(np.shape(rasterdataPres)[0])
        AreapresCrossMax = rasterArea[ind_1, ind_presCrossMax]
        # get mean max for each cross section for pressure
        dCrossMean = np.nansum(rasterdataDepth*rasterArea, axis=1)/np.nansum(rasterArea, axis=1)
        dCrossMax = np.nanmax(rasterdataDepth, 1)
        # also get the Area corresponding to those cells
        ind_dCrossMax = np.nanargmax(rasterdataDepth, 1)
        ind_1 = np.arange(np.shape(rasterdataDepth)[0])
        AreadCrossMax = rasterArea[ind_1, ind_dCrossMax]

        p_cross_all[i] = presCrossMax
        #   Determine runout according to maximum and averaged values
        # search in max values
        lindex = np.nonzero(presCrossMax > p_lim)[0]
        if lindex.any():
            cupper = min(lindex)
            clower = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % p_lim)
            cupper = 0
            clower = 0
        # search in mean values
        lindex = np.nonzero(presCrossMean > p_lim)[0]
        if lindex.any():
            cupper_m = min(lindex)
            clower_m = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % p_lim)
            cupper_m = 0
            clower_m = 0
        # Mean max dpp of Cross-Section
        ampp[i] = np.nansum((presCrossMax*AreapresCrossMax)[cupper:clower+1]) / \
            np.nansum(AreapresCrossMax[cupper:clower+1])
        mmpp[i] = max(presCrossMax[cupper:clower+1])

        amd[i] = np.nansum((dCrossMax*AreadCrossMax)[cupper:clower+1]) / \
            np.nansum(AreadCrossMax[cupper:clower+1])
        mmd[i] = max(dCrossMax[cupper:clower+1])
    #    Runout
        runout[i] = s_coord[clower]
        runout_mean[i] = s_coord[clower_m]

        elevRel[i] = dataDEM[cupper, int(np.floor(n/2)+1)]
        deltaH[i] = dataDEM[cupper, int(np.floor(n/2)+1)] - dataDEM[clower, int(np.floor(n/2)+1)]

        # analyze mass
        try:
            releaseMass[i], entrainedMass[i], grIndex[i], grGrad[i] = read_write(fname_mass[i])
            if (releaseMass == releaseMass[0]).any():
                releaseMass = releaseMass[0]
            else:
                log.warning('Release masses differs between simulations!')
        except IndexError:
            releaseMass[i] = np.NaN
            entrainedMass[i] = np.NaN
            grIndex[i] = np.NaN
            grGrad[i] = np.NaN

        # log.info('%s\t%10.4f\t%10.4f\t%10.4f' % (i+1, runout[i], ampp[i], amd[i]))
        log.info('{: <15} {:<15.4f} {:<15.4f} {:<15.4f}'.format(*[i+1, runout[i], ampp[i], amd[i]]))

    resAnalysis['runout'] = runout
    resAnalysis['runout_mean'] = runout_mean
    resAnalysis['AMPP'] = ampp
    resAnalysis['MMPP'] = mmpp
    resAnalysis['AMD'] = amd
    resAnalysis['MMD'] = mmd
    resAnalysis['elevRel'] = elevRel
    resAnalysis['deltaH'] = deltaH
    resAnalysis['relMass'] = releaseMass
    resAnalysis['entMass'] = entrainedMass
    resAnalysis['growthIndex'] = grIndex
    resAnalysis['growthGrad'] = grGrad

    # prepare for plot
    inputPlot = {}
    inputPlot['p_mean'] = p_cross_all.mean(axis=0)
    inputPlot['p_median'] = np.median(p_cross_all, axis=0)
    inputPlot['p_percentile'] = sp.percentile(p_cross_all, [2.5, 50, 97.5], axis=0)
    inputPlot['dataPressure'] = newRasters['newRasterPressure']
    inputPlot['runout'] = resAnalysis['runout']
    inputPlot['runout_mean'] = resAnalysis['runout_mean']
    inputPlot['pressureLimit'] = p_lim

    outAimec.visu_runout(raster_transfo, inputPlot, cfgPath, cfgFlags)

    return resAnalysis


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
    growthIndex = totMassResults[2]/totMassResults[0]
    growthGrad = (totMassResults[2] - totMassResults[0]) / (timeResults[2] - timeResults[0])
    return relMass, entMass, growthIndex, growthGrad
