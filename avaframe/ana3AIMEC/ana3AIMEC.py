"""
    Main logic for AIMEC post processing
"""

import os
import logging
import glob
import numpy as np
import copy


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.outAIMEC as outAimec

# create local logger
log = logging.getLogger(__name__)

# -----------------------------------------------------------
# Aimec read inputs tools
# -----------------------------------------------------------


def readAIMECinputs(avalancheDir, dirName='com1DFA'):
    """ Read inputs for AIMEC postprocessing

    Reads the requiered files location for AIMEC postpocessing
    given an avalanche directory

    Parameters
    ----------
    avalancheDir : str
        path to directory of avalanche to analyze
    dirName : str
        optional string with name of the module results to analyze

    Returns
    -------
    cfgPath : dict
        dictionary with path to data to analyze
    """
    cfgPath = {}
    pathPressure = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName,
                                'dfa_pressure')
    pathFlowHeight = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName,
                                  'dfa_depth')
    pathMassBalance = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName,
                                   'dfa_mass_balance')
    pathSpeed = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName,
                             'dfa_speed')

    if not os.path.exists(pathMassBalance):
        os.makedirs(pathMassBalance)

    profileLayer = glob.glob(os.path.join(avalancheDir, 'Inputs', 'LINES',
                                          '*aimec*.shp'))
    try:
        message = 'There should be exactly one path_aimec.shp file containing the avalanche path in ' + avalancheDir + '/Inputs/LINES/'
        assert len(profileLayer) == 1, message
    except AssertionError:
        raise
    cfgPath['profileLayer'] = ''.join(profileLayer)

    splitPointLayer = glob.glob(os.path.join(avalancheDir, 'Inputs', 'POINTS',
                                             '*.shp'))
    try:
        message = 'There should be exactly one .shp file containing the split point in ' + avalancheDir + '/Inputs/POINTS/'
        assert len(splitPointLayer) == 1, message
    except AssertionError:
        raise
    cfgPath['splitPointSource'] = ''.join(splitPointLayer)

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
    cfgPath['speedfileList'] = getFileList(pathSpeed)

    pathResult = os.path.join(avalancheDir, 'Outputs', 'ana3AIMEC', dirName)
    cfgPath['pathResult'] = pathResult

    projectName = os.path.basename(avalancheDir)
    cfgPath['projectName'] = projectName
    pathName = os.path.basename(profileLayer[0])
    cfgPath['pathName'] = pathName
    cfgPath['dirName'] = 'com1DFA'

    return cfgPath


def getFileList(path2Folder):
    """ Get sorted list of all files in folder path2Folder"""
    fileList = [path2Folder +
                os.path.sep +
                str(name) for name in
                sorted(os.listdir(path2Folder))
                if os.path.isfile(os.path.join(path2Folder, name))]
    return fileList

# -----------------------------------------------------------
# Aimec main
# -----------------------------------------------------------


def mainAIMEC(cfgPath, cfg):
    """ Main logic for AIMEC postprocessing

    Reads the required files location for ana3AIMEC postpocessing
    given an avalanche directory

    Parameters
    ----------
    cfgPath : dict
        dictionary with path to data to analyze
    cfg : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini

    Returns
    -------
    rasterTransfo: dict
        domain transformation information
    newRasters: dict
        raster data expressed in the new coordinates
    resAnalysis: dict
        results of ana3AIMEC analysis
    """

    # Extract input config parameters
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    pressureLimit = float(cfgSetup['pressureLimit'])
    interpMethod = cfgSetup['interpMethod']

    log.info('Prepare data for post-ptocessing')
    # Make domain transformation
    log.info("Creating new deskewed raster and preparing new raster assignment function")
    rasterTransfo = makeDomainTransfo(cfgPath, cfgSetup, cfgFlags)

    # transform pressure_data, depth_data and speed_data in new raster
    newRasters = {}
    # assign pressure data
    log.info("Assigning pressure data to deskewed raster")
    newRasters['newRasterPressure'] = assignData(cfgPath['pressurefileList'],
                                                 rasterTransfo, interpMethod)
    # assign depth data
    log.info("Assigning depth data to deskewed raster")
    newRasters['newRasterDepth'] = assignData(cfgPath['depthfileList'],
                                              rasterTransfo, interpMethod)
    # assign speed data
    if cfgPath['speedfileList']:
        log.info("Assigning speed data to deskewed raster")
        newRasters['newRasterSpeed'] = assignData(cfgPath['speedfileList'],
                                                  rasterTransfo, interpMethod)

    # assign dem data
    log.info("Assigning dem data to deskewed raster")
    newRasterDEM = assignData([cfgPath['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.info('Analyzing data in path coordinate system')
    resAnalysis = analyzeData(rasterTransfo, pressureLimit, newRasters,
                              cfgPath, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of results')
    outAimec.resultVisu(cfgPath, cfgFlags, rasterTransfo, resAnalysis,
                        pressureLimit)

    # -----------------------------------------------------------
    # write results to file
    # -----------------------------------------------------------
    log.info('Writing results to file')

    outAimec.resultWrite(cfgPath, cfgSetup, resAnalysis)

    return rasterTransfo, newRasters, resAnalysis

# -----------------------------------------------------------
# Aimec processing tools
# -----------------------------------------------------------


def makeDomainTransfo(cfgPath, cfgSetup, cfgFlags):
    """ Make domain transformation

    This function returns the information about the domain transformation
    Data given on a regular grid is projected on a nonuniform grid following
    a polyline to end up with "straightend raster"

    Parameters
    ----------
    cfgPath : dict
        dictionary with path to data to analyze
    cfgSetup : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini
        regarding domain transformation (domain width w, runOutAngle or
        interpolation method)
    cfgFlags : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini
        regarding plotting and writing results flags

    Returns
    -------
    rasterTransfo: dict
        domain transformation information:
            gridx: 2D numpy array
                x coord of the new raster points in old coord system
            gridy: 2D numpy array
                y coord of the new raster points in old coord system
            s: 1D numpy array
                new coord system in the polyline direction
            l: 1D numpy array
                new coord system in the cross direction
            x: 1D numpy array
                coord of the resampled polyline in old coord system
            y: 1D numpy array
                coord of the resampled polyline in old coord system
            rasterArea: 2D numpy array
                real area of the cells of the new raster
            indRunoutPoint: int
                index for start of the runout area (in s)
    """
    # Read input parameters
    rasterSource = cfgPath['pressurefileList'][0]
    demSource = cfgPath['demSource']
    ProfileLayer = cfgPath['profileLayer']
    DefaultName = cfgPath['projectName']

    w = float(cfgSetup['domainWidth'])
    runoutAngle = float(cfgSetup['runoutAngle'])
    interpMethod = cfgSetup['interpMethod']

    log.info('Data-file %s analysed' % rasterSource)
    # read data
    # read raster data
    sourceData = IOf.readRaster(rasterSource)
    dem = IOf.readRaster(demSource)
    header = sourceData['header']
    xllc = header.xllcenter
    yllc = header.yllcenter
    cellsize = header.cellsize
    rasterdata = sourceData['rasterData']
    # Initialize transformation dictionary
    rasterTransfo = {}
    rasterTransfo['domainWidth'] = w
    rasterTransfo['xllc'] = xllc
    rasterTransfo['yllc'] = yllc
    rasterTransfo['cellsize'] = cellsize

    # read avaPath
    Avapath = shpConv.readLine(ProfileLayer, DefaultName, sourceData)
    # read split point
    splitPoint = shpConv.readPoints(cfgPath['splitPointSource'], sourceData)
    # add 'z' coordinate to the avaPath
    Avapath = geoTrans.projectOnRaster(dem, Avapath)
    # reverse avaPath if necessary
    _, Avapath = geoTrans.checkProfile(Avapath, projSplitPoint=None)

    log.info('Creating new raster along polyline: %s' % ProfileLayer)

    # Get new Domain Boundaries DB
    # input: ava path
    # output: Left and right side points for the domain
    rasterTransfo = geoTrans.path2domain(Avapath, rasterTransfo)

    # Make transformation matrix
    rasterTransfo = makeTransfoMat(rasterTransfo)

    # calculate the real area of the new cells as well as the scoord
    rasterTransfo = getSArea(rasterTransfo)

    log.info('Size of rasterdata- old: %d x %d - new: %d x %d' % (
        np.size(rasterdata, 0), np.size(rasterdata, 1),
        np.size(rasterTransfo['gridx'], 0),
        np.size(rasterTransfo['gridx'], 1)))

    ##########################################################################
    # affect values
    rasterTransfo['header'] = header
    # put back scale and origin
    rasterTransfo['s'] = rasterTransfo['s']*cellsize
    rasterTransfo['l'] = rasterTransfo['l']*cellsize
    rasterTransfo['gridx'] = rasterTransfo['gridx']*cellsize + xllc
    rasterTransfo['gridy'] = rasterTransfo['gridy']*cellsize + yllc
    rasterTransfo['rasterArea'] = rasterTransfo['rasterArea']*cellsize*cellsize
    # (x,y) coordinates of the resamples avapth (centerline where l = 0)
    n = np.shape(rasterTransfo['l'])[0]
    indCenter = int(np.floor(n/2))
    rasterTransfo['x'] = rasterTransfo['gridx'][:, indCenter]
    rasterTransfo['y'] = rasterTransfo['gridy'][:, indCenter]


    #################################################################
    # add 'z' coordinate to the centerline
    rasterTransfo = geoTrans.projectOnRaster(dem, rasterTransfo)
    # find projection of split point on the centerline centerline
    projPoint = geoTrans.findSplitPoint(rasterTransfo, splitPoint)
    rasterTransfo['indSplit'] = projPoint['indSplit']
    # prepare find start of runout area points
    angle, tmp, delta_ind = geoTrans.prepareAngleProfile(runoutAngle,
                                                         rasterTransfo)
    # find the runout point: first point under runoutAngle
    indRunoutPoint = geoTrans.findAngleProfile(tmp, delta_ind)
    if (runoutAngle <= angle[indRunoutPoint] and
            runoutAngle >= angle[indRunoutPoint+1]):
        rasterTransfo['indRunoutPoint'] = indRunoutPoint
        rasterTransfo['runoutAngle'] = runoutAngle
        log.info('Measuring run-out length from the %s ° point' % runoutAngle)
    else:
        log.warning('No %s ° point found. Check splitPoint position or runoutAngle value.' % runoutAngle)
        rasterTransfo['indRunoutPoint'] = indRunoutPoint
        rasterTransfo['runoutAngle'] = (angle[indRunoutPoint] +
                                        angle[indRunoutPoint+1])/2
        log.info('Measuring run-out length from the %s ° point' % rasterTransfo['runoutAngle'])

    slRaster = transform(rasterSource, rasterTransfo, interpMethod)

    ###########################################################################
    # visualisation
    inputData = {}
    inputData['slRaster'] = slRaster
    inputData['xyRaster'] = rasterdata
    inputData['Avapath'] = Avapath

    outAimec.visuTransfo(rasterTransfo, inputData, cfgPath, cfgFlags)

    return rasterTransfo


def splitSection(DB, i):
    """ Splits the ith segment of domain boundary DB in the s direction
    (direction of the path)

    Parameters
    ----------
    DB: dict
        domain Boundary dictionary:
            DBXl: 1D numpy array
                x coord of the left boundary
            DBXr: 1D numpy array
                x coord of the right boundary
            DBYl: 1D numpy array
                y coord of the left boundary
            DBYr: 1D numpy array
                y coord of the right boundary
    i: int
        index of the segment of DB to split

    Returns
    -------
    (x,y) coordinates of the ith left and right splited Boundaries
    bxl: 1D numpy array
    byl: 1D numpy array
    bxr: 1D numpy array
    byr: 1D numpy array
    m: int
        number of elements on the new segments
    """
    # left edge
    xl0 = DB['DBXl'][i]
    xl1 = DB['DBXl'][i+1]
    yl0 = DB['DBYl'][i]
    yl1 = DB['DBYl'][i+1]
    dxl = xl1 - xl0
    dyl = yl1 - yl0
    Vl = np.array((dxl, dyl))
    zl = np.linalg.norm(Vl)

    # right edge
    xr0 = DB['DBXr'][i]
    xr1 = DB['DBXr'][i+1]
    yr0 = DB['DBYr'][i]
    yr1 = DB['DBYr'][i+1]
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


def makeTransfoMat(rasterTransfo):
    """ Make transformation matrix.

    Takes a Domain Boundary and finds the (x,y) coordinates of the new
    raster point (the one following the path)

    Parameters
    ----------
    rasterTransfo: dict
        dictionary containing:
            domainWidth: float
            cellsize: float
            DBXl: 1D numpy array
                x coord of the left boundary
            DBXr: 1D numpy array
                x coord of the right boundary
            DBYl: 1D numpy array
                y coord of the left boundary
            DBYr: 1D numpy array
                y coord of the right boundary

    Returns
    -------
    rasterTransfo: dict
        rasterTransfo dictionary updated with
            gridx: 2D numpy array
                x coord of the new raster points in old coord system
            gridy: 2D numpy array
                y coord of the new raster points in old coord system
            l: 1D numpy array
                new coord system in the cross direction
    """
    w = rasterTransfo['domainWidth']
    cellsize = rasterTransfo['cellsize']
    # number of points describing the avaPath
    n_pnt = np.shape(rasterTransfo['DBXr'])[0]
    # Working with no dimentions
    # (the cellsize scaling will be readded at the end)
    # lcoord is the distance from the polyline (cross section)
    # maximum step should be smaller then the cellsize
    nTot = np.ceil(w/cellsize)
    # take the next odd integer. This ensures that the lcoord = 0 exists
    nTot = int(nTot+1) if ((nTot % 2) == 0) else int(nTot)
    n2Tot = int(np.floor(nTot/2))
    lcoord = np.linspace(-n2Tot, n2Tot, nTot)  # this way, 0 is in lcoord

    # initialize new_rasters
    newGridRasterX = np.array([])  # xcoord of the points of the new raster
    newGridRasterY = np.array([])  # ycoord of the points of the new raster
    # loop on each section of the path
    for i in range(n_pnt-1):
        # split edges in segments
        bxl, byl, bxr, byr, m = splitSection(rasterTransfo, i)
        # bxl, byl, bxr, byr reprensent the s direction (olong path)
        # loop on segments of section
        for j in range(m-1):
            # this is the cross section segment (l direction)
            x = np.linspace(bxl[j], bxr[j], nTot)  # line coordinates x
            y = np.linspace(byl[j], byr[j], nTot)  # line coordinates y
            # save x and y coordinates of the new raster points
            if i == 0 and j == 0:
                newGridRasterX = x.reshape(1, nTot)
                newGridRasterY = y.reshape(1, nTot)
            else:
                newGridRasterX = np.append(newGridRasterX, x.reshape(1, nTot),
                                           axis=0)
                newGridRasterY = np.append(newGridRasterY, y.reshape(1, nTot),
                                           axis=0)

    # add last column
    x = np.linspace(bxl[m-1], bxr[m-1], nTot)  # line coordinates x
    y = np.linspace(byl[m-1], byr[m-1], nTot)  # line coordinates y
    newGridRasterX = np.append(newGridRasterX, x.reshape(1, nTot), axis=0)
    newGridRasterY = np.append(newGridRasterY, y.reshape(1, nTot), axis=0)

    rasterTransfo['l'] = lcoord
    rasterTransfo['gridx'] = newGridRasterX
    rasterTransfo['gridy'] = newGridRasterY

    return rasterTransfo


def getSArea(rasterTransfo):
    """ Get the s curvilinear coorsinate and Area on the new raster

    Find the scoord corresponding to the transformation and the Area of
    the cells of the new raster

    Parameters
    ----------
    rasterTransfo: dict
        dictionary containing:
            domainWidth: float
            cellsize: float
            DBXl: 1D numpy array
                x coord of the left boundary
            DBXr: 1D numpy array
                x coord of the right boundary
            DBYl: 1D numpy array
                y coord of the left boundary
            DBYr: 1D numpy array
                y coord of the right boundary
            gridx: 2D numpy array
                x coord of the new raster points in old coord system
            gridy: 2D numpy array
                y coord of the new raster points in old coord system

    Returns
    -------
    rasterTransfo: dict
        rasterTransfo dictionary updated with
            s: 1D numpy array
                new coord system in the polyline direction
            rasterArea: 2D numpy array
                real area of the cells of the new raster
    """
    xcoord = rasterTransfo['gridx']
    ycoord = rasterTransfo['gridy']
    # add ghost lines and columns to the coord matrix
    # in order to perform dx and dy calculation
    n, m = np.shape(xcoord)
    xcoord = np.append(xcoord, xcoord[:, -2].reshape(n, 1), axis=1)
    ycoord = np.append(ycoord, ycoord[:, -2].reshape(n, 1), axis=1)
    n, m = np.shape(xcoord)
    xcoord = np.append(xcoord, xcoord[-2, :].reshape(1, m), axis=0)
    ycoord = np.append(ycoord, ycoord[-2, :].reshape(1, m), axis=0)
    n, m = np.shape(xcoord)
    # calculate dx and dy for each point in the l direction
    dxl = xcoord[0:n-1, 1:m]-xcoord[0:n-1, 0:m-1]
    dyl = ycoord[0:n-1, 1:m]-ycoord[0:n-1, 0:m-1]
    # # deduce the distance in l direction
    # Vl2 = (dxl*dxl + dyl*dyl)
    # Vl = np.sqrt(Vl2)
    # calculate dx and dy for each point in the s direction
    dxs = xcoord[1:n, 0:m-1]-xcoord[0:n-1, 0:m-1]
    dys = ycoord[1:n, 0:m-1]-ycoord[0:n-1, 0:m-1]
    # deduce the distance in s direction
    Vs2 = (dxs*dxs + dys*dys)
    Vs = np.sqrt(Vs2)

    # Method 1
    # calculate area of each cell using Gauss's area formula
    # sum_i{x_(i)*y_(i+1) + y_(i)*x_(i+1)}
    # i = 1 : top left in the matrix : [0:n-1, 0:m-1]
    x1 = xcoord[0:n-1, 0:m-1]
    y1 = ycoord[0:n-1, 0:m-1]
    # i = 2 : top right in the matrix : [1:n, 0:m-1]
    x2 = xcoord[1:n, 0:m-1]
    y2 = ycoord[1:n, 0:m-1]
    # i = 3 : bottom right in the matrix : [1:n, 1:m]
    x3 = xcoord[1:n, 1:m]
    y3 = ycoord[1:n, 1:m]
    # i = 4 : bottom left in the matrix : [0:n-1, 1:m]
    x4 = xcoord[0:n-1, 1:m]
    y4 = ycoord[0:n-1, 1:m]

    Area = np.abs(x1*y2-y1*x2 + x2*y3-y2*x3 + x3*y4-y3*x4 + x4*y1-y4*x1)/2

    # Method 2
    # calculate area of each cell assuming they are parallelogramms
    # (which is wrong)
    # newAreaRaster = np.abs(dxl*dys - dxs*dyl)

    # save Area matrix
    rasterTransfo['rasterArea'] = Area
    # get scoord
    ds = Vs[:, int(np.floor(m/2))-1]
    scoord = np.cumsum(ds)-ds[0]
    rasterTransfo['s'] = scoord

    return rasterTransfo


def transform(fname, rasterTransfo, interpMethod):
    """ Transfer data from old raster to new raster

    Affect value to the points of the new raster (after domain transormation)

    Parameters
    ----------
    fname: str
        path to rasterfile to transform
    rasterTransfo: dict
        transformation information
    interpMethod: str
        interpolation method to chose between 'nearest' and 'bilinear'

    Returns
    -------
    newData: 2D numpy array
        new_data = z, pressure or depth... corresponding to fname on
        the new raster
    """
    name = os.path.basename(fname)
    data = IOf.readRaster(fname)

    # read tranformation info
    newGridRasterX = rasterTransfo['gridx']
    newGridRasterY = rasterTransfo['gridy']

    n, m = np.shape(newGridRasterX)
    xx = newGridRasterX
    yy = newGridRasterY
    Points = {}
    Points['x'] = xx.flatten()
    Points['y'] = yy.flatten()
    iib = len(Points['x'])
    Points, ioob = geoTrans.projectOnRasterVect(data, Points,
                                                     interp=interpMethod)
    newData = Points['z'].reshape(n, m)
    log.info('Data-file: %s - %d raster values transferred - %d out of original raster bounds!' %
             (name, iib-ioob, ioob))

    return newData


def assignData(fnames, rasterTransfo, interpMethod):
    """ Transfer data from old raster to new raster

    Loops through paths in fnames and calls transfom
    Transform affects values to the points of the new raster
    (after domain transormation)

    Parameters
    ----------
    fname: list of str
        list of path to rasterfile to transform
    rasterTransfo: dict
        transformation information
    interpMethod: str
        interpolation method to chose between 'nearest' and 'bilinear'

    Returns
    -------
    newData: 2D numpy array
        new_data = z, pressure or depth... corresponding to fname on
        the new raster
    """

    maxtopo = len(fnames)
    avalData = np.array(([None] * maxtopo))

    log.info('Transfer data of %d file(s) from old to new raster' % maxtopo)
    for i in range(maxtopo):
        fname = fnames[i]
        avalData[i] = transform(fname, rasterTransfo, interpMethod)

    return avalData


# -----------------------------------------------------------
# Aimec analysis tools
# -----------------------------------------------------------

def analyzeData(rasterTransfo, pLim, newRasters, cfgPath, cfgFlags):
    """ Analyse pressure and depth transformed data

    Main function for result analysis
    Plots the desiered figures and saves results as txt files

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    pLim: float
        numerical value of the pressure limit to use
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgPath: dict
        path to data to analyse
    cfgFlags: configparser
        configparser with plot and write flags

    Returns
    -------
    Plots the desiered figures and saves results as txt files
    resAnalysis: dict
        resAnalysis dictionnary containing all results
    """

    resAnalysis = analyzeFields(rasterTransfo, pLim, newRasters, cfgPath)

    resAnalysis = analyzeArea(rasterTransfo, resAnalysis, pLim, newRasters, cfgPath, cfgFlags)

    outAimec.visuSimple(rasterTransfo, resAnalysis, newRasters, cfgPath, cfgFlags)
    outAimec.visuRunout(rasterTransfo, resAnalysis, pLim, newRasters, cfgPath, cfgFlags)

    return resAnalysis


def analyzeFields(rasterTransfo, pLim, newRasters, cfgPath):
    """ Analyse pressure and depth transformed data

    Analyse pressure depth and speed.
    Calculate runout, Max Peak Pressure, Average PP...
    Get mass and entrainement

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    pLim: float
        numerical value of the pressure limit to use
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgPath: dict
        path to data to analyse

    Returns
    -------
    resAnalysis: dict
        resAnalysis dictionnary containing all results:
            -runout: 2D numpy array
                    containing for each simulation analyzed the x and
                    y coord of the runout point as well as the runout distance
                    measured from the begining of the run-out area. run-out
                    calculated with the MAX pressure in each cross section
            -runoutMean: 2D numpy array
                    containing for each simulation analyzed the x
                    and y coord of the runout point as well as the runout
                    distance measured from the begining of the run-out area.
                    run-out calculated with the MEAN pressure in each cross
                    section
            -AMPP: 1D numpy array
                    containing for each simulation analyzed the average
                    max peak pressure
            -MMPP: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak pressure
            -AMD: 1D numpy array
                    containing for each simulation analyzed the average
                    max peak flow depth
            -MMD: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak flow depth
            -AMS: 1D numpy array
                    containing for each simulation analyzed the average
                    max peak speed
            -MMS: 1D numpy array
                    containing for each simulation analyzed the max max
                    peak speed
            -elevRel: 1D numpy array
                    containing for each simulation analyzed the
                    elevation of the release area (based on first point with
                    peak pressure > pLim)
            -deltaH: 1D numpy array
                    containing for each simulation analyzed the
                    elevation fall difference between elevRel and altitude of
                    run-out point
            -relMass: 1D numpy array
                    containing for each simulation analyzed the
                    release mass
            -entMass: 1D numpy array
                    containing for each simulation analyzed the
                    entrained mass
            -finalMass: 1D numpy array
                    containing for each simulation analyzed the
                    final mass
            -relativMassDiff: 1D numpy array
                    containing for each simulation analyzed
                    the final mass diff with ref (in %)
            -growthIndex: 1D numpy array
                    containing for each simulation analyzed the
                    growth index
            -growthGrad: 1D numpy array
                    containing for each simulation analyzed the
                    growth gradient
            -pCrossAll: 2D numpy array
                    containing for each simulation analyzed the
                    max peak pressure in each cross section
            -pressureLimit: float
                    pressure threshold
            -runoutAngle: float
                    angle of the slope at the beginning of the run-out
                    area (given in input)
    """
    # read inputs
    fname = cfgPath['pressurefileList']
    fnameMass = cfgPath['massfileList']

    dataPressure = newRasters['newRasterPressure']
    dataDepth = newRasters['newRasterDepth']
    dataSpeed = newRasters['newRasterSpeed']
    dataDEM = newRasters['newRasterDEM']
    scoord = rasterTransfo['s']
    lcoord = rasterTransfo['l']
    x = rasterTransfo['x']
    y = rasterTransfo['y']
    indRunoutPoint = rasterTransfo['indRunoutPoint']
    sBeta = scoord[indRunoutPoint]

    resAnalysis = {}

    # initialize Arrays
    nTopo = len(fname)
    massDiffers = False
    runout = np.empty((3, nTopo))
    runoutMean = np.empty((3, nTopo))
    ampp = np.empty((nTopo))
    mmpp = np.empty((nTopo))
    amd = np.empty((nTopo))
    mmd = np.empty((nTopo))
    ams = np.empty((nTopo))
    mms = np.empty((nTopo))
    elevRel = np.empty((nTopo))
    deltaH = np.empty((nTopo))
    grIndex = np.empty((nTopo))
    grGrad = np.empty((nTopo))
    releaseMass = np.empty((nTopo))
    entrainedMass = np.empty((nTopo))
    finalMass = np.empty((nTopo))
    relativMassDiff = np.empty((nTopo))

    n = np.shape(lcoord)[0]
    pCrossAll = np.zeros((nTopo, len(scoord)))
    # For each data set
    for i in range(nTopo):
        rasterdataPres = dataPressure[i]
        rasterdataDepth = dataDepth[i]
        rasterdataSpeed = dataSpeed[i]
        rasterArea = rasterTransfo['rasterArea']
        # rasterArea[np.where(np.isnan(rasterdataPres))] = np.nan

        # Mean max in each Cross-Section for each field
        ampp[i], mmpp[i], cInd, pCrossAll[i] = getMaxMeanValues(
            rasterdataPres, rasterArea, pLim, cInd=None)
        amd[i], mmd[i], cInd, _ = getMaxMeanValues(rasterdataDepth, rasterArea,
                                                   pLim, cInd=cInd)
        ams[i], mms[i], cInd, _ = getMaxMeanValues(rasterdataSpeed, rasterArea,
                                                   pLim, cInd=cInd)
        #    Runout
        cupper = cInd['cupper']
        clower = cInd['clower']
        clowerm = cInd['clowerm']
        runout[0, i] = scoord[clower] - sBeta
        runout[1, i] = x[clower]
        runout[2, i] = y[clower]
        runoutMean[0, i] = scoord[clowerm] - sBeta
        runoutMean[1, i] = x[clower]
        runoutMean[2, i] = y[clower]

        elevRel[i] = dataDEM[cupper, int(np.floor(n/2)+1)]
        deltaH[i] = dataDEM[cupper, int(np.floor(n/2)+1)] - dataDEM[clower, int(np.floor(n/2)+1)]

        # analyze mass
        releaseMass[i], entrainedMass[i], finalMass[i], grIndex[i], grGrad[i] = readWrite(
            fnameMass[i])
        relativMassDiff[i] = (finalMass[i]-finalMass[0])/finalMass[0]*100
        if not (releaseMass[i] == releaseMass[0]):
            massDiffers = True
        log.info('{: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10} {: <10}'.format(
            'Sim number ', 'Runout ', 'ampp ', 'mmpp ', 'amd ', 'mmd ', 'GI ', 'GR '))
        log.info('{: <10} {:<10.4f} {:<10.4f} {:<10.4f} {:<10.4f} {:<10.4f} {:<10.4f} {:<10.4f}'.format(
            *[i+1, runout[0, i], ampp[i], mmpp[i], amd[i], mmd[i], grIndex[i], grGrad[i]]))
    if massDiffers:
        log.warning('Release masses differs between simulations!')
    # affect values to output dictionary
    resAnalysis['runout'] = runout
    resAnalysis['runoutMean'] = runoutMean
    resAnalysis['AMPP'] = ampp
    resAnalysis['MMPP'] = mmpp
    resAnalysis['AMD'] = amd
    resAnalysis['MMD'] = mmd
    resAnalysis['AMS'] = ams
    resAnalysis['MMS'] = mms
    resAnalysis['elevRel'] = elevRel
    resAnalysis['deltaH'] = deltaH
    resAnalysis['relMass'] = releaseMass
    resAnalysis['entMass'] = entrainedMass
    resAnalysis['finalMass'] = finalMass
    resAnalysis['relativMassDiff'] = relativMassDiff
    resAnalysis['growthIndex'] = grIndex
    resAnalysis['growthGrad'] = grGrad
    resAnalysis['pCrossAll'] = pCrossAll
    resAnalysis['pressureLimit'] = pLim
    resAnalysis['runoutAngle'] = rasterTransfo['runoutAngle']

    return resAnalysis


def analyzeArea(rasterTransfo, resAnalysis, pLim, newRasters, cfgPath, cfgFlags):
    """Compare results to reference.

    Compute True positive, False negative... areas.

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    resAnalysis: dict
        resAnalysis dictionnary containing all results to update
    pLim: float
        numerical value of the pressure limit to use
    newRasters: dict
        dictionary containing pressure, velocity and flow depth rasters after
        transformation
    cfgPath: dict
        path to data to analyse
    cfgFlags: configparser
        configparser with plot and write flags

    Returns
    -------
    resAnalysis: dict
        updated resAnalysis dictionnary:
            -TP: float
                ref = True sim2 = True
            -FN: float
                ref = False sim2 = True
            -FP: float
                ref = True sim2 = False
            -TN: float
                ref = False sim2 = False
    """
    fname = cfgPath['pressurefileList']

    dataPressure = newRasters['newRasterPressure']
    cellarea = rasterTransfo['rasterArea']
    indRunoutPoint = rasterTransfo['indRunoutPoint']

    # initialize Arrays
    nTopo = len(fname)
    TP = np.empty((nTopo))
    FN = np.empty((nTopo))
    FP = np.empty((nTopo))
    TN = np.empty((nTopo))

    # take first simulation as reference
    newMask = copy.deepcopy(dataPressure[0])
    # prepare mask for area resAnalysis
    newMask = np.where(np.isnan(newMask), 0, newMask)
    newMask = np.where(newMask < pLim, 0, newMask)
    newMask = np.where(newMask >= pLim, 1, newMask)
    # comparison rasterdata with mask
    log.info('{: <15} {: <15} {: <15} {: <15} {: <15}'.format(
        'Sim number ', 'TP ', 'FN ', 'FP ', 'TN'))
    # rasterinfo
    nStart = indRunoutPoint

    for i in range(nTopo):
        rasterdata = dataPressure[i]

        """
        area
        # true positive: reality(mask)=1, model(rasterdata)=1
        # false negative: reality(mask)=1, model(rasterdata)=0
        # false positive: reality(mask)=0, model(rasterdata)=1
        # true negative: reality(mask)=0, model(rasterdata)=0
        """
        # for each pressure-file pLim is introduced (1/3/.. kPa),
        # where the avalanche has stopped
        newRasterData = copy.deepcopy(rasterdata)
        # prepare mask for area resAnalysis
        newRasterData = np.where(np.isnan(newRasterData), 0, newRasterData)
        newRasterData = np.where(newRasterData < pLim, 0, newRasterData)
        newRasterData = np.where(newRasterData >= pLim, 1, newRasterData)

        # inputs for plot
        inputs = {}
        inputs['dataPressure'] = dataPressure
        inputs['refRasterMask'] = newMask
        inputs['newRasterMask'] = newRasterData
        inputs['nStart'] = nStart
        inputs['i'] = i

        outAimec.visuComparison(rasterTransfo, resAnalysis, inputs, cfgPath,
                                cfgFlags)

        tpInd = np.where((newMask[nStart:] == 1) &
                         (newRasterData[nStart:] == 1))
        fpInd = np.where((newMask[nStart:] == 0) &
                         (newRasterData[nStart:] == 1))
        fnInd = np.where((newMask[nStart:] == 1) &
                         (newRasterData[nStart:] == 0))
        tnInd = np.where((newMask[nStart:] == 0) &
                         (newRasterData[nStart:] == 0))

        # subareas
        tp = np.nansum(cellarea[tpInd[0] + nStart, tpInd[1]])
        fp = np.nansum(cellarea[fpInd[0] + nStart, fpInd[1]])
        fn = np.nansum(cellarea[fnInd[0] + nStart, fnInd[1]])
        tn = np.nansum(cellarea[tnInd[0] + nStart, tnInd[1]])

        # take reference (first simulation) as normalizing area
        areaSum = tp + fn

        TP[i] = tp
        FN[i] = fn
        FP[i] = fp
        TN[i] = tn

        log.info('{: <15} {:<15.4f} {:<15.4f} {:<15.4f} {:<15.4f}'.format(
            *[i+1, tp/areaSum, fn/areaSum, fp/areaSum, tn/areaSum]))

    resAnalysis['TP'] = TP
    resAnalysis['FN'] = FN
    resAnalysis['FP'] = FP
    resAnalysis['TN'] = TN

    return resAnalysis


def readWrite(fname_ent):
    """Get mass balance information

    Read mass balance files to get mass properties of the simulation
    (total mass, entrained mass...). Checks for mass conservation

    Parameters
    ----------
    fname_ent: list of str
        list of pah to mass balance files

    Returns
    -------
    relMass: float
        release mass
    entMass: float
        entrained mass
    finalMass: float
        final mass
    growthIndex: float
    growthGrad: float
    """
    #    load data
    #    time, total mass, entrained mass
    massTime = np.loadtxt(fname_ent, delimiter=',', skiprows=1)
    timeResults = [massTime[0, 0], massTime[-1, 0]]
    totMassResults = [massTime[0, 1], massTime[-1, 1]]
    relMass = totMassResults[0]
    entMass = np.sum(massTime[:, 2])
    finalMass = totMassResults[1]
    # check mass balance
    log.info('Total mass change between first and last time step in sim %s is: %.1f kg' %
             (int(os.path.splitext(os.path.basename(fname_ent))[0]), totMassResults[1] - relMass))
    log.info('Total entrained mass in sim %s is: %.1f kg' %
             (int(os.path.splitext(os.path.basename(fname_ent))[0]), entMass))
    if (totMassResults[1] - relMass) == 0:
        diff = np.abs((totMassResults[1] - relMass) - entMass)
        if diff > 0:
            log.warning('Conservation of mass is not satisfied')
            log.warning('Total mass change and total entrained mass differ from %.4f kg' % (diff))
        else:
            log.info('Total mass change and total entrained mass differ from %.4f kg' % (diff))
    else:
        diff = np.abs((totMassResults[1] - relMass) - entMass)/(totMassResults[1] - relMass)
        if diff*100 > 0.05:
            log.warning('Conservation of mass is not satisfied')
            log.warning('Total mass change and total entrained mass differ from %.4f %%' % (diff*100))
        else:
            log.info('Total mass change and total entrained mass differ from %.4f %%' % (diff*100))
#   growth results
    growthIndex = totMassResults[1]/totMassResults[0]
    growthGrad = (totMassResults[1] - totMassResults[0]) / (timeResults[1] - timeResults[0])
    return relMass, entMass, finalMass, growthIndex, growthGrad


def getMaxMeanValues(rasterdataA, rasterArea, pLim, cInd=None):
    """Compute average, max values in each cross section for a given input raster

    Read mass balance files to get mass properties of the simulation
    (total mass, entrained mass...). Checks for mass conservation

    Parameters
    ----------
    rasterdataA: 2D numpy array
        raster data
    rasterArea: 2D numpy array
        raster area corresponding to rasterdataA
    pLim: float
        pressure threshold
    cInd: dict
        indexes of bounds of the avalanche (from release to run-out)

    Returns
    -------
    ama: float
        average maximum of rasterdataA
    mma: float
        maximum maximum of rasterdataA
    cInd: dict
        indexes of bounds of the avalanche (from release to run-out)
    aCrossMax: 1D numpy array
        max of rasterdataA in each cross section
    """
    # get mean max for each cross section for A field
    rasterArea = np.where(np.where(np.isnan(rasterdataA), 0, rasterdataA) > 0, rasterArea, 0)
    AreaSum = np.nansum(rasterArea, axis=1)
    AreaSum = np.where(AreaSum > 0, AreaSum, 1)
    aCrossMean = np.nansum(rasterdataA*rasterArea, axis=1)/AreaSum
    # aCrossMean = np.nanmean(rasterdataA, axis=1)
    aCrossMax = np.nanmax(rasterdataA, 1)
    # also get the Area corresponding to those cells
    indACrossMax = np.nanargmax(rasterdataA, 1)
    ind1 = np.arange(np.shape(rasterdataA)[0])
    AreaACrossMax = rasterArea[ind1, indACrossMax]

    #   Determine runout according to maximum and averaged values
    # search in max values
    if not cInd:
        lindex = np.nonzero(aCrossMax > pLim)[0]
        if lindex.any():
            cupper = min(lindex)
            clower = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %4.2f, too high?' % pLim)
            cupper = 0
            clower = 0
        # search in mean values
        lindex = np.nonzero(aCrossMean > pLim)[0]
        if lindex.any():
            cupperm = min(lindex)
            clowerm = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %4.2f, too high?' % pLim)
            cupperm = 0
            clowerm = 0
        cInd = {}
        cInd['cupper'] = cupper
        cInd['clower'] = clower
        cInd['cupperm'] = cupperm
        cInd['clowerm'] = clowerm
    else:
        cupper = cInd['cupper']
        clower = cInd['clower']
        cupperm = cInd['cupperm']
        clowerm = cInd['clowerm']

    # Mean max of of a in each Cross-Section
    ama = np.nansum((aCrossMax*AreaACrossMax)[cupper:clower+1]) / \
        np.nansum(AreaACrossMax[cupper:clower+1])
    mma = max(aCrossMax[cupper:clower+1])

    return ama, mma, cInd, aCrossMax
