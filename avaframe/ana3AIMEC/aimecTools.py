"""
    Main logic for AIMEC post processing
"""

import logging
import pathlib
import numpy as np
import copy


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.outAIMEC as outAimec
import avaframe.out3Plot.plotUtils as pU

# create local logger
log = logging.getLogger(__name__)

# -----------------------------------------------------------
# Aimec read inputs tools
# -----------------------------------------------------------


def readAIMECinputs(avalancheDir, pathDict, dirName='com1DFA'):
    """ Read inputs for AIMEC postprocessing

    Reads the required geometry files location for AIMEC postpocessing
    given an avalanche directory; avalanche path, split point and DEM

    Parameters
    ----------
    avalancheDir : str
        path to directory of avalanche to analyze
    pathDict: dict
        dictionary with paths to simulation results that shall be analaysed
    dirName: str
        name of results directory (e.g. comModule, simName, ...)

    Returns
    -------
    pathDict : dict
        updated dictionary with path to geometry input data
    """

    refDir = pathlib.Path(avalancheDir, 'Inputs', 'LINES')
    profileLayer = list(refDir.glob('*aimec*.shp'))
    try:
        message = ('There should be exactly one path_aimec.shp file containing the avalanche path in %s/Inputs/LINES/' %
            avalancheDir)
        assert len(profileLayer) == 1, message
    except AssertionError:
        raise
    pathDict['profileLayer'] = profileLayer[0]

    refDir = pathlib.Path(avalancheDir, 'Inputs', 'POINTS')
    splitPointLayer = list(refDir.glob('*.shp'))
    try:
        message = 'There should be exactly one .shp file containing the split point in %s/Inputs/POINTS/' % avalancheDir
        assert len(splitPointLayer) == 1, message
    except AssertionError:
        raise
    pathDict['splitPointSource'] = splitPointLayer[0]

    refDir = pathlib.Path(avalancheDir, 'Inputs')
    demSource = list(refDir.glob('*.asc'))
    try:
        assert len(demSource) == 1, 'There should be exactly one topography .asc file in %s/Inputs/' % avalancheDir
    except AssertionError:
        raise
    pathDict['demSource'] = demSource[0]

    pathResult = pathlib.Path(avalancheDir, 'Outputs', 'ana3AIMEC', dirName)
    pathDict['pathResult'] = pathResult

    projectName = pathlib.Path(avalancheDir).stem
    pathDict['projectName'] = projectName
    pathName = profileLayer[0].stem
    pathDict['pathName'] = pathName

    return pathDict


def fetchReferenceSimNo(pathDict, cfgSetup):
    """ Define reference simulation used for aimec analysis from configuration file,
        if no ordering is performed - set simulation 0 as reference simulation

        Parameters
        -----------
        pathDict: dict
            dictionary wiht paths to result files for given result types (ppr, pfd, ..)
            optional key colorParameter - gives ordering of simulations
        cfgSetup: configParser object
            configuration for aimec - referenceSimValue, varParList used here

        Returns
        --------
        pathDict: dict
            updated pathDict with key: referenceFile to define reference simulation for aimec
    """

    if pathDict['colorParameter'] != [] and cfgSetup['referenceSimValue'] != '':
        typeCP = type(pathDict['colorParameter'][0])
        if typeCP == str:
            colorValues = [x.lower() for x in pathDict['colorParameter']]
            indexRef = colorValues.index(typeCP(cfgSetup['referenceSimValue'].lower()))
        elif typeCP in [float, int]:
            colorValues = np.asarray(pathDict['colorParameter'])
            indexRef = (np.abs(colorValues - typeCP(cfgSetup['referenceSimValue']))).argmin()
        else:
            indexRef = pathDict['colorParameter'].index(typeCP(cfgSetup['referenceSimValue']))
        pathDict['referenceFile'] = indexRef
        log.info('Reference Simulation is based on %s = %s - closest value found is: %s' %
                 (cfgSetup['varParList'].split('|')[0], cfgSetup['referenceSimValue'], str(colorValues[indexRef])))
    else:
        pathDict['referenceFile'] = 0
        log.info('Reference Simulation is based on first simulation in folder')

    return pathDict


def makeDomainTransfo(pathDict, cfgSetup):
    """ Make domain transformation

    This function returns the information about the domain transformation
    Data given on a regular grid is projected on a nonuniform grid following
    a polyline to end up with "straightend raster"

    Parameters
    ----------
    pathDict : dict
        dictionary with path to data to analyze
    cfgSetup : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini
        regarding domain transformation (domain width w, startOfRunoutAreaAngle or
        interpolation method, resType and referenceFile to get header info)

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
            indStartOfRunout: int
                index for start of the runout area (in s)
    """
    # Read input parameters
    demSource = pathDict['demSource']
    ProfileLayer = pathDict['profileLayer']
    splitPointSource = pathDict['splitPointSource']
    DefaultName = pathDict['projectName']

    w = float(cfgSetup['domainWidth'])
    startOfRunoutAreaAngle = float(cfgSetup['startOfRunoutAreaAngle'])

    # get info on reference result file to be analysed - to get info on xllcenter, yllcenter, cellSize
    nRef = pathDict['referenceFile']
    rasterSource = pathDict[cfgSetup['resType']][nRef]
    log.debug('Data-file %s analysed and domain transformation done' % demSource)

    # read data header from result file and dem data
    dem = IOf.readRaster(demSource)
    rasterResult = IOf.readRaster(rasterSource)
    header = rasterResult['header']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    cellSize = header['cellsize']
    rasterdata = rasterResult['rasterData']
    # Initialize transformation dictionary
    rasterTransfo = {}
    rasterTransfo['domainWidth'] = w
    rasterTransfo['xllc'] = xllc
    rasterTransfo['yllc'] = yllc
    rasterTransfo['cellSize'] = cellSize

    # read avaPath
    avaPath = shpConv.readLine(ProfileLayer, DefaultName, dem)
    # read split point
    splitPoint = shpConv.readPoints(splitPointSource, dem)
    # add 'z' coordinate to the avaPath
    avaPath, _ = geoTrans.projectOnRaster(dem, avaPath)
    # reverse avaPath if necessary
    _, avaPath = geoTrans.checkProfile(avaPath, projSplitPoint=None)

    log.debug('Creating new raster along polyline: %s' % ProfileLayer)

    # Get new Domain Boundaries DB
    # input: ava path
    # output: Left and right side points for the domain
    rasterTransfo = geoTrans.path2domain(avaPath, rasterTransfo)

    # Make transformation matrix
    rasterTransfo = makeTransfoMat(rasterTransfo)

    # calculate the real area of the new cells as well as the scoord
    rasterTransfo = getSArea(rasterTransfo)

    log.debug('Size of rasterdata- old: %d x %d - new: %d x %d' % (
        np.size(rasterdata, 0), np.size(rasterdata, 1),
        np.size(rasterTransfo['gridx'], 0),
        np.size(rasterTransfo['gridx'], 1)))

    ##########################################################################
    rasterTransfo['header'] = header
    # put back scale and origin
    rasterTransfo['s'] = rasterTransfo['s']*cellSize
    rasterTransfo['l'] = rasterTransfo['l']*cellSize
    rasterTransfo['gridx'] = rasterTransfo['gridx']*cellSize + xllc
    rasterTransfo['gridy'] = rasterTransfo['gridy']*cellSize + yllc
    rasterTransfo['rasterArea'] = rasterTransfo['rasterArea']*cellSize*cellSize
    # (x,y) coordinates of the resamples avapth (centerline where l = 0)
    n = np.shape(rasterTransfo['l'])[0]
    indCenter = int(np.floor(n/2))
    rasterTransfo['x'] = rasterTransfo['gridx'][:, indCenter]
    rasterTransfo['y'] = rasterTransfo['gridy'][:, indCenter]

    #################################################################
    # add 'z' coordinate to the centerline
    rasterTransfo, _ = geoTrans.projectOnRaster(dem, rasterTransfo)
    # find projection of split point on the centerline centerline
    projPoint = geoTrans.findSplitPoint(rasterTransfo, splitPoint)
    rasterTransfo['indSplit'] = projPoint['indSplit']
    # prepare find start of runout area points
    angle, tmp, ds = geoTrans.prepareAngleProfile(startOfRunoutAreaAngle, rasterTransfo)
    # find the runout point: first point under startOfRunoutAreaAngle
    indStartOfRunout = geoTrans.findAngleProfile(tmp, ds, cfgSetup.getfloat('dsMin'))
    rasterTransfo['indStartOfRunout'] = indStartOfRunout
    rasterTransfo['xBetaPoint'] = rasterTransfo['x'][indStartOfRunout]
    rasterTransfo['yBetaPoint'] = rasterTransfo['y'][indStartOfRunout]
    rasterTransfo['startOfRunoutAreaAngle'] = angle[indStartOfRunout]
    log.info('Start of run-out area at the %.2f ° point of coordinates (%.2f, %.2f)' %
        (rasterTransfo['startOfRunoutAreaAngle'], rasterTransfo['xBetaPoint'], rasterTransfo['yBetaPoint']))

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
            cellSize: float
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
    cellSize = rasterTransfo['cellSize']
    # number of points describing the avaPath
    n_pnt = np.shape(rasterTransfo['DBXr'])[0]
    # Working with no dimentions
    # (the cellsize scaling will be readded at the end)
    # lcoord is the distance from the polyline (cross section)
    # maximum step should be smaller then the cellsize
    nTot = np.ceil(w/cellSize)
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
    """ Get the s curvilinear coordinate and area on the new raster

    Find the scoord corresponding to the transformation and the Area of
    the cells of the new raster

    Parameters
    ----------
    rasterTransfo: dict
        dictionary containing:
            domainWidth: float
            cellSize: float
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

    area = np.abs(x1*y2-y1*x2 + x2*y3-y2*x3 + x3*y4-y3*x4 + x4*y1-y4*x1)/2

    # Method 2
    # calculate area of each cell assuming they are parallelogramms
    # (which is wrong)
    # newAreaRaster = np.abs(dxl*dys - dxs*dyl)

    # save Area matrix
    rasterTransfo['rasterArea'] = area
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
    name = fname.name
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
    Points, ioob = geoTrans.projectOnRaster(data, Points, interp=interpMethod)
    newData = Points['z'].reshape(n, m)
    log.debug('Data-file: %s - %d raster values transferred - %d out of original raster bounds!' %
             (name, iib-ioob, ioob))

    return newData


def assignData(fnames, rasterTransfo, interpMethod):
    """ Transfer data from old raster to new raster

    Loops through paths in fnames and calls transfom
    Transform affects values to the points of the new raster
    (after domain transormation)

    Parameters
    ----------
    fnames: list of str
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

    log.debug('Transfer data of %d file(s) from old to new raster' % maxtopo)
    for i in range(maxtopo):
        fname = fnames[i]
        avalData[i] = transform(fname, rasterTransfo, interpMethod)

    return avalData


def analyzeMass(fnameMass, resAnalysis):
    """ Analyse Mass data

    Parameters
    ----------
    fnameMass: list
        list of path to mass data to analyse
    resAnalysis: dict
        resAnalysis dictionary to update with mass inifo

    Returns
    -------
    resAnalysis: dict
        resAnalysis dictionary to updated with mass inifo:
            relMass: 1D numpy array
                containing for each simulation analyzed the
                release mass
            entMass: 1D numpy array
                containing for each simulation analyzed the
                entrained mass
            finalMass: 1D numpy array
                containing for each simulation analyzed the
                final mass
            relativMassDiff: 1D numpy array
                containing for each simulation analyzed
                the final mass diff with ref (in %)
            growthIndex: 1D numpy array
                containing for each simulation analyzed the
                growth index
            growthGrad: 1D numpy array
                containing for each simulation analyzed the
                growth gradient
    """
    # initialize Arrays
    nTopo = len(fnameMass)
    massDiffers = False
    grIndex = np.empty((nTopo))
    grGrad = np.empty((nTopo))
    releaseMass = np.empty((nTopo))
    entrainedMass = np.empty((nTopo))
    finalMass = np.empty((nTopo))
    relativMassDiff = np.empty((nTopo))
    time = 0

    log.debug('Analyzing mass')
    log.debug('{: <10} {: <10} {: <10}'.format('Sim number ', 'GI ', 'GR '))
    # For each data set
    for i in range(nTopo):
        # analyze mass
        releaseMass[i], entrainedMass[i], finalMass[i], grIndex[i], grGrad[i], entMassFlow, totalMass, time = readWrite(
            fnameMass[i], time)
        if i == 0:
            entMassFlowArray = np.zeros((nTopo, np.size(time)))
            totalMassArray = np.zeros((nTopo, np.size(time)))
        entMassFlowArray[i] = entMassFlow
        totalMassArray[i] = totalMass
        relativMassDiff[i] = (finalMass[i]-finalMass[0])/finalMass[0]*100
        if not (releaseMass[i] == releaseMass[0]):
            massDiffers = True
        log.info('{: <10} {:<10.4f} {:<10.4f}'.format(*[i+1, grIndex[i], grGrad[i]]))
    if massDiffers:
        log.warning('Release masses differs between simulations!')

    resAnalysis['relMass'] = releaseMass
    resAnalysis['entMass'] = entrainedMass
    resAnalysis['entMassFlowArray'] = entMassFlowArray
    resAnalysis['totalMassArray'] = totalMassArray
    resAnalysis['time'] = time
    resAnalysis['finalMass'] = finalMass
    resAnalysis['relativMassDiff'] = relativMassDiff
    resAnalysis['growthIndex'] = grIndex
    resAnalysis['growthGrad'] = grGrad

    return resAnalysis


def computeRunOut(rasterTransfo, thresholdValue, resultsAreaAnalysis, transformedDEMRasters):
    """ Compute runout based on peak field results

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    thresholdValue: float
        numerical value of the threshold limit to use
    resultsAreaAnalysis : dict
        PResCrossMax: 2D numpy array
            containing for each simulation analyzed the
            max of the peak result in each cross section
        PResCrossMean: 2D numpy array
            containing for each simulation analyzed the
            mean of the peak result in each cross section

    Returns
    -------
    runout: 2D numpy array
        containing for each simulation analyzed the x and
        y coord of the runout point as well as the runout distance
        measured from the begining of the path. run-out
        calculated with the MAX result in each cross section
    runoutMean: 2D numpy array
        containing for each simulation analyzed the x
        and y coord of the runout point as well as the runout
        distance measured from the begining of the path.
        run-out calculated with the MEAN result in each cross
        section
    elevRel: 1D numpy array
        containing for each simulation analyzed the
        elevation of the release area (based on first point with
        peak field > thresholdValue)
    deltaH: 1D numpy array
        containing for each simulation analyzed the
        elevation fall difference between elevRel and altitude of
        run-out point
    """
    # read inputs
    scoord = rasterTransfo['s']
    lcoord = rasterTransfo['l']
    n = np.shape(lcoord)[0]
    n = int(np.floor(n/2)+1)
    x = rasterTransfo['x']
    y = rasterTransfo['y']

    resType = resultsAreaAnalysis['resType']
    PResCrossMax = resultsAreaAnalysis[resType]['aCrossMax']
    PResCrossMean = resultsAreaAnalysis[resType]['aCrossMean']

    # initialize Arrays
    nTopo = len(PResCrossMax)
    runout = np.empty((3, nTopo))
    runoutMean = np.empty((3, nTopo))
    elevRel = np.empty((nTopo))
    deltaH = np.empty((nTopo))

    log.debug('Computig runout')
    # For each data set
    for i in range(nTopo):
        lindex = np.nonzero(PResCrossMax[i] > thresholdValue)[0]
        if lindex.any():
            cupper = min(lindex)
            clower = max(lindex)
        else:
            log.error('No max values > threshold found. threshold = %4.2f, too high?' % thresholdValue)
            cupper = 0
            clower = 0
        # search in mean values
        lindex = np.nonzero(PResCrossMean[i] > thresholdValue)[0]
        if lindex.any():
            cupperm = min(lindex)
            clowerm = max(lindex)
        else:
            log.error('No average values > threshold found. threshold = %4.2f, too high?' % thresholdValue)
            cupperm = 0
            clowerm = 0
        cInd = {}
        cInd['cupper'] = cupper
        cInd['clower'] = clower
        cInd['cupperm'] = cupperm
        cInd['clowerm'] = clowerm
        #    Runout
        cupper = cInd['cupper']
        clower = cInd['clower']
        clowerm = cInd['clowerm']
        runout[0, i] = scoord[clower]
        runout[1, i] = x[clower]
        runout[2, i] = y[clower]
        runoutMean[0, i] = scoord[clowerm]
        runoutMean[1, i] = x[clower]
        runoutMean[2, i] = y[clower]
        elevRel[i] = transformedDEMRasters[cupper, n]
        deltaH[i] = transformedDEMRasters[cupper, n] - transformedDEMRasters[clower, n]

    return runout, runoutMean, elevRel, deltaH


def analyzeField(rasterTransfo, transformedRasters, dataType, resultsAreaAnalysis):
    """ Analyse transformed field

    Analyse transformed rasters in transformedRasters and for each one, compute
    the Max and Mean values in each cross section, as well as the
    overall maximum

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    transformedRasters: list
        list containing rasters after transformation
    dataType: str
        type of the data to analyze ('ppr', 'pfd' or 'pfv')
    resultsAreaAnalysis: dict
        result dictionary to be updated

    Returns
    -------
    Updates the resultsAreaAnalysis input dictionary with a sub dictionary
    of the name dataType containing:
        -maxaCrossMax: 1D numpy array
            containing for each simulation analyzed the overall maximum
        -aCrossMax: 2D numpy array
            containing for each simulation analyzed the
            max of the field in each cross section
        -aCrossMean: 2D numpy array
            containing for each simulation analyzed the
            mean of the field in each cross section
    """
    # read inputs
    scoord = rasterTransfo['s']
    rasterArea = rasterTransfo['rasterArea']

    # initialize Arrays
    nTopo = len(transformedRasters)
    maxaCrossMax = np.empty((nTopo))
    aCrossMax = np.zeros((nTopo, len(scoord)))
    aCrossMean = np.zeros((nTopo, len(scoord)))
    log.debug('Analyzing %s' % (dataType))
    log.debug('{: <10} {: <10}'.format('Sim number ', 'maxCrossMax '))
    # For each data set
    for i in range(nTopo):
        rasterData = transformedRasters[i]

        # Max Mean in each Cross-Section for each field
        maxaCrossMax[i], aCrossMax[i], aCrossMean[i] = getMaxMeanValues(rasterData, rasterArea)
        log.debug('{: <10} {:<10.4f}'.format(*[i+1, maxaCrossMax[i]]))

    resultsAreaAnalysis[dataType] = {'transformedRasters': transformedRasters,
                                     'maxaCrossMax': maxaCrossMax,
                                     'aCrossMax': aCrossMax,
                                     'aCrossMean': aCrossMean}

    return resultsAreaAnalysis


def analyzeArea(rasterTransfo, runoutLength, data, cfgSetup, pathDict, cfgFlags):
    """Compare results to reference.

    Compute True positive, False negative... areas.

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    resAnalysis: dict
        resAnalysis dictionary containing all results to update
    data: list
        list of transformed rasters
    cfgSetup: dict
        numerical value of the limit to use for the runout computation
        as well as the levels for the contour line plot
    pathDict: dict
        path to data to analyse
    cfgFlags: configparser
        configparser with plot and write flags

    Returns
    -------
    TP: float
        ref = True sim2 = True
    FN: float
        ref = False sim2 = True
    FP: float
        ref = True sim2 = False
    TN: float
        ref = False sim2 = False
    """
    nRef = pathDict['referenceFile']
    cellarea = rasterTransfo['rasterArea']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    thresholdValue = cfgSetup.getfloat('thresholdValue')
    contourLevels = fU.splitIniValueToArraySteps(cfgSetup['contourLevels'])

    # initialize Arrays
    nTopo = len(data)
    TP = np.empty((nTopo))
    FN = np.empty((nTopo))
    FP = np.empty((nTopo))
    TN = np.empty((nTopo))

    # rasterinfo
    nStart = indStartOfRunout
    # inputs for plot
    inputs = {}
    inputs['runoutLength'] = runoutLength
    inputs['refData'] = data[nRef]
    inputs['nStart'] = nStart
    inputs['resType'] = cfgSetup['resType']

    thresholdArray = contourLevels
    thresholdArray = np.append(thresholdArray, thresholdValue)
    inputs['thresholdArray'] = thresholdArray
    inputs['diffLim'] = cfgSetup.getfloat('diffLim')

    for i in range(nTopo):
        rasterdata = data[i]

        """
        area
        # true positive: result1(mask)=1, result2(rasterdata)=1
        # false negative: result1(mask)=1, result2(rasterdata)=0
        # false positive: result1(mask)=0, result2(rasterdata)=1
        # true negative: result1(mask)=0, result2(rasterdata)=0
        """
        # take first simulation as reference
        refMask = copy.deepcopy(data[nRef])
        # prepare mask for area resAnalysis
        refMask = np.where(np.isnan(refMask), 0, refMask)
        refMask = np.where(refMask < thresholdValue, 0, refMask)
        refMask = np.where(refMask >= thresholdValue, 1, refMask)
        # comparison rasterdata with mask
        log.debug('{: <15} {: <15} {: <15} {: <15} {: <15}'.format(
            'Sim number ', 'TP ', 'FN ', 'FP ', 'TN'))
        newRasterData = copy.deepcopy(rasterdata)
        # prepare mask for area resAnalysis
        newRasterData = np.where(np.isnan(newRasterData), 0, newRasterData)
        newRasterData = np.where(newRasterData < thresholdValue, 0, newRasterData)
        newRasterData = np.where(newRasterData >= thresholdValue, 1, newRasterData)

        tpInd = np.where((refMask[nStart:] == 1) &
                         (newRasterData[nStart:] == 1))
        fpInd = np.where((refMask[nStart:] == 0) &
                         (newRasterData[nStart:] == 1))
        fnInd = np.where((refMask[nStart:] == 1) &
                         (newRasterData[nStart:] == 0))
        tnInd = np.where((refMask[nStart:] == 0) &
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
        if areaSum > 0:
            log.debug('{: <15} {:<15.4f} {:<15.4f} {:<15.4f} {:<15.4f}'.format(
                      *[i+1, tp/areaSum, fn/areaSum, fp/areaSum, tn/areaSum]))
        if tp + fp == 0:
            log.warning('Simulation %s did not reach the run-out area for threshold %.2f %s' %
                        (i, thresholdValue, pU.cfgPlotUtils['unit' + cfgSetup['resType']]))
        # inputs for plot
        inputs['compData'] = rasterdata
        # masked data for the dataThreshold given in the ini file
        inputs['refRasterMask'] = refMask
        inputs['newRasterMask'] = newRasterData
        inputs['i'] = i

        # only plot comparisons of simulations to reference
        if cfgFlags['savePlot']:
            if i != nRef:
                compPlotPath = outAimec.visuComparison(rasterTransfo, inputs, pathDict,
                                        cfgFlags)
            elif pathDict['numSim'] == 1:
                compPlotPath = outAimec.visuComparison(rasterTransfo, inputs, pathDict,
                                        cfgFlags)
                log.warning('only one simulation, area comparison not meaningful')

    return TP, FN, FP, TN, compPlotPath


def readWrite(fname_ent, time):
    """ Get mass balance information
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
        entMassentMassFlow: float
            entrained mass flow (kg/s)
        finalMass: float
            final mass
        growthIndex: float
        growthGrad: float
    """

    #    load data
    #    time, total mass, entrained mass
    massTime = np.loadtxt(fname_ent, delimiter=',', skiprows=1)
    timeSimulation = massTime[:, 0]
    dt = timeSimulation[1:] - timeSimulation[:-1]
    timeResults = [massTime[0, 0], massTime[-1, 0]]
    totMassResults = [massTime[0, 1], massTime[-1, 1]]
    relMass = totMassResults[0]
    if np.size(time) == 1:
        time = np.arange(0, int(timeResults[1]), 0.1)
    entMassFlow = np.interp(time, massTime[1:, 0], massTime[1:, 2]/dt)
    totalMass = np.interp(time, massTime[:, 0], massTime[:, 1])
    entrainedMass = np.sum(massTime[:, 2])
    finalMass = totMassResults[1]

    # check mass balance
    log.info('Total mass change between first and last time step in sim %s is: %.1f kg' %
             (fname_ent.stem, totMassResults[1] - relMass))
    log.info('Total entrained mass in sim %s is: %.1f kg' %
             (fname_ent.stem, entrainedMass))
    if (totMassResults[1] - relMass) == 0:
        diff = np.abs((totMassResults[1] - relMass) - entrainedMass)
        if diff > 0:
            log.warning('Conservation of mass is not satisfied')
            log.warning('Total mass change and total entrained mass differ from %.4f kg' % (diff))
        else:
            log.info('Total mass change and total entrained mass differ from %.4f kg' % (diff))
    else:
        diff = np.abs((totMassResults[1] - relMass) - entrainedMass)/(totMassResults[1] - relMass)
        if diff*100 > 0.05:
            log.warning('Conservation of mass is not satisfied')
            log.warning('Total mass change and total entrained mass differ from %.4f %%' % (diff*100))
        else:
            log.info('Total mass change and total entrained mass differ from %.4f %%' % (diff*100))

    # growth results
    growthIndex = totMassResults[1]/totMassResults[0]
    growthGrad = (totMassResults[1] - totMassResults[0]) / (timeResults[1] - timeResults[0])

    return relMass, entrainedMass, finalMass, growthIndex, growthGrad, entMassFlow, totalMass, time


def getMaxMeanValues(rasterdataA, rasterArea):
    """Compute average, max values in each cross section for a given input raster

    Read mass balance files to get mass properties of the simulation
    (total mass, entrained mass...). Checks for mass conservation

    Parameters
    ----------
    rasterdataA: 2D numpy array
        raster data
    rasterArea: 2D numpy array
        raster area corresponding to rasterdataA

    Returns
    -------
    mma: float
        maximum maximum of rasterdataA
    aCrossMax: 1D numpy array
        max of rasterdataA in each cross section
    aCrossMean: 1D numpy array
        mean of rasterdataA in each cross section (area weighted)
    """
    # get mean max for each cross section for A field
    rasterArea = np.where(np.where(np.isnan(rasterdataA), 0, rasterdataA) > 0, rasterArea, 0)
    areaSum = np.nansum(rasterArea, axis=1)
    areaSum = np.where(areaSum > 0, areaSum, 1)
    aCrossMean = np.nansum(rasterdataA*rasterArea, axis=1)/areaSum
    # aCrossMean = np.nanmean(rasterdataA, axis=1)
    aCrossMax = np.nanmax(rasterdataA, 1)
    # maximum of field a
    maxaCrossMax = np.nanmax(aCrossMax)

    return maxaCrossMax, aCrossMax, aCrossMean
