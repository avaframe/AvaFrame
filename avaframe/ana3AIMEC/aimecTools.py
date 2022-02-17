"""
    Main logic for AIMEC post processing
"""

import logging
import pathlib
import numpy as np
import pandas as pd
import copy


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
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


def fetchReferenceSimNo(avaDir, inputsDF, comModule, cfgSetup):
    """ Define reference simulation used for aimec analysis.

        if com1DFA is used and a varParList is provided, the simulations
        are ordered and the referenceSimValue is used to define the reference.

        otherwise, if a referenceSimName is provided, the reference is based on this name

        otherwise, the first file found is used as reference

        Parameters
        -----------
        avaDir : str
            path to avalanche directory
        inputsDF: dataFrame
            dataFrame with simulations to analyze and path to result files
        comModule: str
            computational module used to produce the results to analyze
        cfgSetup: configParser object
            configuration for aimec - referenceSimValue, varParList used here

        Returns
        --------
        refSimulation: str
            name of the simulation used as reference
        inputsDF:dataFrame
            dataFrame with simulations to analyze and path to result files
            If com1DFA = comModule and a variation parameter was specified, the
            com1DFA configuration is merged to the inputsDF
        colorVariation: boolean
            True if a color variation should be applied in the plots
    """
    inputDir = pathlib.Path(avaDir, 'Outputs', comModule, 'peakFiles')
    if inputDir.is_dir() == False:
        message = 'Input directory %s does not exist - check anaMod' % inputDir
        log.error(message)
        raise FileNotFoundError(message)
    # if the simulations come from com1DFA, it is possible to order the files and define a reference
    # if com1DFA check for configuration files to fetch parameter values for ordering
    if comModule == 'com1DFA' and cfgSetup['varParList'] != '':
        # fetch parameters that shall be used for ordering
        varParList = cfgSetup['varParList'].split('|')

        # create dataFrame with ordered configurations
        configurationDF = cfgUtils.orderSimFiles(avaDir, inputDir, varParList, cfgSetup.getboolean('ascendingOrder'),
                                                 resFiles=False)
        # Merge inputsDF with the configurationDF. Make sure to keep the indexing from inputs and to merge on 'simName'
        inputsDF = inputsDF.reset_index().merge(configurationDF, on='simName').set_index('index')

        # add value of first parameter used for ordering for colorcoding in plots
        sortingParameter = inputsDF[varParList[0]].to_list()
        typeCP = type(sortingParameter[0])
        if cfgSetup['referenceSimValue'] != '':
            typeCP = type(sortingParameter[0])
            if typeCP == str:
                sortingValues = [x.lower() for x in sortingParameter]
                indexRef = sortingValues.index(typeCP(cfgSetup['referenceSimValue'].lower()))
                valRef = sortingParameter[indexRef]
            elif typeCP in [float, int]:
                colorValues = np.asarray(sortingParameter)
                indexRef = (np.abs(colorValues - typeCP(cfgSetup['referenceSimValue']))).argmin()
                valRef = colorValues[indexRef]
            else:
                indexRef = sortingParameter.index(typeCP(cfgSetup['referenceSimValue']))
                valRef = sortingParameter[indexRef]
            log.info('Reference Simulation is based on %s = %s - closest value found is: %s' %
                     (varParList[0], cfgSetup['referenceSimValue'], str(valRef)))
            refSimulation = inputsDF[inputsDF[varParList[0]] == valRef]['simName'].to_list()[0]
            colorVariation = True
        else:
            # reference simulation
            refSimulation = inputsDF.iloc[0]['simName']  # inputsDF.head(1)['simName'].values[0]
            colorVariation = False

    elif cfgSetup['referenceSimName'] != '':
        colorVariation = False
        # if no colorVariation info and refeenceSimValue available but referenceSimName is given - set
        # simulation with referenceSimName in name as referene simulation
        simFound = False
        for inputIndex, inputsDFrow in inputsDF.iterrows():
            if cfgSetup['referenceSimName'] in str(inputsDFrow['simName']):
                refSimulation = inputsDFrow['simName']
                log.info('Reference Simulation is based on provided simName: %s' % cfgSetup['referenceSimName'])
                simFound = True
                break
        if not simFound:
            refSimulation = inputsDF.iloc[0]['simName']
            log.info('Reference Simulation is based on first simulation in folder')
    else:
        colorVariation = False
        refSimulation = inputsDF.iloc[0]['simName']
        log.info('Reference Simulation is based on first simulation in folder')

    return refSimulation, inputsDF, colorVariation


def computeCellSizeSL(cfgSetup, refResultHeader):
    """ Get the new (s, l) coordinate cell size
        read by default the reference result file cell size.
        If a 'cellSizeSL' is specified in cfgSetup then use this one

        Parameters
        -----------
        refResultHeader: dict
            dictionary wiht the raster header
        cfgSetup: configParser object
            configuration for aimec - with field cellSizeSL

        Returns
        --------
        cellSizeSL: float
            cell size to be used for the (s, l) coordinates
    """
    if cfgSetup['cellSizeSL'] == '':
        cellSizeSL = refResultHeader['cellsize']
        log.info('cellSizeSL is read from the reference header and is : %.2f m' % cellSizeSL)
    else:
        try:
            cellSizeSL = cfgSetup.getfloat('cellSizeSL')
            log.info('cellSizeSL is read from the configuration file and is : %.2f m' % cellSizeSL)
        except ValueError:
            message = ('cellSizeSL is read from the configuration file but should be a number, you provided: %s'
                       % cfgSetup['cellSizeSL'])
            log.error(message)
            raise ValueError(message)

    return cellSizeSL


def makeDomainTransfo(pathDict, inputsDF, cfgSetup):
    """ Make domain transformation

    This function returns the information about the domain transformation
    Data given on a regular grid is projected on a nonuniform grid following
    a polyline to end up with "straightend raster"

    Parameters
    ----------
    pathDict : dict
        dictionary with paths to dem and lines for Aimec analysis
    inputsDF : dataFrame
        dataframe with simulations to analyze and associated path to raster data
    cfgSetup : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini
        regarding domain transformation (domain width w, and new cellsize,
        startOfRunoutAreaAngle or interpolation method, resType and
        referenceFile to get header info)

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
    refSimulation = pathDict['refSimulation']
    refResultSource = inputsDF[inputsDF['simName']==refSimulation][cfgSetup['resType']].to_list()[0]
    ProfileLayer = pathDict['profileLayer']
    splitPointSource = pathDict['splitPointSource']
    DefaultName = pathDict['projectName']

    w = cfgSetup.getfloat('domainWidth')
    startOfRunoutAreaAngle = cfgSetup.getfloat('startOfRunoutAreaAngle')

    log.debug('Data-file %s analysed and domain transformation done' % demSource)

    # read dem and reference result file
    dem = IOf.readRaster(demSource)
    refResult = IOf.readRaster(refResultSource)
    # get the cell size for the (s, l) raster
    cellSizeSL = computeCellSizeSL(cfgSetup, refResult['header'])
    # Initialize transformation dictionary
    rasterTransfo = {}
    rasterTransfo['domainWidth'] = w
    rasterTransfo['cellSizeSL'] = cellSizeSL

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

    ##########################################################################
    # put back the scale due to the desired cellsize
    rasterTransfo['s'] = rasterTransfo['s']*cellSizeSL
    rasterTransfo['l'] = rasterTransfo['l']*cellSizeSL
    rasterTransfo['gridx'] = rasterTransfo['gridx']*cellSizeSL
    rasterTransfo['gridy'] = rasterTransfo['gridy']*cellSizeSL
    rasterTransfo['rasterArea'] = rasterTransfo['rasterArea']*cellSizeSL*cellSizeSL
    # (x,y) coordinates of the resampled avapth (centerline where l = 0)
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
    cellSize = rasterTransfo['cellSizeSL']
    # number of points describing the avaPath
    n_pnt = np.shape(rasterTransfo['DBXr'])[0]
    # Working with no dimentions
    # (the cellsize scaling will be readded at the end)
    # lcoord is the distance from the polyline (cross section)
    # the maximum step size should be smaller then the cellsize
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

    Assign value to the points of the new raster (after domain transormation)

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
    log.debug('Data-file: %s - %d raster values transferred - %d out of original raster bounds!' % (name, iib-ioob, ioob))

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


def analyzeMass(fnameMass, simName, refSimName, resAnalysisDF, time=None):
    """ Analyse Mass data

    Parameters
    ----------
    fnameMass: str
        path to mass data to analyse
    simName: str
        simulation Name
    refSimName: str
        reference simulation Name
    resAnalysisDF: dataFrame
        results from Aimec Analysis to update with mass info
    time: None or 1D numpy array
        None at the first call of the function, then a 1D array
        that will be used to analyze all other mass results

    Returns
    -------
    resAnalysisDF: dataFrame
        results from Aimec Analysis updated with mass inifo:
            relMass: float
                release mass
            entMass: float
                entrained mass
            finalMass: float
                final mass
            relativMassDiff: float
                the final mass diff with ref (in %)
            growthIndex: float
                growth index
            growthGrad: float
                growth gradient
            entMassFlowArray: 1D numpy array
                entrained mass function of time
            totalMassArray: 1D numpy array
                total mass function of time
        time: 1D numpy array
            time array for mass analysis
    """
    log.debug('Analyzing mass')
    log.debug('{: <10} {: <10} {: <10}'.format('Sim number ', 'GI ', 'GR '))
    # analyze mass
    releasedMass, entrainedMass, finalMass, grIndex, grGrad, entMassFlow, totalMass, time = readWrite(fnameMass, time)

    resAnalysisDF.loc[simName, 'relMass'] = releasedMass
    resAnalysisDF.loc[simName, 'finalMass'] = finalMass
    resAnalysisDF.loc[simName, 'entMass'] = entrainedMass
    resAnalysisDF.loc[simName, 'growthIndex'] = grIndex
    resAnalysisDF.loc[simName, 'growthGrad'] = grGrad

    releasedMassRef = resAnalysisDF.loc[refSimName, 'relMass']
    finalMassRef = resAnalysisDF.loc[refSimName, 'finalMass']
    relativMassDiff = (finalMass-finalMassRef)/finalMassRef*100
    if not (releasedMass == releasedMassRef):
        log.warning('Release masses differ between simulations!')
    log.info('{: <10} {:<10.4f} {:<10.4f}'.format(*[simName, grIndex, grGrad]))
    resAnalysisDF.loc[simName, 'relativMassDiff'] = relativMassDiff

    if simName == refSimName:
        resAnalysisDF['entMassFlowArray'] = np.nan
        resAnalysisDF['entMassFlowArray'] = resAnalysisDF['entMassFlowArray'].astype(object)
        resAnalysisDF['totalMassArray'] = np.nan
        resAnalysisDF['totalMassArray'] = resAnalysisDF['totalMassArray'].astype(object)

    resAnalysisDF.at[simName, 'entMassFlowArray'] = entMassFlow
    resAnalysisDF.at[simName, 'totalMassArray'] = totalMass
    return resAnalysisDF, time


def computeRunOut(cfgSetup, rasterTransfo, resAnalysisDF, transformedRasters, simName):
    """ Compute runout based on peak field results

    Parameters
    ----------
    cfgSetup: confiParser
        aimec analysis configuration
    rasterTransfo: dict
        transformation information
    resAnalysisDF : dataFrame
        analysis results from aimec containing:
        PResCrossMax: 1D numpy array
            max of the peak result in each cross section
        PResCrossMean: 1D numpy array
            mean of the peak result in each cross section
    transformedRasters: dict
        dict with transformed dem and peak results
    simName: str
        simulation ID

    Returns
    -------
    resAnalysisDF : dataFrame
        result dataFrame updated withfor each simulation the:
            xRunout: float
                x coord of the runout point
                measured from the begining of the path. run-out
                calculated with the MAX result in each cross section
            yRunout: float
                y coord of the runout point
                measured from the begining of the path. run-out
                calculated with the MAX result in each cross section
            sRunout: float
                runout distance measured from the begining of the path.
                run-out calculated with the MAX result in each cross section
            xMeanRunout: float
                x coord of the runout point
                measured from the begining of the path. run-out
                calculated with the MEAN result in each cross section
            yMeanRunout: float
                y coord of the runout point
                measured from the begining of the path. run-out
                calculated with the MEAN result in each cross section
            sMeanRunout: float
                runout distance measured from the begining of the path.
                run-out calculated with the MEAN result in each cross section
            elevRel: float
                elevation of the release area (based on first point with
                peak field > thresholdValue)
            deltaH: float
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
    gridx = rasterTransfo['gridx']
    gridy = rasterTransfo['gridy']

    resType = cfgSetup['resType']
    thresholdValue = cfgSetup.getfloat('thresholdValue')
    transformedDEMRasters = transformedRasters['newRasterDEM']
    PResRasters = transformedRasters['newRaster' + resType.upper()]
    PResCrossMax = resAnalysisDF.loc[simName, resType + 'CrossMax']
    PResCrossMean = resAnalysisDF.loc[simName, resType + 'CrossMean']

    log.debug('Computing runout')
    lindex = np.nonzero(PResCrossMax > thresholdValue)[0]
    if lindex.any():
        cUpper = min(lindex)
        cLower = max(lindex)
    else:
        log.error('No max values > threshold found. threshold = %4.2f, too high?' % thresholdValue)
        cUpper = 0
        cLower = 0
    # search in mean values
    lindex = np.nonzero(PResCrossMean > thresholdValue)[0]
    if lindex.any():
        cUpperm = min(lindex)
        cLowerm = max(lindex)
    else:
        log.error('No average values > threshold found. threshold = %4.2f, too high?' % thresholdValue)
        cUpperm = 0
        cLowerm = 0
    resAnalysisDF.loc[simName, 'sRunout'] = scoord[cLower]
    index = np.nanargmax(PResRasters[cLower, :])
    resAnalysisDF.loc[simName, 'lRunout'] = lcoord[index]
    resAnalysisDF.loc[simName, 'xRunout'] = gridx[cLower, index]
    resAnalysisDF.loc[simName, 'yRunout'] = gridy[cLower, index]
    resAnalysisDF.loc[simName, 'sMeanRunout'] = scoord[cLowerm]
    resAnalysisDF.loc[simName, 'xMeanRunout'] = x[cLowerm]
    resAnalysisDF.loc[simName, 'yMeanRunout'] = y[cLowerm]
    resAnalysisDF.loc[simName, 'elevRel'] = transformedDEMRasters[cUpper, n]
    resAnalysisDF.loc[simName, 'deltaH'] = transformedDEMRasters[cUpper, n] - transformedDEMRasters[cLower, n]

    return resAnalysisDF


def analyzeField(simName, rasterTransfo, transformedRaster, dataType, resAnalysisDF):
    """ Analyse transformed field

    Analyse transformed raster: compute the Max and Mean values in each cross section
    as well as the overall maximum

    Parameters
    ----------
    simName: str
        simulation name
    rasterTransfo: dict
        transformation information
    transformedRaster: 2D numpy array
        raster after transformation
    dataType: str
        type of the data to analyze ('ppr', 'pfd', 'pfv', ...)
    resAnalysisDF: dataFrame
        result dataFrame to be updated

    Returns
    -------
    Updates the resAnalysisDF input dataFrame with:
        -maxaCrossMax: float
            overall maximum
        -aCrossMax: 1D numpy array
            containing the max of the field in each cross section
        -aCrossMean: 1D numpy array
            containing the mean of the field in each cross section
    """
    # read inputs
    rasterArea = rasterTransfo['rasterArea']

    # initialize Arrays
    log.debug('Analyzing %s' % (dataType))
    log.debug('{: <10} {: <10}'.format('Sim number ', 'maxCrossMax '))

    # Max Mean in each Cross-Section for each field
    maxaCrossMax, aCrossMax, aCrossMean = getMaxMeanValues(transformedRaster, rasterArea)
    log.debug('{: <10} {:<10.4f}'.format(*[simName, maxaCrossMax]))

    resAnalysisDF.loc[simName, 'max' + dataType + 'CrossMax'] = maxaCrossMax
    resAnalysisDF.at[simName, dataType + 'CrossMax'] = aCrossMax
    resAnalysisDF.at[simName, dataType + 'CrossMean'] = aCrossMean

    return resAnalysisDF


def analyzeArea(rasterTransfo, resAnalysisDF, simName, newRasters, cfgSetup, pathDict):
    """Compare results to reference.

    Compute True positive, False negative... areas.

    Parameters
    ----------
    rasterTransfo: dict
        transformation information
    resAnalysisDF: dataFrame
        dataFrame containing Aimec results to update
    simName: str
        simulation Name
    newRasters: dict
        dict with tranformed raster for reference and curent simulation
    cfgSetup: confiParser
        numerical value of the limit to use for the runout computation
        as well as the levels for the contour line plot
    pathDict: dict
        path to data dem data and lines for aimec analysis

    Returns
    -------
    resAnalysisDF: dataFrame
        dataFrame containing Aimec results updated with:
            TP: float
                ref = True sim2 = True
            FN: float
                ref = False sim2 = True
            FP: float
                ref = True sim2 = False
            TN: float
                ref = False sim2 = False
    """
    resType = cfgSetup['resType']
    refSimulationName = pathDict['refSimulation']
    cellarea = rasterTransfo['rasterArea']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    thresholdValue = cfgSetup.getfloat('thresholdValue')
    contourLevels = fU.splitIniValueToArraySteps(cfgSetup['contourLevels'])

    # rasterinfo
    nStart = indStartOfRunout
    # inputs for plot
    inputs = {}
    inputs['runoutLength'] = resAnalysisDF.loc[refSimulationName, 'sRunout']
    inputs['refData'] = newRasters['newRefRaster' + resType.upper()]
    inputs['nStart'] = nStart
    inputs['resType'] = resType

    thresholdArray = contourLevels
    thresholdArray = np.append(thresholdArray, thresholdValue)
    inputs['thresholdArray'] = thresholdArray
    inputs['diffLim'] = cfgSetup.getfloat('diffLim')

    rasterdata = newRasters['newRaster' + resType.upper()]

    # take first simulation as reference
    refMask = copy.deepcopy(inputs['refData'])
    # prepare mask for area resAnalysis
    refMask = np.where(np.isnan(refMask), 0, refMask)
    refMask = np.where(refMask <= thresholdValue, 0, refMask)
    refMask = np.where(refMask > thresholdValue, 1, refMask)
    # comparison rasterdata with mask
    log.debug('{: <15} {: <15} {: <15} {: <15} {: <15}'.format('Sim number ', 'TP ', 'FN ', 'FP ', 'TN'))
    compRasterMask = copy.deepcopy(rasterdata)
    # prepare mask for area resAnalysis
    compRasterMask = np.where(np.isnan(compRasterMask), 0, compRasterMask)
    compRasterMask = np.where(compRasterMask <= thresholdValue, 0, compRasterMask)
    compRasterMask = np.where(compRasterMask > thresholdValue, 1, compRasterMask)

    tpInd = np.where((refMask[nStart:] == 1) & (compRasterMask[nStart:] == 1))
    fpInd = np.where((refMask[nStart:] == 0) & (compRasterMask[nStart:] == 1))
    fnInd = np.where((refMask[nStart:] == 1) & (compRasterMask[nStart:] == 0))
    tnInd = np.where((refMask[nStart:] == 0) & (compRasterMask[nStart:] == 0))

    # subareas
    tp = np.nansum(cellarea[tpInd[0] + nStart, tpInd[1]])
    fp = np.nansum(cellarea[fpInd[0] + nStart, fpInd[1]])
    fn = np.nansum(cellarea[fnInd[0] + nStart, fnInd[1]])
    tn = np.nansum(cellarea[tnInd[0] + nStart, tnInd[1]])

    # take reference (first simulation) as normalizing area
    areaSum = tp + fn

    if areaSum > 0:
        log.debug('{: <15} {:<15.4f} {:<15.4f} {:<15.4f} {:<15.4f}'.format(
                  *[simName, tp/areaSum, fn/areaSum, fp/areaSum, tn/areaSum]))
    if tp + fp == 0:
        log.warning('Simulation %s did not reach the run-out area for threshold %.2f %s' %
                    (simName, thresholdValue, pU.cfgPlotUtils['unit' + cfgSetup['resType']]))

    resAnalysisDF.loc[simName, 'TP'] = tp
    resAnalysisDF.loc[simName, 'TN'] = tn
    resAnalysisDF.loc[simName, 'FP'] = fp
    resAnalysisDF.loc[simName, 'FN'] = fn
    # inputs for plot
    inputs['compData'] = rasterdata
    # masked data for the dataThreshold given in the ini file
    inputs['refRasterMask'] = refMask
    inputs['compRasterMask'] = compRasterMask
    inputs['simName'] = simName

    compPlotPath = None
    if simName != refSimulationName:
        # only plot comparisons of simulations to reference
        compPlotPath = outAimec.visuComparison(rasterTransfo, inputs, pathDict)

    return resAnalysisDF, compPlotPath


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
    if time is None:
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
