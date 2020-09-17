"""
    Main logic for AIMEC post processing

    This file is part of Avaframe.
"""

import sys
import os
import time
import logging
import glob
import math
import numpy as np
import scipy as sp
import copy
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.image import NonUniformImage
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in3Utils.ascUtils as IOf
import avaframe.out3SimpPlot.outAIMEC as outAimec
from avaframe.out3SimpPlot.plotSettings import *

# create local logger
log = logging.getLogger(__name__)

# -----------------------------------------------------------
# Aimec read inputs tools
# -----------------------------------------------------------


def readAIMECinputs(avalancheDir, dirName='com1DFA'):
    """
    Reads the requiered files location for AIMEC postpocessing
    given an avalanche directory
    """
    cfgPath = {}
    pathPressure = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName, 'dfa_pressure')
    pathFlowHeight = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName, 'dfa_depth')
    pathMassBalance = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName, 'dfa_mass_balance')
    pathSpeed = os.path.join(avalancheDir, 'Work', 'ana3AIMEC', dirName, 'dfa_speed')

    if not os.path.exists(pathMassBalance):
        os.makedirs(pathMassBalance)

    profileLayer = glob.glob(os.path.join(avalancheDir, 'Inputs', 'LINES', '*aimec*.shp'))
    cfgPath['profileLayer'] = ''.join(profileLayer)

    splitPointLayer = glob.glob(os.path.join(avalancheDir, 'Inputs', 'POINTS', '*.shp'))
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
    """ Get sorted list of all files in folder """
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
    newRasters['newRasterPressure'] = assignData(cfgPath['pressurefileList'], rasterTransfo,
                                   interpMethod)
    # assign depth data
    log.info("Assigning depth data to deskewed raster")
    newRasters['newRasterDepth'] = assignData(cfgPath['depthfileList'], rasterTransfo,
                                interpMethod)
    # assign speed data
    if cfgPath['speedfileList']:
        log.info("Assigning speed data to deskewed raster")
        newRasters['newRasterSpeed'] = assignData(cfgPath['speedfileList'], rasterTransfo,
                                    interpMethod)

    # assign dem data
    log.info("Assigning dem data to deskewed raster")
    newRasterDEM = assignData([cfgPath['demSource']], rasterTransfo,
                              interpMethod)
    newRasters['newRasterDEM'] = newRasterDEM[0]

    # Analyze data
    log.info('Analyzing data in path coordinate system')
    resAnalysis = analyzeData(rasterTransfo, pressureLimit, newRasters, cfgPath, cfgFlags)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    log.info('Visualisation of results')
    outAimec.resultVisu(cfgPath, rasterTransfo, resAnalysis, pressureLimit)

    # -----------------------------------------------------------
    # write results to file
    # -----------------------------------------------------------
    log.info('Writing results to file')

    outAimec.resultWrite(cfgPath, cfgSetup, resAnalysis)

# -----------------------------------------------------------
# Aimec processing tools
# -----------------------------------------------------------


def makeDomainTransfo(cfgPath, cfgSetup, cfgFlags):
    """
    Make domain transformation :
    This function returns the information about the domain transformation
    Data given on a regular grid is projected on a nonuniform grid following
    a polyline

    input: cfgPath, cfgSetup, cfgFlags
    ouput: rasterTransfo as a dictionary
            -gridx: x coord of the new raster points in old coord system (as 2D array)
            -gridy: y coord of the new raster points in old coord system (as 2D array)
            -s: new coord system in the polyline direction (as 1D array)
            -l: new coord system in the cross direction (as 1D array)
            -x: coord of the resampled polyline in old coord system (as 1D array)
            -y: coord of the resampled polyline in old coord system (as 1D array)
            -rasterArea: real area of the cells of the new raster (as 2D array)
            -indRunoutPoint: index for start of the runout area (in s)

    """
    # Read input parameters
    rasterSource = cfgPath['pressurefileList'][0]
    demSource = cfgPath['demSource']
    ProfileLayer = cfgPath['profileLayer']
    outpath = cfgPath['pathResult']
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
    xllc = header.xllcorner
    yllc = header.yllcorner
    cellsize = header.cellsize
    rasterdata = sourceData['rasterData']
    # Initialize transformation dictionary
    rasterTransfo = {}
    rasterTransfo['domainWidth'] = w
    rasterTransfo['xllc'] = xllc
    rasterTransfo['yllc'] = yllc
    rasterTransfo['cellsize'] = cellsize

    # read avaPath
    Avapath = shpConv.readLine(ProfileLayer, DefaultName, sourceData['header'])
    # read split point
    splitPoint = shpConv.readPoints(cfgPath['splitPointSource'], sourceData['header'])
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
        np.size(rasterTransfo['gridx'], 0), np.size(rasterTransfo['gridx'], 1)))

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
    indCenter = int(np.floor(n/2)+1)
    rasterTransfo['x'] = rasterTransfo['gridx'][:, indCenter]
    rasterTransfo['y'] = rasterTransfo['gridy'][:, indCenter]

    #################################################################
    # add 'z' coordinate to the centerline
    rasterTransfo = geoTrans.projectOnRaster(dem, rasterTransfo)
    # find projection of split point on the centerline centerline
    projPoint = geoTrans.findSplitPoint(rasterTransfo, splitPoint)
    rasterTransfo['indSplit'] = projPoint['indSplit']
    # prepare find start of runout area points
    angle, tmp, delta_ind = geoTrans.prepareAngleProfile(runoutAngle, rasterTransfo)
    # find the runout point: first point under runoutAngle
    indRunoutPoint = geoTrans.findAngleProfile(tmp, delta_ind)
    if runoutAngle<angle[indRunoutPoint] and runoutAngle>angle[indRunoutPoint+1]:
        rasterTransfo['indRunoutPoint'] = indRunoutPoint
        rasterTransfo['runoutAngle'] = runoutAngle
        log.info('Measuring run-out length from the %s ° point' % runoutAngle)
    else:
        log.warning('No %s ° point found. Check splitPoint position or runoutAngle value.'% runoutAngle)
        rasterTransfo['indRunoutPoint'] = indRunoutPoint
        rasterTransfo['runoutAngle'] = (angle[indRunoutPoint] + angle[indRunoutPoint+1])/2
        log.info('Measuring run-out length from the %s ° point' % rasterTransfo['runoutAngle'])

    avalData = transform(rasterSource, rasterTransfo, interpMethod)

    ###########################################################################
    # visualisation
    inputData = {}
    inputData['avalData'] = avalData
    inputData['sourceData'] = sourceData
    inputData['Avapath'] = Avapath

    outAimec.visuTransfo(rasterTransfo, inputData, cfgPath, cfgFlags)

    return rasterTransfo


def split_section(DB, i):
    """
    Splits the ith segment of domain boundary DB in the s direction
    (direction of the path)
    input: - DB domain Boundary dictionary
           - i number of the segment of DB to split
    ouput: - (x,y) coordinates of the ith left and right splited Boundaries
            (bxl, byl, bxr, byr). each is a 1D array of size m
            - m number of ellements on the new segments
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
        input: rasterTransfo dictionary containing:
                -domainWidth
                -cellsize
                -DBXl: x coord of the left boundary
                -DBXr: x coord of the right boundary
                -DBYl: y coord of the left boundary
                -DBYr: y coord of the right boundary

        ouput: rasterTransfo dictionary updated with
                -gridx: x coord of the new raster points in old coord system (as 2D array)
                -gridy: y coord of the new raster points in old coord system (as 2D array)
                -l: new coord system in the cross direction (as 1D array)
    """
    w = rasterTransfo['domainWidth']
    cellsize = rasterTransfo['cellsize']
    # number of points describing the avaPath
    n_pnt = np.shape(rasterTransfo['DBXr'])[0]
    # Working with no dimentions (the cellsize scaling will be readded at the end)
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
        bxl, byl, bxr, byr, m = split_section(rasterTransfo, i)
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
                newGridRasterX = np.append(newGridRasterX, x.reshape(1, nTot), axis=0)
                newGridRasterY = np.append(newGridRasterY, y.reshape(1, nTot), axis=0)

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
    """
    Find the scoord corresponding to the transformation and the Area of
    the cells of the new raster
    input: rasterTransfo dictionary to fill in output
            -gridx: x coord of the new raster points in old coord system (as 2D array)
            -gridy: y coord of the new raster points in old coord system (as 2D array)

    ouput: rasterTransfo dictionary updated with
            -s: new coord system in the polyline direction (as 1D array)
            -rasterArea: real area of the cells of the new raster (as 2D array)
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

    # calculate area of each cell
    newAreaRaster = np.abs(dxl*dys - dxs*dyl)
    rasterTransfo['rasterArea'] = newAreaRaster

    # get scoord
    ds = Vs[:, int(np.floor(m/2))-1]
    scoord = np.cumsum(ds)-ds[0]
    rasterTransfo['s'] = scoord

    return rasterTransfo


def transform(fname, rasterTransfo, interpMethod):
    """
    Affect value to the points of the new raster (after domain transormation)
    input:
            -fname = name of rasterfile to transform
            -rasterTransfo = transformation info
            -interpolation method to chose between 'nearest' and 'bilinear'
    ouput:
            -new_data = z, pressure or depth... corresponding to fname on the new raster
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
    Points, iib, ioob = geoTrans.projectOnRasterVect(data, Points, interp=interpMethod)
    newData = Points['z'].reshape(n, m)
    log.info('Data-file: %s - %d raster values transferred - %d out of original raster bounds!' %
             (name, iib-ioob, ioob))

    return newData


def assignData(fnames, rasterTransfo, interpMethod):
    """
    Affect value to the points of the new raster (after domain transormation)
    input:
            -fnames = list of names of rasterfiles to transform
            -rasterTransfo = transformation info
            -interpolation method to chose between 'nearest' and 'bilinear'
    ouput: avalData = z, pressure or depth... corresponding to fnames on the new rasters
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
    """
    Analyse pressure and depth deskewed data
    """

    resAnalysis = analyzeFields(rasterTransfo, pLim, newRasters, cfgPath)

    resAnalysis = analyzeArea(rasterTransfo, resAnalysis, pLim, newRasters, cfgPath, cfgFlags)


    outAimec.visuSimple(rasterTransfo, resAnalysis, newRasters, cfgPath, cfgFlags)
    outAimec.visuRunout(rasterTransfo, resAnalysis, pLim, newRasters, cfgPath, cfgFlags)

    return resAnalysis


def analyzeFields(rasterTransfo, pLim, newRasters, cfgPath):
    """
    Analyse pressure depth and speed.
    Calculate runout, Max Peak Pressure, Average PP...
    Get mass and entrainement
    input:
            -rasterTransfo = transformation info
            -pressure threshold value
            -newRasters dictionnary with the data in the new coord system
            -cfgPath
    output: resAnalysis dictionnary containing all results
            -runout: 2D array containing for each simulation analyzed the x and y coord
                     of the runout point as well as the runout distance measured from
                     the begining of the run-out area. run-out calculated with the MAX pressure in each cross section
            -runoutMean: 2D array containing for each simulation analyzed the x and y coord
                     of the runout point as well as the runout distance measured from
                     the begining of the run-out area. run-out calculated with the MEAN pressure in each cross section
            -AMPP: 1D array containing for each simulation analyzed the average max peak pressure
            -MMPP: 1D array containing for each simulation analyzed the max max peak pressure
            -AMD: 1D array containing for each simulation analyzed the average max peak flow depth
            -MMD: 1D array containing for each simulation analyzed the max max peak flow depth
            -AMS: 1D array containing for each simulation analyzed the average max peak speed
            -MMS: 1D array containing for each simulation analyzed the max max peak speed
            -elevRel: 1D array containing for each simulation analyzed the elevation of
                      the release area (based on first point with peak pressure > pLim)
            -deltaH: 1D array containing for each simulation analyzed the elevation fall
                    difference between elevRel and altitude of run-out point
            -relMass: 1D array containing for each simulation analyzed the release mass
            -entMass: 1D array containing for each simulation analyzed the entrained mass
            -finalMass: 1D array containing for each simulation analyzed the final mass
            -relativMassDiff: 1D array containing for each simulation analyzed the final mass diff with ref (in %)
            -growthIndex: 1D array containing for each simulation analyzed the growth index
            -growthGrad: 1D array containing for each simulation analyzed the growth gradient
            -pCrossAll: 2D array containing for each simulation analyzed the max peak pressure in each cross section
            -pressureLimit: pressure threshold pLim
            -runoutAngle: angle of the slope at the beginning of the run-out area (given in input)

    """
    # read inputs
    fname = cfgPath['pressurefileList']
    fnameMass = cfgPath['massfileList']
    outpath = cfgPath['pathResult']

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
        rasterArea[np.where(np.isnan(rasterdataPres))] = np.nan

        # Mean max in each Cross-Section for each field
        ampp[i], mmpp[i], cInd, pCrossAll[i] = getMaxMeanValues(rasterdataPres, rasterArea, pLim, cInd=None)
        amd[i], mmd[i], cInd, _ = getMaxMeanValues(rasterdataDepth, rasterArea, pLim, cInd=cInd)
        ams[i], mms[i], cInd, _ = getMaxMeanValues(rasterdataSpeed, rasterArea, pLim, cInd=cInd)
        #    Runout
        cupper = cInd['cupper']
        clower = cInd['clower']
        cupperm = cInd['cupperm']
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
        releaseMass[i], entrainedMass[i], finalMass[i], grIndex[i], grGrad[i] = readWrite(fnameMass[i])
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
    """
    Compare results to reference.
    Compute True positive, False negative... areas.
    input:
            -rasterTransfo:transformation info
            -pLim: pressure threshold value
            -newRasters dictionnary with the data in the new coord system
            -cfgPath
            -cfgFlags
    output: resAnalysis dictionnary completed with area information
            -TP: ref = True sim2 = True
            -FN: ref = False sim2 = True
            -FP: ref = True sim2 = False
            -TN: ref = False sim2 = False
    """
    fname = cfgPath['pressurefileList']

    dataPressure = newRasters['newRasterPressure']
    scoord = rasterTransfo['s']
    lcoord = rasterTransfo['l']
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
    newMask = np.where(newMask< pLim, 0, newMask)
    newMask = np.where(newMask>= pLim, 1, newMask)
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
        # for each pressure-file pLim is introduced (1/3/.. kPa), where the avalanche has stopped
        newRasterData = copy.deepcopy(rasterdata)
        # prepare mask for area resAnalysis
        newRasterData = np.where(np.isnan(newRasterData), 0, newRasterData)
        newRasterData = np.where(newRasterData<pLim, 0, newRasterData)
        newRasterData = np.where(newRasterData>=pLim, 1, newRasterData)

        if cfgFlags.getboolean('savePlot') and i > 0:
            # read paths
            pathResult = cfgPath['pathResult']
            projectName = cfgPath['dirName']
            outFileName = projectName + '_' +  str(i) +  '_compToRef'
            outname = os.path.join(pathResult, 'pics', outFileName)
            if not os.path.exists(os.path.dirname(outname)):
                os.makedirs(os.path.dirname(outname))
            fig = plt.figure(figsize=(figW*2, figH))
            y_lim = scoord[indRunoutPoint+20]+np.nanmax(resAnalysis['runout'][0])
        #    for figure: referenz-simulation bei pLim=1
            ax1 = plt.subplot(121)
            ax1.set_title('Reference Peak Presseure in the RunOut area' + '\n' +  'Pressure threshold: %.1f kPa' % pLim)
            cmap = cmapPres
            cmap.set_under(color='w')
            im = NonUniformImage(ax1, extent=[lcoord.min(), lcoord.max(),
                                              scoord.min(), scoord.max()], cmap=cmap)
            im.set_clim(vmin=pLim, vmax=np.nanmax((dataPressure[0])[nStart:]))
            im.set_data(lcoord, scoord, dataPressure[0])
            ref0 = ax1.images.append(im)
            cbar = ax1.figure.colorbar(im, extend='both', ax=ax1, use_gridspec=True)
            cbar.ax.set_ylabel('peak pressure [kPa]')
            ax1.set_xlim([lcoord.min(), lcoord.max()])
            ax1.set_ylim([scoord[indRunoutPoint-20], y_lim])
            ax1.set_xlabel('l [m]')
            ax1.set_ylabel('s [m]')

            ax2 = plt.subplot(122)
            ax2.set_title(
                'Difference between current and reference in the RunOut area' + '\n' +  'Blue = FN, Red = FP')
            colorsList = [[0, 0, 1], [1, 1, 1], [1, 0, 0]]
            cmap = matplotlib.colors.ListedColormap(colorsList)
            cmap.set_under(color='b')
            cmap.set_over(color='r')
            cmap.set_bad(color='k')
            im = NonUniformImage(ax2, extent=[lcoord.min(), lcoord.max(),
                                              scoord.min(), scoord.max()], cmap=cmap)
            im.set_clim(vmin=-0.000000001, vmax=0.000000001)
            im.set_data(lcoord, scoord, newRasterData-newMask)
            ref0 = ax2.images.append(im)
            # cbar = ax2.figure.colorbar(im, ax=ax2, extend='both', use_gridspec=True)
            # cbar.ax.set_ylabel('peak pressure [kPa]')
            ax2.set_xlim([lcoord.min(), lcoord.max()])
            ax2.set_ylim([scoord[indRunoutPoint-20], y_lim])
            ax2.set_xlabel('l [m]')
            ax2.set_ylabel('s [m]')
            plt.subplots_adjust(wspace=0.3)
            # plt.show()

            fig.savefig(outname)
            plt.close(fig)

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
    """
    Read mass balance files to get mass properties of the simulation
    (total mass, entrained mass...). Checks for mass conservation
    input:
            -fname_ent: list of mass balance files
    output:
            -relMass: release mass
            -entMass: entrained mass
            -growthIndex
            -growthGrad
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
             (int(os.path.splitext(os.path.basename(fname_ent))[0]), totMassResults[1]- relMass))
    log.info('Total entrained mass in sim %s is: %.1f kg' % (int(os.path.splitext(os.path.basename(fname_ent))[0]), entMass))
    if (totMassResults[1]- relMass)==0:
        diff = np.abs((totMassResults[1]- relMass) - entMass)
        if diff>0:
            log.warning('Conservation of mass is not satisfied')
            log.warning('Total mass change and total entrained mass differ from %.4f kg' % (diff))
        else:
            log.info('Total mass change and total entrained mass differ from %.4f kg' % (diff))
    else:
        diff = np.abs((totMassResults[1]- relMass) - entMass)/(totMassResults[1]- relMass)
        if diff*100>0.05:
            log.warning('Conservation of mass is not satisfied' )
            log.warning('Total mass change and total entrained mass differ from %.4f %%' % (diff*100))
        else:
            log.info('Total mass change and total entrained mass differ from %.4f %%' % (diff*100))
#   growth results
    growthIndex = totMassResults[1]/totMassResults[0]
    growthGrad = (totMassResults[1] - totMassResults[0]) / (timeResults[1] - timeResults[0])
    return relMass, entMass, finalMass, growthIndex, growthGrad


def getMaxMeanValues(rasterdataA, rasterArea, pLim, cInd=None):
    """
    Compute average and max of in each cross section for a given input raster
    input:
            -rasterdataA: raster data
            -rasterArea: raster area corresponding to rasterdataA
            -pLim: pressure threshold
            -cInd: index of bounds of the avalanche (from release to run-out)
    output:
            -ama: average maximum of rasterdataA
            -mma: maximum maximum of rasterdataA
            -cInd: index of bounds of the avalanche (from release to run-out)
            -aCrossMax: 1D Aarray containing max of rasterdataA in each cross section
    """
    # get mean max for each cross section for A field
    aCrossMean = np.nansum(rasterdataA*rasterArea, axis=1)/np.nansum(rasterArea, axis=1)
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
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % pLim)
            cupper = 0
            clower = 0
            # search in mean values
            lindex = np.nonzero(aCrossMean > pLim)[0]
        if lindex.any():
            cupperm = min(lindex)
            clowerm = max(lindex)
        else:
            log.error('No average pressure values > threshold found. threshold = %10.4f, too high?' % pLim)
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
