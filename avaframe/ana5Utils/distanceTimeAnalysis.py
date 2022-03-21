""" functions to convert to a different coordinate system and
    produce a range-time diagram from simulation results
    options: 1) range-time diagram from radar's field of view capped to this
             2) thalweg-time diagram from beta point of view along avalanche path
"""


import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import logging
import pathlib

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.outTopo as oP
import avaframe.in3Utils.geoTrans as gT
import avaframe.out3Plot.outDistanceTimeAnalysis as dtAnaPlots
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.out3Plot.plotUtils as pU
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.in3Utils.fileHandlerUtils as fU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getRadarLocation(cfg):
    """ fetch radar location and direction of of field of view coordinates

        Parameters
        ------------
        cfg: configparser object
            configuration settings

        Returns
        --------
        radarFov: numpy array
            array with x and y coordinates of radar location and point in direction of field of view

    """

    # fetch coordinate info
    locationInfo = cfg['GENERAL']['radarLocation'].split('|')

    # set radar location and point in direction of field of view as list
    radarFov = np.asarray([[float(locationInfo[0]), float(locationInfo[1])], [float(locationInfo[2]),
        float(locationInfo[3])]])

    return radarFov


def setDemOrigin(demSims):
    """ reset original origin to simulation DEM - required if for ava sim is set to something different

        Parameters
        ------------
        demSims: dict
            dictionary with header and data of DEM used to run ava simulations
            must contain a key called originalHeader - where the original origin is given

        Returns
        --------
        demOriginal: dict
            updated DEM with xllcenter and yllcenter set to original origin
            cellsize and nrows and columns kept from simulation DEM
    """
    if 'key' not in demSims:
        message = 'DEM dictionary does not have originalHeader key - required to get original origin'
        log.error(message)
        raise KeyError(message)

    # fetch dem data and header - use simulation dem data but xllcenter und yllcenter from original dem
    demOriginal = {'header': demSims['originalHeader'], 'rasterData': demSims['rasterData']}

    return demOriginal


def setupRangeTimeDiagram(dem, cfgRangeTime):
    """ create setup for creating range-time diagrams from avalanche simulation results
        with respect to a radar's field of view

        Parameters
        -----------
        dem: dict
            DEM dictionary with header and raster data
        cfgRangeTime: configparser object
            configuration settings for range-time diagrams

        Returns
        --------
        rangeMasked: numpy masked array
            radar range retrieved from DEM and masked with radar's field of view
    """

    # load parameters from configuration
    rgWidth = cfgRangeTime['GENERAL'].getfloat('rgWidth')

    # fetch radar field of view coordinates
    radarFov = getRadarLocation(cfgRangeTime)

    # create masked array of radar range with DEM extent using radar field of view and range gates
    rangeMasked, rangeGates = radarMask(dem, radarFov, cfgRangeTime['GENERAL'].getfloat('aperture'),
        cfgRangeTime)

    # create range gates and initialise mti array
    rArray = np.ma.getdata(rangeMasked)
    nRangeGate = len(rangeGates)
    mti = np.zeros((nRangeGate, 1))

    # fetch field of view mask
    bmaskRadar = np.ma.getmask(rangeMasked)

    # add all info to mti dict
    mtiInfo = {'mti': mti, 'rangeGates': rangeGates, 'rArray': rArray, 'bmaskRadar': bmaskRadar,
        'rangeList': [], 'timeList': [], 'rangeMasked': rangeMasked,
        'demOriginal': dem, 'type': 'rangeTime', 'referencePointName': 'radar'}

    return mtiInfo


def radarMask(demOriginal, radarFov, aperture, cfgRangeTime):
    """ function to create masked array of radar range using dem and radar's field of view

        Parameters
        -----------
        dem: dict
            DEM dictionary with header info and rasterData
        radarFov: numpy array
            radars field of view min and max points
        aperture: float
            aperature angle
        cfgRangeTime: configparser object
            configuration settings here used: avalancheDir info

        Returns
        --------
        radarRange: numpy masked array
            masked array of radar range covering full DEM extent
        rangeGates: numpy array
            range gates
    """

    # fetch header info - required for creating coordinate grid
    xllc = demOriginal['header']['xllcenter']
    yllc = demOriginal['header']['yllcenter']
    cellSize = demOriginal['header']['cellsize']
    rasterdata = demOriginal['rasterData']

    # Set coordinate grid with given origin
    X, Y = oP._setCoordinateGrid(xllc, yllc, cellSize, rasterdata)

    if (np.any(radarFov[0] < np.amin(X)) or np.any(radarFov[0] > np.amax(X))  or
        np.any(radarFov[1] < np.amin(Y))  or np.any(radarFov[1] > np.amax(Y))):
        message = 'Radar location outside of DEM, x= %.2f, %.2f, y= %.2f, %.2f' % (radarFov[0][0],
            radarFov[0][1], radarFov[1][0], radarFov[1][1])
        log.error(message)
        raise AssertionError(message)

    # define radar field of view as dict
    radarpath = {'x': np.asarray(radarFov[0]),
                 'y': np.asarray(radarFov[1])}
    # add z-coordinate from dem data
    radarProfile = gT.projectOnRaster(demOriginal, radarpath, interp='bilinear')

    # radar location - position of radar [0] and point in direction of field of view [1]
    rX = radarProfile[0]['x'][0]
    rY = radarProfile[0]['y'][0]
    rZ = radarProfile[0]['z'][0]
    # a point that gives direction - point in direction of field of view
    rDx = radarProfile[0]['x'][1]
    rDy = radarProfile[0]['y'][1]
    rDz = radarProfile[0]['z'][1]

    # convert centerline of radar field of view to spherical coordinates
    rD, phiD, thetaD = gT.cartToSpherical(X=rDx-rX, Y=rDy-rY, Z=rDz-rZ)

    # domain in spherical coordinates with center beeing radar position
    r, phi, theta = gT.cartToSpherical(X=X-rX, Y=Y-rY, Z=rasterdata-rZ)

    # create mask where in spherical coordinates azimuth angle is bigger than field of view of radar
    # this is defined by the radar position where it is looking and the aperature angle
    maskPhi = ~np.logical_and(phi > phiD - aperture, phi < phiD + aperture)

    # compute radar range field (distance from radar -masked with field of view)
    # TODO: check if we need to mask also where there are nans in dem
    radarRange = np.ma.array(r, mask=maskPhi)

    # create range gates
    minR = int(np.floor(np.nanmin(radarRange)))
    maxR = int(np.ceil(np.nanmax(radarRange)))
    rgWidth = cfgRangeTime['GENERAL'].getfloat('rgWidth')
    rangeGates= np.arange(minR+rgWidth/2, maxR+rgWidth/2, rgWidth)

    # create plot of range distance already masked with radar field of view
    dtAnaPlots.radarFieldOfViewPlot(radarFov, aperture, radarRange, X, Y, cfgRangeTime['GENERAL'],
        rangeGates, demOriginal)

    return radarRange, rangeGates


def minRangeSimulation(flowF, rangeMasked, threshold):
    """ fetch info on the min range of simulation with respect to masked range array
        the masked range array gives info what the visible range is

        Parameters
        -------------
        flowF: numpy nd array
            flow parameter result field
        rangeMasked: masked numpy array
            masked array of radar range
        threshold: float
            used to create mask for getting line of sight range,
            data below threshold not taken into account

        Returns
        --------
        losRange: float
            minimum distance in line of sight to avalanche front
            defined by the flowF

    """

    # fetch mask where result field smaller than threshold
    maskAva = np.ma.masked_where(flowF < threshold, flowF)
    bmaskAva = np.ma.getmask(maskAva)
    # fetch mask of radar field of view
    maskPhi = np.ma.getmask(rangeMasked)

    # merge masks of result field smaller threshold and outside radar field of view
    maskFull = ~np.logical_and(~maskPhi, ~bmaskAva)

    # first get full data of the masked radar range array
    r = np.ma.getdata(rangeMasked)
    # mask this with full mask (result below threshold and outside field of view)
    radarRange = np.ma.array(r, mask=maskFull)

    # line of sight min of masked array
    losRange = radarRange.min()

    return losRange


def extractMeanValuesAtRangeGates(cfgRangeTime, flowF, mtiInfo):
    """ extract average values of flow parameter result at certain distance intervals (range gates)
        add info to be used for colorcoding range-time diagram

        Parameters
        -----------
        cfgRangeTime: configparser object
            configuration settings for creating range-time diagram
        flowF: numpy array
            flow parameter result field
        mtiInfo: dict
            info here used: rArray, bmaskRadar rangeMasked, rangeGates, mti, timeList

        Returns
        --------
        mtiInfo: dict
            updated mtiInfo dict with info on mti values
    """

    # load parameters from configuration file and mtiInfo
    rgWidth = cfgRangeTime['GENERAL'].getfloat('rgWidth')
    # if threshold = 0.0 mean over all cells - if threshold > 0.0 mean only over all cells where there
    # is avalanche flow
    threshold = cfgRangeTime['GENERAL'].getfloat('thresholdResult')
    rangeGates = mtiInfo['rangeGates']
    rArray = mtiInfo['rArray']
    bmaskRadar = mtiInfo['bmaskRadar']
    rangeMasked = mtiInfo['rangeMasked']
    mti = mtiInfo['mti']
    mtiNew = np.zeros((len(rangeGates), 1))

    # get mask of results below threshold
    maskAva = np.ma.masked_where(flowF < threshold, flowF)
    bmaskAva = np.ma.getmask(maskAva)
    # and join to a mask with "visible part of avalanche" - mask where not field of view and not
    # ava result above threshold
    bmaskAvaRadar = ~np.logical_and(~bmaskRadar, ~bmaskAva)
    # mask range with this mask to obtain 'visible part of avalanche' ranges
    rmaskAva_radar= np.ma.array(np.ma.getdata(rangeMasked), mask=bmaskAvaRadar) # masked ranges

    # min and max radar range of masked radar range
    # if no data in masked array
    if not rmaskAva_radar.mask.all():
        #TODO why int?
        minRAva = int(np.floor(rmaskAva_radar.min()))
        maxRAva = int(np.ceil(rmaskAva_radar.max()))
    else:
        minRAva = np.amin(rangeGates)
        maxRAva = np.amax(rangeGates)

    # create array of capped range gates
    smallRangeGatesIndex = np.where((rangeGates > minRAva) & (rangeGates < maxRAva))[0]

    # loop over each range gate within visible part to populate mti array
    for indexRI in smallRangeGatesIndex:
        # create mask for range slice within range gate
        bmaskRange = ~np.logical_and(rArray > rangeGates[indexRI]-rgWidth/2,
            rArray < rangeGates[indexRI]+rgWidth/2)
        bmaskAvaRadarRangeslice = ~np.logical_and(~bmaskRange, ~bmaskAvaRadar)
        if np.any(~bmaskAvaRadarRangeslice):
            # only update if range_slice mask is not empty
            mtiValue = np.mean(np.ma.array(np.ma.getdata(maskAva), mask=bmaskAvaRadarRangeslice))
            mtiNew[indexRI] = mtiValue
            if cfgRangeTime['GENERAL'].getboolean('debugPlot'):
                dtAnaPlots.plotMaskForMTI(cfgRangeTime['GENERAL'], bmaskRange, bmaskAvaRadar,
                bmaskAvaRadarRangeslice, mtiInfo)

    # add average values of this time step to full mti array
    if mtiInfo['timeList'] == []:
        mti = mtiNew
    else:
        mti = np.hstack((mti, mtiNew))

    # update mtiInfo dict
    mtiInfo['mti'] = mti

    return mtiInfo


def extractFrontInSim(flowF, mtiInfo, threshold):
    """ extract front in simulation result raster with respect to radar and add range info

        Parameters
        -----------
        flowF: numpy array
            simulation result of flow parameter
        mtiInfo: dict
            info on here used: rangeList, demOriginal, rangeMasked
        threshold: float
            used to create mask for getting line of sight range,
            data below threshold not taken into account

        Returns
        --------
        mtiInfo: dict
            updated dictionary new info on rangeList - list of front distance to radar
            for respective time step
    """

    # get line of sight distance to identify front, use threshold to mask flowF
    los = minRangeSimulation(flowF, mtiInfo['rangeMasked'], threshold)

    # update lists of time step and front location
    mtiInfo['rangeList'].append(los)

    return mtiInfo


def setupThalwegTimeDiagram(dem, cfgRangeTime):
    """ initialize parameters and prerequisites for generating thalweg-time diagrams

        Parameters
        -----------
        dem: dict
            dictionary with info on dem header and data
        cfgRangeTime: confiparser object
            configuration settings of range-time diagram

        Returns
        --------
        mtiInfo: dict
            dictionary with arrays, lists for range-time diagram:
            mti for colorcoding (average values of range gates/ crossprofiles)
            timeList time step info
            rangeList range info
            rangeGates info on cross profile locations
            rangeRaster - info on distance from beta point in path following coordinate system
            betaPoint - coordinates of start of runout zone point
            betaPointAngle - angle of beta point used
    """

    # set required parameters from/to configuration
    avaDir = cfgRangeTime['GENERAL']['avalancheDir']
    cfgSetup = cfgRangeTime['GENERAL']
    resType = cfgRangeTime['GENERAL']['rangeTimeResType']
    cfgSetup['resType'] = resType

    # fetch info on avalanche path, split point, path to dem
    pathDict = {}
    pathDict = aT.readAIMECinputs(avaDir, pathDict, dirName='com1DFA')

    # generate data required to perform domain tranformation
    rasterTransfo = aT.makeDomainTransfoOnTheGo(avaDir, dem, cfgSetup, pathDict)

    # tranform DEM data
    log.debug("Assigning dem data to deskewed raster")
    newRasterDEM = aT.transform(dem, rasterTransfo, cfgSetup['interpMethod'])

    # create raster with info on distance from beta point
    # fetch info on distance of runout area start point index in s array
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    # create range raster with respect to runout area start point
    # minus in upstream direction of runout area start point
    rangeGates = rasterTransfo['s'][:] - rasterTransfo['s'][indStartOfRunout]
    rangeRaster = np.repeat([rangeGates], newRasterDEM.shape[1], axis=0).T

    # create dictionary with info on front distance, mean values, range distance array
    mtiInfo = {'mti': '', 'rangeList': [], 'timeList': []}
    # mtiInfo['rangeGates'] = rasterTransfo['s'][indStartOfRunout] - rasterTransfo['s'][:]
    mtiInfo['rangeGates'] = rangeGates
    mtiInfo['betaPoint'] = [rasterTransfo['x'][indStartOfRunout], rasterTransfo['y'][indStartOfRunout]]
    mtiInfo['betaPointAngle'] = rasterTransfo['startOfRunoutAreaAngle']
    mtiInfo['rangeRaster'] = rangeRaster
    mtiInfo['rasterTransfo'] = rasterTransfo
    mtiInfo['type'] = 'thalwegTime'
    mtiInfo['z'] = rasterTransfo['z']
    mtiInfo['s'] = rasterTransfo['s']
    mtiInfo['referencePointName'] = 'beta point'

    return mtiInfo


def extractFrontAndMeanValues(cfgRangeTime, flowF, demHeader, mtiInfo):
    """ extract avalanche front and mean values of flow parameter result field
        used for colorcoding range-time diagram

        Parameters
        ------------
        cfgRangeTime: configparser object
            configuration settings
        flowF: numpy nd array
            flow parameter result field
        demHeader: dict
            dictionary with info on DEM header used for locating fields
        mtiInfo: dict
            dictionary with arrays, lists for range-time diagram:
            mti for colorcoding (average values of range gates/ crossprofiles)
            timeList time step info
            rangeList range info
            rangeGates info on cross profile locations
            rangeRaster - info on distance from beta point in path following coordinate system
            rasterTransfo - information on domain transformation to path following coordinate system

        Returns
        --------
        mtiInfo: dict
            updated dictionary with info on average value along crossprofiles (mti)
            and distance to front from reference point along path (rangeList)
    """

    # setup input parameters
    rasterTransfo = mtiInfo['rasterTransfo']

    # add header to result field - required for coordinate transformation
    fieldFull = {'header': demHeader, 'rasterData': flowF}

    # transorm result field into avalanche path following coordinate system
    slRaster = aT.transform(fieldFull, rasterTransfo, cfgRangeTime['GENERAL']['interpMethod'])

    # fetch raster area and compute mean, max values for each cross-profile
    # TODO: average over cells – weighted will cell area but only projected area (aimec function)
    rasterArea = rasterTransfo['rasterArea']
    maxaCrossMax, aCrossMax, aCrossMean = aT.getMaxMeanValues(slRaster, rasterArea)

    # add mean values for each cross-section for actual time step to mean values array
    # reshape is required as rows- mean values for each crossprofile and cols: time steps
    if mtiInfo['timeList'] == []:
        mtiInfo['mti'] = aCrossMean.reshape(-1,1)
    else:
        mtiInfo['mti'] = np.hstack((mtiInfo['mti'], aCrossMean.reshape(-1,1)))

    # extract avalanche front distance to reference point in path following coordinate system
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    lindex = np.nonzero(slRaster > cfgRangeTime['GENERAL'].getfloat('thresholdResult'))[0]
    if lindex.any():
        cUpper = min(lindex)
        cLower = max(lindex)
        rangeValue = rasterTransfo['s'][cLower] - rasterTransfo['s'][indStartOfRunout]
        mtiInfo['rangeList'].append(rangeValue)
    else:
        mtiInfo['rangeList'].append(np.nan)

    # plot avalanche front and transformed raster field
    if cfgRangeTime['GENERAL'].getboolean('debugPlot'):
        dtAnaPlots.plotRangeRaster(slRaster, rasterTransfo, cfgRangeTime['GENERAL'],
            mtiInfo['rangeRaster'], cLower)

    return mtiInfo


def initializeRangeTime(modName, cfg, dem, simHash):
    """ initialize generation of range-time diagram for visualizing simulation data

        Parameters
        -----------
        modName:
            module name
        cfg: configparser object
            configuration settings from computational module
        dem: dict
            dictionary with DEM header and data
        simHash: str
            unique simulation ID

        Returns
        --------
        mtiInfo: dict
            info on time steps (timeList) and distance (rangeList) of avalanche front
            averaged values of result parameter (mti) for each range gate (rangeGates) for colorcoding
            x, y origin of DEM, cellSize and plotTitle
            optional info on masks for radar field of view (rangeMasked), etc.
    """

    # fetch configuration and add info
    cfgRangeTime = cfgUtils.getModuleConfig(modName)
    cfgRangeTime['GENERAL']['tEnd'] = cfg['GENERAL']['tEnd']
    cfgRangeTime['GENERAL']['avalancheDir'] = cfg['GENERAL']['avalancheDir']
    cfgRangeTime['GENERAL']['simHash'] = simHash
    # fetch time steps for creating range time diagram
    dtRangeTime = fU.splitTimeValueToArrayInterval(cfgRangeTime['GENERAL'])

    if cfg['VISUALISATION'].getboolean('TTdiagram'):
        mtiInfo = setupThalwegTimeDiagram(dem, cfgRangeTime)
        mtiInfo['plotTitle'] = 'thalweg-time diagram %s' % simHash
        mtiInfo['textbox'] = 'beta point: %.2f, %.2f' % (mtiInfo['betaPoint'][0],
            mtiInfo['betaPoint'][1])
    else:
        mtiInfo = setupRangeTimeDiagram(dem, cfgRangeTime)
        mtiInfo['plotTitle'] = 'range-time diagram %s' % simHash

    mtiInfo['xOrigin'] = dem['header']['xllcenter']
    mtiInfo['yOrigin'] = dem['header']['yllcenter']
    mtiInfo['cellSize'] = dem['header']['cellsize']

    return mtiInfo, dtRangeTime, cfgRangeTime


def fetchRangeTimeInfo(cfgRangeTime, cfg, dtRangeTime, t, demHeader, fields, mtiInfo):
    """ determine avalanche front and average values of flow parameter along path
        update mtiInfo dictionary with this information

        Parameters
        -----------
        cfgRangeTime: configparser
            configuration settings for range-time diagram
        cfg: configParser object
            configuration settings of computational module
        dtRangeTime: list
            list of time steps where avalanche front shall be exported
        t: float
            actual simulation time step
        demHeader: dict
            dictionary with DEM header
        fields: dict
            dictionary with flow parameter result fields
        mtiInfo: dict
            info on time steps (timeList) and distance (rangeList) of avalanche front_
            averaged values of result parameter (mti) for each range gate (rangeGates) for colorcodig
            optional info on masks for radar field of view (rangeMasked), etc.

        Returns
        --------
        mtiInfo: dict
            updated dictionary
        dtRangeTime: list
            updated list of time steps
    """

    # load result type for fields
    rangeTimeResType = cfgRangeTime['GENERAL']['rangeTimeResType']

    # log info
    log.info('Processing frame %d of time step: %.2f for range-time plot' % (len(mtiInfo['timeList'])+1, t))

    # extract front and average values of result parameter
    if cfg['VISUALISATION'].getboolean('TTdiagram'):
        mtiInfo = extractFrontAndMeanValues(cfgRangeTime, fields[rangeTimeResType], demHeader, mtiInfo)
    else:
        # extract front in simulation result for each time step
        mtiInfo = extractFrontInSim(fields[rangeTimeResType], mtiInfo,
             cfgRangeTime['GENERAL'].getfloat('thresholdResult'))
        # extract average values along range gates for coloring range-time plot
        mtiInfo = extractMeanValuesAtRangeGates(cfgRangeTime, fields[rangeTimeResType], mtiInfo)

    # append time step info
    mtiInfo['timeList'].append(t)

    # remove saving time steps that have already been saved
    dtRangeTime = com1DFA.updateSavingTimeStep(dtRangeTime, cfgRangeTime['GENERAL'], t)

    return mtiInfo, dtRangeTime


def fetchTimeStepFromName(pathNames):
    """ split path name to fetch info on time step

        Parameters
        -----------
        pathName: pathlib path or list of paths
            path to file including time step info as _t%.2f_

        Returns
        -------
        timeStep: list
            actual time steps
        indexTime: numpy array
            index of timeStep list to order them in ascending order
    """

    if isinstance(pathNames, list) is False:
        pathNames = [pathNames]

    timeSteps = []
    for fName in pathNames:
        timeSteps.append(float(fName.stem.split('_t')[1]))

    indexTime = np.argsort(np.array(timeSteps))

    return timeSteps, indexTime


def approachVelocity(mtiInfo):
    """ compute approach velocity based on front location and time step

        Parameters
        -----------
        mtiInfo: dict
            info on distance to front and time steps

        Returns
        --------
        maxVel: float
            max value of approach velocity
    """

    # load lists
    timeList = mtiInfo['timeList']
    rangeList = mtiInfo['rangeList']

    # as these might be not in ascending order timestep wise - order first
    timeListSorted = np.asarray(timeList)[np.argsort(timeList)]
    # convert to only positive distance values (move reference point to end of ava)
    # to get correct distance increments for velocity computation
    # nanmin required as if FV is zero range is nan
    rangeListSorted = np.asarray(rangeList)[np.argsort(timeList)] + abs(np.nanmin(np.asarray(rangeList)))
    maxVel = 0.0
    rangeVel = 0.0
    # compute deltaDistance/deltaTime to get approach velocity
    # use interval of 2*dtSave for little smmoothing
    for index, range in enumerate(rangeListSorted[1:-1]):
        vel = (rangeListSorted[index+1] - rangeListSorted[index-1]) / (timeListSorted[index+1] - timeListSorted[index-1])
        if abs(vel) > maxVel:
            maxVel = abs(vel)
            rangeVel = rangeListSorted[index] - abs(np.nanmin(np.asarray(rangeList)))
            timeVel = timeListSorted[index]

    return maxVel, rangeVel, timeVel
