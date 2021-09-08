"""
    Main functions for python DFA kernel
"""

import logging
import time
import pathlib
import numpy as np
import copy
import pickle
from datetime import datetime
import pandas as pd
import matplotlib.path as mpltPath
import pandas as pds
from itertools import product

# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.com1DFA.timeDiscretizations as tD
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.com1DFA.DFAfunctionsCython as DFAfunC
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
import avaframe.out3Plot.outDebugPlots as debPlot

#######################################
# Set flags here
#######################################
# create local logger
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')
# set feature leapfrog time stepping
featLF = False


def com1DFAMain(cfg, avaDir, cuSimName, inputSimFiles, outDir, relThField):
    """ Run main com1DFA model

    This will compute a dense flow avalanche

    Parameters
    ----------
    cfg : dict
        configuration read from ini file
    cuSimName: str
        name of simulation
    inputSimFiles: dict
        dictionary with input files
    avaDir : str or pathlib object
        path to avalanche directory
    outDir: str or pathlib object
        path to Outputs
    relThField: 2D array
        release thickness field with varying release thickness if '', release thickness is taken from
        (a) shapefile or (b) configuration file

    Returns
    -------
    reportDictList : list
        list of dictionaries that contain information on simulations that can be used for report generation
    """

    # Setup configuration
    cfgGen = cfg['GENERAL']

    # create required input from files
    demOri, inputSimLines = prepareInputData(inputSimFiles)

    # find out which simulations to perform
    relName, inputSimLines, badName = prepareRelase(cfg, inputSimFiles['releaseScenario'], inputSimLines)

    log.info('Perform %s simulation' % cuSimName)

    # +++++++++PERFORM SIMULAITON++++++++++++++++++++++
    # for timing the sims
    startTime = time.time()
    particles, fields, dem, reportAreaInfo = initializeSimulation(cfg, demOri, inputSimLines, cuSimName, relThField)

    # ------------------------
    #  Start time step computation
    Tsave, particlesList, fieldsList, infoDict = DFAIterate(cfg, particles, fields, dem)

    # write mass balance to File
    writeMBFile(infoDict, avaDir, cuSimName)

    tcpuDFA = '%.2f' % (time.time() - startTime)
    log.info(('cpu time DFA = %s s' % (tcpuDFA)))

    if 'particles' in cfgGen['resType']:
        # export particles dictionaries of saving time steps
        outDirData = outDir / 'particles'
        fU.makeADir(outDirData)
        savePartToPickle(particlesList, outDirData, cuSimName)

    # Result parameters to be exported
    exportFields(cfg, Tsave, fieldsList, demOri, outDir, cuSimName)

    # write report dictionary
    reportDict = createReportDict(avaDir, cuSimName, relName, inputSimLines, cfgGen, reportAreaInfo)
    # add time and mass info to report
    reportDict = reportAddTimeMassInfo(reportDict, tcpuDFA, infoDict)

    return particlesList, fieldsList, Tsave, dem, reportDict, cfg


def prepareRelase(cfg, rel, inputSimLines):
    """ get Simulation to run for a given release


    Parameters
    ----------
    cfg : dict
        configuration parameters - keys: relTh, secRelArea, secondaryRelTh
    rel : str
        path to release file

    Returns
    -------
    relName : str
        release name
    relDict : list
        release dictionary
    badName : boolean
        changed release name
    """

    # load info
    entResInfo = inputSimLines['entResInfo']

    # Set release areas and simulation name
    relName = rel.stem
    simName = relName
    badName = False
    if '_' in relName:
        badName = True
        log.warning('Release area scenario file name includes an underscore \
        the suffix _AF will be added for the simulation name')
    releaseLine = inputSimLines['releaseLine']
    for k in range(len(releaseLine['d0'])):
        if releaseLine['d0'][k] == 'None':
            releaseLine['d0'][k] = cfg['GENERAL'].getfloat('relTh')
        else:
            releaseLine['d0'][k] = float(releaseLine['d0'][k])
    inputSimLines['releaseLine'] = releaseLine
    log.debug('Release area scenario: %s - perform simulations' % (relName))

    if cfg.getboolean('GENERAL', 'secRelArea'):
        if entResInfo['flagSecondaryRelease'] == 'No':
            message = 'No secondary release file found'
            log.error(message)
            raise FileNotFoundError(message)
        secondaryReleaseLine = inputSimLines['secondaryReleaseLine']
        for k in range(len(secondaryReleaseLine['d0'])):
            if secondaryReleaseLine['d0'][k] == 'None':
                secondaryReleaseLine['d0'][k] = cfg['GENERAL'].getfloat('secondaryRelTh')
            else:
                secondaryReleaseLine['d0'][k] = float(secondaryReleaseLine['d0'][k])
    else:
        inputSimLines['entResInfo']['flagSecondaryRelease'] = 'No'
        secondaryReleaseLine = None

    inputSimLines['secondaryReleaseLine'] = secondaryReleaseLine

    return relName, inputSimLines, badName


def prepareInputData(inputSimFiles):
    """ Fetch input data

    Parameters
    ----------
    relFiles : str
        path to release file
    inputSimFiles : dict
        demFile : str
            path to dem file
        secondaryReleaseFile : str
            path to secondaryRelease file
        entFiles : str
            path to entrainment file
        resFile : str
            path to resistance file
        entResInfo : flag dict
            flag if Yes entrainment and/or resistance areas found and used for simulation
            flag True if a Secondary Release file found and activated

    Returns
    -------
    demOri : dict
        dictionary with original dem
    inputSimLines : dict
        releaseLine : dict
            release line dictionary
        secondaryReleaseLine : dict
            secondaryRelease line dictionary
        entLine : dict
            entrainment line dictionary
        resLine : dict
            resistance line dictionary
        entrainmentArea : str
            entrainment file name
        resistanceArea : str
            resistance file name
        entResInfo : flag dict
            flag if Yes entrainment and/or resistance areas found and used for simulation
            flag True if a Secondary Release file found and activated
    """

    # load data
    entResInfo = inputSimFiles['entResInfo']
    relFile = inputSimFiles['releaseScenario']

    # get dem information
    demOri = IOf.readRaster(inputSimFiles['demFile'], noDataToNan=True)

    # get line from release area polygon
    releaseLine = shpConv.readLine(relFile, 'release1', demOri)
    releaseLine['file'] = relFile

    # get line from secondary release area polygon
    if entResInfo['flagSecondaryRelease'] == 'Yes':
        secondaryReleaseFile = inputSimFiles['secondaryReleaseFile']
        secondaryReleaseLine = shpConv.readLine(secondaryReleaseFile, '', demOri)
        secondaryReleaseLine['fileName'] = [secondaryReleaseFile]
    else:
        secondaryReleaseLine = None

    # get line from entrainement area polygon
    if entResInfo['flagEnt'] == 'Yes':
        entFile = inputSimFiles['entFile']
        entLine = shpConv.readLine(entFile, '', demOri)
        entrainmentArea = entFile.name
        entLine['fileName'] = entFile
    else:
        entLine = None
        entrainmentArea = ''

    # get line from resistance area polygon
    if entResInfo['flagRes'] == 'Yes':
        resFile = inputSimFiles['resFile']
        resLine = shpConv.readLine(resFile, '', demOri)
        resistanceArea = resFile.name
        resLine['fileName'] = resFile
    else:
        resLine = None
        resistanceArea = ''

    inputSimLines = {'releaseLine': releaseLine, 'secondaryReleaseLine': secondaryReleaseLine,
                     'entLine': entLine, 'resLine': resLine, 'entrainmentArea': entrainmentArea,
                     'resistanceArea': resistanceArea, 'entResInfo': entResInfo}

    return demOri, inputSimLines


def createReportDict(avaDir, logName, relName, inputSimLines, cfgGen, reportAreaInfo):
    """ create simulaton report dictionary

    Parameters
    ----------
    logName : str
        simulation scenario name
    relName : str
        release name
    relDict : dict
        release dictionary
    cfgGen : configparser
        general configuration file
    entrainmentArea : str
        entrainment file name
    resistanceArea : str
        resistance file name

    Returns
    -------
    reportST : dict
        simulation scenario dictionary
    """

    # load parameters
    dateTimeInfo = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    entInfo = reportAreaInfo['entrainment']
    resInfo = reportAreaInfo['resistance']
    entrainmentArea = inputSimLines['entrainmentArea']
    resistanceArea = inputSimLines['resistanceArea']
    relDict = inputSimLines['releaseLine']

    # Create dictionary
    reportST = {}
    reportST = {}
    reportST = {'headerLine': {'type': 'title', 'title': 'com1DFA Simulation'},
                'avaName': {'type': 'avaName', 'name': avaDir},
                'simName': {'type': 'simName', 'name': logName},
                'time': {'type': 'time', 'time': dateTimeInfo},
                'Simulation Parameters': {
                'type': 'list',
                'Program version': 'development',
                'Parameter set': '',
                'Release Area Scenario': relName,
                'Entrainment': entInfo,
                'Resistance': resInfo,
                'Parameter variation on': '',
                'Parameter value': '',
                'Mu': cfgGen['mu'],
                'Density [kgm-3]': cfgGen['rho'],
                'Friction model': cfgGen['frictModel']},
                'Release Area': {'type': 'columns', 'Release area scenario': relName, 'Release Area': relDict['Name'],
                                 'Release thickness [m]': relDict['d0']}}

    if entInfo == 'Yes':
        reportST.update({'Entrainment area':
                            {'type': 'columns',
                            'Entrainment area scenario': entrainmentArea,
                            'Entrainment thickness [m]': cfgGen.getfloat('hEnt'),
                            'Entrainment density [kgm-3]': cfgGen['rhoEnt']}})
    if resInfo == 'Yes':
        reportST.update({'Resistance area': {'type': 'columns', 'Resistance area scenario': resistanceArea}})

    reportST['Release Area'].update(reportAreaInfo['Release area info'])

    return reportST


def reportAddTimeMassInfo(reportDict, tcpuDFA, infoDict):
    """ Add time and mass info to report """

    # add mass info
    reportDict['Simulation Parameters'].update({'Initial mass [kg]': ('%.2f' % infoDict['initial mass'])})
    reportDict['Simulation Parameters'].update({'Final mass [kg]': ('%.2f' % infoDict['final mass'])})
    reportDict['Simulation Parameters'].update({'Entrained mass [kg]': ('%.2f' % infoDict['entrained mass'])})
    reportDict['Simulation Parameters'].update({'Entrained volume [m3]': ('%.2f' % infoDict['entrained volume'])})

    # add stop info
    reportDict['Simulation Parameters'].update(infoDict['stopInfo'])

    # add computation time to report dict
    reportDict['Simulation Parameters'].update({'Computation time [s]': tcpuDFA})

    return reportDict


def initializeMesh(cfg, demOri, num):
    """ Create rectangular mesh

    Reads the DEM information, computes the normal vector field and
    boundries to the DEM. Also generates the grid for the neighbour search

    Parameters
    ----------
    demOri : dict
        dictionary with initial dem information
    num : int
        chose between 4, 6 or 8 (using then 4, 6 or 8 triangles) or
        1 to use the simple cross product method

    Returns
    -------
    dem : dict
        dictionary relocated in (0,0) and completed with normal field and
        boundaries as well as neighbour search grid information
    """

    demOri = geoTrans.remeshDEM(cfg, demOri)
    dem = setDEMoriginToZero(demOri)
    dem['originOri'] = {'xllcenter': demOri['header']['xllcenter'], 'yllcenter': demOri['header']['yllcenter']}

    # read dem header
    headerDEM = dem['header']
    nColsDEM = headerDEM['ncols']
    nRowsDEM = headerDEM['nrows']
    cszDEM = headerDEM['cellsize']

    # get normal vector of the grid mesh
    Nx, Ny, Nz = DFAtls.getNormalMesh(dem, num)
    # if no normal available, put 0 for Nx and Ny and 1 for Nz
    dem['Nx'] = np.where(np.isnan(Nx), 0., Nx)
    dem['Ny'] = np.where(np.isnan(Ny), 0., Ny)
    # build no data mask (used to find out of dem particles)
    outOfDEM = np.where(np.isnan(dem['rasterData']), 1, 0).astype(np.bool).flatten()
    dem['Nz'] = Nz
    dem['outOfDEM'] = outOfDEM

    # Prepare SPH grid
    headerNeighbourGrid = {}
    cszNeighbourGrid = cfg.getfloat('sphKernelRadius')
    headerNeighbourGrid['cellsize'] = cszNeighbourGrid
    headerNeighbourGrid['ncols'] = np.ceil(nColsDEM * cszDEM / cszNeighbourGrid)
    headerNeighbourGrid['nrows'] = np.ceil(nRowsDEM * cszDEM / cszNeighbourGrid)
    headerNeighbourGrid['xllcenter'] = 0
    headerNeighbourGrid['yllcenter'] = 0
    dem['headerNeighbourGrid'] = headerNeighbourGrid

    # get real Area
    areaRaster = DFAtls.getAreaMesh(Nx, Ny, Nz, cszDEM, num)
    dem['areaRaster'] = areaRaster
    projArea = nColsDEM * nRowsDEM * cszDEM * cszDEM
    actualArea = np.nansum(areaRaster)
    log.debug('Largest cell area: %.2f m²' % (np.nanmax(areaRaster)))
    log.debug('Projected Area : %.2f' % projArea)
    log.debug('Total Area : %.2f' % actualArea)

    return demOri, dem


def setDEMoriginToZero(demOri):
    """ set origin of DEM to 0,0 """

    dem = copy.deepcopy(demOri)
    dem['header']['xllcenter'] = 0
    dem['header']['yllcenter'] = 0

    return dem


def initializeSimulation(cfg, demOri, inputSimLines, logName, relThField):
    """ create simulaton report dictionary

    Parameters
    ----------
    cfg : str
        simulation scenario name
    demOri : dict
        dictionary with original dem
    inputSimLines : dict
        releaseLine : dict
            release line dictionary
        secondaryReleaseLine : dict
            secondary release line dictionary
        entLine : dict
            entrainment line dictionary
        resLine : dict
            resistance line dictionary
    logName : str
        simulation scenario name
    relThField : 2D numpy array
        inhomogeneous release thickness if wanted (relThField='' by default  - in this case
        release thickness from (a) shapefile or if not provided (b) configuration file is used)

    Returns
    -------
    particles : dict
        particles dictionary at initial time step
        list of secondary release particles to be used
    fields : dict
        fields dictionary at initial time step
    dem : dict
        dictionary with new dem (lower left center at origin)
    """
    cfgGen = cfg['GENERAL']
    methodMeshNormal = cfg.getfloat('GENERAL', 'methodMeshNormal')
    thresholdPointInPoly = cfgGen.getfloat('thresholdPointInPoly')

    # -----------------------
    # Initialize mesh
    log.debug('Initializing Mesh')
    demOri, dem = initializeMesh(cfgGen, demOri, methodMeshNormal)

    # ------------------------
    log.debug('Initializing main release area')
    # process release info to get it as a raster
    releaseLine = inputSimLines['releaseLine']
    # check if release features overlap between features
    prepareArea(releaseLine, demOri, thresholdPointInPoly, combine=True, checkOverlap=True)
    if len(relThField) == 0:
        # if no release thickness field or function - set release according to shapefile or ini file
        # this is a list of release rasters that we want to combine
        releaseLine = prepareArea(releaseLine, demOri, np.sqrt(2), thList=releaseLine['d0'], combine=True, checkOverlap=False)
    else:
        # if relTh provided - set release thickness with field or function
        releaseLine = prepareArea(releaseLine, demOri, np.sqrt(2), combine=True, checkOverlap=False)
        releaseLine['rasterData'] = releaseLine['rasterData'] * relThField

    # compute release area
    header = dem['header']
    csz = header['cellsize']
    relRaster = releaseLine['rasterData']
    relRasterOnes = np.where(relRaster > 0, 1., 0.)
    relAreaActual = np.nansum(relRasterOnes*dem['areaRaster'])
    relAreaProjected = np.sum(csz*csz*relRasterOnes)
    reportAreaInfo = {'Release area info': {'Projected Area [m2]': '%.2f' % (relAreaProjected),
                      'Actual Area [m2]': '%.2f' % (relAreaActual)}}

    # ------------------------
    # initialize simulation
    # create primary release area particles and fields
    releaseLine['header'] = demOri['header']
    particles = initializeParticles(cfgGen, releaseLine, dem, logName=logName)
    particles, fields = initializeFields(cfgGen, dem, particles)

    # ------------------------
    # process secondary release info to get it as a list of rasters
    if inputSimLines['entResInfo']['flagSecondaryRelease'] == 'Yes':
        log.info('Initializing secondary release area')
        secondaryReleaseInfo = inputSimLines['secondaryReleaseLine']
        secondaryReleaseInfo['header'] = demOri['header']

        # fetch secondary release areas
        secondaryReleaseInfo = prepareArea(secondaryReleaseInfo, demOri, np.sqrt(2), thList=secondaryReleaseInfo['d0'], combine=False)
        # remove overlap with main release areas
        noOverlaprasterList = []
        for secRelRatser, secRelName in zip(secondaryReleaseInfo['rasterData'], secondaryReleaseInfo['Name']):
            noOverlaprasterList.append(geoTrans.checkOverlap(secRelRatser, relRaster, 'Secondary release ' + secRelName, 'Release', crop=True))

        secondaryReleaseInfo['flagSecondaryRelease'] = 'Yes'
        secondaryReleaseInfo['rasterList'] = noOverlaprasterList
    else:
        secondaryReleaseInfo = {}
        secondaryReleaseInfo['flagSecondaryRelease'] = 'No'

    particles['secondaryReleaseInfo'] = secondaryReleaseInfo

    # initialize entrainment and resistance
    # get info of simType and whether or not to initialize resistance and entrainment
    simTypeActual = cfgGen['simTypeActual']
    rhoEnt = cfgGen.getfloat('rhoEnt')
    hEnt = cfgGen.getfloat('hEnt')
    entrMassRaster, reportAreaInfo = initializeMassEnt(demOri, simTypeActual, inputSimLines['entLine'], reportAreaInfo, thresholdPointInPoly)
    # check if entrainment and release overlap
    entrMassRaster = geoTrans.checkOverlap(entrMassRaster, relRaster, 'Entrainment', 'Release', crop=True)
    # check for overlap with the secondary release area
    if secondaryReleaseInfo['flagSecondaryRelease'] == 'Yes':
        for secRelRaster in secondaryReleaseInfo['rasterList']:
            entrMassRaster = geoTrans.checkOverlap(entrMassRaster, secRelRaster, 'Entrainment', 'Secondary release ', crop=True)
    # surfacic entrainment mass available (unit kg/m²)
    fields['entrMassRaster'] = entrMassRaster*rhoEnt*hEnt
    entreainableMass = np.nansum(fields['entrMassRaster']*dem['areaRaster'])
    log.info('Mass available for entrainment: %.2f kg' % (entreainableMass))

    log.debug('Initializing resistance area')
    cResRaster, reportAreaInfo = initializeResistance(cfgGen, demOri, simTypeActual, inputSimLines['resLine'], reportAreaInfo, thresholdPointInPoly)
    fields['cResRaster'] = cResRaster

    return particles, fields, dem, reportAreaInfo


def initializeParticles(cfg, releaseLine, dem, logName=''):
    """ Initialize DFA simulation

    Create particles and fields dictionary according to config parameters
    release raster and dem

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    relRaster: 2D numpy array
        release depth raster
    dem : dict
        dictionary with dem information

    Returns
    -------
    particles : dict
        particles dictionary at initial time step
    fields : dict
        fields dictionary at initial time step
    """

    # get simulation parameters
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    avaDir = cfg['avalancheDir']
    massPerParticleDeterminationMethod = cfg['massPerParticleDeterminationMethod']
    interpOption = cfg.getfloat('interpOption')

    # read dem header
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    csz = header['cellsize']
    relRaster = releaseLine['rasterData']
    areaRaster = dem['areaRaster']

    # initialize arrays
    partPerCell = np.zeros(np.shape(relRaster), dtype=np.int64)
    # find all non empty cells (meaning release area)
    indRelY, indRelX = np.nonzero(relRaster)

    # derive mass per particle to define number of particles per cell:
    if massPerParticleDeterminationMethod == 'MPPDIR':
        massPerPart = cfg.getfloat('massPerPart')
        log.debug('Number of particles defined by: mass per particle %s' % cfg['massPerPart'])
    elif massPerParticleDeterminationMethod == 'MPPDH':
        deltaTh = cfg.getfloat('deltaTh')
        ds = min(csz, cfg.getfloat('sphKernelRadius'))
        massPerPart = rho * ds * ds * deltaTh
        log.debug('Number of particles defined by: release thickness per particle: %s' % cfg['deltaTh'])
        log.debug('mass per particle is %.2f' % massPerPart)

    # make option available to read initial particle distribution from file
    if cfg.getboolean('initialiseParticlesFromFile'):
        # TODO: this is for development purposes, change or remove in the future
        # If initialisation from file
        if cfg['particleFile']:
            inDirPart = pathlib.Path(cfg['particleFile'])
        else:
            inDirPart = pathlib.Path(avaDir, 'Outputs', 'com1DFAOrig')

        searchDir = inDirPart / 'particles'
        inDirPart = list(searchDir.glob(('*' + cfg['releaseScenario'] + '_' + '*' + cfg['simTypeActual'] + '*')))
        if inDirPart == []:
            messagePart = 'Initialise particles from file - no particles file found for releaseScenario: %s and simType: %s' % \
                          (cfg['releaseScenario'], cfg['simTypeActual'])
            log.error(messagePart)
            raise FileNotFoundError(messagePart)
        elif len(inDirPart) > 1:
            log.warning('More than one file found for Initialise particle from file: took %s' % inDirPart[0])
            inDirPart = inDirPart[0]
        else:
            inDirPart = inDirPart[0]

        log.info('Initial particle distribution read from file!! %s' % (inDirPart))
        Particles, TimeStepInfo = readPartFromPickle(inDirPart)
        particles = Particles[0]
        Xpart = particles['x']
        Mpart = particles['m']
        Hpart = np.ones(len(Xpart))
        NPPC = np.ones(len(Xpart))
        particles['Npart'] = len(Xpart)
        particles['s'] = np.zeros(np.shape(Xpart))
        particles['l'] = np.zeros(np.shape(Xpart))
    else:
        # initialize random generator
        rng = np.random.default_rng(int(cfg['seed']))

        Npart = 0
        NPPC = np.empty(0)
        Xpart = np.empty(0)
        Ypart = np.empty(0)
        Mpart = np.empty(0)
        Hpart = np.empty(0)
        # loop on non empty cells
        for indRelx, indRely in zip(indRelX, indRelY):
            # compute number of particles for this cell
            hCell = relRaster[indRely, indRelx]
            volCell = areaRaster[indRely, indRelx] * hCell
            massCell = volCell * rho
            xpart, ypart, mPart, nPart = placeParticles(massCell, indRelx, indRely, csz,
                                         massPerPart, rng, cfg['initPartDistType'].lower())
            Npart = Npart + nPart
            partPerCell[indRely, indRelx] = nPart
            # initialize particles position, mass, height...
            NPPC = np.append(NPPC, nPart*np.ones(nPart))
            Xpart = np.append(Xpart, xpart)
            Ypart = np.append(Ypart, ypart)
            Mpart = np.append(Mpart, mPart * np.ones(nPart))

        Hpart = DFAfunC.projOnRaster(Xpart, Ypart, relRaster, csz, ncols, nrows, interpOption)
        # create dictionnary to store particles properties
        particles = {}
        particles['Npart'] = Npart
        particles['x'] = Xpart
        particles['y'] = Ypart
        particles['s'] = np.zeros(np.shape(Xpart))
        particles['l'] = np.zeros(np.shape(Xpart))
        # adding z component
        particles, _ = geoTrans.projectOnRaster(dem, particles, interp='bilinear')
        particles['m'] = Mpart

    particles['massPerPart'] = massPerPart
    particles['mTot'] = np.sum(particles['m'])
    particles['h'] = Hpart
    particles['NPPC'] = NPPC
    particles['ux'] = np.zeros(np.shape(Xpart))
    particles['uy'] = np.zeros(np.shape(Xpart))
    particles['uz'] = np.zeros(np.shape(Xpart))
    particles['stoppCriteria'] = False
    kineticEne = np.sum(0.5 * Mpart * DFAtls.norm2(particles['ux'], particles['uy'], particles['uz']))
    particles['kineticEne'] = kineticEne
    particles['potentialEne'] = np.sum(gravAcc * Mpart * particles['z'])
    particles['peakKinEne'] = kineticEne
    particles['peakMassFlowing'] = 0
    particles['simName'] = logName
    particles['xllcenter'] = dem['originOri']['xllcenter']
    particles['yllcenter'] = dem['originOri']['yllcenter']

    if not cfg.getboolean('initialiseParticlesFromFile'):
        particles = checkParticlesInRelease(particles, releaseLine, cfg.getfloat('thresholdPointInPoly'))

    # initialize time
    t = 0
    particles['t'] = t

    relCells = np.size(indRelY)
    partPerCell = particles['Npart']/relCells

    log.info('Initialized particles. MTot = %.2f kg, %s particles in %.2f cells.' %
             (particles['mTot'], particles['Npart'], relCells))
    log.info('Mass per particle = %.2f kg and particles per cell = %.2f.' %
             (particles['mTot']/particles['Npart'], partPerCell))

    if debugPlot:
        debPlot.plotPartIni(particles, dem)

    return particles


def initializeFields(cfg, dem, particles):
    """Initialize fields and update particles flow depth

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    dem : dict
        dictionary with dem information
    particles : dict
        particles dictionary at initial time step

    Returns
    -------
    particles : dict
        particles dictionary at initial time step updated with the flow depth
    fields : dict
        fields dictionary at initial time step
    """
    # read dem header
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    PFV = np.zeros((nrows, ncols))
    PP = np.zeros((nrows, ncols))
    FD = np.zeros((nrows, ncols))
    fields = {}
    fields['pfv'] = PFV
    fields['ppr'] = PP
    fields['pfd'] = FD
    fields['FV'] = PFV
    fields['P'] = PP
    fields['FD'] = FD
    fields['Vx'] = PFV
    fields['Vy'] = PFV
    fields['Vz'] = PFV

    particles = DFAfunC.getNeighborsC(particles, dem)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)

    return particles, fields


def placeParticles(massCell, indx, indy, csz, massPerPart, rng, initPartDistType):
    """ Create particles in given cell

    Compute number of particles to create in a given cell.
    Place particles in cell according to the chosen pattern (random semirandom
    or ordered)

    Parameters
    ----------
    massCell: float
        mass of snow in cell
    indx: int
        column index of the cell
    indy: int
        row index of the cell
    csz : float
        cellsize
    massPerPart : float
        maximum mass per particle
    initPartDistType: str
        initial particle distribution: random, semiRandom, uniform

    Returns
    -------
    xpart : 1D numpy array
        x position of particles
    ypart : 1D numpy array
        y position of particles
    mPart : 1D numpy array
        mass of particles
    nPart : int
        number of particles created
    """

    initPartDistType = initPartDistType.lower()

    if initPartDistType == 'random':
        # number of particles needed (floating number)
        nFloat = massCell / massPerPart
        nPart = int(np.floor(nFloat))
        # adding 1 with a probability of the residual proba
        proba = nFloat - nPart
        if rng.random(1) < proba:
            nPart = nPart + 1
        nPart = np.maximum(nPart, 1)
    else:
        n = (np.ceil(np.sqrt(massCell / massPerPart))).astype('int')
        nPart = n*n
        d = csz/n
        pos = np.linspace(0., csz-d, n) + d/2.
        x, y = np.meshgrid(pos, pos)
        x = x.flatten()
        y = y.flatten()

    mPart = massCell / nPart

    # TODO make this an independent function
    #######################
    # start ###############
    if initPartDistType == 'semirandom':
        # place particles equaly distributed with a small variation
        xpart = csz * (- 0.5 + indx) + x + (rng.random(nPart) - 0.5) * d
        ypart = csz * (- 0.5 + indy) + y + (rng.random(nPart) - 0.5) * d
    elif initPartDistType == 'random':
        # place particles randomly in the cell
        xpart = csz * (rng.random(nPart) - 0.5 + indx)
        ypart = csz * (rng.random(nPart) - 0.5 + indy)
    else:
        # place particles equaly distributed
        xpart = csz * (- 0.5 + indx) + x
        ypart = csz * (- 0.5 + indy) + y
        if initPartDistType != 'uniform':
            log.warning('Chosen value for initial particle distribution type not available: %s uniform is used instead' %
                        initPartDistType)

    return xpart, ypart, mPart, nPart


def initializeMassEnt(dem, simTypeActual, entLine, reportAreaInfo, thresholdPointInPoly):
    """ Initialize mass for entrainment

    Parameters
    ----------
    dem: dict
        dem dictionary
    simTypeActual: str
        simulation type
    entLine: dict
        entrainment line dictionary
    reportAreaInfo: dict
        simulation area information dictionary
    thresholdPointInPoly: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside

    Returns
    -------
    entrMassRaster : 2D numpy array
        raster of available mass for entrainment
    reportAreaInfo: dict
        simulation area information dictionary completed with entrainment area info
    """
    # read dem header
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    if 'ent' in simTypeActual:
        entrainmentArea = entLine['fileName']
        log.info('Initializing entrainment area: %s' % (entrainmentArea))
        log.info('Entrainment area features: %s' % (entLine['Name']))
        entLine = prepareArea(entLine, dem, thresholdPointInPoly)
        entrMassRaster = entLine['rasterData']
        reportAreaInfo['entrainment'] = 'Yes'
    else:
        entrMassRaster = np.zeros((nrows, ncols))
        reportAreaInfo['entrainment'] = 'No'

    return entrMassRaster, reportAreaInfo


def initializeResistance(cfg, dem, simTypeActual, resLine, reportAreaInfo, thresholdPointInPoly):
    """ Initialize resistance matrix

    Parameters
    ----------
    dem: dict
        dem dictionary
    simTypeActual: str
        simulation type
    resLine: dict
        resistance line dictionary
    reportAreaInfo: dict
        simulation area information dictionary
    thresholdPointInPoly: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside

    Returns
    -------
    cResRaster : 2D numpy array
        raster of resistance coefficients
    reportAreaInfo: dict
        simulation area information dictionary completed with entrainment area info
    """
    d = cfg.getfloat('dRes')
    cw = cfg.getfloat('cw')
    sres = cfg.getfloat('sres')
    # read dem header
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    if simTypeActual in ['entres', 'res']:
        resistanceArea = resLine['fileName']
        log.info('Initializing resistance area: %s' % (resistanceArea))
        log.info('Resistance area features: %s' % (resLine['Name']))
        resLine = prepareArea(resLine, dem, thresholdPointInPoly)
        mask = resLine['rasterData']
        cResRaster = 0.5 * d * cw / (sres*sres) * mask
        reportAreaInfo['resistance'] = 'Yes'
    else:
        cResRaster = np.zeros((nrows, ncols))
        reportAreaInfo['resistance'] = 'No'

    return cResRaster, reportAreaInfo


def DFAIterate(cfg, particles, fields, dem):
    """ Perform time loop for DFA simulation

     Save results at desired intervals

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at initial time step
        secondaryReleaseParticles : list
            list of secondary release area particles dictionaries at initial time step
    fields : dict
        fields dictionary at initial time step
    dem : dict
        dictionary with dem information

    Returns
    -------
    particlesList : list
        list of particles dictionary
    fieldsList : list
        list of fields dictionary (for each time step saved)
    Tcpu : dict
        computation time dictionary
    infoDict : dict
        Dictionary of all simulations carried out
    """

    cfgGen = cfg['GENERAL']
    # Initialise cpu timing
    Tcpu = {}
    Tcpu['TimeLoop'] = 0
    Tcpu['Force'] = 0.
    Tcpu['ForceSPH'] = 0.
    Tcpu['Pos'] = 0.
    Tcpu['Neigh'] = 0.
    Tcpu['Field'] = 0.

    # Load configuration settings
    tEnd = cfgGen.getfloat('tEnd')
    dtSave = fU.splitTimeValueToArrayInterval(cfgGen)
    sphOption = cfgGen.getint('sphOption')
    log.debug('using sphOption %s:' % sphOption)
    # desired output fields
    resTypes = fU.splitIniValueToArraySteps(cfgGen['resType'])
    # make sure to save all desiered resuts for first and last time step for
    # the report
    resTypesReport = fU.splitIniValueToArraySteps(cfg['REPORT']['plotFields'])
    resTypesLast = list(set(resTypes + resTypesReport))
    # derive friction type
    # turn friction model into integer
    frictModelsList = ['samosat', 'coulomb', 'voellmy']
    frictModel = cfgGen['frictModel'].lower()
    frictType = frictModelsList.index(frictModel) + 1
    log.debug('Friction Model used: %s, %s' % (frictModelsList[frictType-1], frictType))

    # Initialise Lists to save fields and add initial time step
    particlesList = []
    fieldsList = []
    timeM = []
    massEntrained = []
    massTotal = []

    # time stepping scheme info
    if featLF:
        log.debug('Use LeapFrog time stepping')
    else:
        log.debug('Use standard time stepping')
    # Initialize time and counters
    nSave = 1
    Tcpu['nSave'] = nSave
    nIter = 1
    nIter0 = 1
    iterate = True
    particles['iterate'] = iterate
    t = particles['t']
    log.debug('Saving results for time step t = %f s', t)
    fieldsList, particlesList = appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypes)
    # add initial time step to Tsave array
    Tsave = [0]
    # derive time step for first iteration
    if cfgGen.getboolean('cflTimeStepping') and nIter > cfgGen.getfloat('cflIterConstant'):
        # overwrite the dt value in the cfg
        dt = tD.getcflTimeStep(particles, dem, cfgGen)
    else:
        # get time step
        dt = cfgGen.getfloat('dt')

    t = t + dt

    # Start time step computation
    while t <= tEnd*(1.+1.e-13) and iterate:
        startTime = time.time()
        log.debug('Computing time step t = %f s, dt = %f s' % (t, dt))
        # Perform computations
        if featLF:
            particles, fields, Tcpu, dt = computeLeapFrogTimeStep(
                cfgGen, particles, fields, dt, dem, Tcpu)
        else:
            particles, fields, Tcpu = computeEulerTimeStep(
                cfgGen, particles, fields, dt, dem, Tcpu, frictType)

        Tcpu['nSave'] = nSave
        particles['t'] = t
        iterate = particles['iterate']

        # write mass balance info
        massEntrained.append(particles['massEntrained'])
        massTotal.append(particles['mTot'])
        timeM.append(t)
        # print progress to terminal
        print("time step t = %f s\r" % t, end="")
        # make sure the array is not empty
        if t >= dtSave[0]:
            Tsave.append(t)
            log.debug('Saving results for time step t = %f s', t)
            log.debug('MTot = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
            log.debug(('cpu time Force = %s s' % (Tcpu['Force'] / nIter)))
            log.debug(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / nIter)))
            log.debug(('cpu time Position = %s s' % (Tcpu['Pos'] / nIter)))
            log.debug(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / nIter)))
            log.debug(('cpu time Fields = %s s' % (Tcpu['Field'] / nIter)))
            fieldsList, particlesList = appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypes)
            if dtSave.size == 1:
                dtSave = [2*cfgGen.getfloat('tEnd')]
            else:
                indSave = np.where(dtSave > t)
                dtSave = dtSave[indSave]

        # derive time step
        if cfgGen.getboolean('cflTimeStepping') and nIter > cfgGen.getfloat('cflIterConstant'):
            # overwrite the dt value in the cfg
            dt = tD.getcflTimeStep(particles, dem, cfgGen)
        else:
            # get time step
            dt = cfgGen.getfloat('dt')

        t = t + dt
        nIter = nIter + 1
        nIter0 = nIter0 + 1
        tcpuTimeLoop = time.time() - startTime
        Tcpu['TimeLoop'] = Tcpu['TimeLoop'] + tcpuTimeLoop

    Tcpu['nIter'] = nIter
    log.info('Ending computation at time t = %f s', t-dt)
    log.debug('Saving results for time step t = %f s', t-dt)
    log.info('MTot = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
    log.info('Computational performances:')
    log.info(('cpu time Force = %s s' % (Tcpu['Force'] / nIter)))
    log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / nIter)))
    log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / nIter)))
    log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / nIter)))
    log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / nIter)))
    log.info(('cpu time TimeLoop = %s s' % (Tcpu['TimeLoop'] / nIter)))
    log.info(('cpu time total other = %s s' %
              ((Tcpu['Force'] + Tcpu['ForceSPH'] + Tcpu['Pos'] + Tcpu['Neigh'] + Tcpu['Field']) /
               nIter)))
    Tsave.append(t-dt)
    fieldsList, particlesList = appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypesLast)

    # create infoDict for report and mass log file
    infoDict = {'massEntrained': massEntrained, 'timeStep': timeM, 'massTotal': massTotal, 'Tcpu': Tcpu,
                'final mass': massTotal[-1], 'initial mass': massTotal[0], 'entrained mass': np.sum(massEntrained),
                'entrained volume': (np.sum(massEntrained)/cfgGen.getfloat('rhoEnt'))}

    # determine if stop criterion is reached or end time
    stopCritNotReached = particles['iterate']
    avaTime = particles['t']
    stopCritPer = cfgGen.getfloat('stopCrit') * 100.
    # update info dict with stopping info for report
    if stopCritNotReached:
        infoDict.update({'stopInfo': {'Stop criterion': 'end Time reached: %.2f' % avaTime, 'Avalanche run time [s]': '%.2f' % avaTime}})
    else:
        infoDict.update({'stopInfo': {'Stop criterion': '< %.2f percent of PKE' % stopCritPer, 'Avalanche run time [s]': '%.2f' % avaTime}})

    return Tsave, particlesList, fieldsList, infoDict


def appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypes):
    """ append fields and optionally particle dictionaries to list for export

        Parameters
        ------------
        particles: dict
            dictionary with particle properties
        fields: dict
            dictionary with all result type fields
        resTypes: list
            list with all result types that shall be exported

        Returns
        -------
        Fields: list
            updated list with desired result type fields dictionary
        Particles: list
            updated list with particles dicionaries

    """

    fieldAppend = {}
    for resType in resTypes:
        if resType == 'particles':
            particlesList.append(copy.deepcopy(particles))
        elif resType != '':
            fieldAppend[resType] = copy.deepcopy(fields[resType])
    fieldsList.append(fieldAppend)

    return fieldsList, particlesList


def writeMBFile(infoDict, avaDir, logName):
    """ write mass balance info to file

        Parameters
        -----------
        infoDict: dict
            info on mass
        avaDir: str or pathlib path
            path to avalanche directory
        logName: str
            simulation name
    """

    t = infoDict['timeStep']
    massEntrained = infoDict['massEntrained']
    massTotal = infoDict['massTotal']

    # write mass balance info to log file
    massDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA')
    fU.makeADir(massDir)
    with open(massDir / ('mass_%s.txt' % logName), 'w') as mFile:
        mFile.write('time, current, entrained\n')
        for m in range(len(t)):
            mFile.write('%.02f,    %.06f,    %.06f\n' %
                    (t[m], massTotal[m], massEntrained[m]))


def computeEulerTimeStep(cfg, particles, fields, dt, dem, Tcpu, frictType):
    """ compute next time step using an euler forward scheme

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at t
    fields : dict
        fields dictionary at t
    dt : float
        time step
    dem : dict
        dictionary with dem information
    Tcpu : dict
        computation time dictionary
    frictType: int
        indicator for chosen type of friction model

    Returns
    -------
    particles : dict
        particles dictionary at t + dt
    fields : dict
        fields dictionary at t + dt
    Tcpu : dict
        computation time dictionary
    """
    # get forces
    startTime = time.time()

    # loop version of the compute force
    log.debug('Compute Force C')
    particles, force, fields = DFAfunC.computeForceC(cfg, particles, fields, dem, dt, frictType)
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce

    # compute lateral force (SPH component of the calculation)
    startTime = time.time()
    if cfg.getint('sphOption') == 0:
        force['forceSPHX'] = np.zeros(np.shape(force['forceX']))
        force['forceSPHY'] = np.zeros(np.shape(force['forceY']))
        force['forceSPHZ'] = np.zeros(np.shape(force['forceZ']))
    else:
        log.debug('Compute Force SPH C')
        particles, force = DFAfunC.computeForceSPHC(cfg, particles, force, dem, gradient=0)
    tcpuForceSPH = time.time() - startTime
    Tcpu['ForceSPH'] = Tcpu['ForceSPH'] + tcpuForceSPH

    # update velocity and particle position
    startTime = time.time()
    # particles = updatePosition(cfg, particles, dem, force)
    log.debug('Update position C')
    particles = DFAfunC.updatePositionC(cfg, particles, dem, force, dt)
    tcpuPos = time.time() - startTime
    Tcpu['Pos'] = Tcpu['Pos'] + tcpuPos

    # release secondary release area?
    if particles['secondaryReleaseInfo']['flagSecondaryRelease'] == 'Yes':
        particles = releaseSecRelArea(cfg, particles, fields, dem)

    # get particles location (neighbours for sph)
    startTime = time.time()
    log.debug('get Neighbours C')
    particles = DFAfunC.getNeighborsC(particles, dem)

    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

    # update fields (compute grid values)
    startTime = time.time()
    log.debug('update Fields C')
    # particles, fields = updateFields(cfg, particles, force, dem, fields)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)
    tcpuField = time.time() - startTime
    Tcpu['Field'] = Tcpu['Field'] + tcpuField

    return particles, fields, Tcpu


def computeLeapFrogTimeStep(cfg, particles, fields, dt, dem, Tcpu):
    """ compute next time step using a Leap Frog scheme


    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at t
    fields : dict
        fields dictionary at t
    dt : float
        time step
    dem : dict
        dictionary with dem information
    Ment : 2D numpy array
        entrained mass raster
    Cres : 2D numpy array
        resistance raster
    Tcpu : dict
        computation time dictionary

    Returns
    -------
    particles : dict
        particles dictionary at t + dt
    fields : dict
        fields dictionary at t + dt
    Tcpu : dict
        computation time dictionary
    dt : float
        time step
    """

    # start timing
    startTime = time.time()
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce

    # dtK5 is half time step
    dtK5 = 0.5 * dt
    # cfg['dt'] = str(dtK5)
    log.debug('dt used now is %f' % dt)

    # load required DEM and mesh info
    csz = dem['header']['cellsize']
    Nx = dem['Nx']
    Ny = dem['Ny']
    Nz = dem['Nz']

    # particle properties
    mass = particles['m']
    xK = particles['x']
    yK = particles['y']
    zK = particles['z']
    uxK = particles['ux']
    uyK = particles['uy']
    uzK = particles['uz']

    # +++++++++++++Time integration using leapfrog 'Drift-Kick-Drif' scheme+++++
    # first predict position at time t_(k+0.5)
    # 'DRIFT'
    xK5 = xK + dt * 0.5 * uxK
    yK5 = yK + dt * 0.5 * uyK
    zK5 = zK + dt * 0.5 * uzK
    # update position from particles
    particles['x'] = xK5
    particles['y'] = yK5
    # For now z-position is taken from DEM - no detachment enforces...
    particles, _ = geoTrans.projectOnRaster(dem, particles, interp='bilinear')
    # TODO: do we need to update also h from particles?? I think yes! also mass, ent, res
    # particles['h'] = ?

    # 'KICK'
    # compute velocity at t_(k+0.5)
    # first compute force at t_(k+0.5)
    startTime = time.time()
    # TODO check  effect of artificial viscosity - update of velocity works here too
    particles, force = DFAfunC.computeForceC(cfg, particles, fields, dem, dtK5)
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce
    startTime = time.time()
    particles, force = DFAfunC.computeForceSPHC(cfg, particles, force, dem)
    tcpuForceSPH = time.time() - startTime
    Tcpu['ForceSPH'] = Tcpu['ForceSPH'] + tcpuForceSPH
    # particles, force = computeForceSPH(cfg, particles, force, dem)
    mass = particles['m']
    uxNew = uxK + (force['forceX'] + force['forceSPHX']) * dt / mass
    uyNew = uyK + (force['forceY'] + force['forceSPHY']) * dt / mass
    uzNew = uzK + (force['forceZ'] + force['forceSPHZ']) * dt / mass

    # 'DRIF'
    # now update position at t_(k+ 1)
    xNew = xK5 + dtK5 * uxNew
    yNew = yK5 + dtK5 * uyNew
    zNew = zK5 + dtK5 * uzNew

    # ++++++++++++++UPDATE Particle Properties
    # update mass required if entrainment
    massNew = mass + force['dM']
    particles['mTot'] = np.sum(massNew)
    particles['x'] = xNew
    particles['y'] = yNew
    particles['s'] = particles['s'] + np.sqrt((xNew-xK)*(xNew-xK) + (yNew-yK)*(yNew-yK))
    # make sure particle is on the mesh (recompute the z component)
    particles, _ = geoTrans.projectOnRaster(dem, particles, interp='bilinear')

    nx, ny, nz = DFAtls.getNormalArray(xNew, yNew, Nx, Ny, Nz, csz)
    nx, ny, nz = DFAtls.normalize(nx, ny, nz)
    particles['m'] = massNew
    # normal component of the velocity
    uN = uxNew*nx + uyNew*ny + uzNew*nz
    # remove normal component of the velocity
    particles['ux'] = uxNew - uN * nx
    particles['uy'] = uyNew - uN * ny
    particles['uz'] = uzNew - uN * nz

    #################################################################
    # this is dangerous!!!!!!!!!!!!!!
    ###############################################################
    # remove particles that are not located on the mesh any more
    particles = DFAtls.removeOutPart(particles, dem, dt)

    # ++++++++++++++GET particles location (neighbours for sph)
    startTime = time.time()
    particles = DFAfunC.getNeighborsC(particles, dem)
    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

    # ++++++++++++++UPDATE FIELDS (compute grid values)
    # update fields (compute grid values)
    startTime = time.time()
    # particles, fields = updateFields(cfg, particles, force, dem, fields)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)
    tcpuField = time.time() - startTime
    Tcpu['Field'] = Tcpu['Field'] + tcpuField

    return particles, fields, Tcpu, dt


def prepareArea(line, dem, radius, thList='', combine=True, checkOverlap=True):
    """ convert shape file polygon to raster

    Parameters
    ----------
    line: dict
        line dictionary
    dem : dict
        dictionary with dem information
    radius : float
        include all cells which center is in the polygon or close enough
    thList: list
        thickness values for all features in the line dictionary
    combine : Boolean
        if True sum up the rasters in the area list to return only 1 raster
        if False return the list of distinct area rasters
        this option works only if thList is not empty
    checkOverlap : Boolean
        if True check if features are overlaping and return an error if it is the case
        if False check if features are overlaping and average the value for overlaping areas

    Returns
    -------
    updates the line dictionary with the rasterData: Either
        Raster : 2D numpy array
            raster of the area (returned if relRHlist is empty OR if combine is set
            to True)
        or
        RasterList : list
            list of 2D numpy array rasters (returned if relRHlist is not empty AND
            if combine is set to True)
    """
    NameRel = line['Name']
    StartRel = line['Start']
    LengthRel = line['Length']
    RasterList = []

    for i in range(len(NameRel)):
        name = NameRel[i]
        start = StartRel[i]
        end = start + LengthRel[i]
        avapath = {}
        avapath['x'] = line['x'][int(start):int(end)]
        avapath['y'] = line['y'][int(start):int(end)]
        avapath['Name'] = name
        # if relTh is given - set relTh
        if thList != '':
            log.info('Release feature %s, relTh= %.2f' % (name, thList[i]))
            Raster = polygon2Raster(dem['header'], avapath, radius, th=thList[i])
        else:
            Raster = polygon2Raster(dem['header'], avapath, radius)
        RasterList.append(Raster)

    # if RasterList not empty check for overlap between features
    Raster = np.zeros(np.shape(dem['rasterData']))
    for rast in RasterList:
        ind1 = Raster > 0
        ind2 = rast > 0
        indMatch = np.logical_and(ind1, ind2)
        if indMatch.any():
            # if there is an overlap, raise error
            if checkOverlap:
                message = 'Features are overlaping - this is not allowed'
                log.error(message)
                raise AssertionError(message)
            else:
                # if there is an overlap, take average of values for the overlapping cells
                Raster = np.where(((Raster > 0) & (rast > 0)),  (Raster + rast)/2, Raster + rast)
        else:
            Raster = Raster + rast
    if debugPlot:
        debPlot.plotAreaDebug(dem, avapath, Raster)
    if combine:
        line['rasterData'] = Raster
        return line
    else:
        line['rasterData'] = RasterList
        return line


def polygon2Raster(demHeader, Line, radius, th=''):
    """ convert line to raster

    Parameters
    ----------
    demHeader: dict
        dem header dictionary
    Line : dict
        line dictionary
    radius : float
        include all cells which center is in the polygon or close enough
    th: float
        thickness value ot the line feature
    Returns
    -------

    Mask : 2D numpy array
        updated raster
    """
    # adim and center dem and polygon
    ncols = demHeader['ncols']
    nrows = demHeader['nrows']
    xllc = demHeader['xllcenter']
    yllc = demHeader['yllcenter']
    csz = demHeader['cellsize']
    xCoord0 = (Line['x'] - xllc) / csz
    yCoord0 = (Line['y'] - yllc) / csz
    if (xCoord0[0] == xCoord0[-1]) and (yCoord0[0] == yCoord0[-1]):
        xCoord = np.delete(xCoord0, -1)
        yCoord = np.delete(yCoord0, -1)
    else:
        xCoord = copy.deepcopy(xCoord0)
        yCoord = copy.deepcopy(yCoord0)
        xCoord0 = np.append(xCoord0, xCoord0[0])
        yCoord0 = np.append(yCoord0, yCoord0[0])

    # get the raster corresponding to the polygon
    polygon = np.stack((xCoord, yCoord), axis=-1)
    path = mpltPath.Path(polygon)
    # add a tolerance to include cells for which the center is on the lines
    # for this we need to know if the path is clockwise or counter clockwise
    # to decide if the radius should be positif or negatif in contains_points
    is_ccw = geoTrans.isCounterClockWise(path)
    r = (radius*is_ccw - radius*(1-is_ccw))
    x = np.linspace(0, ncols-1, ncols)
    y = np.linspace(0, nrows-1, nrows)
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    points = np.stack((X, Y), axis=-1)
    mask = path.contains_points(points, radius=r)
    Mask = mask.reshape((nrows, ncols)).astype(int)
    # thickness field is provided, then return array with ones
    if th != '':
        log.debug('REL set from dict, %.2f' % th)
        Mask = np.where(Mask > 0, th, 0.)
    else:
        Mask = np.where(Mask > 0, 1., 0.)

    if debugPlot:
        debPlot.plotRemovePart(xCoord0, yCoord0, demHeader, X, Y, Mask, mask)

    return Mask


def checkParticlesInRelease(particles, line, radius):
    """ remove particles laying outside the polygon

    Parameters
    ----------
    particles : dict
        particles dictionary
    line: dict
        line dictionary
    radius: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside

    Returns
    -------
    particles : dict
        particles dictionary where particles outside of the polygon have been removed

    """
    NameRel = line['Name']
    StartRel = line['Start']
    LengthRel = line['Length']
    Mask = np.full(np.size(particles['x']), False)
    for i in range(len(NameRel)):
        name = NameRel[i]
        start = StartRel[i]
        end = start + LengthRel[i]
        avapath = {}
        avapath['x'] = line['x'][int(start):int(end)]
        avapath['y'] = line['y'][int(start):int(end)]
        avapath['Name'] = name
        mask = pointInPolygon(line['header'], particles, avapath, radius)
        Mask = np.logical_or(Mask, mask)

    nRemove = len(Mask)-np.sum(Mask)
    if nRemove > 0:
        particles = DFAtls.removePart(particles, Mask, nRemove, 'because they are not within the release polygon')
        log.debug('removed %s particles because they are not within the release polygon' % (nRemove))

    return particles


def pointInPolygon(demHeader, points, Line, radius):
    """ find particles within a polygon

    Parameters
    ----------
    demHeader: dict
        dem header dictionary
    points: dict
        points to check
    Line : dict
        line dictionary
    radius: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside
    Returns
    -------

    Mask : 1D numpy array
        Mask of particles to keep
    """
    xllc = demHeader['xllcenter']
    yllc = demHeader['yllcenter']
    xCoord0 = (Line['x'] - xllc)
    yCoord0 = (Line['y'] - yllc)
    if (xCoord0[0] == xCoord0[-1]) and (yCoord0[0] == yCoord0[-1]):
        xCoord = np.delete(xCoord0, -1)
        yCoord = np.delete(yCoord0, -1)
    else:
        xCoord = copy.deepcopy(xCoord0)
        yCoord = copy.deepcopy(yCoord0)
        xCoord0 = np.append(xCoord0, xCoord0[0])
        yCoord0 = np.append(yCoord0, yCoord0[0])

    # get the raster corresponding to the polygon
    polygon = np.stack((xCoord, yCoord), axis=-1)
    path = mpltPath.Path(polygon)
    # add a tolerance to include cells for which the center is on the lines
    # for this we need to know if the path is clockwise or counter clockwise
    # to decide if the radius should be positif or negatif in contains_points
    is_ccw = geoTrans.isCounterClockWise(path)
    r = (radius*is_ccw - radius*(1-is_ccw))
    points2Check = np.stack((points['x'], points['y']), axis=-1)
    mask = path.contains_points(points2Check, radius=r)
    mask = np.where(mask > 0, True, False)

    if debugPlot:
        debPlot.plotPartAfterRemove(points, xCoord0, yCoord0, mask)

    return mask


def releaseSecRelArea(cfg, particles, fields, dem):
    """ Release secondary release area if trigered

    Initialize particles of the trigured secondary release area and add them
    to the simulation (particles dictionary)
    """

    secondaryReleaseInfo = particles['secondaryReleaseInfo']
    flowDepthField = fields['FD']
    secRelRasterList = secondaryReleaseInfo['rasterData']
    secRelRasterNameList = secondaryReleaseInfo['Name']
    count = 0
    indexRel = []
    for secRelRaster, secRelRasterName in zip(secRelRasterList, secRelRasterNameList):
        # do the two arrays intersect (meaning a flowing particle entered the
        # secondary release area)
        mask = (secRelRaster > 0) & (flowDepthField > 0)
        if mask.any():
            # create secondary release area particles
            log.info('Initializing secondary release area feature %s' % secRelRasterName)
            secRelInfo = shpConv.extractFeature(secondaryReleaseInfo, count)
            secRelInfo['rasterData'] = secRelRaster
            secRelParticles = initializeParticles(cfg, secRelInfo, dem)
            # release secondary release area by just appending the particles
            log.info('Releasing secondary release area %s at t = %.2f s' % (secRelRasterName, particles['t']))
            particles = DFAtls.mergeParticleDict(particles, secRelParticles)
            # save index of secRel feature
            indexRel.append(secRelRasterName)
        count = count + 1

    secondaryReleaseInfo['rasterData'] = secRelRasterList
    particles['secondaryReleaseInfo'] = secondaryReleaseInfo
    for item in indexRel:
        iR = secRelRasterNameList.index(item)
        # remove it from the secondary release area list
        secRelRasterList.pop(iR)
        secondaryReleaseInfo = shpConv.removeFeature(secondaryReleaseInfo, iR)
        secRelRasterNameList.pop(iR)

    # update secondaryReleaseInfo
    secondaryReleaseInfo['rasterData'] = secRelRasterList
    particles['secondaryReleaseInfo'] = secondaryReleaseInfo

    return particles


def savePartToPickle(dictList, outDir, logName):
    """ Save each dictionary from a list to a pickle in outDir; works also for one dictionary instead of list

        Parameters
        ---------
        dictList: list or dict
            list of dictionaries or single dictionary
        outDir: str
            path to output directory
        logName : str
            simulation Id
    """

    if isinstance(dictList, list):
        for dict in dictList:
            pickle.dump(dict, open(outDir / ("particles_%s_%09.4f.p" % (logName, dict['t'])), "wb"))
    else:
        pickle.dump(dictList, open(outDir / ("particles_%s_%09.4f.p" % (logName, dictList['t'])), "wb"))


def readPartFromPickle(inDir, flagAvaDir=False):
    """ Read pickles within a directory and return List of dicionaries read from pickle

        Parameters
        -----------
        inDir: str
            path to input directory
        flagAvaDir: bool
            if True inDir corresponds to an avalanche directory and pickles are
            read from avaDir/Outputs/com1DFA/particles
    """

    if flagAvaDir:
        inDir = pathlib.Path(inDir, 'Outputs', 'com1DFA', 'particles')

    # search for all pickles within directory
    PartDicts = sorted(list(inDir.glob('*.p')))

    # initialise list of particle dictionaries
    Particles = []
    TimeStepInfo = []
    for particles in PartDicts:
        particles = pickle.load(open(particles, "rb"))
        Particles.append(particles)
        TimeStepInfo.append(particles['t'])

    return Particles, TimeStepInfo


def savePartToCsv(particleProperties, dictList, outDir):
    """ Save each particle dictionary from a list to a csv file; works also for one dictionary instead of list

        Parameters
        ---------
        particleProperties: str
            all particle properties that shall be saved to csv file (e.g.: m, velocityMagnitude, ux,..)
        dictList: list or dict
            list of dictionaries or single dictionary
        outDir: str
            path to output directory; particlesCSV will be created in this outDir
    """

    # set output directory
    outDir = outDir / 'particlesCSV'
    fU.makeADir(outDir)

    # read particle properties to be saved
    particleProperties = particleProperties.split('|')

    # write particles locations and properties to csv file
    nParticles = len(dictList)
    count = 0
    for m in range(nParticles):
        particles = dictList[count]
        simName = particles['simName']
        csvData = {}
        csvData['X'] = particles['x'] + particles['xllcenter']
        csvData['Y'] = particles['y'] + particles['yllcenter']
        csvData['Z'] = particles['z']

        for partProp in particleProperties:
            if partProp == 'velocityMagnitude':
                ux = particles['ux']
                uy = particles['uy']
                uz = particles['uz']
                csvData[partProp] = DFAtls.norm(ux, uy, uz)
            else:
                if partProp in particles:
                    csvData[partProp] = particles[partProp]
                else:
                    log.warning('Particle property %s does not exist' % (partProp))
        csvData['time'] = particles['t']

        # create pandas dataFrame and save to csv
        outFile = outDir / ('particles%s.csv.%d' % (simName, count))
        particlesData = pds.DataFrame(data=csvData)
        particlesData.to_csv(outFile, index=False)
        count = count + 1


def exportFields(cfg, Tsave, fieldsList, demOri, outDir, logName):
    """ export result fields to Outputs directory according to result parameters and time step
        that can be specified in the configuration file

        Parameters
        -----------
        cfg: dict
            configurations
        Tsave: list
            list of time step that corresponds to each dict in Fields
        Fields: list
            list of Fields for each dtSave
        outDir: str
            outputs Directory


        Returns
        --------
        exported peak fields are saved in Outputs/com1DFA/peakFiles

    """

    resTypesGen = fU.splitIniValueToArraySteps(cfg['GENERAL']['resType'])
    resTypesReport = fU.splitIniValueToArraySteps(cfg['REPORT']['plotFields'])
    if resTypesGen == []:
        resTypesGen = resTypesReport
    if 'particles' in resTypesGen:
        resTypesGen.remove('particles')
    numberTimes = len(Tsave)-1
    countTime = 0
    for timeStep in Tsave:
        if (countTime == numberTimes):
            # for last time step we need to add the report fields
            resTypes = list(set(resTypesGen + resTypesReport))
        else:
            resTypes = resTypesGen
        for resType in resTypes:
            resField = fieldsList[countTime][resType]
            if resType == 'ppr':
                resField = resField * 0.001
            dataName = logName + '_' + resType + '_' + 't%.2f' % (Tsave[countTime]) + '.asc'
            # create directory
            outDirPeak = outDir / 'peakFiles' / 'timeSteps'
            fU.makeADir(outDirPeak)
            outFile = outDirPeak / dataName
            IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
            if countTime == numberTimes:
                log.debug('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f - FINAL time step ' % (resType, Tsave[countTime]))
                dataName = logName + '_' + resType + '.asc'
                # create directory
                outDirPeakAll = outDir / 'peakFiles'
                fU.makeADir(outDirPeakAll)
                outFile = outDirPeakAll / dataName
                IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
            else:
                log.debug('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f ' % (resType, Tsave[countTime]))
        countTime = countTime + 1


def prepareVarSimDict(standardCfg, inputSimFiles, variationDict):
    """ Prepare a dictionary with simulations that shall be run with varying parameters following the variation dict

        Parameters
        -----------
        standardCfg : configParser object
            default configuration or local configuration
        inputSimFiles: dict
            info dict on available input data
        variationDict: dict
            dictionary with parameter to be varied as key and list of it's values


        Returns
        -------
        simDict: dict
            dicionary with info on simHash, releaseScenario, release area file path,
            simType and contains full configuration configparser object for simulation run

    """

    # get list of simulation types that are desired
    if 'simTypeList' in variationDict:
        simTypeList = variationDict['simTypeList']
        del variationDict['simTypeList']
    else:
        simTypeList = standardCfg['GENERAL']['simTypeList'].split('|')
    # get a list of simulation types that are desired AND available
    simTypeList = getSimTypeList(simTypeList, inputSimFiles)

    # set simTypeList (that has been checked if available) as parameter in variationDict
    variationDict['simTypeList'] = simTypeList
    # create a dataFrame with all possible combinations of the variationDict values
    variationDF = pd.DataFrame(product(*variationDict.values()), columns=variationDict.keys())

    # generate a list of full simulation info for all release area scenarios and simTypes
    # simulation info must contain: simName, releaseScenario, relFile, configuration as dictionary
    simDict = {}
    for rel in inputSimFiles['relFiles']:
        relName = rel.stem
        if '_' in relName:
            relNameSim = relName + '_AF'
        else:
            relNameSim = relName
        cfgSim = cfgUtils.convertConfigParserToDict(standardCfg)
        for row in variationDF.itertuples():
            for parameter in variationDict:
                cfgSim['GENERAL'][parameter] = row._asdict()[parameter]
            cfgSim['GENERAL']['simTypeActual'] = row._asdict()['simTypeList']
            cfgSim['GENERAL']['releaseScenario'] = relName
            # convert back to configParser object
            cfgSimObject = cfgUtils.convertDictToConfigParser(cfgSim)
            # create unique hash for simulation configuration
            simHash = cfgUtils.cfgHash(cfgSimObject)
            simName = relNameSim + '_' + row._asdict()['simTypeList'] + '_' + cfgSim['GENERAL']['modelType'] + '_' + simHash
            simDict[simName] = {'simHash': simHash, 'releaseScenario': relName,
                                'simType': row._asdict()['simTypeList'], 'relFile': rel,
                                'cfgSim': cfgSimObject}

    return simDict


def getSimTypeList(simTypeList, inputSimFiles):
    """ Define available simulation types of requested types

        Parameters
        -----------
        standardCfg : configParser object
            default configuration or local configuration
        inputSimFiles: dict
            info dict on available input data

        Returns
        --------
        simTypeList: list
            list of requested simTypes where also the required input data is available

    """

    # read entrainment resistance info
    entResInfo = inputSimFiles['entResInfo']

    # define simulation type
    if 'available' in simTypeList:
        if entResInfo['flagEnt'] == 'Yes' and entResInfo['flagRes'] == 'Yes':
            simTypeList.append('entres')
        elif entResInfo['flagEnt'] == 'Yes' and entResInfo['flagRes'] == 'No':
            simTypeList.append('ent')
        elif entResInfo['flagEnt'] == 'No' and entResInfo['flagRes'] == 'Yes':
            simTypeList.append('res')
        # always add null simulation
        simTypeList.append('null')
        simTypeList.remove('available')

    # remove duplicate entries
    simTypeList = set(simTypeList)
    simTypeList = sorted(list(simTypeList), reverse=False)

    if 'ent' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagEnt'] == 'No':
            message = 'No entrainment file found'
            log.error(message)
            raise FileNotFoundError(message)
    if 'res' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagRes'] == 'No':
            message = 'No resistance file found'
            log.error(message)
            raise FileNotFoundError(message)

    return simTypeList
