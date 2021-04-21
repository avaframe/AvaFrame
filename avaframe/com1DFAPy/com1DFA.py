"""
    Main functions for python DFA kernel
"""

import logging
import time
import os
import numpy as np
import glob
import copy
import pickle
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local imports
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.makePalette as makePalette
import avaframe.com1DFAPy.timeDiscretizations as tD
import avaframe.com1DFAPy.DFAtools as DFAtls
import avaframe.com1DFAPy.DFAfunctionsCython as DFAfunC
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils

#######################################
# Set flags here
#######################################
# create local logger
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')
# set feature flag for initial particle distribution
# particles are homegeneosly distributed with a little random variation
flagSemiRand = False
# particles are randomly distributed
flagRand = True
# set feature leapfrog time stepping
featLF = False
featCFL = False
featCFLConstrain = False


def com1DFAMain(cfg, avaDir, relThField):
    """ Run main model

    This will compute a dense flow avalanche

    Parameters
    ----------
    cfg : dict
        configuration read from ini file
    avaDir : str
        path to avalanche directory
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

    # add avalanache directory to configuration
    cfgGen['avalancheDir'] = avaDir

     # turn friction models into integers
    frictModelsList = ['samosAT', 'Coulomb']
    cfgGen['frictType'] = str(frictModelsList.index(cfgGen['frictModel']) + 1)
    log.info('Friction Model used: %s, %s' % (cfgGen['frictModel'], cfgGen['frictType']))

    # Log current avalanche directory
    log.debug('Your current avalanche name: %s' % avaDir)

    # Create output and work directories
    # set module name, reqiured as long we are in dev phase
    # - because need to create e.g. Output folder for com1DFAPy to distinguish from
    modName = 'com1DFAPy'
    workDir, outDir = inDirs.initialiseRunDirs(avaDir, modName)

    # Load input data
    flagDev = cfg['FLAGS'].getboolean('flagDev')

    # fetch input data - dem, release-, entrainment- and resistance areas
    demFile, relFiles, entFiles, resFile, entResInfo = gI.getInputDataCom1DFAPy(
        avaDir, cfg['FLAGS'], flagDev)

    # Counter for release area loop
    countRel = 0

    # Setup simulation dictionaries for report genereation
    reportDictList = []

    # Loop through release areas
    for rel in relFiles:

        demOri, releaseLine, entLine, resLine, entrainmentArea, resistanceArea = prepareInputData(demFile, rel, entFiles, resFile)

        # find out which simulations to perform
        relName, cuSim, relDict, badName = getSimulation(cfg, rel, entResInfo)
        releaseLine['d0'] = relDict['d0']

        for sim in cuSim:
            #logName = sim + '_' + cfgGen['mu']
            logName = relName + '_' + sim + '_dfa_' + cfgGen['mu']
            log.info('Perform %s simulation' % logName)
            # add simType to configuration
            cfgGen['simTypeActual'] = sim

            # +++++++++PERFORM SIMULAITON++++++++++++++++++++++
            # for timing the sims
            startTime = time.time()
            particles, fields, dem, reportAreaInfo = initializeSimulation(cfg, demOri, releaseLine, entLine, resLine, logName, relThField)
            relFiles = releaseLine['file']
            # ------------------------
            #  Start time step computation
            Tsave, Particles, Fields, infoDict = DFAIterate(cfgGen, particles, fields, dem)

            # write mass balance to File
            writeMBFile(infoDict, avaDir, logName)

            tcpuDFA = '%.2f' % (time.time() - startTime)
            log.info(('cpu time DFA = %s s' % (tcpuDFA)))

            if 'particles' in cfgGen['resType']:
                # export particles dictionaries of saving time steps
                outDirData = os.path.join(outDir, 'particles')
                fU.makeADir(outDirData)
                savePartToPickle(Particles, outDirData)

            # Result parameters to be exported
            exportFields(cfgGen, Tsave, Fields, rel, demOri, outDir, logName)

            # write report dictionary
            reportDict = createReportDict(avaDir, logName, relName, relDict, cfgGen, entrainmentArea, resistanceArea, reportAreaInfo)
            # add time and mass info to report
            reportDict = reportAddTimeMassInfo(reportDict, tcpuDFA, cfgGen, infoDict)

            # Add to report dictionary list
            reportDictList.append(reportDict)

            # Count total number of simulations
            countRel = countRel + 1
    log.debug('Avalanche Simulations performed')

    return Particles, Fields, Tsave, dem, reportDictList


def getSimulation(cfg, rel, entResInfo):
    """ get Simulation to run for a given release


    Parameters
    ----------
    cfgFlags : dict
        Flags configuration parameters
    rel : str
        path to release file

    Returns
    -------
    relName : str
        release name
    cuSim : list
        list of simulations to run
    relDict : list
        release dictionary
    badName : boolean
        changed release name
    """

    cfgFlags = cfg['FLAGS']

    # read list of desired simulation types
    simTypeList = cfg['GENERAL']['simTypeList'].split('|')

    # Set release areas and simulation name
    relName = os.path.splitext(os.path.basename(rel))[0]
    simName = relName
    badName = False
    if '_' in relName:
        badName = True
        log.warning('Release area scenario file name includes an underscore \
        the suffix _AF will be added')
        simName = relName + '_AF'
    relDict = shpConv.SHP2Array(rel)
    for k in range(len(relDict['d0'])):
        if relDict['d0'][k] == 'None':
            relDict['d0'][k] = cfg['GENERAL'].getfloat('relTh')
        else:
            relDict['d0'][k] = float(relDict['d0'][k])

    log.info('Release area scenario: %s - perform simulations' % (relName))

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

    if 'ent' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagEnt'] == 'No':
            log.error('No entrainment file found')
            raise FileNotFoundError
    if 'res' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagRes'] == 'No':
            log.error('No resistance file found')
            raise FileNotFoundError

    return relName, simTypeList, relDict, badName


def prepareInputData(demFile, relFile, entFiles, resFile):
    """ Fetch input data

    Parameters
    ----------
    demFile : str
        path to dem file
    relFiles : str
        path to release file
    entFiles : str
        path to entrainment file
    resFile : str
        path to resistance file

    Returns
    -------
    demOri : dict
        dictionary with original dem
    releaseLine : dict
        release line dictionary
    entLine : dict
        entrainment line dictionary
    resLine : dict
        resistance line dictionary
    entrainmentArea : str
        entrainment file name
    resistanceArea : str
        resistance file name
    """
    # get dem information
    demOri = IOf.readRaster(demFile)
    # get line from release area polygon
    releaseLine = shpConv.readLine(relFile, 'release1', demOri)
    releaseLine['file'] = relFile
    # get line from entrainement area polygon
    if entFiles:
        entLine = shpConv.readLine(entFiles, '', demOri)
        entrainmentArea = os.path.splitext(os.path.basename(entFiles))[0]
        entLine['Name'] = [entrainmentArea]
    else:
        entLine = None
        entrainmentArea = ''
    # get line from resistance area polygon
    if resFile:
        resLine = shpConv.readLine(resFile, '', demOri)
        resistanceArea = os.path.splitext(os.path.basename(resFile))[0]
        resLine['Name'] = [resistanceArea]
    else:
        resLine = None
        resistanceArea = ''

    return demOri, releaseLine, entLine, resLine, entrainmentArea, resistanceArea


def createReportDict(avaDir, logName, relName, relDict, cfgGen, entrainmentArea, resistanceArea, reportAreaInfo):
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

    # load parameters set in configuration file
    dateTimeInfo = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    entInfo = reportAreaInfo['entrainment']
    resInfo = reportAreaInfo['resistance']

    # Create dictionary
    reportST = {}
    reportST = {}
    reportST = {'headerLine': {'type': 'title', 'title': 'com1DFA Simulation'},
                'avaName': {'type': 'avaName', 'name': avaDir},
                'simName': {'type': 'simName', 'name': logName},
                'time': {'type': 'time', 'time': dateTimeInfo},
                'Simulation Parameters': {
                'type': 'list',
                'Program version': 'developement',
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


def reportAddTimeMassInfo(reportDict, tcpuDFA, cfgGen, infoDict):
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

    dem = setDEMoriginToZero(demOri)

    # read dem header
    headerDEM = dem['header']
    nColsDEM = headerDEM.ncols
    nRowsDEM = headerDEM.nrows
    cszDEM = headerDEM.cellsize

    # get normal vector of the grid mesh
    Nx, Ny, Nz = DFAtls.getNormalMesh(dem, num)
    dem['Nx'] = np.where(np.isnan(Nx), 0., Nx)
    dem['Ny'] = np.where(np.isnan(Ny), 0., Ny)
    # build no data mask (used to find out of dem particles)
    bad = np.where(np.isnan(Nx), True, False)
    dem['Nz'] = Nz
    dem['Bad'] = bad

    # Prepare SPH grid
    headerNeighbourGrid = IOf.cASCheader()
    cszNeighbourGrid = cfg.getfloat('sphKernelRadius')
    headerNeighbourGrid.cellsize = cszNeighbourGrid
    headerNeighbourGrid.ncols = np.floor(nColsDEM * cszDEM / cszNeighbourGrid) + 1
    headerNeighbourGrid.nrows = np.floor(nRowsDEM * cszDEM / cszNeighbourGrid) + 1
    dem['headerNeighbourGrid'] = headerNeighbourGrid

    # get real Area
    areaRaster = DFAtls.getAreaMesh(Nx, Ny, Nz, cszDEM, num)
    dem['areaRaster'] = areaRaster
    projArea = nColsDEM * nRowsDEM * cszDEM * cszDEM
    actualArea = np.nansum(areaRaster)
    log.info('Largest cell area: %.2f m²' % (np.nanmax(areaRaster)))
    log.debug('Projected Area : %.2f' % projArea)
    log.debug('Total Area : %.2f' % actualArea)

    return dem


def setDEMoriginToZero(demOri):
    """ set origin of DEM to 0,0 """

    dem = copy.deepcopy(demOri)
    dem['header'].xllcenter = 0
    dem['header'].yllcenter = 0
    dem['header'].xllcorner = -dem['header'].cellsize/2
    dem['header'].yllcorner = -dem['header'].cellsize/2

    return dem


def initializeSimulation(cfg, demOri, releaseLine, entLine, resLine, logName, relThField):
    """ create simulaton report dictionary

    Parameters
    ----------
    cfg : str
        simulation scenario name
    demOri : dict
        dictionary with original dem
    releaseLine : dict
        release line dictionary
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
    fields : dict
        fields dictionary at initial time step
    dem : dict
        dictionary with new dem (lower left center at origin)
    """
    cfgGen = cfg['GENERAL']
    methodMeshNormal = cfg.getfloat('GENERAL', 'methodMeshNormal')

    # -----------------------
    # Initialize mesh
    dem = initializeMesh(cfgGen, demOri, methodMeshNormal)
    # ------------------------
    # process release info to get it as a raster
    if len(relThField) == 0:
        # if no release thickness field or function - set release according to shapefile or ini file
        relRaster = prepareArea(releaseLine, demOri, relThList=releaseLine['d0'])
    else:
        # if relTh provided - set release thickness with field or function
        relRaster = prepareArea(releaseLine, demOri)
        relRaster = relRaster * relThField

    # compute release area
    header = dem['header']
    csz = header.cellsize
    relRasterOnes = np.where(relRaster > 0, 1., 0.)
    relAreaActual = np.nansum(relRasterOnes*dem['areaRaster'])
    relAreaProjected = np.sum(csz*csz*relRasterOnes)
    reportAreaInfo = {'Release area info': {'Projected Area [m2]':  '%.2f' % (relAreaProjected),
                'Actual Area [m2]': '%.2f' % (relAreaActual)}}

    # ------------------------
    # initialize simulation

    # create particles, create resistance and
    # entrainment matrix, initialize fields, get normals and neighbours
    particles, fields = initializeParticles(cfgGen, relRaster, dem, logName=logName)

    # initialize entrainment and resistance
    # get info of simType and whether or not to initialize resistance and entrainment
    simTypeActual = cfgGen['simTypeActual']
    rhoEnt = cfgGen.getfloat('rhoEnt')
    hEnt = cfgGen.getfloat('hEnt')
    entrMassRaster, reportAreaInfo = initializeMassEnt(demOri, simTypeActual, entLine, reportAreaInfo)
    cResRaster, reportAreaInfo = initializeResistance(cfgGen, demOri, simTypeActual, resLine, reportAreaInfo)
    # surfacic entrainment mass available (unit kg/m²)
    fields['entrMassRaster'] = entrMassRaster*rhoEnt*hEnt
    entreainableMass = np.nansum(fields['entrMassRaster']*dem['areaRaster'])
    log.info('Mass available for entrainment: %.2f kg' % (entreainableMass))
    fields['cResRaster'] = cResRaster

    return particles, fields, dem, reportAreaInfo


def initializeParticles(cfg, relRaster, dem, logName=''):
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

    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    csz = header.cellsize
    areaRaster = dem['areaRaster']
    totalMassRaster = np.nansum(areaRaster*relRaster*rho)

    # initialize arrays
    partPerCell = np.zeros(np.shape(relRaster), dtype=np.int64)
    FD = np.zeros((nrows, ncols))
    # find all non empty cells (meaning release area)
    indRelY, indRelX = np.nonzero(relRaster)

    # derive mass per particle to define number of particles per cell:
    if massPerParticleDeterminationMethod == 'MPPDIR':
        massPerPart = cfg.getfloat('massPerPart')
        log.info('Number of particles defined by: mass per particle %s' % cfg['massPerPart'])
    elif massPerParticleDeterminationMethod == 'MPPDH':
        deltaTh = cfg.getfloat('deltaTh')
        ds = min(csz, cfg.getfloat('sphKernelRadius'))
        massPerPart = rho * ds * ds * deltaTh
        log.info('Number of particles defined by: release thickness per particle: %s' % cfg['deltaTh'])
        log.info('mass per particle is %.2f' % massPerPart)

    # make option available to read initial particle distribution from file
    if cfg.getboolean('initialiseParticlesFromFile'):
        # TODO: this is for development purposes, change or remove in the future
        # If initialisation from file
        if cfg['particleFile']:
            inDirPart = cfg['particleFile']
        else:
            inDirPart = os.path.join(avaDir, 'Outputs', 'com1DFA')

        partDirName = logName
        inDirPart = os.path.join(inDirPart, 'particles', partDirName)
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
        Apart = np.empty(0)
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
            xpart, ypart, mPart, nPart = placeParticles(massCell, indRelx, indRely, csz, massPerPart, rng)
            Npart = Npart + nPart
            partPerCell[indRely, indRelx] = nPart
            # initialize particles position, mass, height...
            NPPC = np.append(NPPC, nPart*np.ones(nPart))
            Apart = np.append(Apart, areaRaster[indRely, indRelx]*np.ones(nPart)/nPart)
            Xpart = np.append(Xpart, xpart)
            Ypart = np.append(Ypart, ypart)
            Mpart = np.append(Mpart, mPart * np.ones(nPart))
            Hpart = np.append(Hpart, hCell * np.ones(nPart))

        Hpart, _ = geoTrans.projectOnGrid(Xpart, Ypart, relRaster, csz=csz, interp='bilinear')
        Mpart = rho * Hpart * Apart
        # create dictionnary to store particles properties
        particles = {}
        particles['Npart'] = Npart
        particles['x'] = Xpart
        particles['y'] = Ypart
        particles['s'] = np.zeros(np.shape(Xpart))
        particles['l'] = np.zeros(np.shape(Xpart))
        # adding z component
        particles, _ = geoTrans.projectOnRaster(dem, particles, interp='bilinear')
        # readjust mass
        mTot = np.sum(Mpart)
        particles['m'] = Mpart*totalMassRaster/mTot

    particles['mTot'] = np.sum(particles['m'])
    particles['h'] = Hpart
    particles['NPPC'] = NPPC
    particles['hNearestNearest'] = Hpart
    particles['hNearestBilinear'] = Hpart
    particles['hBilinearNearest'] = Hpart
    particles['hBilinearBilinear'] = Hpart
    particles['hSPH'] = Hpart
    particles['GHX'] = np.zeros(np.shape(Xpart))
    particles['GHY'] = np.zeros(np.shape(Xpart))
    particles['GHZ'] = np.zeros(np.shape(Xpart))
    particles['ux'] = np.zeros(np.shape(Xpart))
    particles['uy'] = np.zeros(np.shape(Xpart))
    particles['uz'] = np.zeros(np.shape(Xpart))
    particles['stoppCriteria'] = False
    kineticEne = np.sum(0.5 * Mpart * DFAtls.norm2(particles['ux'], particles['uy'], particles['uz']))
    particles['kineticEne'] = kineticEne
    particles['potentialEne'] = np.sum(gravAcc * Mpart * particles['z'])
    particles['peakKinEne'] = kineticEne

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

    particles = DFAfunC.getNeighboursC(particles, dem)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)

    # initialize time
    t = 0
    particles['t'] = t

    relCells = np.size(indRelY)
    partPerCell = particles['Npart']/relCells
    log.info('Expeced mass. Mexpected = %.2f kg.' % (totalMassRaster))
    log.info('Initializted simulation. MTot = %.2f kg, %s particles in %.2f cells.' %
             (particles['mTot'], particles['Npart'], relCells))
    log.info('Mass per particle = %.2f kg and particles per cell = %.2f.' %
             (particles['mTot']/particles['Npart'], partPerCell))

    if debugPlot:
        x = np.arange(ncols) * csz
        y = np.arange(nrows) * csz
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        ref0, im = pU.NonUnifIm(ax, x, y, areaRaster, 'x [m]', 'y [m]',
                                extent=[x.min(), x.max(), y.min(), y.max()],
                                cmap=cmap, norm=None)

        ax.plot(Xpart, Ypart, 'or', linestyle='None')
        pU.addColorBar(im, ax, None, 'm²')
        plt.show()

    return particles, fields


def placeParticles(massCell, indx, indy, csz, massPerPart, rng):
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
    if flagRand:
        # number of particles needed (floating number)
        nFloat = massCell / massPerPart
        nPart = int(np.floor(nFloat))
        # adding 1 with a probability of the residual proba
        proba = nFloat - nPart
        if rng.random(1) < proba:
            nPart = nPart + 1
        #nPart = (nFloor + rng.binomial(1, proba)).astype('int')
        # TODO: ensure that there is at last one particle
        nPart = np.maximum(nPart, 1)
    else:
        n = (np.floor(np.sqrt(massCell / massPerPart)) + 1).astype('int')
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
    if flagSemiRand:
        # place particles equaly distributed with a small variation
        xpart = csz * (- 0.5 + indx) + x + (rng.random(nPart) - 0.5) * d
        ypart = csz * (- 0.5 + indy) + y + (rng.random(nPart) - 0.5) * d
    elif flagRand:
        # place particles randomly in the cell
        xpart = csz * (rng.random(nPart) - 0.5 + indx)
        ypart = csz * (rng.random(nPart) - 0.5 + indy)
    else:
        # place particles equaly distributed
        xpart = csz * (- 0.5 + indx) + x
        ypart = csz * (- 0.5 + indy) + y
    return xpart, ypart, mPart, nPart


def initializeMassEnt(dem, simTypeActual, entLine, reportAreaInfo):
    """ Initialize mass for entrainment

    Parameters
    ----------
    dem: dict
        dem dictionary

    Returns
    -------
    entrMassRaster : 2D numpy array
        raster of available mass for entrainment
    """
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    if simTypeActual in ['entres', 'ent']:
        # entrainmentArea = os.path.splitext(os.path.basename(entFiles))[0]
        entrainmentArea = entLine['Name']
        log.info('Entrainment area: %s' % (entrainmentArea))
        entrMassRaster = prepareArea(entLine, dem)
        reportAreaInfo['entrainment'] = 'Yes'
    else:
        entrMassRaster = np.zeros((nrows, ncols))
        reportAreaInfo['entrainment'] = 'No'

    return entrMassRaster, reportAreaInfo


def initializeResistance(cfg, dem, simTypeActual, resLine, reportAreaInfo):
    """ Initialize resistance matrix

    Parameters
    ----------
    dem: dict
        dem dictionary

    Returns
    -------
    cResRaster : 2D numpy array
        raster of resistance coefficients
    """
    d = cfg.getfloat('dRes')
    cw = cfg.getfloat('cw')
    sres = cfg.getfloat('sres')
    # read dem header
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    if simTypeActual in ['entres', 'res']:
        resistanceArea = resLine['Name']
        log.info('Resistance area: %s' % (resistanceArea))
        mask = prepareArea(resLine, dem)
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
    fields : dict
        fields dictionary at initial time step
    dem : dict
        dictionary with dem information

    Returns
    -------
    Particles : list
        list of particles dictionary
    Fields : list
        list of fields dictionary (for each time step saved)
    Tcpu : dict
        computation time dictionary
    infoDict : dict
        Dictionary of all simulations carried out
    """

    # Initialise cpu timing
    Tcpu = {}
    Tcpu['Force'] = 0.
    Tcpu['ForceVect'] = 0.
    Tcpu['ForceSPH'] = 0.
    Tcpu['Pos'] = 0.
    Tcpu['Neigh'] = 0.
    Tcpu['Field'] = 0.

    # Load configuration settings
    tEnd = cfg.getfloat('tEnd')
    dtSave = fU.splitIniValueToArrayInterval(cfg['tSteps'], cfg)
    # remove time step o as this is anyway saved
    if dtSave[0] == 0.0:
        dtSave = dtSave[1:]
    sphOption = cfg.getint('sphOption')
    log.info('using sphOption %s:' % sphOption)
    # desired output fields
    resTypesString = cfg['resType']
    resTypes = resTypesString.split('_')

    # Initialise Lists to save fields
    if 'particles' in resTypes:
        Particles = [copy.deepcopy(particles)]
    else:
        Particles = ''
    fieldAppend = {}
    for resType in resTypes:
        fieldAppend[resType] = fields[resType]
    Fields = [copy.deepcopy(fieldAppend)]
    #Fields = [copy.deepcopy(fields)]
    Tsave = [0]

    if featLF:
        log.info('Use LeapFrog time stepping')
    else:
        log.info('Use standard time stepping')
    # Initialize time and counters
    nSave = 1
    Tcpu['nSave'] = nSave
    nIter = 1
    nIter0 = 1
    iterate = True
    particles['iterate'] = iterate
    t = particles['t']
    timeM = []
    massEntrained = []
    massTotal = []
    # ++++++++++++++++if you want to use cfl time step+++++++++++++++++++
    # CALL TIME STEP:
    # to play around with the courant number
    if featCFL:
        dtStable = tD.getcflTimeStep(particles, dem, cfg)
    elif featCFLConstrain:
        dtStable = tD.getcfldTwithConstraints(particles, dem, cfg)

    # dt overwrites dt in .ini file, so comment this block if you dont want to use cfl
    # ++++++++++++++++++++++++++++++++++++++++++++++
    # get time step
    dt = cfg.getfloat('dt')
    t = t + dt

    # Start time step computation
    while t <= tEnd*(1.+1.e-13) and iterate:
        log.debug('Computing time step t = %f s', t)

        # Perform computations
        if featLF:
            particles, fields, Tcpu, dt = computeLeapFrogTimeStep(
                cfg, particles, fields, dt, dem, Tcpu)
        else:
            particles, fields, Tcpu = computeEulerTimeStep(
                cfg, particles, fields, dt, dem, Tcpu)

        Tcpu['nSave'] = nSave
        particles['t'] = t
        iterate = particles['iterate']

        # write mass balance info
        massEntrained.append(particles['massEntrained'])
        massTotal.append(particles['mTot'])
        timeM.append(t)
        if t >= dtSave[0]:
            Tsave.append(t)
            log.info('Saving results for time step t = %f s', t)
            log.info('MTot = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
            log.info(('cpu time Force = %s s' % (Tcpu['Force'] / nIter)))
            log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / nIter)))
            log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / nIter)))
            log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / nIter)))
            log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / nIter)))
            log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / nIter)))
            if 'particles' in resTypes:
                Particles.append(copy.deepcopy(particles))
            fieldAppend = {}
            for resType in resTypes:
                fieldAppend[resType] = fields[resType]
            Fields.append(copy.deepcopy(fieldAppend))
            if len(dtSave) > 1:
                dtSave = dtSave[1:]
            else:
                dtSave = [cfg.getfloat('tEnd')]

        # ++++++++++++++++if you want to use cfl time step+++++++++++++++++++
        # CALL TIME STEP:
        # to play around with the courant number
        if featCFL:
            dtStable = tD.getcflTimeStep(particles, dem, cfg)
        elif featCFLConstrain:
            dtStable = tD.getcfldTwithConstraints(particles, dem, cfg)

        # dt overwrites dt in .ini file, so comment this block if you dont want to use cfl
        # ++++++++++++++++++++++++++++++++++++++++++++++
        # get time step
        dt = cfg.getfloat('dt')
        t = t + dt
        nIter = nIter + 1
        nIter0 = nIter0 + 1

    Tcpu['nIter'] = nIter
    log.info('Ending computation at time t = %f s', t-dt)
    log.info('Saving results for time step t = %f s', t-dt)
    log.info('MTot = %f kg, %s particles' % (particles['mTot'], particles['Npart']))
    log.info(('cpu time Force = %s s' % (Tcpu['Force'] / nIter)))
    log.info(('cpu time ForceVect = %s s' % (Tcpu['ForceVect'] / nIter)))
    log.info(('cpu time ForceSPH = %s s' % (Tcpu['ForceSPH'] / nIter)))
    log.info(('cpu time Position = %s s' % (Tcpu['Pos'] / nIter)))
    log.info(('cpu time Neighbour = %s s' % (Tcpu['Neigh'] / nIter)))
    log.info(('cpu time Fields = %s s' % (Tcpu['Field'] / nIter)))
    Tsave.append(t)
    if 'particles' in resTypes:
        Particles.append(copy.deepcopy(particles))
    for resType in resTypes:
        fieldAppend[resType] = fields[resType]
    Fields.append(copy.deepcopy(fieldAppend))

    infoDict = {'massEntrained': massEntrained, 'timeStep': timeM, 'massTotal': massTotal, 'Tcpu': Tcpu,
                'final mass': massTotal[-1], 'initial mass': massTotal[0], 'entrained mass': np.sum(massEntrained),
                'entrained volume': (np.sum(massEntrained)/cfg.getfloat('rhoEnt'))}

    # determine if stop criterion is reached or end time
    stopCritNotReached = particles['iterate']
    avaTime = particles['t']
    stopCritPer = cfg.getfloat('stopCrit') *100.
    if stopCritNotReached:
        infoDict.update({'stopInfo': {'Stop criterion': 'end Time reached: %.2f' % avaTime, 'Avalanche run time [s]': '%.2f' % avaTime}})
    else:
        infoDict.update({'stopInfo': {'Stop criterion': '< %.2f percent of PKE' % stopCritPer, 'Avalanche run time [s]': '%.2f' % avaTime}})

    return Tsave, Particles, Fields, infoDict


def writeMBFile(infoDict, avaDir, logName):
    """ write mass balance info to file """

    t = infoDict['timeStep']
    massEntrained = infoDict['massEntrained']
    massTotal = infoDict['massTotal']

    # write mass balance info to log file
    massDir = os.path.join(avaDir, 'Outputs', 'com1DFAPy')
    fU.makeADir(massDir)
    with open(os.path.join(massDir, 'mass_%s.txt' % logName), 'w') as mFile:
        mFile.write('time, current, entrained\n')
        for m in range(len(t)):
            mFile.write('%.02f,    %.06f,    %.06f\n' %
                    (t[m], massTotal[m], massEntrained[m]))


def computeEulerTimeStep(cfg, particles, fields, dt, dem, Tcpu):
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
    particles, force, fields = DFAfunC.computeForceC(cfg, particles, fields, dem, dt)
    tcpuForce = time.time() - startTime
    Tcpu['Force'] = Tcpu['Force'] + tcpuForce

    # compute lateral force (SPH component of the calculation)
    startTime = time.time()
    if cfg.getint('sphOption') == 0:
        force['forceSPHX'] = np.zeros(np.shape(force['forceX']))
        force['forceSPHY'] = np.zeros(np.shape(force['forceY']))
        force['forceSPHZ'] = np.zeros(np.shape(force['forceZ']))
    else:
        particles, force = DFAfunC.computeForceSPHC(cfg, particles, force, dem, gradient=0)
    tcpuForceSPH = time.time() - startTime
    Tcpu['ForceSPH'] = Tcpu['ForceSPH'] + tcpuForceSPH

    hmin = cfg.getfloat('hmin')

    # update velocity and particle position
    startTime = time.time()
    # particles = updatePosition(cfg, particles, dem, force)
    particles = DFAfunC.updatePositionC(cfg, particles, dem, force)
    tcpuPos = time.time() - startTime
    Tcpu['Pos'] = Tcpu['Pos'] + tcpuPos

    # get particles location (neighbours for sph)
    startTime = time.time()
    particles = DFAfunC.getNeighboursC(particles, dem)

    tcpuNeigh = time.time() - startTime
    Tcpu['Neigh'] = Tcpu['Neigh'] + tcpuNeigh

    # update fields (compute grid values)
    startTime = time.time()
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
    csz = dem['header'].cellsize
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
    # force = computeForceVect(cfg, particles, dem, Ment, Cres, dtK5)
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
    particles = removeOutPart(cfg, particles, dem)

    # ++++++++++++++GET particles location (neighbours for sph)
    startTime = time.time()
    particles = DFAfunC.getNeighboursC(particles, dem)
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


def prepareArea(releaseLine, dem, relThList=''):
    """ convert shape file polygon to raster

    Parameters
    ----------
    releaseLine: dict
        line dictionary
    dem : dict
        dictionary with dem information
    relThList: list
        release thickness values for all release features

    Returns
    -------

    Raster : 2D numpy array
        raster
    """
    NameRel = releaseLine['Name']
    StartRel = releaseLine['Start']
    LengthRel = releaseLine['Length']
    Raster = np.zeros(np.shape(dem['rasterData']))
    RasterList = []

    for i in range(len(NameRel)):
        name = NameRel[i]
        start = StartRel[i]
        end = start + LengthRel[i]
        avapath = {}
        avapath['x'] = releaseLine['x'][int(start):int(end)]
        avapath['y'] = releaseLine['y'][int(start):int(end)]
        avapath['Name'] = name
        # if relTh is given - set relTh
        if relThList != '':
            log.info('Release feature %s, relTh= %.2f' % (name, relThList[i]))
            Raster = np.zeros(np.shape(dem['rasterData']))
            Raster = polygon2Raster(dem['header'], avapath, Raster, relTh=relThList[i])
            RasterList.append(Raster)
        else:
            Raster = polygon2Raster(dem['header'], avapath, Raster, relTh='')

    # if relTh provided by a field or function - create release Raster with ones
    if relThList != '':
        Raster = np.zeros(np.shape(dem['rasterData']))
        for rast in RasterList:
            ind1 = Raster > 0
            ind2 = rast > 0
            indMatch = np.logical_and(ind1, ind2)
            try:
                message = 'Release area features are overlaping - this is not allowed'
                assert indMatch.any() == False, message
                Raster = Raster + rast
            except AssertionError:
                log.error('%s' % message)
                raise

    return Raster


def polygon2Raster(demHeader, Line, Mask, relTh=''):
    """ convert line to raster

    Parameters
    ----------
    demHeader: dict
        dem header dictionary
    Line : dict
        line dictionary
    Mask : 2D numpy array
        raster to update
    Returns
    -------

    Mask : 2D numpy array
        updated raster
    """
    # adim and center dem and polygon
    ncols = demHeader.ncols
    nrows = demHeader.nrows
    xllc = demHeader.xllcenter
    yllc = demHeader.yllcenter
    csz = demHeader.cellsize
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
    x = np.linspace(0, ncols-1, ncols)
    y = np.linspace(0, nrows-1, nrows)
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    points = np.stack((X, Y), axis=-1)
    mask = path.contains_points(points)
    mask = mask.reshape((nrows, ncols)).astype(int)
    # mask = geoTrans.poly2maskSimple(xCoord, yCoord, ncols, nrows)
    Mask = Mask + mask
    # set release thickness unless release thickness field is provided, then return array with ones
    if relTh != '':
        log.info('REL set from dict, %.2f' % relTh)
        Mask = np.where(Mask > 0, relTh, 0.)
    else:
        Mask = np.where(Mask > 0, 1., 0.)

    if debugPlot:
        x = np.arange(ncols) * csz
        y = np.arange(nrows) * csz
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        ax.set_title('Release area')
        cmap = copy.copy(mpl.cm.get_cmap("Greys"))
        ref0, im = pU.NonUnifIm(ax, x, y, Mask, 'x [m]', 'y [m]',
                                extent=[x.min(), x.max(), y.min(), y.max()],
                                cmap=cmap, norm=None)
        ax.plot(xCoord0 * csz, yCoord0 * csz, 'r', label='release polyline')
        plt.legend()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(im, cax=cax)
        plt.show()

    return Mask


def plotPosition(fig, ax, particles, dem, data, Cmap, unit, plotPart=False, last=False):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    X = PointsX[0, :]
    Y = PointsY[:, 0]
    Z = dem['rasterData']
    x = particles['x'] + xllc
    y = particles['y'] + yllc
    xx = np.arange(ncols) * csz + xllc
    yy = np.arange(nrows) * csz + yllc
    try:
        # Get the images on an axis
        cb = ax.images[-1].colorbar
        if cb:
            cb.remove()
    except IndexError:
        pass

    ax.clear()
    ax.set_title('t=%.2f s' % particles['t'])
    cmap, _, lev, norm, ticks = makePalette.makeColorMap(
        Cmap, 0.0, np.nanmax(data), continuous=True)
    cmap.set_under(color='w')
    ref0, im = pU.NonUnifIm(ax, xx, yy, data, 'x [m]', 'y [m]',
                         extent=[x.min(), x.max(), y.min(), y.max()],
                         cmap=cmap, norm=norm)

    Cp1 = ax.contour(X, Y, Z, levels=10, colors='k')
    pU.addColorBar(im, ax, ticks, unit)
    if plotPart:
        # ax.plot(x, y, '.b', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 1], y[NPPC == 1], '.c', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 4], y[NPPC == 4], '.b', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 9], y[NPPC == 9], '.r', linestyle='None', markersize=1)
        # ax.plot(x[NPPC == 16], y[NPPC == 16], '.m', linestyle='None', markersize=1)
        # load variation colormap
        variable = particles['h']
        cmap, _, _, norm, ticks = makePalette.makeColorMap(
            pU.cmapDepth, np.amin(variable), np.amax(variable), continuous=True)
        # set range and steps of colormap
        cc = variable
        sc = ax.scatter(x, y, c=cc, cmap=cmap, marker='.')

        if last:
            pU.addColorBar(sc, ax, ticks, 'm', 'Flow Depth')

    plt.pause(0.1)
    return fig, ax


def plotContours(fig, ax, particles, dem, data, Cmap, unit, last=False):
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    X = PointsX[0, :]
    Y = PointsY[:, 0]
    Z = dem['rasterData']
    try:
        # Get the images on an axis
        cb = ax.images[-1].colorbar
        if cb:
            cb.remove()
    except IndexError:
        pass

    ax.clear()
    ax.set_title('t=%.2f s' % particles['t'])
    cmap, _, lev, norm, ticks = makePalette.makeColorMap(
        Cmap, 0.0, np.nanmax(data), continuous=True)
    cmap.set_under(color='w')

    CS = ax.contour(X, Y, data, levels=8, origin='lower', cmap=cmap,
                    linewidths=2)
    lev = CS.levels

    if last:
        # pU.addColorBar(im, ax, ticks, unit, 'Flow Depth')
        CB = fig.colorbar(CS)
        ax.clabel(CS, inline=1, fontsize=8)

    plt.pause(0.1)
    return fig, ax, cmap, lev


def removeOutPart(cfg, particles, dem):
    """ find and remove out of raster particles

    Parameters
    ----------
    cfg : configparser
        DFA parameters
    particles : dict
        particles dictionary
    dem : dict
        dem dictionary

    Returns
    -------
    particles : dict
        particles dictionary
    """
    dt = cfg.getfloat('dt')
    header = dem['header']
    nrows = header.nrows
    ncols = header.ncols
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    Bad = dem['Bad']

    x = particles['x']
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    x = x + ux*dt
    y = y + uy*dt

    # find coordinates in normalized ref (origin (0,0) and cellsize 1)
    Lx = (x - xllc) / csz
    Ly = (y - yllc) / csz
    mask = np.ones(len(x), dtype=bool)
    indOut = np.where(Lx <= 1.5)
    mask[indOut] = False
    indOut = np.where(Ly <= 1.5)
    mask[indOut] = False
    indOut = np.where(Lx >= ncols-1.5)
    mask[indOut] = False
    indOut = np.where(Ly >= nrows-1.5)
    mask[indOut] = False

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they exited the domain' % (nRemove))

    x = particles['x']
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    mask = np.ones(len(x), dtype=bool)
    y = particles['y']
    ux = particles['ux']
    uy = particles['uy']
    indXDEM = particles['indXDEM']
    indYDEM = particles['indYDEM']
    indOut = np.where(Bad[indYDEM, indXDEM], False, True)
    mask = np.logical_and(mask, indOut)
    indOut = np.where(Bad[indYDEM+np.sign(uy).astype('int'), indXDEM], False, True)
    mask = np.logical_and(mask, indOut)
    indOut = np.where(Bad[indYDEM, indXDEM+np.sign(ux).astype('int')], False, True)
    mask = np.logical_and(mask, indOut)
    indOut = np.where(Bad[indYDEM+np.sign(uy).astype('int'), indXDEM+np.sign(ux).astype('int')], False, True)
    mask = np.logical_and(mask, indOut)

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they exited the domain' % (nRemove))

    return particles


def removeSmallPart(hmin, particles, dem):
    """ find and remove too small particles

    Parameters
    ----------
    hmin : float
        minimum depth
    particles : dict
        particles dictionary
    dem : dict
        dem dictionary

    Returns
    -------
    particles : dict
        particles dictionary
    """
    h = particles['h']

    indOut = np.where(h < hmin)
    mask = np.ones(len(h), dtype=bool)
    mask[indOut] = False

    nRemove = len(mask)-np.sum(mask)
    if nRemove > 0:
        particles = removePart(particles, mask, nRemove)
        log.info('removed %s particles because they were too thin' % (nRemove))
        particles = DFAfunC.getNeighboursC(particles, dem)

    return particles


def removePart(particles, mask, nRemove):
    """ remove given particles

    Parameters
    ----------
    particles : dict
        particles dictionary
    mask : 1D numpy array
        particles to keep
    nRemove : int
        number of particles removed

    Returns
    -------
    particles : dict
        particles dictionary
    """
    particles['Npart'] = particles['Npart'] - nRemove
    particles['NPPC'] = particles['NPPC'][mask]
    particles['x'] = particles['x'][mask]
    particles['y'] = particles['y'][mask]
    particles['z'] = particles['z'][mask]
    particles['s'] = particles['s'][mask]
    particles['l'] = particles['l'][mask]
    particles['ux'] = particles['ux'][mask]
    particles['uy'] = particles['uy'][mask]
    particles['uz'] = particles['uz'][mask]
    particles['m'] = particles['m'][mask]
    particles['h'] = particles['h'][mask]
    particles['inCellDEM'] = particles['inCellDEM'][mask]
    particles['indXDEM'] = particles['indXDEM'][mask]
    particles['indYDEM'] = particles['indYDEM'][mask]
    particles['partInCell'] = particles['partInCell'][mask]

    return particles


def splitPart(cfg, particles, dem):
    massPerPart = cfg.getfloat('massPerPart')
    m = particles['m']
    nSplit = np.round(m/massPerPart)
    Ind = np.where(nSplit > 1)[0]
    if np.size(Ind) > 0:
        for ind in Ind:
            mNew = m[ind] / nSplit[ind]
            nAdd = (nSplit[ind]-1).astype('int')
            log.debug('Spliting particle %s in %s' % (ind, nAdd+1))
            particles['Npart'] = particles['Npart'] + nAdd
            particles['NPPC'] = np.append(particles['NPPC'], particles['NPPC'][ind]*np.ones((nAdd)))
            particles['x'] = np.append(particles['x'], particles['x'][ind]*np.ones((nAdd)))
            particles['y'] = np.append(particles['y'], particles['y'][ind]*np.ones((nAdd)))
            particles['z'] = np.append(particles['z'], particles['z'][ind]*np.ones((nAdd)))
            particles['s'] = np.append(particles['s'], particles['s'][ind]*np.ones((nAdd)))
            particles['l'] = np.append(particles['l'], particles['l'][ind]*np.ones((nAdd)))
            particles['ux'] = np.append(particles['ux'], particles['ux'][ind]*np.ones((nAdd)))
            particles['uy'] = np.append(particles['uy'], particles['uy'][ind]*np.ones((nAdd)))
            particles['uz'] = np.append(particles['uz'], particles['uz'][ind]*np.ones((nAdd)))
            particles['m'] = np.append(particles['m'], mNew*np.ones((nAdd)))
            particles['m'][ind] = mNew
            particles['h'] = np.append(particles['h'], particles['h'][ind]*np.ones((nAdd)))
            particles['inCellDEM'] = np.append(particles['inCellDEM'], particles['inCellDEM'][ind]*np.ones((nAdd)))
            particles['indXDEM'] = np.append(particles['indXDEM'], particles['indXDEM'][ind]*np.ones((nAdd)))
            particles['indYDEM'] = np.append(particles['indYDEM'], particles['indYDEM'][ind]*np.ones((nAdd)))
            particles['partInCell'] = np.append(particles['partInCell'], particles['partInCell'][ind]*np.ones((nAdd)))

    return particles


def savePartToPickle(dictList, outDir):
    """ Save each dictionary from a list to a pickle in outDir; works also for one dictionary instead of list

        Parameters
        ---------
        dictList: list or dict
            list of dictionaries or single dictionary
        outDir: str
            path to output directory
    """

    if isinstance(dictList, list):
        for dict in dictList:
            pickle.dump(dict, open(os.path.join(outDir, "particles%09.4f.p" % (dict['t'])), "wb"))
    else:
        pickle.dump(dictList, open(os.path.join(outDir, "particles%09.4f.p" % (dictList['t'])), "wb"))


def readPartFromPickle(inDir, flagAvaDir=False):
    """ Read pickles within a directory and return List of dicionaries read from pickle

        Parameters
        -----------
        inDir: str
            path to input directory
        flagAvaDir: bool
            if True inDir corresponds to an avalanche directory and pickles are
            read from avaDir/Outputs/com1DFAPy/particles
    """

    if flagAvaDir:
        inDir = os.path.join(inDir, 'Outputs', 'com1DFAPy', 'particles')

    # search for all pickles within directory
    PartDicts = sorted(glob.glob(os.path.join(inDir, '*.p')))

    # initialise list of particle dictionaries
    Particles = []
    TimeStepInfo = []
    for particles in PartDicts:
        particles = pickle.load(open(particles, "rb"))
        Particles.append(particles)
        TimeStepInfo.append(particles['t'])

    return Particles, TimeStepInfo


def exportFields(cfgGen, Tsave, Fields, relFile, demOri, outDir, logName):
    """ export result fields to Outputs directory according to result parameters and time step
        that can be specified in the configuration file

        Parameters
        -----------
        cfgGen: dict
            configurations
        Tsave: list
            list of time step that corresponds to each dict in Fields
        Fields: list
            list of Fields for each dtSave
        relFile: str
            path to release area shapefile
        outDir: str
            outputs Directory


        Returns
        --------
        exported peak fields are saved in Outputs/com1DFAPy/peakFiles

    """

    resTypesString = cfgGen['resType']
    resTypes = resTypesString.split('_')
    numberTimes = len(Tsave)-1
    countTime = 0
    for timeStep in Tsave:
        for resType in resTypes:
            resField = Fields[countTime][resType]
            if resType == 'ppr':
                resField = resField * 0.001
            dataName = logName + '_' + resType + '_'  + 't%.2f' % (Tsave[countTime]) +'.asc'
            # create directory
            outDirPeak = os.path.join(outDir, 'peakFiles', 'timeSteps')
            fU.makeADir(outDirPeak)
            outFile = os.path.join(outDirPeak, dataName)
            IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
            if countTime == numberTimes:
                log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f - FINAL time step ' % (resType,Tsave[countTime]))
                dataName = logName + '_' + resType + '_' +'.asc'
                # create directory
                outDirPeakAll = os.path.join(outDir, 'peakFiles')
                fU.makeADir(outDirPeakAll)
                outFile = os.path.join(outDirPeakAll, dataName)
                IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
            else:
                log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f ' % (resType,Tsave[countTime]))
        countTime = countTime + 1


def exportFieldsOld(cfgGen, Tsave, Fields, relFile, demOri, outDir, logName):
    """ export result fields to Outputs directory according to result parameters and time step
        that can be specified in the configuration file

        Parameters
        -----------
        cfgGen: dict
            configurations
        Tsave: list
            list of time step that corresponds to each dict in Fields
        Fields: list
            list of Fields for each dtSave
        relFile: str
            path to release area shapefile
        outDir: str
            outputs Directory


        Returns
        --------
        exported peak fields are saved in Outputs/com1DFAPy/peakFiles

    """

    resTypesString = cfgGen['resType']
    resTypes = resTypesString.split('_')
    tSteps = fU.getTimeIndex(cfgGen, Fields)
    if -1 not in tSteps:
        tSteps.append(-1)
        log.info('-1 added to tStep')
    for tStep in tSteps:
        finalFields = Fields[tStep]
        for resType in resTypes:
            resField = finalFields[resType]
            if resType == 'ppr':
                resField = resField * 0.001
            dataName = logName + '_' + resType + '_'  + 't%.2f' % (Tsave[tStep]) +'.asc'
            # create directory
            outDirPeak = os.path.join(outDir, 'peakFiles', 'timeSteps')
            fU.makeADir(outDirPeak)
            outFile = os.path.join(outDirPeak, dataName)
            IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
            if tStep == -1:
                log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f - FINAL time step ' % (resType,Tsave[tStep]))
                dataName = logName + '_' + resType + '_' +'.asc'
                # create directory
                outDirPeakAll = os.path.join(outDir, 'peakFiles')
                fU.makeADir(outDirPeakAll)
                outFile = os.path.join(outDirPeakAll, dataName)
                IOf.writeResultToAsc(demOri['header'], resField, outFile, flip=True)
            else:
                log.info('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f ' % (resType,Tsave[tStep]))


def analysisPlots(Particles, Fields, cfg, demOri, dem, outDir):
    """ create analysis plots during simulation run """

    cfgGen = cfg['GENERAL']
    partRef = Particles[0]
    Z0 = partRef['z'][0]
    rho = cfgGen.getfloat('rho')
    gravAcc = cfgGen.getfloat('gravAcc')
    mu = cfgGen.getfloat('mu')
    repeat = True
    while repeat:
        fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        T = np.array([0])
        Z = np.array([0])
        U = np.array([0])
        S = np.array([0])
        for part, field in zip(Particles, Fields):
            T = np.append(T, part['t'])
            S = np.append(S, part['s'][0])
            Z = np.append(Z, part['z'][0])
            U = np.append(U, DFAtls.norm(part['ux'][0], part['uy'][0], part['uz'][0]))
            fig, ax = plotPosition(
                fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True)
            fig.savefig(os.path.join(outDir, 'particlest%f.%s' % (part['t'], pU.outputFormat)))

        fig, ax = plotPosition(
                fig, ax, part, demOri, dem['Nz'], pU.cmapDEM2, '', plotPart=True, last=True)
        fig.savefig(os.path.join(outDir, 'particlesFinal.%s' % (pU.outputFormat)))
        value = input("[y] to repeat:\n")
        if value != 'y':
            repeat = False

    fieldEnd = Fields[-1]
    partEnd = Particles[-1]
    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig3, ax3 = plt.subplots(figsize=(pU.figW, pU.figH))
    fig1, ax1 = plotPosition(
        fig1, ax1, partEnd, demOri, fieldEnd['FD'], pU.cmapPres, 'm', plotPart=False)
    fig2, ax2 = plotPosition(
        fig2, ax2, partEnd, demOri, fieldEnd['FV'], pU.cmapPres, 'm/s', plotPart=False)
    fig3, ax3 = plotPosition(
        fig3, ax3, partEnd, demOri, fieldEnd['P']/1000, pU.cmapPres, 'kPa', plotPart=False)
    plt.show()
