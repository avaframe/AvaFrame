"""
    Fetch input data for avalanche simulations
"""

# Load modules
import os
import glob
import pathlib
import logging
import numpy as np
import pandas as pd
import shapely as shp

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.com1DFA.deriveParameterSet as dP
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.com1DFA import com1DFA


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def readDEM(avaDir):
    """ read the ascii DEM file from a provided avalanche directory

    Parameters
    ----------
    avaDir : str
        path to avalanche directory

    Returns
    -------
    dem : dict
        dict with header and raster data
    """

    # get dem file name
    demSource = getDEMPath(avaDir)

    log.debug('Read DEM: %s' % demSource)

    dem = IOf.readRaster(demSource)

    return(dem)


def getDEMPath(avaDir):
    """ get the DEM file path from a provided avalanche directory

    Parameters
    ----------
    avaDir : str
        path to avalanche directory

    Returns
    -------
    demFile : str (first element of list)
        full path to DEM .asc file
    """

    # if more than one .asc file found throw error
    inputDir = pathlib.Path(avaDir, 'Inputs')
    demFile = list(inputDir.glob('*.asc'))
    if len(demFile) > 1:
        message = 'There should be exactly one topography .asc file in %s/Inputs/' % (avaDir)
        log.error(message)
        raise AssertionError(message)

    # if is no .asc file found - throw error
    filesFound = list(inputDir.glob('*.*'))
    if len(demFile) == 0 and len(filesFound):
        for fileF in filesFound:
            message = 'DEM file format not correct in %s/Inputs/ - only .asc is allowed but %s is provided' % (avaDir, fileF.name)
            log.error(message)
            raise AssertionError(message)
    elif len(demFile) == 0:
        message = 'No topography .asc file in %s/Inputs/' % (avaDir)
        log.error(message)
        raise FileNotFoundError(message)

    return demFile[0]


def getDEMFromConfig(avaDir, fileName=''):
    """ get dem file path in avaDir/Inputs or if fileName in avaDir/Inputs/fileName

        Parameters
        -----------
        avaDir: str or pathlib path
            path to avalancheDir
        fileName: pathlib path
            path to dem file with filename in avaDir/Inputs

        Returns
        --------
        demFile: pathlib path
            path to dem file
    """

    if fileName == '':
        demFilePath = getDEMPath(avaDir)
    else:
        inputDir = pathlib.Path(avaDir, 'Inputs')
        demFile = inputDir / fileName
        if demFile.is_file() is False:
            message = 'Dem file: %s does not exist' % (str(demFile))
            log.error(message)
            raise FileNotFoundError(message)

    return demFile


def getInputData(avaDir, cfg):
    """ Fetch input datasets required for simulation, duplicated function because
        simulation type set differently in com1DFAOrig compared to com1DFA:
        TODO: remove duplicate once it is not required anymore

    Parameters
    ----------
    avaDir : str
        path to avalanche directory
    cfg : dict
        configuration read from com1DFA simulation ini file

    Returns
    -------
    demFile[0] : str (first element of list)
        list of full path to DEM .asc file
    relFiles : list
        list of full path to release area scenario .shp files
    entFile : str
        full path to entrainment area .shp file
    resFile : str
        full path to resistance area .shp file
    wallFile: str
        full path to wall line .shp file
    entResInfo : flag dict
        flag if Yes entrainment and/or resistance areas found and used for simulation
    """

    # Set directories for inputs, outputs and current work
    inputDir = os.path.join(avaDir, 'Inputs')

    # Set flag if there is an entrainment or resistance area
    entResInfo= {'flagEnt': 'No', 'flagRes': 'No'}

    # Initialise release areas, default is to look for shapefiles
    if cfg['releaseScenario'] != '':
        releaseDir = 'REL'
        relFiles = []
        releaseFiles = cfg['releaseScenario'].split('|')
        for rel in releaseFiles:
            if '.shp' in rel:
                relf = os.path.join(inputDir, releaseDir, rel)
            else:
                relf = os.path.join(inputDir, releaseDir, '%s.shp' % (rel))
            if not os.path.isfile(relf):
                message = 'No release scenario called: %s' % (relf)
                log.error(message)
                raise FileNotFoundError(message)
            relFiles.append(relf)
        log.debug('Release area file is specified to be: %s' % relFiles)
    else:
        releaseDir = 'REL'
        relFiles = sorted(glob.glob(inputDir+os.sep + releaseDir+os.sep + '*.shp'))
    log.info('Release area files are: %s' % relFiles)

    # Initialise resistance areas
    resFile, entResInfo['flagRes'] = getAndCheckInputFiles(inputDir, 'RES', 'Resistance', fileExt='shp')
    if resFile is None:
        resFile = ''
    # Initialise entrainment areas
    entFile, entResInfo['flagEnt'] = getAndCheckInputFiles(inputDir, 'ENT', 'Entrainment', fileExt='shp')
    if entFile is None:
        entFile = ''
    # Initialise dam line
    wallFile, entResInfo['flagWall'] = getAndCheckInputFiles(inputDir, 'DAM', 'Dam', fileExt='shp')
    # Initialise DEM
    demFile = getDEMPath(avaDir)

    return demFile, relFiles, entFile, resFile, wallFile, entResInfo


def getInputDataCom1DFA(avaDir):
    """ Fetch input datasets required for simulation, duplicated function because
    now fetch all available files simulation type set differently in com1DFA compared
    to com1DFAOrig: TODO: remove duplicate once it is not required anymore

    Parameters
    ----------
    avaDir : str or pathlib object
        path to avalanche directory
    cfg : dict
        configuration read from com1DFA simulation ini file

    Returns
    -------
    inputSimFiles: dict
        dictionary with all the input files

        - demFile : str (first element of list), list of full path to DEM .asc file
        - relFiles : list, list of full path to release area scenario .shp files
        - secondaryReleaseFile : str, full path to secondary release area .shp file
        - entFile : str, full path to entrainment area .shp file
        - resFile : str, full path to resistance area .shp file
        - entResInfo : flag dict
        flag if Yes entrainment and/or resistance areas found and used for simulation
        flag True if a Secondary Release file found and activated

    """

    # Set directories for inputs, outputs and current work
    inputDir = pathlib.Path(avaDir, 'Inputs')

    # Set flag if there is an entrainment or resistance area
    entResInfo = {}

    releaseDir = 'REL'
    releaseDir = inputDir / 'REL'
    relFiles = sorted(list(releaseDir.glob('*.shp')))
    log.info('Release area files are: %s' % [str(relFilestr) for relFilestr in relFiles])

    # check if relThFile is available
    relThFile, entResInfo['releaseThicknessFile'] = getAndCheckInputFiles(inputDir, 'RELTH',
                                                                          'release thickness data', fileExt='asc')

    # Initialise secondary release areas
    secondaryReleaseFile, entResInfo['flagSecondaryRelease'] = getAndCheckInputFiles(inputDir, 'SECREL',
                                                                                     'Secondary release', fileExt='shp')

    # Initialise resistance areas
    resFile, entResInfo['flagRes'] = getAndCheckInputFiles(inputDir, 'RES', 'Resistance', fileExt='shp')

    # Initialise entrainment areas
    entFile, entResInfo['flagEnt'] = getAndCheckInputFiles(inputDir, 'ENT', 'Entrainment', fileExt='shp')

    # Initialise dam line
    damFile, entResInfo['dam'] = getAndCheckInputFiles(inputDir, 'DAM', 'Dam', fileExt='shp')

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    # return DEM, first item of release, entrainment and resistance areas
    inputSimFiles = {'demFile': demFile, 'relFiles': relFiles, 'secondaryReleaseFile': secondaryReleaseFile,
                     'entFile': entFile, 'resFile': resFile, 'damFile': damFile, 'entResInfo': entResInfo,
                     'relThFile': relThFile}

    return inputSimFiles


def getAndCheckInputFiles(inputDir, folder, inputType, fileExt='shp'):
    """Fetch fileExt files and check if they exist and if it is not more than one

    Raises error if there is more than one fileExt file.

    Parameters
    ----------
    inputDir : pathlib object or str
        path to avalanche input directory
    folder : str
        subfolder name where the shape file should be located (SECREL, ENT or RES)
    inputType : str
        type of input (used for the logging messages). Secondary release or Entrainment or Resistance
    fileExt: str
        file extension e.g. shp, asc - optional default is shp

    Returns
    -------
    OutputFile: str
        path to file checked
    available: str
        Yes or No depending of if there is a shape file available (if No, OutputFile is None)
    """
    available = 'No'
    # Initialise secondary release areas
    dir = pathlib.Path(inputDir, folder)
    OutputFile = list(dir.glob('*.%s' % fileExt))
    if len(OutputFile) < 1:
        OutputFile = None
    elif len(OutputFile) > 1:
        message = 'More than one %s .%s file in %s/%s/ not allowed' % (inputType, fileExt, inputDir, folder)
        log.error(message)
        raise AssertionError(message)
    else:
        available = 'Yes'
        OutputFile = OutputFile[0]

    return OutputFile, available


def getThicknessInputSimFiles(inputSimFiles, avaDir):
    """ add thickness of shapefiles to dictionary

        Parameters
        -----------
        inputSimFiles: dict
            dictionary with info on release and entrainment file paths
        avaDir: str or pathlib path
            path to avalanche directory

        Returns
        --------
        inputSimFiles: dict
            updated dictionary with thickness info read from shapefile attributes
            now includes one separate dictionary for each release, entrainment or secondary release
            scenario with a thickness and id value for each feature (given as list)
    """

    # create pathlib Path
    avaDir = pathlib.Path(avaDir)

    # check if thickness info is required from entrainment and secondary release according to simType
    thTypeList = ['entFile', 'secondaryReleaseFile']

    # fetch thickness attribute of entrainment area and secondary release
    for thType in ['entFile', 'secondaryReleaseFile']:
        if inputSimFiles[thType] != None:
            thicknessList, idList, ci95List = shpConv.readThickness(inputSimFiles[thType])
            inputSimFiles[inputSimFiles[thType].stem] = {'thickness': thicknessList, 'id': idList,
                'ci95': ci95List}

    # initialize release scenario list
    releaseScenarioList = []

    # fetch thickness attribute of release areas and add info to input dict
    for releaseA in inputSimFiles['relFiles']:
        # fetch thickness and id info from input data
        thicknessList, idList, ci95List = shpConv.readThickness(releaseA)
        inputSimFiles[releaseA.stem] = {'thickness': thicknessList, 'id': idList, 'ci95': ci95List}
        # append release scenario name to list
        releaseScenarioList.append(releaseA.stem)

    # append release scenario names
    inputSimFiles['releaseScenarioList'] = releaseScenarioList

    return inputSimFiles


def updateThicknessCfg(inputSimFiles, avaDir, modName, cfgInitial):
    """ add available release scenarios to ini file and
        set thickness values in ini files

        Parameters
        -----------
        inputSimFiles: dict
            dictionary with info on release and entrainment file paths
        avaDir: str or pathlib path
            path to avalanche directory
        modName : computational module
            computational module
        cfgInitial: configParser object
            configParser object with the current (and possibly overridden) configuration

        Returns
        --------
        inputSimFiles: dict
            updated dictionary with thickness info read from shapefile attributes
            now includes one separate dictionary for each release, entrainment or secondary release
            scenario with a thickness and id value for each feature (given as list)
        cfgInitial: configparser object
            updated config object with release scenario, thickness info, etc.

    """

    # create pathlib Path
    avaDir = pathlib.Path(avaDir)

    # get name of module as string
    modNameString = str(pathlib.Path(modName.__file__).stem)

    # check if thickness info is required from entrainment and secondary release according to simType
    simTypeList = cfgInitial['GENERAL']['simTypeList'].split('|')
    thTypeList = []
    if any(simType in ['ent', 'entres', 'available'] for simType in simTypeList):
        thTypeList.append('entFile')
    if cfgInitial['GENERAL'].getboolean('secRelArea'):
        thTypeList.append('secondaryReleaseFile')

    # initialize release scenario list
    releaseScenarioIni = cfgInitial['INPUT']['releaseScenario']
    if releaseScenarioIni == '':
        releaseScenarioList = inputSimFiles['releaseScenarioList']
    else:
        releaseScenarioList = cfgInitial['INPUT']['releaseScenario'].split('|')

    # add input data info to cfg object
    # fetch thickness attribute of release areas and add info to input dict
    for releaseA in releaseScenarioList:
        # update configuration with thickness value to be used for simulations
        cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, releaseA, 'relTh')
        if cfgInitial['GENERAL'].getboolean('relThFromFile'):
            if inputSimFiles['relThFile'] == None:
                message = 'relThFromFile set to True but no relTh file found'
                log.error(message)
                raise FileNotFoundError(message)
            else:
                cfgInitial['INPUT']['relThFile'] = str(pathlib.Path('RELTH', inputSimFiles['relThFile'].name))

    # add entrainment and secondary release thickness in input data info and in cfg object
    if inputSimFiles['entFile'] != None and 'entFile' in thTypeList:
        cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, inputSimFiles['entFile'].stem, 'entTh')
        cfgInitial['INPUT']['entrainmentScenario'] = inputSimFiles['entFile'].stem
    if inputSimFiles['secondaryReleaseFile'] != None and 'secondaryReleaseFile' in thTypeList:
        cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles,
            inputSimFiles['secondaryReleaseFile'].stem, 'secondaryRelTh')
        cfgInitial['INPUT']['secondaryReleaseScenario'] = inputSimFiles['secondaryReleaseFile'].stem

    # create cfg string from release scenario list and add to cfg object
    releaseScenarioName = cfgUtils.convertToCfgList(releaseScenarioList)
    if cfgInitial['INPUT']['releaseScenario'] == '':
        cfgInitial['INPUT']['releaseScenario'] = releaseScenarioName
    else:
        for relIniFileName in cfgInitial['INPUT']['releaseScenario'].split('|'):
            if relIniFileName not in releaseScenarioList:
                message = ('Chosen release scenario: %s not available' % relIniFileName)
                log.error(message)
                raise FileNotFoundError(message)
        else:
            log.info('Chosen release scenarios: %s' % cfgInitial['INPUT']['releaseScenario'])

    return cfgInitial


def initializeDEM(avaDir, demPath=''):
    """ check for dem and load to dict

        Parameters
        -----------
        avaDir: str or pathlib path
            path to avalanche directory
        demPath: str or pathlib Path
            path to dem relative to Inputs - optional if not provided read DEM from Inputs

        Returns
        --------
        demOri: dict
            dem dictionary with header and data
    """

    if demPath == '':
        dem = readDEM(avaDir)
    else:
        # build full path and load data to dict
        demFile = pathlib.Path(avaDir, 'Inputs', demPath)
        dem = IOf.readRaster(demFile, noDataToNan=True)

    return dem


def selectReleaseFile(inputSimFiles, releaseScenario):
    """ select release scenario

        Parameters
        -----------
        inputSimFiles: dict
            dictionary with info on input data
        releaseScenario: str
            name of release scenario


        Returns
        -------
        inputSimFiles: dict
            dictionary with info on input data updated with releaseScenario
    """


    # fetch release file path for scenario
    relFiles = inputSimFiles['relFiles']
    for relF in relFiles:
        if relF.stem == releaseScenario:
            releaseScenarioPath = relF

    inputSimFiles['releaseScenario'] =  releaseScenarioPath

    return inputSimFiles


def fetchReleaseFile(inputSimFiles, releaseScenario, cfgSim, releaseList):
    """ select release scenario, update configuration to only include thickness info
        of current scenario and return file path

        Parameters
        -----------
        inputSimFiles: dict
            dictionary with info on input data
        releaseScenario: str
            name of release scenario
        cfgSim: conigparser object
            configuration of simulation
        releaseList: list
            list of available release scenarios

        Returns
        -------
        releaseScenarioPath: pathlib path
            file path to release scenario shp file
        cfgSim: configparser object
            updated cfg object, removed thickness info from not other release scenarios than used
            one and rename thickness values of chosen scenario to relThThickness, relThId, ...
    """

    # fetch release files paths
    relFiles = inputSimFiles['relFiles']

    foundScenario = False
    for relF in relFiles:
        if relF.stem == releaseScenario:
            releaseScenarioPath = relF
            foundScenario = True

    if foundScenario is False:
        message = 'Release area scenario %s not found - check input data' % (releaseScenario)
        log.error(message)
        raise FileNotFoundError(message)

    # update config entry for release scenario, thickness and id
    cfgSim['INPUT']['releaseScenario'] = str(releaseScenario)
    if cfgSim['GENERAL']['relThFromShp'] == 'True':
        for scenario in releaseList:
            if scenario == releaseScenario:
                cfgSim['INPUT']['relThId'] = cfgSim['INPUT'][scenario + '_' + 'relThId']
                cfgSim['INPUT']['relThThickness'] = cfgSim['INPUT'][scenario + '_' + 'relThThickness']
                cfgSim['INPUT']['relThCi95'] = cfgSim['INPUT'][scenario + '_' + 'relThCi95']
            # remove thickness, id and ci95 values specified by releaseScenario
            cfgSim['INPUT'].pop(scenario + '_' + 'relThId')
            cfgSim['INPUT'].pop(scenario + '_' + 'relThThickness')
            cfgSim['INPUT'].pop(scenario + '_' + 'relThCi95')

    return releaseScenarioPath, cfgSim


def createReleaseStats(avaDir, cfg):
    """ create a csv file with info on release shp file on:
        max, mean and min elevation, slope and projected and real area

        Parameters
        -------------

        Returns
        --------
        fPath: pathlib Path
            file path
    """

    # fetch input data (dem and release area shp file)
    inputData = getInputDataCom1DFA(avaDir)

    # get line from release area polygon
    relFiles = inputData['relFiles']

    # initialize dem and get real area of dem
    dem = IOf.readRaster(inputData['demFile'], noDataToNan=True)
    dem['originalHeader'] = dem['header'].copy()
    methodMeshNormal = cfg.getfloat('GENERAL', 'methodMeshNormal')
    # get normal vector of the grid mesh
    dem = DFAtls.getNormalMesh(dem, methodMeshNormal)
    dem = DFAtls.getAreaMesh(dem, methodMeshNormal)
    header = dem['header']
    csz = header['cellsize']

    # loop over all relFiles and compute information saved to dictionary of dataframes
    relInfo = {}
    relPath = pathlib.Path(avaDir, 'Outputs', 'com1DFA', 'releaseInfoFiles')
    fU.makeADir(relPath)
    relDFDict = {}
    # loop over all relFiles
    for relFile in relFiles:
        # for relF in relFiles:
        releaseLine = {}
        releaseLine = shpConv.readLine(relFile, 'release1', dem)
        releaseLine['file'] = relFile
        releaseLine['type'] = 'Release'
        releaseLine['thicknessSource'] = ['artificial'] * len(releaseLine['id'])
        releaseLine['thickness'] = ['1.'] * len(releaseLine['id'])
        # convert release line to a raster with 1 set inside of the release line
        # and compute projected and actual areas
        areaActualList, areaProjectedList, releaseLine = computeAreasFromRasterAndLine(releaseLine, dem)

        # compute max, min elevation, mean slope and save all to dict
        relInfo = computeRelStats(releaseLine, dem)

        # create dict and dataFrame
        relD = {'release feature': relInfo['featureNames'], 'MaxZ [m]': relInfo['zMax'],
            'MinZ [m]': relInfo['zMin'], 'slope [deg]': relInfo['meanSlope'],
            'projected area [ha]':[areaP/10000. for areaP in areaProjectedList],
            'actual area [ha]': [area/10000. for area in areaActualList]}
        relDF = pd.DataFrame(data=relD, index=np.arange(len(relInfo['featureNames'])))
        relInfo = relPath / ('%s.csv' % relFile.stem)
        relDF.to_csv(relInfo, index=False)
        log.info('Written relInfo file for %s to %s' % (relFile.stem, str(relPath)))
        relDFDict[relFile.stem] = relDF

    return relDFDict


def computeAreasFromRasterAndLine(line, dem):
    """ compute the area covered by a polygon by creating a raster from polygon
        projected area and actual area using a dem info

        Parameters
        -----------
        line: dict
            dictionary with info on line
            x, y coordinates start, end of each line feature
        dem: dict
            dictionary with dem data, header and areaRaster

        Returns
        --------
        areaActual: float
            actual area taking slope into account by using dem area
        areaProjected: float
            projected area in xy plane
    """

    line = com1DFA.prepareArea(line, dem, 0.01, combine=False, checkOverlap=False)

    csz = dem['header']['cellsize']
    # create dict for raster data for each feature
    areaProjectedList = []
    areaActualList = []
    for index, lineRaster in enumerate(line['rasterData']):
        lineRasterOnes = np.where(lineRaster > 0, 1., 0.)
        areaActualList.append(np.nansum(lineRasterOnes*dem['areaRaster']))
        areaProjectedList.append(np.sum(csz*csz*lineRasterOnes))

    return areaActualList, areaProjectedList, line


def computeRelStats(line, dem):
    """ compute stats of a polygon and a dem
        actual area (taking slope into account), projected area,
        max, mean, min elevation, mean slope

        Parameters
        -----------
        line: dict
            dictionary with info on line (x, y, Start, Length, rasterData, ...)
        dem: dict
            dictionary with info on dem (header, rasterData, normals, areaRaster)

        Returns
        --------
        lineDict: dict
            dictionary with stats info:
            featureNames, zMax, zMin, Slope, Area, AreaP
    """

    # compute projected areas using shapely
    projectedAreas = computeAreasFromLines(line)

    # create dict for raster data for each feature
    lineDict = {}
    for index, relRaster in enumerate(line['rasterData']):
        zArray = np.where(relRaster > 0, dem['rasterData'], np.nan)

        # compute slope of release area
        _, _, NzNormed = DFAtls.normalize(dem['Nx'], dem['Ny'], dem['Nz'])
        nzArray = np.where(relRaster > 0, NzNormed, np.nan)
        slopeArray = np.rad2deg(np.arccos(nzArray))
        lineDict.setdefault('meanSlope', []).append(np.nanmean(slopeArray))

        # compute elevation stats
        lineDict.setdefault('zMax', []).append(np.nanmax(zArray))
        lineDict.setdefault('zMean', []).append(np.nanmean(zArray))
        lineDict.setdefault('zMin', []).append(np.nanmin(zArray))

        # create a dataframe
        lineDict.setdefault('featureNames', []).append(line['Name'][index])

        # print info
        log.debug('++++ Release feature %s ++++++++++' % line['Name'][index])
        log.debug('maximum elevation: %.2f' % np.nanmax(zArray))
        log.debug('mean elevation: %.2f' % np.nanmean(zArray))
        log.debug('minimum elevation: %.2f' % np.nanmin(zArray))
        log.debug('mean slope: %.2f' % np.nanmean(slopeArray))

    return lineDict


def computeAreasFromLines(line):
    """ compute the area of a polygon in xy using shapely

        Parameters
        -----------
        line: dict
            dictionary with info on line
            x, y coordinates start, end of each line feature

        Returns
        --------
        projectedAreas: list
            list of projected area for each polygon in line dict
    """

    projectedAreas = []
    name = line['Name']
    start = line['Start']
    Length = line['Length']
    # fetch individual polygons
    for i in range(len(name)):
        end = start[i] + Length[i]
        x = line['x'][int(start[i]):int(end)]
        y = line['y'][int(start[i]):int(end)]

        # create shapely polygon
        for m in range(len(x)):
            avaPoly = shp.Polygon(list(zip(x, y)))
        projectedAreas.append(shp.area(avaPoly))

    return projectedAreas


def getInputPaths(avaDir):
    """ Fetch paths to dem and first release area shp file found

    Parameters
    ----------
    avaDir : str or pathlib object
        path to avalanche directory


    Returns
    -------
    demFile : pathlib path
        full path to DEM .asc file
    relFiles : list
        list of full paths to release area scenario .shp files found in avaDir/Inputs/REL
    relFieldFiles : list
        list of full paths to release area thickness .asc files found in avaDir/Inputs/RELTH

    """

    # Set directories for inputs, outputs and current work
    inputDir = pathlib.Path(avaDir, 'Inputs')

    # fetch release area shp files
    releaseShpDir = inputDir / 'REL'
    relFiles = sorted(list(releaseShpDir.glob('*.shp')))
    log.info('Release area files are: %s' % [str(relFilestr) for relFilestr in relFiles])

    # fetch release thickness fields
    releaseFieldDir = inputDir / 'RELTH'
    relFieldFiles = sorted(list(releaseFieldDir.glob('*.asc')))
    if len(relFieldFiles) > 0:
        log.info('Release area files are: %s' % [str(relFFilestr) for relFFilestr in relFieldFiles])
    else:
        relFieldFiles = None

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    return demFile, relFiles, relFieldFiles