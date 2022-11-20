"""
    Fetch input data for avalanche simulations
"""

# Load modules
import os
import glob
import pathlib
import logging
import shutil

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.com1DFA.deriveParameterSet as dP
import avaframe.in2Trans.shp_to_ascii as shp_to_ascii
import avaframe.in2Trans.cellsize_change as cellsize_change


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
    if resFile == None:
        resFile = ''
    # Initialise entrainment areas
    entFile, entResInfo['flagEnt'] = getAndCheckInputFiles(inputDir, 'ENT', 'Entrainment', fileExt='shp')
    if entFile == None:
        entFile = ''
    # Initialise DEM
    demFile = getDEMPath(avaDir)

    return demFile, relFiles, entFile, resFile, entResInfo


def getInputDataCom1DFA(avaDir, cfg):
    """ Fetch input datasets required for simulation, duplicated function because
        simulation type set differently in com1DFA compared to com1DFAOrig:
        TODO: remove duplicate once it is not required anymore

    Parameters
    ----------
    avaDir : str or pathlib object
        path to avalanche directory
    cfg : dict
        configuration read from com1DFA simulation ini file

    Returns
    -------
    inputSimFiles: dict
        dictionary with all the input files:
        demFile : str (first element of list)
            list of full path to DEM .asc file
        relFiles : list
            list of full path to release area scenario .shp files
        secondaryReleaseFile : str
            full path to secondary release area .shp file
        entFile : str
            full path to entrainment area .shp file
        resFile : str
            full path to resistance area .shp file
        entResInfo : flag dict
            flag if Yes entrainment and/or resistance areas found and used for simulation
            flag True if a Secondary Release file found and activated
    """

    # Set directories for inputs, outputs and current work
    inputDir = pathlib.Path(avaDir, 'Inputs')

    # Set flag if there is an entrainment or resistance area
    entResInfo = {}

    # Initialise release areas, default is to look for shapefiles
    if cfg['INPUT']['releaseScenario'] != '':
        releaseDir = 'REL'
        relFiles = []
        releaseFiles = cfg['INPUT']['releaseScenario'].split('|')
        for rel in releaseFiles:
            if '.shp' in rel:
                relf = inputDir / releaseDir / rel
            else:
                relf = inputDir / releaseDir / ('%s.shp' % (rel))
            if not relf.is_file():
                message = 'No release scenario called: %s' % (relf)
                log.error(message)
                raise FileNotFoundError(message)
            relFiles.append(relf)
        log.debug('Release area file is specified to be: %s' % relFiles)
    else:
        releaseDir = 'REL'
        releaseDir = inputDir / 'REL'
        relFiles = sorted(list(releaseDir.glob('*.shp')))
    log.info('Release area files are: %s' % [str(relFilestr) for relFilestr in relFiles])

    # check for release thickness file if relThFromFile
    if cfg['GENERAL'].getboolean('relThFromFile'):
        relThFile, entResInfo['releaseThicknessFile'] = getAndCheckInputFiles(inputDir, 'RELTH',
                                                                              'release thickness data', fileExt='asc')
    else:
        relThFile = ''

    # Initialise secondary release areas
    secondaryReleaseFile, entResInfo['flagSecondaryRelease'] = getAndCheckInputFiles(inputDir, 'SECREL',
                                                                                     'Secondary release', fileExt='shp')

    # Initialise resistance areas
    resFile, entResInfo['flagRes'] = getAndCheckInputFiles(inputDir, 'RES', 'Resistance', fileExt='shp')

    # Initialise entrainment areas
    entFile, entResInfo['flagEnt'] = getAndCheckInputFiles(inputDir, 'ENT', 'Entrainment', fileExt='shp')

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    # return DEM, first item of release, entrainment and resistance areas
    inputSimFiles = {'demFile': demFile, 'relFiles': relFiles, 'secondaryReleaseFile': secondaryReleaseFile,
                     'entFile': entFile, 'resFile': resFile, 'entResInfo': entResInfo, 'relThFile': relThFile}

    return inputSimFiles

def getCellsize(avaDir):
    
    modName = 'r_avaflow'
    modNameCom1DFA = 'com1DFA'
    cfgAvaflow = cfgUtils.getModuleConfig_avaflow(modName, fileOverride='', toPrint=False)
    cfgCom1DFA = cfgUtils.getModuleConfig_avaflow2(modNameCom1DFA, fileOverride='', toPrint=False)
    
    cfgGen = cfgAvaflow["GENERAL"]
    cfgCom1DFA = cfgCom1DFA["GENERAL"]
    
    paramflag = cfgGen["paramflag"]
    cellsize_ravaflow = cfgGen["cellsize"]
    cellsize_com1DFA = cfgCom1DFA["meshCellSize"]
    
    if paramflag == "1":
        cellsize_return = cellsize_com1DFA
    else:
        cellsize_return = cellsize_ravaflow
        
    return cellsize_return, cellsize_ravaflow
    
    

def getInputDataravaflow(avaDir, cfg):
    
    # Transfer DEM file from com1DFA directory to avaflow input
    inputDir = pathlib.Path(avaDir, 'Inputs')
    demFile = list(inputDir.glob('*.asc'))
    if len(demFile) > 1:
        message = 'There should be exactly one topography .asc file in %s/Inputs/' % (avaDir)
        log.error(message)
        raise AssertionError(message)
    
    if len(demFile) == 0:
        print('No topography .asc-file supplied in com1DFA directory -> Topography file in avaflow directory will be used' )

    if len(demFile) == 1:
        demFp = demFile[0]
        
        orig = demFp
        target = os.path.join(avaDir, 'Avaflow_Input')
        
        shutil.copy(orig, target)
    
    # Transfer release file from com1DFA directory to avaflow input
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/REL')
    shpFile = list(inputDir_rel.glob('*.shp'))
    
    if len(shpFile) == 0:
        print('No release .shp-file supplied in com1DFA directory -> Release file in avaflow directory will be used')

    if len(shpFile) >= 1:
        # Transfer file also to avaflow input
        src_dir = os.path.join(avaDir, 'Inputs', 'REL')
        dest_dir = os.path.join(avaDir, 'Avaflow_Input')
        # fetching all the files in the source directory
        files=os.listdir(src_dir)
     
        for fname in files:
             
            shutil.copy2(os.path.join(src_dir,fname), dest_dir)
       
    # Transfer secondary release file from com1DFA directory to avaflow input
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/SECREL')
    shpFile = list(inputDir_rel.glob('*.shp'))
    
    if len(shpFile) == 0:
        print('No secondary release .shp-file')

    if len(shpFile) == 1:
        # Transfer file also to avaflow input
        src_dir = os.path.join(avaDir, 'Inputs', 'SECREL')
        dest_dir = os.path.join(avaDir, 'Avaflow_Input')
        # fetching all the files in the source directory
        files=os.listdir(src_dir)
     
        for fname in files:
             
            shutil.copy2(os.path.join(src_dir,fname), dest_dir)
            
    # Transfer entrainment area file from com1DFA directory to avaflow input
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/ENT')
    shpFile = list(inputDir_rel.glob('*.shp'))
    
    if len(shpFile) == 0:
        print('No entrainment area .shp-file')

    if len(shpFile) == 1:
        # Transfer file also to avaflow input
        src_dir = os.path.join(avaDir, 'Inputs', 'ENT')
        dest_dir = os.path.join(avaDir, 'Avaflow_Input')
        # fetching all the files in the source directory
        files=os.listdir(src_dir)
     
        for fname in files:
             
            shutil.copy2(os.path.join(src_dir,fname), dest_dir)
            
    # Transform those files to ascii
    shp_to_ascii.shptoasc(avaDir)
    
    ### Change the cellsize of the ascii grids to defined value in cfg-file 
    # cellsize = getCellsize(avaDir)
    # cellsize_change.Cellsize_change(avaDir, int(cellsize[0]))
    
    #Internal storage in r-avaflow needs to be increased, sofar only larger grids are possible
    cellsize_change.Cellsize_change(avaDir, 20)
    


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


def getThickness(inputSimFiles, avaDir, modName, cfgInitial):
    """ add thickness of shapefiles to dictionary, create one ini file per releaseScenario and
        set thickness values in ini files

        Parameters
        -----------
        inputSimFiles: dict
            dictionary with info on release and entrainment file paths
        avaDir: str or pathlib path
            path to avalanche directory
        modName : computational module
            computational module
        cfg: configParser object
            configParser object with the current (and possibly overridden) configuration

        Returns
        --------
        inputSimFiles: dict
            updated dictionary with thickness info read from shapefile attributes
            now includes one separate dictionary for each release, entrainment or secondary release
            scenario with a thickness and id value for each feature (given as list)
        cfgFilesRels: list
            list of updated ini files - one for each release Scenario

    """

    # create pathlib Path
    avaDir = pathlib.Path(avaDir)

    # get name of module as string
    modNameString = str(pathlib.Path(modName.__file__).stem)

    # initialise list for cfgFiles
    cfgFilesRels = []

    # check if thickness info is required from entrainment and secondary release according to simType
    simTypeList = cfgInitial['GENERAL']['simTypeList'].split('|')
    thTypeList = []
    if any(simType in ['ent', 'entres', 'available'] for simType in simTypeList):
        thTypeList.append('entFile')
    if cfgInitial['GENERAL'].getboolean('secRelArea'):
        thTypeList.append('secondaryReleaseFile')

    # fetch thickness attribute of entrainment area and secondary release
    for thType in ['entFile', 'secondaryReleaseFile']:
        if inputSimFiles[thType] != None:
            thicknessList, idList, ci95List = shpConv.readThickness(inputSimFiles[thType])
            inputSimFiles[inputSimFiles[thType].stem] = {'thickness': thicknessList, 'id': idList,
                'ci95': ci95List}

    # fetch thickness attribute of release areas and create cfg file for each release scenario
    for releaseA in inputSimFiles['relFiles']:
        # fetch thickness and id info from input data
        thicknessList, idList, ci95List = shpConv.readThickness(releaseA)
        inputSimFiles[releaseA.stem] = {'thickness': thicknessList, 'id': idList, 'ci95': ci95List}

        # add input data info
        cfgInitial['INPUT'] = {'DEM': inputSimFiles['demFile'].stem, 'releaseScenario': releaseA.stem}
        # update configuration with thickness value to be used for simulations
        cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, releaseA.stem, 'relTh')
        if cfgInitial['GENERAL'].getboolean('relThFromFile'):
            cfgInitial['INPUT']['relThFile'] = str(pathlib.Path('RELTH', inputSimFiles['relThFile'].name))

        # add entrainment and secondary release thickness in input data info
        if inputSimFiles['entFile'] != None and 'entFile' in thTypeList:
            cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, inputSimFiles['entFile'].stem, 'entTh')
            cfgInitial['INPUT']['entrainmentScenario'] = inputSimFiles['entFile'].stem

        if inputSimFiles['secondaryReleaseFile'] != None and 'secondaryReleaseFile' in thTypeList:
            cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles,
                inputSimFiles['secondaryReleaseFile'].stem, 'secondaryRelTh')
            cfgInitial['INPUT']['secondaryReleaseScenario'] = inputSimFiles['secondaryReleaseFile'].stem

        # create new ini file for each release scenario with updated info on thickness values (and parameter variation)
        cfgFileRelease = avaDir / 'Outputs' / modNameString / ('%s_com1DFACfg.ini' % releaseA.stem)
        with open(cfgFileRelease, 'w') as configfile:
            cfgInitial.write(configfile)
        cfgFilesRels.append(cfgFileRelease)

    return inputSimFiles, cfgFilesRels


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


def selectReleaseScenario(inputSimFiles, cfg):
    """ select release scenario and remove other release files in inputSimFiles dictionary

        Parameters
        -----------
        inputSimFiles: dict
            dictionary with info on input data
        cfg: conigparser object
            configuration, here Flag releaseScenario is used

        Returns
        -------
        inputSimFiles: dict
            updated dictionary with only one releaseScenario
    """

    relFiles = inputSimFiles['relFiles']

    foundScenario = False
    for relF in relFiles:
        if relF.stem == cfg['releaseScenario']:
            releaseScenario = relF
            foundScenario = True

    if foundScenario is False:
        message = 'Release area scenario %s not found - check input data' % (cfg['releaseScenario'])
        log.error(message)
        raise FileNotFoundError(message)

    inputSimFiles['relFiles'] = [releaseScenario]
    # add release area scenario
    inputSimFiles['releaseScenario'] = releaseScenario

    return inputSimFiles
