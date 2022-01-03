"""
    Fetch input data for avalanche simulations
"""

# Load modules
import os
import glob
import pathlib
import logging

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.com1DFA.deriveParameterSet as dP



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

    inputDir = pathlib.Path(avaDir, 'Inputs')
    demFile = list(inputDir.glob('*.asc'))
    if not len(demFile) == 1:
        message = 'There should be exactly one topography .asc file in %s/Inputs/' % (avaDir)
        log.error(message)
        raise AssertionError(message)

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
    resFile, entResInfo['flagRes'] = getAndCheckInputFiles(inputDir, 'RES', 'Resistance')
    if resFile==None:
        resFile = ''
    # Initialise entrainment areas
    entFile, entResInfo['flagEnt'] = getAndCheckInputFiles(inputDir, 'ENT', 'Entrainment')
    if entFile==None:
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
    entResInfo= {'flagEnt': 'No', 'flagRes': 'No'}

    # Initialise release areas, default is to look for shapefiles
    if cfg['releaseScenario'] != '':
        releaseDir = 'REL'
        relFiles = []
        releaseFiles = cfg['releaseScenario'].split('|')
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

    # Initialise secondary release areas
    secondaryReleaseFile, entResInfo['flagSecondaryRelease'] = getAndCheckInputFiles(inputDir, 'SECREL', 'Secondary release')

    # Initialise resistance areas
    resFile, entResInfo['flagRes'] = getAndCheckInputFiles(inputDir, 'RES', 'Resistance')

    # Initialise entrainment areas
    entFile, entResInfo['flagEnt'] = getAndCheckInputFiles(inputDir, 'ENT', 'Entrainment')

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    # return DEM, first item of release, entrainment and resistance areas
    inputSimFiles = {'demFile': demFile, 'relFiles': relFiles, 'secondaryReleaseFile': secondaryReleaseFile,
                     'entFile': entFile, 'resFile': resFile, 'entResInfo': entResInfo}


    return inputSimFiles


def getAndCheckInputFiles(inputDir, folder, inputType):
    """Fetch shape files and check if they exist and if it is not more than one

    Raises error if there is more than one shape file.

    Parameters
    ----------
    inputDir : pathlib object or str
        path to avalanche input directory
    folder : str
        subfolder name where the shape file should be located (SECREL, ENT or RES)
    inputType : str
        type of input (used for the logging messages). Secondary release or Entrainment or Resistance

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
    OutputFile = list(dir.glob('*.shp'))
    if len(OutputFile) < 1:
        OutputFile = None
    elif len(OutputFile) > 1:
        message = 'More than one %s .shp file in %s/%s/ not allowed' % (inputType, inputDir, folder)
        log.error(message)
        raise AssertionError(message)
    else:
        available = 'Yes'
        OutputFile = OutputFile[0]

    return OutputFile, available


def getThickness(inputSimFiles, avaDir, modName, cfgFile):
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
        cfgFile: str
            path to cfgFile for reading overall config of comModule - if empty read local or default ini file

        Returns
        --------
        inputSimFiles: dict
            updated dictionary with thickness info read from shapefile attributes
        cfgFilesRels: list
            list of updated ini files - one for each release Scenario

    """

    # create pathlib Path
    avaDir = pathlib.Path(avaDir)

    # get name of module as string
    modNameString = str(pathlib.Path(modName.__file__).stem)

    # initialise list for cfgFiles
    cfgFilesRels = []

    # fetch thickness attribute of entrainment area and secondary release
    for thType in ['entFile', 'secondaryReleaseFile']:
        if inputSimFiles[thType] != None:
            thicknessList, idList = shpConv.readThickness(inputSimFiles[thType])
            inputSimFiles[inputSimFiles[thType].stem] = {'thickness': thicknessList, 'id': idList}

    # fetch thickness attribute of release areas
    for releaseA in inputSimFiles['relFiles']:
        cfgInitial = cfgUtils.getModuleConfig(modName, fileOverride=cfgFile, toPrint=False)
        thicknessList, idList = shpConv.readThickness(releaseA)
        inputSimFiles[releaseA.stem] = {'thickness': thicknessList, 'id': idList}

        # update configuration with thickness value to be used for simulations
        cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, releaseA.stem, 'relTh')

        # set input data info
        cfgInitial['INPUT'] = {'DEM': inputSimFiles['demFile'].stem , 'releaseScenario': releaseA.stem}
        # set entrainment and secondary release thickness  and input data info
        if inputSimFiles['entFile'] != None:
            cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, inputSimFiles['entFile'].stem, 'entTh')
            cfgInitial['INPUT']['entrainmentScenario'] =inputSimFiles['entFile'].stem
        if inputSimFiles['secondaryReleaseFile'] != None:
            cfgInitial = dP.getThicknessValue(cfgInitial, inputSimFiles, inputSimFiles['secondaryReleaseFile'].stem, 'secondaryRelTh')
            cfgInitial['INPUT']['secondaryReleaseScenario'] =inputSimFiles['secondaryReleaseFile'].stem

        # create new ini file for each release scenario with updated info on thickness values (and parameter variation)
        cfgFileRelease = avaDir / 'Outputs' / modNameString / ('%s_com1DFACfg.ini' % releaseA.stem)
        with open(cfgFileRelease, 'w') as configfile:
            cfgInitial.write(configfile)
        cfgFilesRels.append(cfgFileRelease)


    return inputSimFiles, cfgFilesRels


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

    for relF in relFiles:
        if relF.stem == cfg['releaseScenario']:
            releaseScenario = relF

    inputSimFiles['relFiles'] = [releaseScenario]

    return inputSimFiles
