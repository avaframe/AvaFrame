"""
    Fetch input data for avalanche simulations
"""

# Load modules
import os
import glob
import logging

# Local imports
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils


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

    demFile = glob.glob(os.path.join(avaDir, 'Inputs', '*.asc'))
    if not len(demFile) == 1:
        message = 'There should be exactly one topography .asc file in %s/Inputs/' % (avaDir)
        log.error(message)
        raise AssertionError(message)

    return demFile[0]


def getInputData(avaDir, cfg, flagDev=False):
    """ Fetch input datasets required for simulation

    Parameters
    ----------
    avaDir : str
        path to avalanche directory
    cfg : dict
        configuration read from com1DFA simulation ini file
    flagDev : bool
        optional - if True: use for devREL folder to get release area scenarios

    Returns
    -------
    demFile[0] : str (first element of list)
        list of full path to DEM .asc file
    relFiles : list
        list of full path to release area scenario .shp files
    entFiles[0] : str (fist element of list)
        list of full path to entrainment area .shp files
    resFiles[0] : str (first element of list)
        list of full path to resistance area .shp files
    flagEntRes : bool
        flag if True entrainment and/or resistance areas found and used for simulation
    """

    # Set directories for inputs, outputs and current work
    inputDir = os.path.join(avaDir, 'Inputs')

    # Set flag if there is an entrainment or resistance area
    entResInfo= {'flagEntRes': False, 'flagEnt': 'No', 'flagRes': 'No'}

    # Initialise release areas, default is to look for shapefiles
    if flagDev is True:
        releaseDir = 'devREL'
        relFiles = glob.glob(inputDir+os.sep + releaseDir+os.sep + '*.shp')
    elif cfg['releaseScenario'] != '':
        releaseDir = 'REL'
        relFiles = []
        releaseFiles = cfg['releaseScenario'].split('|')
        for rel in releaseFiles:
            if '.shp' in rel:
                relf = os.path.join(inputDir, releaseDir, rel)
            else:
                relf = os.path.join(inputDir, releaseDir, '%s.shp' % (rel))
            if os.path.isfile(relf) is True:
                relFiles.append(relf)
            else:
                log.error('No release scenario called: %s' % (relf))
        log.debug('Release area file is specified to be: %s' % relFiles)
    else:
        releaseDir = 'REL'
        relFiles = sorted(glob.glob(inputDir+os.sep + releaseDir+os.sep + '*.shp'))
    log.info('Release area files are: %s' % relFiles)

    # Initialise resistance areas
    if cfg.getboolean('noResistance'):
        resFiles = []
        resFiles.append('')
    else:
        resFiles = glob.glob(inputDir+os.sep + 'RES' + os.sep+'*.shp')
        if len(resFiles) < 1:
            log.debug('No resistance file found')
            resFiles.append('')  # Kept this for future enhancements
        else:
            try:
                message = 'There shouldn\'t be more than one resistance .shp file in ' + inputDir + '/RES/'
                assert len(resFiles) < 2, message
            except AssertionError:
                raise
            entResInfo['flagRes'] = 'Yes'
            entResInfo['flagEntRes'] = True

    # Initialise entrainment areas
    if cfg.getboolean('noEntrainment'):
        entFiles = []
        entFiles.append('')
    else:
        entFiles = glob.glob(inputDir+os.sep + 'ENT' + os.sep+'*.shp')
        if len(entFiles) < 1:
            log.debug('No entrainment file found')
            entFiles.append('')  # Kept this for future enhancements
        else:
            try:
                message = 'There shouldn\'t be more than one entrainment .shp file in ' + inputDir + '/ENT/'
                assert len(entFiles) < 2, message
            except AssertionError:
                raise
            entResInfo['flagEnt'] = 'Yes'
            entResInfo['flagEntRes'] = True

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    # return DEM, first item of release, entrainment and resistance areas
    return demFile, relFiles, entFiles[0], resFiles[0], entResInfo


def getInputDataCom1DFAPy(avaDir, cfg, flagDev=False):
    """ Fetch input datasets required for simulation, duplicated function because
        simulation type set differently in com1DFAPy compared to com1DFA:
        TODO: remove duplicate once it is not required anymore

    Parameters
    ----------
    avaDir : str
        path to avalanche directory
    cfg : dict
        configuration read from com1DFA simulation ini file
    flagDev : bool
        optional - if True: use for devREL folder to get release area scenarios

    Returns
    -------
    demFile[0] : str (first element of list)
        list of full path to DEM .asc file
    inputSimFiles: dict
        dictionary with all the input files:
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
    """

    # Set directories for inputs, outputs and current work
    inputDir = os.path.join(avaDir, 'Inputs')

    # Set flag if there is an entrainment or resistance area
    entResInfo= {'flagEnt': 'No', 'flagRes': 'No'}

    # Initialise release areas, default is to look for shapefiles
    if flagDev is True:
        releaseDir = 'devREL'
        relFiles = glob.glob(inputDir+os.sep + releaseDir+os.sep + '*.shp')
    elif cfg['FLAGS']['releaseScenario'] != '':
        releaseDir = 'REL'
        relFiles = []
        releaseFiles = cfg['FLAGS']['releaseScenario'].split('|')
        for rel in releaseFiles:
            if '.shp' in rel:
                relf = os.path.join(inputDir, releaseDir, rel)
            else:
                relf = os.path.join(inputDir, releaseDir, '%s.shp' % (rel))
            if not os.path.isfile(relf):
                message = 'No release scenario called: %s' % (relf)
                log.error(message)
                raise FileNotFoundError(message)
        log.debug('Release area file is specified to be: %s' % relFiles)
    else:
        releaseDir = 'REL'
        relFiles = sorted(glob.glob(inputDir+os.sep + releaseDir+os.sep + '*.shp'))
    log.info('Release area files are: %s' % relFiles)

    # Initialise secondary release areas
    secondaryReleaseFile, entResInfo['flagSecondaryRelease'] = getAndCheckInputFiles(inputDir, 'SECREL', 'Secondary release', flag=cfg.getboolean('GENERAL', 'secRelArea'))

    # Initialise resistance areas
    resFile, entResInfo['flagRes'] = getAndCheckInputFiles(inputDir, 'RES', 'Resistance')

    # Initialise entrainment areas
    entFile, entResInfo['flagEnt'] = getAndCheckInputFiles(inputDir, 'ENT', 'Entrainment')

    # Initialise DEM
    demFile = getDEMPath(avaDir)

    # return DEM, first item of release, entrainment and resistance areas
    inputSimFiles = {'relFiles': relFiles, 'secondaryReleaseFile': secondaryReleaseFile,
                     'entFile': entFile, 'resFile': resFile, 'entResInfo': entResInfo}
    return demFile, inputSimFiles


def getAndCheckInputFiles(inputDir, folder, inputType, flag=False):
    """Fetch input shape files and check them

    Raises error if there is more than one shape file.
    If the flag is set to true, raise an error if there is no shape file

    Parameters
    ----------
    inputDir : str
        path to avalanche input directory
    folder : str
        subfolder name where the shape file should be located (SECREL, ENT or RES)
    inputType : str
        type of input (used for the logging messages). Secondary release or Entrainment or Resistance
    flag : bool
        optional (false as default) - if True: raise error if there is no shape file

    Returns
    -------
    OutputFile: str
        path to file checked
    available: str
        Yes or No depending of if there is a shape file available (if No, OutputFile is None)
    """
    available = 'No'
    # Initialise secondary release areas
    OutputFile = glob.glob(inputDir+os.sep + folder + os.sep+'*.shp')
    if len(OutputFile) < 1:
        if flag:
            message = 'No %s area file found' % (inputType)
            log.error(message)
            raise FileNotFoundError(message)
        else:
            OutputFile = None
    elif len(OutputFile) > 1:
        message = 'There shouldn\'t be more than one %s .shp file in %s/%s/' % (inputType, inputDir, folder)
        log.error(message)
        raise AssertionError(message)
    else:
        available = 'Yes'
        OutputFile = OutputFile[0]
        log.info('%s area file is: %s' % (inputType, OutputFile))

    return OutputFile, available
