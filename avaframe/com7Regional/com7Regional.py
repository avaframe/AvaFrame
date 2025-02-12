"""Module for handling regional avalanche simulations."""

import pathlib
import shutil
import logging
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

import avaframe.in3Utils.initializeProject as initProj
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils, cfgHandling
from avaframe.in3Utils import logUtils
from avaframe.in2Trans import rasterUtils

# create local logger
log = logging.getLogger(__name__)

def com7RegionalMain(cfgMain, cfg):
    """Run com7Regional with given configuration.
    
    This function processes multiple avalanche directories in parallel, running simulations
    for each directory.
    
    Parameters
    ----------
    cfgMain : configparser.ConfigParser
        Main configuration settings
    cfg : configparser.ConfigParser
        Regional configuration settings with potential overrides
        
    Returns
    -------
    allPeakFilesDir : pathlib.Path or None
        Path to the directory containing all peak files, if copyPeakFiles is True
    allTimeStepsDir : pathlib.Path or None
        Path to the directory containing all time step files, if copyPeakFiles is True
    mergedRastersDir : pathlib.Path or None
        Path to the directory containing merged rasters, if mergeOutput is True

    Notes
    -----
    The function expects the following directory structure:
    avalancheDir/
    └── regionalDir/
        ├── avalanche1/
        ├── avalanche2/
        └── ...

    Where:
    - avalancheDir: Main directory specified in cfgMain
    - regionalDir: Subdirectory specified in cfg['GENERAL']['regionalDir']
    """
    # Define the regional directory in relation to the avalanche directory
    regionalDirFromCfg = str(cfg['GENERAL']['regionalDir'])
    regionalDir = pathlib.Path(cfgMain['MAIN']['avalancheDir']) / regionalDirFromCfg

    # List valid avalanche directories within the regional directory
    avaDirs = findAvaDirs(regionalDir)

    # Get total number of simulations
    log.info(f"Getting total number of simulations to perform...")
    with logUtils.silentLogger():
        totalSims = getTotalNumberOfSims(avaDirs, cfgMain, cfg)
    log.info(f"Found {totalSims} (new) simulations to perform across {len(avaDirs)} directories")

    # Get number of processes based on number of avaDirs
    nProcesses = cfgUtils.getNumberOfProcesses(cfgMain, len(avaDirs))

    # Set nCPU for com1 to 1 to avoid nested parallelization
    cfgMain['MAIN']['nCPU'] = '1'

    # Track progress and results
    completed = 0
    nSuccesses = 0

    # Process avalanche directories within the regional folder concurrently
    with ProcessPoolExecutor(max_workers=nProcesses) as executor:
        # Submit each avalanche directory to the executor
        futures = {executor.submit(processAvaDirCom1Regional, cfgMain, cfg, avaDir):
                       avaDir for avaDir in avaDirs}
        # Log results as each future completes
        for future in as_completed(futures):
            avaDir = futures[future]
            try:
                resultDir, status = future.result()
                completed += 1
                
                if status == "Success":
                    nSuccesses += 1

                log.info(f"{status} in directory: {resultDir.relative_to(pathlib.Path(regionalDir))} "
                         f"- Overall progress: {completed}/{len(avaDirs)}")
            except Exception as e:
                log.error(f"Error processing {avaDir}: {e}")

    log.info(f"Processing complete. Success in {nSuccesses} out of {len(avaDirs)} directories.")

    # Copy/move peak files if configured
    allPeakFilesDir = None
    allTimeStepsDir = None
    if cfg['GENERAL'].getboolean('copyPeakFiles'):
        allPeakFilesDir, allTimeStepsDir = moveOrCopyPeakFiles(cfg, regionalDir, avaDirs)

    # Merge output rasters if configured
    mergedRastersDir = None
    if cfg['GENERAL'].getboolean('mergeOutput'):
        mergedRastersDir = mergeOutputRasters(cfg, regionalDir, allPeakFilesDir)

    return allPeakFilesDir, allTimeStepsDir, mergedRastersDir

def findAvaDirs(Dir):
    """Find all valid avalanche directories within a given directory.

    A directory is considered a valid avalanche directory if it contains an "Inputs" folder.

    Parameters
    ----------
    Dir : pathlib.Path or str
        Path to the directory to search in

    Returns
    -------
    avaDirs : list
        List of pathlib.Path objects pointing to valid avalanche directories
    """
    avaDirs = [pathlib.Path(p).parent for p in pathlib.Path(Dir).glob("*/Inputs")]
    log.info(f"Found a total of '{len(avaDirs)}' avalanche directories in: {Dir}:")
    for avaDir in avaDirs:
        log.info(f"'{avaDir.name}'")

    return avaDirs

def getTotalNumberOfSims(avaDirs, cfgMain, cfgCom7):
    """Get total number of simulations across all avalanche directories.
    
    Parameters
    ----------
    avaDirs : list
        List of avalanche directories
    cfgMain : configparser.ConfigParser
        Main configuration
    cfgCom7 : configparser.ConfigParser
        Regional configuration with potential overrides
        
    Returns
    -------
    int
        Total number of simulations
    """
    totalSims = 0
    for avaDir in avaDirs:
        # Create copies of configs to avoid modifying originals
        cfgMainCopy = cfgUtils.convertDictToConfigParser(cfgUtils.convertConfigParserToDict(cfgMain))
        cfgMainCopy['MAIN']['avalancheDir'] = str(avaDir)
        
        # Get com1DFA config with regional overrides (same as in processAvaDirCom1Regional)
        cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False,
                                           onlyDefault=cfgCom7['com1DFA_com1DFA_override'].getboolean('defaultConfig'))
        cfgCom1DFA, _ = cfgHandling.applyCfgOverride(cfgCom1DFA, cfgCom7, com1DFA, addModValues=False)
        
        # Get simulations for this directory
        try:
            simDict, _, _, _ = com1DFA.com1DFAPreprocess(cfgMainCopy, "cfgFromObject", cfgCom1DFA)
            totalSims += len(simDict)
        except Exception as e:
            log.warning(f"Could not get simulations for {avaDir}: {e}")
            continue
            
    return totalSims

def processAvaDirCom1Regional(cfgMain, cfgCom7, avalancheDir):
    """Run com1DFA simulation in a specific avalanche directory with regional settings.

    Note: This function calls com1DFA within each avalanche directory within the input directory. 
    If wanted it may be used as a template to call another operation within each directory, such as com2AB, ana5Utils, etc.

    Parameters
    ----------
    cfgMain : configparser.ConfigParser
        Main configuration settings
    cfgCom7 : configparser.ConfigParser
        Regional configuration settings with potential overrides
    avalancheDir : pathlib.Path or str
        Path to the avalanche directory to process

    Returns
    -------
    avalancheDir : pathlib.Path or str
        Path to the avalanche directory that was processed
    status : str
        Status of the simulation, "Success" if completed
    """
    # Initialize log for each process
    log = logUtils.initiateLogger(avalancheDir, logName='runCom1DFA')
    log.info('COM1DFA PROCESS CALLED BY COM7REGIONAL RUN')
    log.info('Current avalanche: %s', avalancheDir)

    # Update cfgMain setting to reflect the current avalancheDir
    cfgMain['MAIN']['avalancheDir'] = str(avalancheDir)

    # Clean input directory of old work and output files from module
    initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)

    # Create com1DFA configuration for the current avalanche directory and override with regional settings
    cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False,
                                          onlyDefault=cfgCom7['com1DFA_com1DFA_override'].getboolean('defaultConfig'))
    cfgCom1DFA, cfgCom7 = cfgHandling.applyCfgOverride(cfgCom1DFA, cfgCom7, com1DFA, addModValues=False)

    # Run com1DFA in the current avalanche directory
    com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgCom1DFA)

    return avalancheDir, "Success"

def findPeakFiles(directory, pattern):
    """Find peak files with given pattern and extension.
    
    Parameters
    ----------
    directory : pathlib.Path
        Directory to search in
    pattern : str
        Base pattern to match, e.g., 'Outputs/**/peakFiles/*' or '*_ppr'
        
    Returns
    -------
    list
        List of found files with both .asc and .tif extensions
    """
    ascFiles = list(directory.glob(pattern + '.asc'))
    tifFiles = list(directory.glob(pattern + '.tif'))
    return ascFiles + tifFiles

def moveOrCopyPeakFiles(cfg, avalancheDir, avaDirs):
    """Consolidate peak files from multiple avalanche directories.

    Creates two directories:
    1. allPeakFiles: Contains peak files from all avalanche directories
    2. allPeakFiles/allTimeSteps: Contains time step files from all avalanche directories

    Parameters
    ----------
    cfg : configparser.ConfigParser
        Configuration containing GENERAL settings:
        - copyPeakFiles: If True, copy/move files; if False, do nothing
        - moveInsteadOfCopy: If True, move files instead of copying
    avalancheDir : pathlib.Path or str
        Base directory where allPeakFiles will be created
    avaDirs : list
        List of avalanche directories

    Returns
    -------
    allPeakFilesDir : pathlib.Path or None
        Path to the created allPeakFiles directory or None if copyPeakFiles is False
    allTimeStepsDir : pathlib.Path or None
        Path to the created allTimeSteps directory or None if copyPeakFiles is False
    """
    if not cfg['GENERAL'].getboolean('copyPeakFiles'):
        log.info("copyPeakFiles is False - no files will be copied or moved")
        return None, None

    # Set up outdirs
    allPeakFilesDir = pathlib.Path(avalancheDir, 'Outputs', 'com7Regional', 'allPeakFiles')
    allTimeStepsDir = allPeakFilesDir / 'allTimeSteps'
    for dirPath in [allPeakFilesDir, allTimeStepsDir]:
        if dirPath.exists():
            shutil.rmtree(str(dirPath))
        dirPath.mkdir(parents=True, exist_ok=True)

    # Get operation type
    operation = shutil.move if cfg['GENERAL'].getboolean('moveInsteadOfCopy') else shutil.copy2
    operationType = 'Moving' if operation == shutil.move else 'Copying'

    # Process each avalanche directory
    for avaDir in avaDirs:
        log.info(f"{operationType} files from: {avaDir.name}")
        
        # Process peak files
        peakFiles = findPeakFiles(avaDir, 'Outputs/**/peakFiles/*')
        for peakFile in peakFiles:
            targetPath = allPeakFilesDir / peakFile.name
            operation(str(peakFile), str(targetPath))

        # Process time step files
        timeStepFiles = findPeakFiles(avaDir, 'Outputs/**/peakFiles/timeSteps/*')
        for timeStepFile in timeStepFiles:
            targetPath = allTimeStepsDir / timeStepFile.name
            operation(str(timeStepFile), str(targetPath))

    return allPeakFilesDir, allTimeStepsDir

def getRasterBounds(rasterFiles):
    """Get the union bounds and validate cell sizes of multiple rasters.

    Parameters
    ----------
    rasterFiles : list
        List of paths to raster files

    Returns
    -------
    bounds : dict
        Dictionary containing xmin, xmax, ymin, ymax
    cellSize : float
        Cell size of the rasters

    Raises
    ------
    ValueError
        If cell sizes of rasters differ
    """
    # Read first raster to get cellSize and initialize bounds
    firstRaster = rasterUtils.readRaster(rasterFiles[0])
    header = firstRaster['header']
    cellSize = float(header['cellsize'])

    # Initialize bounds with first raster
    bounds = {
        'xmin': header['xllcenter'] - cellSize/2,
        'xmax': header['xllcenter'] + (header['ncols'] - 0.5) * cellSize,
        'ymin': header['yllcenter'] - cellSize/2,
        'ymax': header['yllcenter'] + (header['nrows'] - 0.5) * cellSize
    }

    # Update bounds with remaining rasters
    for rasterFile in rasterFiles:
        raster = rasterUtils.readRaster(rasterFile)
        header = raster['header']

        # Validate cell size
        if abs(float(header['cellsize']) - cellSize) > 1e-6:
            raise ValueError("Cell sizes of rasters differ")

        # Update bounds
        xmin = header['xllcenter'] - cellSize/2
        xmax = header['xllcenter'] + (header['ncols'] - 0.5) * cellSize
        ymin = header['yllcenter'] - cellSize/2
        ymax = header['yllcenter'] + (header['nrows'] - 0.5) * cellSize

        bounds['xmin'] = min(bounds['xmin'], xmin)
        bounds['xmax'] = max(bounds['xmax'], xmax)
        bounds['ymin'] = min(bounds['ymin'], ymin)
        bounds['ymax'] = max(bounds['ymax'], ymax)

    return bounds, cellSize


def mergeRasters(rasterFiles, bounds, cellSize, noDataValue=0, mergeMethod='max'):
    """Merge multiple rasters into a single raster.

    Parameters
    ----------
    rasterFiles : list
        List of paths to raster files
    bounds : dict
        Dictionary containing xmin, xmax, ymin, ymax
    cellSize : float
        Cell size of the rasters
    noDataValue : float, optional
        Value to use for no data cells, default: 0
    mergeMethod : str, optional
        Method to use for merging overlapping cells, default: 'max'
        Valid options: 'max', 'min', 'mean', 'sum', 'count'

    Returns
    -------
    mergedHeader : dict
        Header information for the merged raster
    mergedData : numpy.ndarray
        2D array containing the merged raster data
    """
    # Calculate dimensions for merged raster
    nrows = int((bounds['ymax'] - bounds['ymin']) / cellSize)
    ncols = int((bounds['xmax'] - bounds['xmin']) / cellSize)

    # Create merged raster header
    mergedHeader = {
        'ncols': ncols,
        'nrows': nrows,
        'xllcenter': bounds['xmin'] + cellSize / 2,
        'yllcenter': bounds['ymin'] + cellSize / 2,
        'cellsize': cellSize,
        'nodata_value': noDataValue,
        'transform': rasterUtils.transformFromASCHeader({
            'ncols': ncols,
            'nrows': nrows,
            'xllcenter': bounds['xmin'] + cellSize / 2,
            'yllcenter': bounds['ymin'] + cellSize / 2,
            'cellsize': cellSize
        })
    }

    # Get driver and CRS from first raster
    firstRaster = rasterUtils.readRaster(rasterFiles[0])
    mergedHeader['driver'] = firstRaster['header']['driver']
    mergedHeader['crs'] = firstRaster['header']['crs']

    # Initialize merged data array
    if mergeMethod == 'count':
        mergedData = np.zeros((nrows, ncols), dtype=np.int32)
    else:
        mergedData = np.full((nrows, ncols), noDataValue, dtype=np.float32)
        if mergeMethod == 'min':
            mergedData.fill(np.inf)
        elif mergeMethod == 'max':
            mergedData.fill(-np.inf)

    # Merge rasters
    for rasterFile in rasterFiles:
        raster = rasterUtils.readRaster(rasterFile)
        header = raster['header']
        data = raster['rasterData']

        # Calculate indices in merged raster
        xstart = int((header['xllcenter'] - cellSize/2 - bounds['xmin']) / cellSize)
        ystart = int((header['yllcenter'] - cellSize/2 - bounds['ymin']) / cellSize)
        
        # Handle data based on merge method
        validMask = ~np.isnan(data) & (data != header['nodata_value']) & (data != 0)
        if mergeMethod == 'max':
            np.maximum.at(mergedData, (slice(ystart, ystart + header['nrows']), slice(xstart, xstart + header['ncols'])), 
                         np.where(validMask, data, -np.inf))
        elif mergeMethod == 'min':
            np.minimum.at(mergedData, (slice(ystart, ystart + header['nrows']), slice(xstart, xstart + header['ncols'])), 
                         np.where(validMask, data, np.inf))
        elif mergeMethod == 'sum':
            mergedData[ystart:ystart + header['nrows'], xstart:xstart + header['ncols']] += np.where(validMask, data, 0)
        elif mergeMethod == 'mean':
            validData = mergedData[ystart:ystart + header['nrows'], xstart:xstart + header['ncols']]
            validData = np.where(validMask, data, validData)
            mergedData[ystart:ystart + header['nrows'], xstart:xstart + header['ncols']] = validData
        elif mergeMethod == 'count':
            mergedData[ystart:ystart + header['nrows'], xstart:xstart + header['ncols']] += validMask.astype(np.int32)

    # Post-process based on merge method
    if mergeMethod in ['max', 'min']:
        # Replace inf values with nodata
        mergedData = np.where(np.isinf(mergedData), noDataValue, mergedData)
    elif mergeMethod == 'count':
        # Replace zero counts with nodata
        mergedData = np.where(mergedData > 0, mergedData, noDataValue)

    return mergedHeader, mergedData


def mergeOutputRasters(cfg, avalancheDir, allPeakFilesDir=None):
    """Merge output rasters (peakFiles) from all avalanche simulations.

    Parameters
    ----------
    cfg : configparser object
        Configuration settings containing:
        - GENERAL.mergeOutput: If True, merge rasters
        - GENERAL.mergeTypes: Types of rasters to merge (e.g., 'ppr|pfv|pft')
        - GENERAL.mergeMethods: Methods to use for merging (e.g., 'max')
    avalancheDir : pathlib.Path or str
        Base directory where merged rasters will be saved
    allPeakFilesDir : pathlib.Path, optional
        Directory containing all peak files. If None, will search in individual avalanche directories.

    Returns
    -------
    mergedRastersDir : pathlib.Path
        Path to the directory containing merged rasters or None if mergeOutput is False
    """
    if not cfg['GENERAL'].getboolean('mergeOutput'):
        log.info("mergeOutput is False - no rasters will be merged")
        return None

    # Set up merged rasters directory
    mergedRastersDir = pathlib.Path(avalancheDir, 'Outputs', 'com7Regional', 'mergedRasters')
    mergedRastersDir.mkdir(parents=True, exist_ok=True)

    # Get merge types and methods from config
    mergeTypes = cfg['GENERAL']['mergeTypes'].lower().split('|')
    mergeMethods = cfg['GENERAL']['mergeMethods'].lower().split('|')
    log.info(f"Merging raster types: {mergeTypes}")
    log.info(f"Using merge methods: {mergeMethods}")

    # Get all peak files directories
    allPeakFilesDirs = []
    if allPeakFilesDir is not None and allPeakFilesDir.exists():
        # If we have a consolidated directory, use that
        allPeakFilesDirs = [allPeakFilesDir]
    else:
        # Otherwise search in individual avalanche directories
        for avaDir in findAvaDirs(avalancheDir):
            peakFilesDir = avaDir / 'Outputs' / 'com1DFA' / 'peakFiles'
            if peakFilesDir.is_dir():
                allPeakFilesDirs.append(peakFilesDir)

    # Process each raster type
    for rasterType in mergeTypes:
        # Collect all raster files for this type
        rasterFiles = []
        for peakFilesDir in allPeakFilesDirs:
            rasterFiles.extend(findPeakFiles(peakFilesDir, f'*_{rasterType}'))

        if not rasterFiles:
            log.warning(f"No {rasterType} rasters found to merge")
            continue

        log.info(f"Merging {len(rasterFiles)} {rasterType} rasters")
        
        # Get bounds and validate cell sizes
        bounds, cellSize = getRasterBounds(rasterFiles)

        # Merge and save rasters
        for mergeMethod in mergeMethods:
            mergedHeader, mergedData = mergeRasters(rasterFiles, bounds, cellSize, noDataValue=0, mergeMethod=mergeMethod)
            outputPath = mergedRastersDir / f'merged_{rasterType}_{mergeMethod}'
            rasterUtils.writeResultToRaster(mergedHeader, mergedData, outputPath, flip=True)
            log.info(f"Saved merged {rasterType} raster (method: {mergeMethod}) to: {outputPath}")

    return mergedRastersDir
