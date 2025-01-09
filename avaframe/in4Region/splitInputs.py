"""Module for splitting and organizing avalanche input data."""

import logging
import shapefile  # pyshp
from shapely.geometry import shape, box
from avaframe.in2Trans import ascUtils
import pathlib

from avaframe.in3Utils.initializeProject import initializeFolderStruct

# create local logger
log = logging.getLogger(__name__)

def splitInputsMain(inputDir, outputDir, cfg):
    """Process and organize avalanche input data into individual folders.

    Parameters
    ----------
    inputDir : pathlib.Path object
        Path to input directory containing release areas (REL) and DEM files
    outputDir : pathlib.Path object
        Path to output directory where organized folders will be created
    cfg : dict
        Configuration settings containing:
        - GENERAL.splitDEM : bool, whether to split DEM
        - GENERAL.bufferSize : float, buffer size for DEM clipping

    Returns
    -------
    none

    Notes
    -----
    Expected input directory structure:
    inputDir/
    ├── REL/
    │   └── release_areas.shp (with optional .prj)
    ├── ENT/ (optional)
    │   └── entrainment.shp (with optional .prj)
    ├── RES/ (optional)
    │   └── resistance.shp (with optional .prj)
    └── dem_file.asc
    """
    # Fetch the necessary input
    inputRELDir = inputDir / 'REL'
    inputShp = next(inputRELDir.glob("*.shp"), None)
    if not inputShp:
        log.error(f"No shapefile found in {inputRELDir}.")
        return
    inputDEM = next(inputDir.glob("*.asc"), None) #ToDo: needs adjustment once we include tif files as well
    if not inputDEM:
        log.error(f"No DEM file found in {inputDir}.")
        return

    # Create the output directory
    outputDir.mkdir(parents=True, exist_ok=True)

    # Step 1: Create the central list
    log.info("Creating folder list...")
    folderList = createFolderList(inputShp)
    # Group folders based on "name" attribute - with identical "name" attributes before the first underscore and update the list
    folderListGrouped = groupFoldersByName(folderList)
    log.info("Finished creating folder list")

    # Step 2: Set up ava directories
    log.info("Running folder initialization for each entry...")
    for entry in folderListGrouped:
        folderName = entry['folderName']
        initializeFolderStruct(str(outputDir / folderName), removeExisting=True)
        log.debug(f"Created folder structure for '{folderName}'.")
    log.info("Finished folder initialization")

    # Step 3: Split and move release areas to each directory
    log.info("Splitting and moving release areas...")
    splitAndMoveReleaseAreas(folderListGrouped, inputShp, outputDir)
    log.info("Finished splitting and moving release areas")

    # Step 4: Clip and move DEM
    log.info("Clipping and moving DEM...")
    demExtent = clipDEMByCentroidAndMove(folderListGrouped, inputDEM, outputDir, cfg)
    log.info("Finished clipping and moving of DEM")

    # Step 5: Clip and move ENT and RES
    log.info("Clipping and moving ENT and RES...")
    clipAndMoveEntRes(inputDir, outputDir, demExtent)
    log.info("Finished clipping and moving ENT and RES")

    # Step 6: Divide release areas into scenarios based on "scenario" attribute
    log.info("Separating release area by scenario attribute...")
    splitByScenarios(folderListGrouped, outputDir)
    log.info("Finished separating by scenarios")

def readShapefile(inputShp):
    """Read the fields, properties, geometries, and spatial reference of an input shapefile.
    To be used in combination with shapefile.Reader. Could be expanded upon to get e.g.
    shapeTypes, bounds, numFeatures and metadata if needed

    # ToDo: maybe move to some other module since its generally useful, e.g. in1Data, in2Trans -> shapeUtils.py

    Parameters
    ----------
    inputShp: pathlib.Path object
        the input shapefile

    Returns
    -------
    fields: list
        the fields of the input shapefile
    fieldNames: list
        the names of the fields
    properties: list
        a list of dictionaries containing the properties of each feature
    geometries: list
        a list of geometry objects
    srs: str
        the spatial reference system fetched from eventual .prj file
    """
    with shapefile.Reader(str(inputShp)) as src:
        fields = src.fields[1:]  # Skip deletion flag
        fieldNames = [field[0] for field in fields]
        properties = []
        geometries = []
        for shapeRecord in src.iterShapeRecords():
            properties.append(dict(zip(fieldNames, shapeRecord.record)))
            geometries.append(shape(shapeRecord.shape.__geo_interface__))

        srs = None
        # Check if .prj file exists and read it
        srsfile = inputShp.with_suffix('.prj')
        if srsfile.is_file():
            with open(srsfile, 'r') as f:
                srs = f.read().strip()
            log.debug(f"Found and read .prj file: {srsfile}")
        else:
            log.debug(f"No .prj file found at: {srsfile}")

    return fields, fieldNames, properties, geometries, srs

def writeShapefile(outputPath, fields, fieldNames, features, srs=None):
    """Write features to a shapefile with given fields and properties.

    Parameters
    ----------
    outputPath: pathlib.Path object
        path where the shapefile will be written
    fields: list
        the fields of the shapefile
    fieldNames: list
        the names of the fields
    features: list
        list of tuples containing (properties, geometry) for each feature
    srs: str, optional
        the spatial reference system for the .prj file

    Returns
    -------
    none
    """
    with shapefile.Writer(str(outputPath)) as dst:
        for field in fields:
            dst.field(*field)
        for properties, geometry in features:
            dst.shape(geometry)
            record = [properties.get(fieldName, '') for fieldName in fieldNames]
            dst.record(*record)
    
    if srs is not None:
        prjOutPath = outputPath.with_suffix('.prj')
        with open(prjOutPath, 'w') as prjFile:
            prjFile.write(srs)
        log.debug(f"Wrote projection file to '{prjOutPath}'")
    
    log.debug(f"Saved shapefile to '{outputPath}'.")

def createFolderList(inputShp):
    """Create a list of entries from each feature in the input shapefile.

    Parameters
    ----------
    inputShp: pathlib.Path object
        path to input shapefile

    Returns
    -------
    folderList: list
        list of dictionaries containing folderName, properties, and geometry for each feature
    """
    folderList = []
    unnamedCount = 1
    fields, fieldNames, properties, geometries, srs = readShapefile(inputShp)

    for i, (properties, geometry) in enumerate(zip(properties, geometries)):
        folderName = properties.get('name', '').strip() or f"unnamedAvalanche{str(unnamedCount).zfill(5)}"
        if not properties.get('name', '').strip():
            unnamedCount += 1

        folderList.append({
            'folderName': folderName,
            'properties': properties,
            'geometry': geometry,
        })

    log.info(f"Created initial folder list with '{len(folderList)}' entries.")
    return folderList

def groupFoldersByName(folderList):
    """Group entries in the folderList by name before the first underscore.

    Parameters
    ----------
    folderList: list
        list of dictionaries containing folderName, properties, and geometry for each feature

    Returns
    -------
    folderListByName: list
        list of dictionaries with entries grouped by name, containing folderName, properties list, and geometries list
    """
    folderListByName = []

    for entry in folderList:
        name = entry['folderName'].split('_')[0] # Get release area name before first underscore
        if not any(e['folderName'] == name for e in folderListByName):
            folderListByName.append({
                'folderName': name,
                'properties': [entry['properties']],
                'geometries': [entry['geometry']]
            })
        else:
            for e in folderListByName:
                if e['folderName'] == name:
                    e['properties'].append(entry['properties'])
                    e['geometries'].append(entry['geometry'])

    log.info(f"Grouped '{len(folderList)}' avalanche directories with identical names before first underscore. "
             f"Updated folder list contains '{len(folderListByName)}' entries")
    return folderListByName

def splitAndMoveReleaseAreas(folderList, inputShp, outputDir):
    """Split release areas into individual shapefiles and write them to their respective folders.

    Parameters
    ----------
    folderList: list
        list of dictionaries containing folderName, properties list, and geometries list
    inputShp: pathlib.Path object
        path to input shapefile
    outputDir: pathlib.Path object
        path to output directory where folders will be created

    Returns
    -------
    none
    """
    # Read the input shapefile
    fields, fieldNames, properties, geometries, srs = readShapefile(inputShp)

    featuresByName = {}
    for entry in folderList:
        name = entry['folderName'].split('_')[0] # Get release area name before first underscore
        # Group entries with the same name
        if name not in featuresByName:
            featuresByName[name] = []
        # add corresponding properties and geometries
        for i, properties in enumerate(entry['properties']):
            featuresByName[name].append((properties, entry['geometries'][i]))

    # Write shapefiles to their respective folders
    for name, features in featuresByName.items():
        shpOutPath = outputDir / name / "Inputs" / "REL" / name
        writeShapefile(shpOutPath, fields, fieldNames, features, srs)
        log.debug(f"Saved release area to '{shpOutPath}'.")

def clipDEMByCentroidAndMove(folderList, inputDEM, outputDir, cfg):
    """Clip the DEM around each release scenario's centerpoint and move to respective folders.
    
    Currently splits around centerpoint of the first feature in a group. ToDo: improve splitting logic

    Parameters
    ----------
    folderList: list
        list of dictionaries containing folderName, properties list, and geometries list
    inputDEM: pathlib.Path object
        path to input DEM file (.asc format)
    outputDir: pathlib.Path object
        path to output directory where clipped DEMs will be saved
    cfg: dict
        configuration settings containing GENERAL.bufferSize for clipping extent

    Returns
    -------
    tuple
        (xmin, ymin, xmax, ymax) representing the DEM extent after clipping
    """
    # Read input DEM
    demData = ascUtils.readRaster(inputDEM)
    header = demData['header']
    raster = demData['rasterData']

    cellsize = header["cellsize"]
    xOrigin = header["xllcenter"]
    yOrigin = header["yllcenter"]

    # Initialize extent variables
    xMin, yMin = float('inf'), float('inf')
    xMax, yMax = float('-inf'), float('-inf')

    # Find centerpoint
    for entry in folderList:
        folderName = entry['folderName']
        for properties in entry['properties']:
            geometry = entry['geometries'][0] #first geometry in list
            center = geometry.centroid

            # Calculate bounding box indices for clipping
            bufferSize = float(cfg['GENERAL']['bufferSize'])
            currXMin = center.x - bufferSize
            currYMin = center.y - bufferSize
            currXMax = center.x + bufferSize
            currYMax = center.y + bufferSize

            # Update overall extent
            xMin = min(xMin, currXMin)
            yMin = min(yMin, currYMin)
            xMax = max(xMax, currXMax)
            yMax = max(yMax, currYMax)

        # Convert bounding box to grid indices
        colStart = max(0, int((currXMin - xOrigin) / cellsize))
        colEnd = min(header["ncols"], int((currXMax - xOrigin) / cellsize))
        # Calculate row indices
        rowStart = max(0, int((yOrigin + header["nrows"] * cellsize - currYMax) / cellsize))
        rowEnd = min(header["nrows"], int((yOrigin + header["nrows"] * cellsize - currYMin) / cellsize))

        # Flip rows for bottom-left origin
        rowStart, rowEnd = header["nrows"] - rowEnd, header["nrows"] - rowStart

        # Validate grid indices
        if colStart >= colEnd or rowStart >= rowEnd: # e.g. if the bounding box falls completely outside the DEM
            log.warning(f"Invalid clipping bounds for '{folderName}'. Skipping.")
            continue

        # Clip DEM section
        clippedDEM = raster[rowStart:rowEnd, colStart:colEnd]
        log.debug(f"Clipped DEM for '{folderName}' with '{len(clippedDEM)}' rows and '{len(clippedDEM[0])}' columns")

        # Prepare header for clipped DEM #ToDo: needs adjustment once we include tif files
        clippedDEMHeader = {
            "ncols": colEnd - colStart,
            "nrows": rowEnd - rowStart,
            "xllcenter": xOrigin + colStart * cellsize,
            "yllcenter": yOrigin + rowStart * cellsize,
            "cellsize": cellsize,
            "nodata_value": header["nodata_value"],
        }

        # Define output path
        demOutPath = outputDir / folderName / "Inputs" / f"{folderName}_DEM.asc"

        # Save clipped DEM
        try:
            ascUtils.writeResultToAsc(clippedDEMHeader, clippedDEM, demOutPath, flip=True)
            log.debug(f"Clipped DEM saved to: {demOutPath}.")
        except Exception as e:
            log.error(f"Failed to save clipped DEM for: {folderName}: {e}")

    # Return the overall extent
    return (xMin, yMin, xMax, yMax)

def clipAndMoveEntRes(inputDir, outputDir, demExtent):
    """Clip and move ENT and RES files based on DEM extent.
    
    Parameters
    ----------
    inputDir : pathlib.Path object
        Path to input directory containing ENT and RES folders
    outputDir : pathlib.Path object
        Path to output directory where clipped files will be saved
    demExtent : tuple
        Overall DEM extent as (xmin, ymin, xmax, ymax), used only for logging

    Returns
    -------
    none
    """
    # Process both ENT and RES directories
    for dirType in ['ENT', 'RES']:
        typeDir = inputDir / dirType
        if not typeDir.exists():
            log.warning(f"No {dirType} directory found in {inputDir}")
            continue

        # Find shapefile in directory
        shpFile = next(typeDir.glob("*.shp"), None)
        if not shpFile:
            log.warning(f"No shapefile found in {typeDir}")
            continue

        # Read shapefile
        fields, fieldNames, properties, geometries, srs = readShapefile(shpFile)

        # Process each output directory (avalanche scenario)
        for entry in outputDir.iterdir():
            if not entry.is_dir():
                continue

            # Read the clipped DEM for this scenario to get its extent
            demFile = entry / 'Inputs' / f"{entry.name}_DEM.asc"
            if not demFile.exists():
                log.warning(f"No DEM file found for {entry.name}, skipping {dirType} processing")
                continue

            # Get DEM extent for this scenario
            demData = ascUtils.readRaster(demFile)
            header = demData['header']
            xMin = header['xllcenter']
            yMin = header['yllcenter']
            xMax = xMin + header['cellsize'] * header['ncols']
            yMax = yMin + header['cellsize'] * header['nrows']
            scenarioBbox = box(xMin, yMin, xMax, yMax)

            # Clip geometries with scenario's DEM extent
            clippedFeatures = []
            for prop, geom in zip(properties, geometries):
                if geom.intersects(scenarioBbox):
                    clippedGeom = geom.intersection(scenarioBbox)
                    if not clippedGeom.is_empty:
                        clippedFeatures.append((prop, clippedGeom))

            if not clippedFeatures:
                log.debug(f"No {dirType} features intersect with DEM extent for {entry.name}")
                continue

            # Create output directory and save clipped shapefile
            targetDir = entry / 'Inputs' / dirType
            targetDir.mkdir(parents=True, exist_ok=True)
            
            outputPath = targetDir / f"{entry.name}_{dirType}.shp"
            writeShapefile(outputPath, fields, fieldNames, clippedFeatures, srs)
            log.debug(f"Clipped {dirType} shapefile saved to: {outputPath}")

def splitByScenarios(folderList, outputDir):
    """Split release areas into separate shapefiles based on their scenario attribute.

    Parameters
    ----------
    folderList: list
        list of dictionaries containing folderName, properties list, and geometries list
    outputDir: pathlib.Path object
        path to output directory where scenario shapefiles will be created

    Returns
    -------
    none

    Notes
    -----
    - If a feature has no scenario attribute or it's empty, it will be marked as 'NULL' and grouped together with other 'NULL' features
    - Intermediate shapefiles are deleted or renamed after scenario splitting
    """
    totalInputFiles = 0
    totalScenarioFiles = 0

    # Loop through each folder
    for folder in folderList:
        inputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / folder['folderName']
        fields, fieldNames, properties, geometries, srs = readShapefile(inputShp)
        totalInputFiles += 1
        
        # Get the scenario attribute values
        if 'scenario' in fieldNames: #Check if scenario attribute exists
            # Create scenario dictionary
            scenarios = {}
            for shapeRecord in shapefile.Reader(str(inputShp)).iterShapeRecords():
                properties = dict(zip(fieldNames, shapeRecord.record))
                scenarioValues = properties.get('scenario', '').split(',')
                for scenario in scenarioValues:
                    # Check if scenario value is empty and set flag
                    if scenario.strip() == '':
                        scenario = 'NULL'
                    # If scenario is not in scenarios dict, add it
                    if scenario not in scenarios:
                        scenarios[scenario] = []
                    scenarios[scenario].append(shapeRecord)
            
            # Write the scenario shapefiles
            for scenario, records in scenarios.items():
                if all(scenario == 'NULL' for scenario in scenarios):
                    outputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / f"{folder['folderName']}_REL"
                elif scenario == 'NULL':
                    outputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / f"{folder['folderName']}_NULL"
                else:
                    outputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / f"{folder['folderName']}_{scenario}"
                
                # Filter out the scenario attribute and remove it
                shapeFeatures = [(dict(zip(fieldNames, record.record)), record.shape) for record in records]
                filteredFields = [field for field in fields if field[0] != 'scenario']
                filteredFieldNames = [name for name in fieldNames if name != 'scenario']
                
                # Write the shapefile
                writeShapefile(outputShp, filteredFields, filteredFieldNames, shapeFeatures, srs)
                totalScenarioFiles += 1

            # Delete the intermediate shapefile
            for ext in ['.shp', '.shx', '.dbf', '.prj']:
                if (inputShp.with_suffix(ext)).exists():
                    (inputShp.with_suffix(ext)).unlink()
        else:
            # If no scenario attribute exists, just rename the file
            outputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / f"{folder['folderName']}_REL"
            shapeFeatures = [(dict(zip(fieldNames, record.record)), record.shape) 
                           for record in shapefile.Reader(str(inputShp)).iterShapeRecords()]
            writeShapefile(outputShp, fields, fieldNames, shapeFeatures, srs)
            totalScenarioFiles += 1
            
            # Delete the intermediate shapefile
            for ext in ['.shp', '.shx', '.dbf', '.prj']:
                if (inputShp.with_suffix(ext)).exists():
                    (inputShp.with_suffix(ext)).unlink()

    log.info(f"Split '{totalInputFiles}' release area shapefiles into '{totalScenarioFiles}' scenarios")
