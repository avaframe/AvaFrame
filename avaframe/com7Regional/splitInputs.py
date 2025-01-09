"""Module for splitting and organizing regional avalanche input data."""

import logging
import shapefile  # pyshp
from shapely.geometry import shape, box, MultiPolygon
import pathlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import matplotlib as mpl
import time

from avaframe.out3Plot import plotUtils as pU
from avaframe.in2Trans import ascUtils
from avaframe.in3Utils.initializeProject import initializeFolderStruct

# create local logger
log = logging.getLogger(__name__)

def splitInputsMain(inputDir, outputDir, cfg, cfgMain):
    """Process and organize avalanche input data into individual avalanche directories based 
    on release area's "group" and "scenario" attributes.

    Parameters
    ----------
    inputDir : pathlib.Path object
        Path to input directory containing release areas (REL) and DEM files
    outputDir : pathlib.Path object
        Path to output directory where organized folders will be created
    cfg : dict
        Configuration settings containing:
        - GENERAL.bufferSize : float, buffer size for DEM clipping
    cfgMain : dict
        Configuration settings containing:
        - FLAGS.createReport : bool, whether to write report
        - FLAGS.savePlot : bool, whether to save plots

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

    # Step 1: Create the directory list
    log.info("Creating folder list based on REL 'group' attribute...")
    dirListGrouped = createDirList(inputShp)
    log.info("Finished creating folder list")

    # Step 2: Set up avalanche directories
    log.info("Initializing folder structure for each entry...")
    for entry in dirListGrouped:
        dirName = entry['dirName']
        initializeFolderStruct(str(outputDir / dirName), removeExisting=True)
        log.debug(f"Created folder structure for '{dirName}'.")
    log.info("Finished folder initialization")

    # Step 3: Split and move release areas to each directory
    log.info("Splitting and moving release areas...")
    splitAndMoveReleaseAreas(dirListGrouped, inputShp, outputDir)
    log.info("Finished splitting and moving release areas")

    # Step 4: Clip and move DEM
    log.info("Clipping and moving DEM...")
    groupExtents = clipDEMByReleaseGroup(dirListGrouped, inputDEM, outputDir, cfg)
    log.info("Finished clipping and moving of DEM")

    # Step 5: Clip and move optional input (currently only ENT and RES)
    log.info("Clipping and moving optional input...")
    groupFeatures = clipAndMoveOptionalInput(inputDir, outputDir, groupExtents)
    log.info("Finished clipping and moving optional input")

    # Step 6: Divide release areas into scenarios
    log.info("Separating release areas by scenarios...")
    splitByScenarios(dirListGrouped, outputDir)
    log.info("Finished separating by scenarios")

    # Step 7: Write reports
    if cfgMain['FLAGS'].getboolean('createReport'):
        log.info("Writing reports...")
        writeScenarioReport(dirListGrouped, outputDir)
        if cfgMain['FLAGS'].getboolean('savePlot'):
            writeVisualReport(dirListGrouped, inputDEM, outputDir, groupExtents, groupFeatures)
        log.info("Finished writing reports")

def readShapefile(inputShp):
    """Read the fields, properties, geometries, and spatial reference of an input shapefile.
    To be used in combination with shapefile.Reader. Could be expanded upon to get e.g.
    shapeTypes, bounds, numFeatures and metadata if needed

    # ToDo: maybe move to some other module e.g. in1Data, in2Trans -> shapeUtils.py or update to use pre-existing function from shpConversion.py

    Parameters
    ----------
    inputShp : pathlib.Path
        Path to input shapefile

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
            if shapeRecord.shape.shapeType == shapefile.NULL:
                log.warning(f"Skipping NULL shape in {inputShp}")
                continue
            properties.append(dict(zip(fieldNames, shapeRecord.record)))
            geometries.append(shape(shapeRecord.shape.__geo_interface__))

        srs = None
        # Check if .prj file exists and read it
        srsfile = inputShp.with_suffix('.prj')
        if srsfile.is_file():
            with open(srsfile, 'r') as f:
                srs = f.read().strip()
            
    return fields, fieldNames, properties, geometries, srs

def writeShapefile(outputPath, fields, fieldNames, features, srs=None):
    """Write features to a shapefile with given fields and properties.

    Parameters
    ----------
    outputPath: pathlib.Path
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

def createDirList(inputShp):
    """Create a list of entries from each feature in the input shapefile, grouped by the 'group' attribute.

    Parameters
    ----------
    inputShp: pathlib.Path object
        path to input shapefile

    Returns
    -------
    dirListGrouped: list
        list of dictionaries containing dirName (group name), properties list, and geometries list,
        where features are grouped by their 'group' attribute
    """
    fields, fieldNames, properties, geometries, srs = readShapefile(inputShp)
    
    # Create dictionary to store groups
    groups = {}
    unnamedCount = 1
    
    for props, geom in zip(properties, geometries):
        propsLower = {key.lower(): value for key, value in props.items()} # Handle case sensitivity
        
        # Get group name from 'group' attribute, fallback to unnamed if not present
        groupName = propsLower.get('group', '').strip() or f"{str(unnamedCount).zfill(5)}"
        if not propsLower.get('group', '').strip():
            unnamedCount += 1
            log.info(f"No 'group' field or empty group found in {inputShp}, using '{groupName}'")
        
        # Initialize group if not exists
        if groupName not in groups:
            groups[groupName] = {
                'dirName': groupName,
                'properties': [],
                'geometries': []
            }
        
        # Add feature to group
        groups[groupName]['properties'].append(props)
        groups[groupName]['geometries'].append(geom)
    
    # Convert dictionary to list and sort by dirName
    dirListGrouped = list(groups.values())
    dirListGrouped.sort(key=lambda x: x['dirName'].lower())

    # Log total number of features
    totalFeatures = sum(len(group['geometries']) for group in dirListGrouped)
    log.info(f"Found '{totalFeatures}' features that were organized into '{len(dirListGrouped)}' groups")
    
    return dirListGrouped

def splitAndMoveReleaseAreas(dirList, inputShp, outputDir):
    """Split release areas into individual shapefiles and write them to their respective folders.

    Parameters
    ----------
    dirList: list
        list of dictionaries containing dirName, properties list, and geometries list
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
    for entry in dirList:
        name = entry['dirName'] # Get release area name
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

def checkFeatureIsolation(geometries, properties, bufferSize, groupName):
    """Check if any feature in the group is isolated from all others.
    
    A feature is considered isolated if its buffered bounding box does not overlap
    with any other feature's buffered bounding box in the group.
    
    Parameters
    ----------
    geometries: list
        List of geometry objects to check
    properties: list
        List of dictionaries containing properties for each geometry
    bufferSize: float
        Buffer size to use when creating bounding boxes
    groupName: str
        Name of the group, used for error messages
        
    Raises
    ------
    ValueError
        If any feature is isolated from all others in the group
    """
    # Skip check if only one feature
    if len(geometries) <= 1:
        log.debug(f"Group '{groupName}' has only one feature, proceeding without isolation check.")
        return

    # Create buffered bounding boxes for each feature
    boundingBoxes = []
    for geom in geometries:
        center = geom.centroid
        
        # Calculate bounding box for this feature
        currXMin = center.x - bufferSize
        currYMin = center.y - bufferSize
        currXMax = center.x + bufferSize
        currYMax = center.y + bufferSize
        
        # Update group extent
        boundingBoxes.append(box(currXMin, currYMin, currXMax, currYMax))

    # Check each feature's bounding box against all others
    for i, bbox in enumerate(boundingBoxes):
        hasOverlap = False
        for j, otherBbox in enumerate(boundingBoxes):
            if i != j and bbox.intersects(otherBbox):
                hasOverlap = True
                break
        
        if not hasOverlap:
            # Find feature name regardless of case (NAME, name, Name etc.)
            featureProps = {key.lower(): value for key, value in properties[i].items()}
            featureName = featureProps.get('name', f'unnamed feature {i+1}').strip()
            
            message = f"Feature '{featureName}' in group '{groupName}' is isolated from all other features - consider splitting into a separate group"
            log.error(message)
            raise ValueError(message)

def clipDEMByReleaseGroup(dirList, inputDEM, outputDir, cfg):
    """Clip the DEM to include all features in each release group.

    Parameters
    ----------
    dirList : list
        List of dictionaries containing dirName, and geometries list
    inputDEM : pathlib.Path
        Path to input DEM file
    outputDir : pathlib.Path
        Path to output directory where DEM will be saved
    cfg : configparser object
        Configuration settings

    Returns
    -------
    dict
        Dictionary with dirName as key and (xMin, xMax, yMin, yMax) as value,
        containing the DEM extents for each group
    """    
    # Read input DEM
    demData = ascUtils.readRaster(inputDEM)
    header = demData['header']
    raster = demData['rasterData']
    cellSize = header['cellsize']
    xOrigin = header['xllcenter']
    yOrigin = header['yllcenter']
    nRows = header['nrows']
    nCols = header['ncols']

    # Process each group
    groupExtents = {}
    for entry in dirList:
        dirName = entry['dirName']
        geometries = entry['geometries']
        
        if not geometries:
            message = f"No geometries found for {dirName}"
            log.error(message)
            raise ValueError(message)
            
        # Get extent of all geometries in group
        bounds = [geom.bounds for geom in geometries]
        xMins, yMins, xMaxs, yMaxs = zip(*bounds)
        
        # Calculate extent with buffer
        bufferSize = cfg['GENERAL'].getfloat('bufferSize')
        xMin = min(xMins) - bufferSize
        xMax = max(xMaxs) + bufferSize
        yMin = min(yMins) - bufferSize
        yMax = max(yMaxs) + bufferSize
        groupExtents[dirName] = (xMin, xMax, yMin, yMax) # Store extent for this group
        
        # Convert extent to grid indices
        colStart = max(0, int((xMin - xOrigin) / cellSize))
        colEnd = min(nCols, int((xMax - xOrigin) / cellSize) + 1)
        
        # Convert y-coordinates to row indices (flipped for bottom-left origin)
        rowStart = max(0, int((yOrigin + nRows * cellSize - yMax) / cellSize))
        rowEnd = min(nRows, int((yOrigin + nRows * cellSize - yMin) / cellSize) + 1)
        
        # Flip row indices for bottom-left origin
        rowStart, rowEnd = nRows - rowEnd, nRows - rowStart
        
        # Clip the DEM data
        clippedData = raster[rowStart:rowEnd, colStart:colEnd]
        
        # Create header for clipped DEM
        clippedHeader = header.copy()
        clippedHeader['ncols'] = colEnd - colStart
        clippedHeader['nrows'] = rowEnd - rowStart
        clippedHeader['xllcenter'] = xOrigin + colStart * cellSize
        clippedHeader['yllcenter'] = yOrigin + rowStart * cellSize
        
        # Write clipped DEM
        outputDEM = outputDir / dirName / 'Inputs' / f"{dirName}_DEM.asc"
        ascUtils.writeResultToAsc(clippedHeader, clippedData, outputDEM, flip=True)
        log.debug(f"Clipped DEM saved to: {outputDEM}")

        # Store DEM extents (reduced by one pixel on each side to ensure DEM > clip extents)
        xMinDEM = clippedHeader['xllcenter'] - (cellSize/2) + cellSize
        yMinDEM = clippedHeader['yllcenter'] - (cellSize/2) + cellSize
        xMaxDEM = clippedHeader['xllcenter'] + (clippedHeader['ncols'] * cellSize) - (cellSize/2) - cellSize
        yMaxDEM = clippedHeader['yllcenter'] + (clippedHeader['nrows'] * cellSize) - (cellSize/2) - cellSize
        groupExtents[dirName] = (xMinDEM, xMaxDEM, yMinDEM, yMaxDEM)
    
    return groupExtents

def clipAndMoveOptionalInput(inputDir, outputDir, groupExtents):
    """Clip and move ENT and RES files based on group DEM extent.

    #ToDo: extend to include other input types
    
    Parameters
    ----------
    inputDir : pathlib.Path
        Path to input directory containing ENT and RES folders
    outputDir : pathlib.Path
        Path to output directory where clipped files will be saved
    groupExtents : dict
        Dictionary with dirName as key and (xMin, xMax, yMin, yMax) as value,
        containing the DEM clipping extents for each group

    Returns
    -------
    dict
        Dictionary containing clipped features for each group and type
        {dirName: {'ENT': [...], 'RES': [...]}}
    """
    groupFeatures = {}
    # Process ENT and RES directories
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

        # Process each output directory that has extents
        for entry in outputDir.iterdir():
            if not entry.is_dir() or entry.name not in groupExtents:
                continue

            # Get extent
            xMin, xMax, yMin, yMax = groupExtents[entry.name]
            scenarioBbox = box(xMin, yMin, xMax, yMax)

            # Initialize group in dictionary if not exists
            if entry.name not in groupFeatures:
                groupFeatures[entry.name] = {'ENT': [], 'RES': []}

            # Clip geometries with groups DEM extent
            clippedFeatures = []
            for prop, geom in zip(properties, geometries):
                if geom.intersects(scenarioBbox):
                    clippedGeom = geom.intersection(scenarioBbox)
                    if not clippedGeom.is_empty:
                        clippedFeatures.append((prop, clippedGeom))
                        groupFeatures[entry.name][dirType].append(clippedGeom)

            if not clippedFeatures:
                log.debug(f"No {dirType} features intersect with DEM extent for {entry.name}")
                continue

            # Create output directory and save clipped shp
            targetDir = entry / 'Inputs' / dirType
            targetDir.mkdir(parents=True, exist_ok=True)
            outputPath = targetDir / f"{entry.name}_{dirType}.shp"
            writeShapefile(outputPath, fields, fieldNames, clippedFeatures, srs)
            log.debug(f"Clipped {dirType} shapefile saved to: {outputPath}")

    return groupFeatures

def getScenarioGroups(inputShp, fieldNames):
    """Group shapefile records by their scenario attribute.
    
    Parameters
    ----------
    inputShp : pathlib.Path
        Path to input shapefile
    fieldNames : list
        List of field names in the shapefile
        
    Returns
    -------
    dict
        Dictionary mapping scenario names to lists of shape records
    """
    scenarios = {}
    for shapeRecord in shapefile.Reader(str(inputShp)).iterShapeRecords():
        properties = {k.lower(): v for k, v in zip(fieldNames, shapeRecord.record)}
        scenarioValues = properties.get('scenario', '').split(',')
        for scenario in scenarioValues:
            # Check if scenario value is empty and set flag
            if scenario.strip() == '':
                scenario = 'NULL'
            # If scenario is not in scenarios dict, add it
            if scenario not in scenarios:
                scenarios[scenario] = []
            scenarios[scenario].append(shapeRecord)
    return scenarios

def writeScenarioShapefile(outputShp, records, fields, fieldNames, srs):
    """Write a shapefile for a specific scenario.
    
    Parameters
    ----------
    outputShp : pathlib.Path
        Path where to write the shapefile
    records : list
        List of shape records for this scenario
    fields : list
        List of field definitions
    fieldNames : list
        List of field names
    srs : str
        Spatial reference system string
    """
    # Filter out the scenario attribute
    shapeFeatures = [(dict(zip(fieldNames, record.record)), record.shape) for record in records]
    filteredFields = [field for field in fields if field[0].lower() != 'scenario']
    filteredFieldNames = [name for name in fieldNames if name.lower() != 'scenario']
    
    # Write the shapefile
    writeShapefile(outputShp, filteredFields, filteredFieldNames, shapeFeatures, srs)

def splitByScenarios(dirList, outputDir):
    """Split release areas into separate shapefiles based on their scenario attribute.

    Parameters
    ----------
    dirList: list
        list of dictionaries containing dirName and list of geometries
    outputDir: pathlib.Path object
        path to output directory

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
    for folder in dirList:
        inputShp = pathlib.Path(outputDir) / folder['dirName'] / "Inputs" / "REL" / folder['dirName']
        fields, fieldNames, properties, geometries, srs = readShapefile(inputShp)
        totalInputFiles += 1
        
        # Get the scenario attribute values
        if 'scenario' in map(str.lower, fieldNames):
            # Group records by scenario
            scenarios = getScenarioGroups(inputShp, fieldNames)
            
            # Write a shapefile for each scenario
            for scenario, records in scenarios.items():
                if all(scenario == 'NULL' for scenario in scenarios):
                    outputShp = pathlib.Path(outputDir) / folder['dirName'] / "Inputs" / "REL" / f"{folder['dirName']}_REL"
                elif scenario == 'NULL':
                    outputShp = pathlib.Path(outputDir) / folder['dirName'] / "Inputs" / "REL" / f"{folder['dirName']}_NULL"
                else:
                    outputShp = pathlib.Path(outputDir) / folder['dirName'] / "Inputs" / "REL" / f"{folder['dirName']}_{scenario}"
                
                writeScenarioShapefile(outputShp, records, fields, fieldNames, srs)
                totalScenarioFiles += 1

            # Delete the intermediate shapefile
            for ext in ['.shp', '.shx', '.dbf', '.prj']:
                if (inputShp.with_suffix(ext)).exists():
                    (inputShp.with_suffix(ext)).unlink()
        else:
            # If no scenario attribute exists, rename the file (necessary for further processing)
            outputShp = pathlib.Path(outputDir) / folder['dirName'] / "Inputs" / "REL" / f"{folder['dirName']}_REL"
            shapeFeatures = [(dict(zip(fieldNames, record.record)), record.shape) 
                           for record in shapefile.Reader(str(inputShp)).iterShapeRecords()]
            writeShapefile(outputShp, fields, fieldNames, shapeFeatures, srs)
            for ext in ['.shp', '.shx', '.dbf', '.prj']:
                if (inputShp.with_suffix(ext)).exists():
                    (inputShp.with_suffix(ext)).unlink()
            
    if totalScenarioFiles == 0:
        log.info("No 'scenario' attribute or only 'NULL' found in release area shapefiles, continuing")
    else:
        log.info(f"Split '{totalInputFiles}' release area shapefiles into '{totalScenarioFiles}' scenarios")

def writeScenarioReport(dirListGrouped, outputDir):
    """Create a report in txt format listing all scenarios and their associated features.
    
    Parameters
    ----------
    dirListGrouped : list
        list of dictionaries containing dirName and list of geometries
    outputDir : pathlib.Path
        Path to output directory where the report will be saved
    
    Returns
    -------
    none
    """
    reportPath = outputDir / 'splitInputs_scenarioReport.txt'
    
    with open(reportPath, 'w') as f:
        f.write("SCENARIO REPORT\n")
        f.write("==============\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Process each group and their scenarios
        for group in sorted(dirListGrouped, key=lambda x: x['dirName'].lower()):
            dirName = group['dirName']
            f.write(f"Group: {dirName}\n")
            f.write("-" * (len(dirName) + 7) + "\n\n")
            
            relPath = pathlib.Path(outputDir) / dirName / "Inputs" / "REL"
            scenarioFiles = sorted(relPath.glob(f"{dirName}_*.shp"), 
                                 key=lambda x: x.stem.split('_')[-1])
            
            if not scenarioFiles:
                f.write("No scenarios found\n\n")
                continue
            
            # Write release areas for each scenario
            for scenFile in scenarioFiles:
                fields, fieldNames, properties, geometries, _ = readShapefile(scenFile)
                scenName = scenFile.stem.split('_')[-1]
                
                f.write(f"Scenario: {scenName}\n")
                f.write(f"No. of release areas: {len(geometries)}\n")
                
                if 'name' in map(str.lower, fieldNames): # Handle case sensitivity
                    nameIdx = [i for i, name in enumerate(fieldNames) if name.lower() == 'name'][0]
                    with shapefile.Reader(str(scenFile)) as shp:
                        records = sorted(shp.records(), key=lambda x: x[nameIdx].lower())
                        for record in records:
                            f.write(f"- {record[nameIdx]}\n")
                else:
                    for i in range(len(geometries)):
                        f.write(f"- Release Area {i+1}\n")
                f.write("\n")
            
            # Write total entrainment and resistance areas for the group
            entPath = pathlib.Path(outputDir) / dirName / "Inputs" / "ENT"
            resPath = pathlib.Path(outputDir) / dirName / "Inputs" / "RES"
            
            entFiles = list(entPath.glob(f"{dirName}_*.shp"))
            if entFiles:
                totalEnt = sum(len(readShapefile(ef)[3]) for ef in entFiles)
                if totalEnt > 0:
                    f.write(f"No. of entrainment areas: {totalEnt}\n")
            
            resFiles = list(resPath.glob(f"{dirName}_*.shp"))
            if resFiles:
                totalRes = sum(len(readShapefile(rf)[3]) for rf in resFiles)
                if totalRes > 0:
                    f.write(f"No. of resistance areas: {totalRes}\n")
            
            f.write("\n")
    
    log.info(f"Scenario report written to: {reportPath}")

def getExteriorCoords(geom):
    """Get the exterior coordinates of a shapely geometry to handle both single and multi-polygon geometries.

    Parameters
    ----------
    geom : shapely.geometry
        The shapely geometry to get the exterior coordinates from.

    Returns
    -------
    list
        A list of tuples containing the x and y coordinates of the geometry exterior.
    """
    if isinstance(geom, MultiPolygon):
        return [poly.exterior.xy for poly in geom.geoms]
    else:
        return [geom.exterior.xy]

def makeVisualReport(dirListGrouped, inputDEM, outputDir, groupExtents, groupFeatures, reportType):
    """Create a visual report showing the DEM extent with either basic or optional inputs.
    
    Parameters
    ----------
    dirListGrouped : list
        List of dictionaries containing dirName and list of geometries
    inputDEM : pathlib.Path
        Path to input DEM file
    outputDir : pathlib.Path
        Path to output directory where the report will be saved
    groupExtents : dict
        Dictionary with dirName as key and (xMin, xMax, yMin, yMax) as value,
        containing the DEM clipping extents for each group
    groupFeatures : dict
        Dictionary containing clipped features for each group and type
    reportType : str
        Type of report to create, either 'basic' or 'optional'
        - 'basic': Shows DEM extent and release areas only
        - 'optional': Shows DEM extent with entrainment and resistance areas
    
    Returns
    -------
    pathlib.Path
        Path to the generated report image
    """
    # Set up figure
    plt.figure(figsize=(10, 8))
    ax = plt.subplot(1, 1, 1)
    
    # Read and plot DEM
    demData = ascUtils.readRaster(inputDEM)
    header = demData['header']
    cellSize = header['cellsize']
    xMin = header['xllcenter']
    yMin = header['yllcenter']
    xMax = xMin + cellSize * header['ncols']
    yMax = yMin + cellSize * header['nrows']
    im = ax.imshow(demData['rasterData'], extent=[xMin, xMax, yMin, yMax],
                  cmap=pU.cmapDEM.reversed(), alpha=1, origin='lower', zorder=1)
    
    # Custom color scheme for groups
    colors = ['#0EF8EA', '#FFA500', '#C71585', '#00FF00', '#FF4500', '#800080',
              '#ADFF2F', '#FF6347', '#8A2BE2', '#FFFF00', '#FF0000']

    # Plot groups
    for idx, group in enumerate(dirListGrouped):
        dirName = group['dirName']
        color = colors[idx % len(colors)]
        # Plot DEM extent using groupExtents
        if dirName in groupExtents:
            xMin, xMax, yMin, yMax = groupExtents[dirName]
            rect = Rectangle((xMin, yMin), xMax - xMin, yMax - yMin,
                           fill=False, linestyle='--', color=color)
            ax.add_patch(rect)
            
            if reportType == 'basic':
                # Plot release areas
                for geom in group['geometries']:
                    for x, y in getExteriorCoords(geom):
                        plt.fill(x, y, alpha=1.0, color=color)
            else:
                # Plot optional features
                if dirName in groupFeatures:
                    for geom in groupFeatures[dirName].get('ENT', []):
                        for x, y in getExteriorCoords(geom):
                            plt.fill(x, y, alpha=0.3, color=color, edgecolor='none')
                    for geom in groupFeatures[dirName].get('RES', []):
                        for x, y in getExteriorCoords(geom):
                            plt.fill(x, y, alpha=0.5, color=color, hatch='xxxx', fill=False, edgecolor=color, linewidth=0.5)
            
            # Place group label using groupExtents
            plt.text(xMin, yMax, dirName, color=color, fontsize=8,
                          transform=mpl.transforms.offset_copy(ax.transData, x=1, y=-7, units='points', fig=ax.figure))

    # Create legend
    mapElements = [Rectangle((0, 0), 1, 1, fill=False, linestyle='--', color='black',
                           label='Clipped DEM Extent')]
    if reportType == 'basic':
        mapElements.append(Patch(facecolor='black', alpha=1.0, label='Release Areas'))
    else:  # optional
        mapElements.extend([
            Patch(facecolor='black', alpha=0.3, label='Entrainment Areas'),
            Patch(facecolor='none', alpha=0.3, hatch='xxxx', label='Resistance Areas', edgecolor='black', linewidth=0.5)
        ])

    plt.legend(handles=mapElements,
                    title='Legend',
                    loc='center left',
                    bbox_to_anchor=(1, 0.5))
    
    # Add DEM colorbar
    cax = ax.inset_axes([1.015, 0, 0.375, 0.02])  # [x, y, width, height]
    plt.colorbar(im, cax=cax, orientation='horizontal', label='Elevation [m]')
    
    # Format plot; add title and labels
    ax.set_aspect('equal')
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useOffset=False))
    regionalDir = inputDEM.parent.parent.name
    reportName = 'Basic' if reportType == 'basic' else 'Optional'
    plt.title(f'Split Inputs Report - {reportName} - {regionalDir}')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Set axis limits based on DEM extents with a small margin
    xMins = [ext[0] for ext in groupExtents.values()]
    xMaxs = [ext[1] for ext in groupExtents.values()]
    yMins = [ext[2] for ext in groupExtents.values()]
    yMaxs = [ext[3] for ext in groupExtents.values()]
    xMin, xMax = min(xMins), max(xMaxs)
    yMin, yMax = min(yMins), max(yMaxs)
    margin = 0.01
    dx = (xMax - xMin) * margin
    dy = (yMax - yMin) * margin
    ax.set_xlim(xMin - dx, xMax + dx)
    ax.set_ylim(yMin - dy, yMax + dy)
    
    reportPath = outputDir / f'splitInputs_visualReport_{reportType}.png'
    plt.savefig(reportPath, dpi=200, bbox_inches='tight')
    plt.close()
    
    return reportPath

def writeVisualReport(dirListGrouped, inputDEM, outputDir, groupExtents, groupFeatures):
    """Write visual reports summarizing the split inputs operation.
    
    Creates two visual reports in PNG format:
    1. Basic report showing DEM extent and release areas
    2. Optional features report showing DEM extent with entrainment and resistance areas
    
    Parameters
    ----------
    dirListGrouped : list
        List of dictionaries containing dirName and list of geometries
    inputDEM : pathlib.Path
        Path to input DEM file
    outputDir : pathlib.Path
        Path to output directory where reports will be saved
    groupExtents : dict
        Dictionary with dirName as key and (xMin, xMax, yMin, yMax) as value,
        containing the DEM clipping extents for each group
    groupFeatures : dict
        Dictionary containing clipped features for each group and type
    
    Returns
    -------
    none
    """
    # Create basic features report
    basicPath = makeVisualReport(dirListGrouped, inputDEM, outputDir, 
                                 groupExtents, groupFeatures, 'basic')
    log.info(f"Visual report (basic) written to: {basicPath}")
    
    # Create optional features report
    optionalPath = makeVisualReport(dirListGrouped, inputDEM, outputDir, 
                                    groupExtents, groupFeatures, 'optional')
    log.info(f"Visual report (optional) written to: {optionalPath}")
