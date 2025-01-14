"""Module for splitting and organizing regional avalanche input data."""

import logging
import shapefile  # pyshp
from shapely.geometry import shape, box, Polygon
from avaframe.in2Trans import ascUtils
import pathlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Patch
import matplotlib as mpl

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

    # Step 2: Set up avalanche directories
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
    groupExtents = clipDEMByReleaseGroup(folderListGrouped, inputDEM, outputDir, cfg)
    log.info("Finished clipping and moving of DEM")

    # Step 5: Clip and move optional input (currently only ENT and RES)
    log.info("Clipping and moving optional input...")
    groupFeatures = clipAndMoveOptionalInput(inputDir, outputDir, groupExtents)
    log.info("Finished clipping and moving optional input")

    # Step 6: Divide release areas into scenarios based on "scenario" attribute
    log.info("Separating release area by scenario attribute...")
    splitByScenarios(folderListGrouped, outputDir)
    log.info("Finished separating by scenarios")

    # Step 7: Write report
    log.info("Writing report...")
    writeReport(folderListGrouped, inputDEM, outputDir, groupExtents, groupFeatures)
    log.info("Finished writing report")

def readShapefile(inputShp):
    """Read the fields, properties, geometries, and spatial reference of an input shapefile.
    To be used in combination with shapefile.Reader. Could be expanded upon to get e.g.
    shapeTypes, bounds, numFeatures and metadata if needed

    # ToDo: maybe move to some other module since its generally useful, e.g. in1Data, in2Trans -> shapeUtils.py

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
        folderName = properties.get('name', '').strip() or f"noName{str(unnamedCount).zfill(5)}"
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

def clipDEMByReleaseGroup(folderList, inputDEM, outputDir, cfg):
    """Clip the DEM to include all features in each release group.

    Parameters
    ----------
    folderList : list
        List of dictionaries containing folderName, and geometries list
    inputDEM : pathlib.Path
        Path to input DEM file
    outputDir : pathlib.Path
        Path to output directory where DEM will be saved
    cfg : configparser object
        Configuration settings

    Returns
    -------
    dict
        Dictionary with folderName as key and (xMin, xMax, yMin, yMax) as value,
        containing the DEM clipping extents for each group
    """
    # Get buffer size from config
    bufferSize = cfg['GENERAL'].getfloat('bufferSize')
    log.info(f"Using buffer size of '{bufferSize}' m")
    
    groupExtents = {}

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
    for entry in folderList:
        folderName = entry['folderName']
        geometries = entry['geometries']

        if not geometries:
            message = f"No geometries found for {folderName}"
            log.error(message)
            raise ValueError(message)

        # Get extent of all geometries in group
        bounds = [geom.bounds for geom in geometries]
        xMins, yMins, xMaxs, yMaxs = zip(*bounds)
        
        # Calculate extent with buffer
        xMin = min(xMins) - bufferSize
        xMax = max(xMaxs) + bufferSize
        yMin = min(yMins) - bufferSize
        yMax = max(yMaxs) + bufferSize

        # Store extent for this group
        groupExtents[folderName] = (xMin, xMax, yMin, yMax)
        
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
        outputDEM = outputDir / folderName / 'Inputs' / f"{folderName}_DEM.asc"
        ascUtils.writeResultToAsc(clippedHeader, clippedData, outputDEM, flip=True)
        log.debug(f"Clipped DEM saved to: {outputDEM}")

    return groupExtents

def clipAndMoveOptionalInput(inputDir, outputDir, groupExtents):
    """Clip and move ENT and RES files based on group DEM extent.
    
    Parameters
    ----------
    inputDir : pathlib.Path
        Path to input directory containing ENT and RES folders
    outputDir : pathlib.Path
        Path to output directory where clipped files will be saved
    groupExtents : dict
        Dictionary with folderName as key and (xMin, xMax, yMin, yMax) as value,
        containing the DEM clipping extents for each group

    Returns
    -------
    dict
        Dictionary containing clipped features for each group and type
        {groupName: {'ENT': [...], 'RES': [...]}}
    """
    # Initialize dictionary to store clipped features by group
    groupFeatures = {}

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

        # Process each output directory
        for entry in outputDir.iterdir():
            if not entry.is_dir():
                continue

            # Get extent for this group
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

            # Create output directory and save clipped shapefile
            targetDir = entry / 'Inputs' / dirType
            targetDir.mkdir(parents=True, exist_ok=True)
            
            outputPath = targetDir / f"{entry.name}_{dirType}.shp"
            writeShapefile(outputPath, fields, fieldNames, clippedFeatures, srs)
            log.debug(f"Clipped {dirType} shapefile saved to: {outputPath}")

    return groupFeatures

def splitByScenarios(folderList, outputDir):
    """Split release areas into separate shapefiles based on their scenario attribute.

    Parameters
    ----------
    folderList: list
        list of dictionaries containing folderName and list of geometries
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

def writeReport(folderListGrouped, inputDEM, outputDir, groupExtents, groupFeatures):
    """Write a visual report summarizing the split inputs operation.

    Parameters
    ----------
    folderListGrouped : list
        list of dictionaries containing folderName and list of geometries
    inputDEM : pathlib.Path
        path to input DEM file
    outputDir : pathlib.Path
        Path to output directory
    groupExtents : dict
        Dictionary containing the clipping extents for each group as (xMin, xMax, yMin, yMax)
    groupFeatures : dict
        Dictionary containing clipped features for each group and type

    Returns
    -------
    none
    """
    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Read and plot DEM
    demData = ascUtils.readRaster(inputDEM)
    header = demData['header']
    cellSize = header['cellsize']
    xMin = header['xllcenter']
    yMin = header['yllcenter']
    xMax = xMin + cellSize * header['ncols']
    yMax = yMin + cellSize * header['nrows']

    log.debug(f"DEM header info: xMin={xMin}, xMax={xMax}, yMin={yMin}, yMax={yMax}")
    log.debug(f"DEM shape: {demData['rasterData'].shape}")
    log.debug(f"DEM min/max values: {np.nanmin(demData['rasterData']):.2f}, {np.nanmax(demData['rasterData']):.2f}")

    # Plot DEM with correct extent
    im = ax.imshow(demData['rasterData'], extent=[xMin, xMax, yMin, yMax],
                  cmap='gray', alpha=0.5, origin='lower')
    
    # Create legend elements for map features
    mapElements = [
        Rectangle((0, 0), 1, 1, fill=False, linestyle='--', color='black',
                 label='Clipping Extent'),
        Patch(facecolor='black', alpha=1.0, label='Release Areas'),
        Patch(facecolor='gray', alpha=0.3, hatch='///', label='Entrainment Areas'),
        Patch(facecolor='gray', alpha=0.3, label='Resistance Areas')
    ]

    # Create legend elements for groups
    groupElements = []

    # Create colormap of distinctive colors
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', 
             '#a65628', '#f781bf', '#66c2a5', '#fc8d62', '#8da0cb',
             '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3',
             '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99', '#b15928',
             '#67001f', '#d6604d', '#4393c3', '#053061', '#66c2a4',
             '#5aae61', '#9970ab', '#762a83', '#fee090', '#bf812d']

    # Process groups and plot them
    for idx, group in enumerate(folderListGrouped):
        folderName = group['folderName']
        geometries = group['geometries']
        color = colors[idx % len(colors)]  # Cycle through colors if more groups than colors
        
        # Get group extent
        xMin, xMax, yMin, yMax = groupExtents[folderName]
        log.debug(f"Group {folderName} extent: x[{xMin:.2f}, {xMax:.2f}], y[{yMin:.2f}, {yMax:.2f}]")
        
        # Plot release areas
        for geom in geometries:
            x, y = geom.exterior.xy
            poly = plt.fill(x, y, alpha=1.0, color=color)

        # Plot ENT areas if they exist
        if folderName in groupFeatures and groupFeatures[folderName]['ENT']:
            for geom in groupFeatures[folderName]['ENT']:
                x, y = geom.exterior.xy
                plt.fill(x, y, alpha=0.3, color=color, hatch='///')

        # Plot RES areas if they exist
        if folderName in groupFeatures and groupFeatures[folderName]['RES']:
            for geom in groupFeatures[folderName]['RES']:
                x, y = geom.exterior.xy
                plt.fill(x, y, alpha=0.3, color=color)

        # Plot clipping extent box
        rect = Rectangle((xMin, yMin), xMax - xMin, yMax - yMin,
                        fill=False, linestyle='--', color=color)
        ax.add_patch(rect)

        # Add group to legend
        groupElements.append(Patch(facecolor=color, alpha=1.0,
                                 label=folderName))

    # Set equal aspect ratio
    ax.set_aspect('equal')

    # Format axis
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useOffset=False))

    # Add title and labels
    regionalDir = inputDEM.parent.parent.name
    plt.title(f'Split Inputs Report - {regionalDir}')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    
    # Create two-part legend
    firstLegend = ax.legend(handles=mapElements, title='Map Features',
                          bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.add_artist(firstLegend)
    ax.legend(handles=groupElements, title='Groups',
             bbox_to_anchor=(1.05, 0.6), loc='upper left')

    plt.grid(True, linestyle='--', alpha=0.3)

    # Set axis limits using global extent from groupExtents
    globalExtent = (
        min(e[0] for e in groupExtents.values()),
        max(e[1] for e in groupExtents.values()),
        min(e[2] for e in groupExtents.values()),
        max(e[3] for e in groupExtents.values())
    )
    ax.set_xlim(globalExtent[0], globalExtent[1])
    ax.set_ylim(globalExtent[2], globalExtent[3])
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(outputDir / f'splitInputsReport_{regionalDir}.png', dpi=300, bbox_inches='tight')
    plt.close()
