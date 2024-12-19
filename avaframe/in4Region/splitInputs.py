"""Helper functions for splitting avalanche input data into individual folders"""
import logging
import shapefile  # pyshp
from shapely.geometry import shape, Point, box, MultiPolygon
from avaframe.in2Trans import ascUtils
import pathlib

from avaframe.in3Utils.initializeProject import initializeFolderStruct

# create local logger
log = logging.getLogger(__name__)

def splitInputsMain(inputDir, outputDir, cfg):
    """Main function for splitting inputs

    Parameters
    ----------
    inputDir: pathlib.Path object
        path to input directory
    outputDir: pathlib.Path object
        path to output directory
    cfg: dict
        configuration settings

    Returns
    -------
    none
    """
    # Fetch the necessary input
    inputRELDir = inputDir / 'REL'
    inputShp = next(inputRELDir.glob("*.shp"), None)
    if not inputShp:
        log.error(f"No shapefile found in {inputRELDir}.")
        return
    inputDEM = next(inputDir.glob("*.asc"), None)
    if not inputDEM:
        log.error(f"No DEM file found in {inputDir}.")
        return

    # Create the output directory
    outputDir.mkdir(parents=True, exist_ok=True)

    # Step 1: Create the central list
    log.info("Creating folder list...")
    folderList = createFolderList(inputShp)
    # Group folders with identical "name" attributes before the first underscore and update the list
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
    if cfg['GENERAL'].getboolean('splitDEM'):
        log.info("Clipping and moving DEM...")
        clipDEMByCentroidAndMove(folderListGrouped, inputDEM, outputDir, cfg)
        log.info("Finished clipping and moving of DEM")

    # (Step 5: Clip and move ENT and RES... other input too?)

    # Step 6: Divide release areas into scenarios based on "scenario" attribute, if it exists
    log.info("Separating release area by scenario attribute...")
    splitByScenarios(folderListGrouped, outputDir)
    log.info("Finished separating by scenarios")

def readShapefile(src):
    """Read the fields, properties, and geometries of an input shapefile.
    To be used in combination with shapefile.Reader. Could be expanded upon to get e.g.
    the srs, shapeTypes, bounds, numFeatures and metadata if needed

    Parameters
    ----------
    src: shapefile.Reader
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
    """

    fields = src.fields[1:]  # Skip deletion flag
    fieldNames = [field[0] for field in fields]
    properties = []
    geometries = []
    for shapeRecord in src.iterShapeRecords():
        properties.append(dict(zip(fieldNames, shapeRecord.record)))
        geometries.append(shape(shapeRecord.shape.__geo_interface__))
    return fields, fieldNames, properties, geometries

def createFolderList(inputShp):
    """Create a list of entries of each feature in the input shapefile."""
    folderList = []
    unnamedCount = 1
    with shapefile.Reader(inputShp) as src:
        fields, fieldNames, properties, geometries = readShapefile(src)

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
    """Group entries in the folderList by name. Return an updated folderList"""
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
    """Split release areas into individual shapefiles and write them to their respective folders."""
    # Read the input shapefile
    with shapefile.Reader(inputShp) as src:
        fields, fieldNames, properties, geometries = readShapefile(src)

        featuresByName = {}
        for entry in folderList:
            name = entry['folderName'].split('_')[0] # Get release area name before first underscore
            # Group entries with the same name
            if name not in featuresByName:
                featuresByName[name] = []
            # add corresponding properties and geometries
            for i, properties in enumerate(entry['properties']):
                featuresByName[name].append((properties, geometries[i]))

        # Write shapefiles to their respective folders
        for name, features in featuresByName.items():
            shpOutPath = outputDir / name / "Inputs" / "REL" / name
            with shapefile.Writer(str(shpOutPath)) as dst:
                for field in fields: #write fields
                    dst.field(*field)
                for i, feature in enumerate(features):
                    properties = feature[0]
                    geometry = [e['geometries'][i] for e in folderList if e['folderName'] == name][0]
                    dst.shape(geometry)
                    record = [properties.get(fieldName, '') for fieldName in fieldNames]
                    dst.record(*record)

            log.debug(f"Saved release area to '{shpOutPath}'.")

def clipDEMByCentroidAndMove(folderList, inputDEM, outputDir, cfg):
    #Todo: Currently splits around centerpoint of the first feature in a group - improve
    """Clip the DEM around each release scenario's centerpoint and move it to respective folders."""
    # Read input DEM
    demData = ascUtils.readRaster(inputDEM)
    header = demData['header']
    raster = demData['rasterData']

    cellsize = header["cellsize"]
    xOrigin = header["xllcenter"]
    yOrigin = header["yllcenter"]

    # Find centerpoint
    for entry in folderList:
        folderName = entry['folderName']
        for properties in entry['properties']:
            geometry = entry['geometries'][0] #first geometry in list
            center = geometry.centroid

            # Calculate bounding box indices for clipping
            bufferSize = float(cfg['GENERAL']['bufferSize'])
            xmin, ymin, xmax, ymax = (
                center.x - bufferSize,
                center.y - bufferSize,
                center.x + bufferSize,
                center.y + bufferSize,
            )

        # Convert bounding box to grid indices
        colStart = max(0, int((xmin - xOrigin) / cellsize))
        colEnd = min(header["ncols"], int((xmax - xOrigin) / cellsize))
        # Calculate row indices
        rowStart = max(0, int((yOrigin + header["nrows"] * cellsize - ymax) / cellsize))
        rowEnd = min(header["nrows"], int((yOrigin + header["nrows"] * cellsize - ymin) / cellsize))

        # Flip rows for bottom-left origin
        rowStart, rowEnd = header["nrows"] - rowEnd, header["nrows"] - rowStart

        # Validate grid indices
        if colStart >= colEnd or rowStart >= rowEnd: # e.g. if the bounding box falls completely outside the DEM
            log.warning(f"Invalid clipping bounds for '{folderName}'. Skipping.")
            continue

        # Clip DEM section
        clippedDEM = raster[rowStart:rowEnd, colStart:colEnd]
        log.debug(f"Clipped DEM for '{folderName}' with '{len(clippedDEM)}' rows and '{len(clippedDEM[0])}' columns")

        # Prepare header for clipped DEM
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
            log.info(f"Clipped DEM saved to '{demOutPath}'.")
        except Exception as e:
            log.error(f"Failed to save clipped DEM for '{folderName}': {e}")

    return demOutPath

def splitByScenarios(folderList, outputDir):
    """Use the scenario attribute to separate by scenarios"""
    # Read the split shapefiles
    for folder in folderList:
        inputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / folder['folderName']
        with shapefile.Reader(str(inputShp)) as src:
            fields, fieldNames, properties, geometries = readShapefile(src)
            # Get the scenario attribute values
            if 'scenario' in fieldNames: #Check if scenario attribute exists
                # Create scenario dictionary
                scenarios = {}
                for shapeRecord in src.iterShapeRecords():
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
                src.close()
                # Check if all scenarios are empty
                all_null = all(scenario == 'NULL' for scenario in scenarios)
                # Write the scenario shapefiles
                for scenario, records in scenarios.items():
                    if all_null:
                        outputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / f"{folder['folderName']}_REL"
                    elif scenario == 'NULL':
                        outputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / f"{folder['folderName']}_NULL"
                    else:
                        outputShp = pathlib.Path(outputDir) / folder['folderName'] / "Inputs" / "REL" / f"{folder['folderName']}_{scenario}"
                    with shapefile.Writer(str(outputShp)) as dst:
                        for field in fields:  # write fields
                            if field[0] != 'scenario':  # skip scenario attribute
                                dst.field(*field)
                        for record in records:  # write geometry
                            dst.shape(record.shape)
                            recordValues = [record.record[i] for i, field in enumerate(fieldNames) if
                                            field != 'scenario']
                            dst.record(*recordValues)
                for file in inputShp.parent.rglob(inputShp.stem + '.*'):  # Delete original shapefile
                        file.unlink()
            else:
                src.close()
                log.info(f"No 'scenario' attribute found in '{inputShp}'. Skipping.")
