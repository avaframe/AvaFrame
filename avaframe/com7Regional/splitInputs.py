"""Module for splitting and organizing regional avalanche input data."""

import logging
import shapefile  # pyshp
from shapely.geometry import box
import pathlib
import time

from avaframe.in2Trans import rasterUtils
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in1Data import getInput
from avaframe.in3Utils.initializeProject import initializeFolderStruct
from avaframe.in2Trans import shpConversion as shpConv
from avaframe.out3Plot.outCom7Regional import createReportPlot

# create local logger
log = logging.getLogger(__name__)


def splitInputsMain(avalancheDir, outputDir, cfg, cfgMain):
    """Process and organize avalanche input data into individual avalanche directories based
    on release area's "group" and "scenario" attributes provided in the release area file. If no
    "group" attribute is provided, one avalanche directory per feature will be created (scenario is
    ignored in this case).

    Parameters
    ----------
    avalancheDir : pathlib.Path object
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
    avalancheDir/
    └── Inputs/
        ├── REL/
        │   └── *.shp         # all release areas
        ├── ENT/              # all entrainment areas (optional)
        │   └── *.shp
        ├── RES/              # all resistance areas (optional)
        │   └── *.shp
        └── *.asc or *.tif    # digital elevation model (DEM)
    """
    # Fetch the necessary input
    inputSimFilesAll = getInput.getInputDataCom1DFA(avalancheDir)

    # extract release shapefile, make sure only one exists
    if len(inputSimFilesAll["relFiles"]) == 1:
        inputShp = inputSimFilesAll["relFiles"][0]
    else:
        log.error(f"Expected only one release area file, found {len(inputSimFilesAll['relFiles'])}.")
        return

    # Get the input DEM
    inputDEM = getInput.getDEMPath(avalancheDir)

    # Create the output directory
    fU.makeADir(outputDir)

    # Step 1: Create the directory list
    log.info("Initializing folder structure for each group...")
    dirListGrouped = createDirList(inputShp)
    log.info("Finished creating folder list")

    # Step 2: Set up avalanche directories
    log.info("Initializing folder structure for each entry...")
    for entry in dirListGrouped:
        dirName = entry["dirName"]
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
    groupFeatures = clipAndMoveOptionalInput(inputSimFilesAll, outputDir, groupExtents)
    log.info("Finished clipping and moving optional input")

    # Step 6: Divide release areas into scenarios
    log.info("Separating release areas by scenarios...")
    splitByScenarios(dirListGrouped, outputDir)
    log.info("Finished separating by scenarios")

    # Step 7: Write reports
    if cfgMain["FLAGS"].getboolean("createReport"):
        log.info("Writing reports...")
        writeScenarioReport(dirListGrouped, outputDir)
        if cfgMain["FLAGS"].getboolean("savePlot"):
            createReportPlots(dirListGrouped, inputDEM, outputDir, groupExtents, groupFeatures)
        log.info("Finished writing reports")


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
    fields, fieldNames, properties, geometries, srs = shpConv.readShapefile(inputShp)

    # Create dictionary to store groups
    groups = {}
    unnamedCount = 1

    for props, geom in zip(properties, geometries):
        propsLower = {key.lower(): value for key, value in props.items()}  # Handle case sensitivity

        # Get group name from 'group' attribute, fallback to unnamed if not present
        groupName = propsLower.get("group", "").strip() or f"{str(unnamedCount).zfill(5)}"
        if not propsLower.get("group", "").strip():
            unnamedCount += 1
            log.info(f"No 'group' field or empty group found in {inputShp}, using '{groupName}'")

        # Initialize group if not exists
        if groupName not in groups:
            groups[groupName] = {
                "dirName": groupName,
                "properties": [],
                "geometries": [],
            }

        # Add feature to group
        groups[groupName]["properties"].append(props)
        groups[groupName]["geometries"].append(geom)

    # Convert dictionary to list and sort by dirName
    dirListGrouped = list(groups.values())
    dirListGrouped.sort(key=lambda x: x["dirName"].lower())

    # Log total number of features
    totalFeatures = sum(len(group["geometries"]) for group in dirListGrouped)
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
    fields, fieldNames, properties, geometries, srs = shpConv.readShapefile(inputShp)

    featuresByName = {}
    for entry in dirList:
        name = entry["dirName"]  # Get release area name
        # Group entries with the same name
        if name not in featuresByName:
            featuresByName[name] = []
        # add corresponding properties and geometries
        for i, properties in enumerate(entry["properties"]):
            featuresByName[name].append((properties, entry["geometries"][i]))

    # Write shapefiles to their respective folders
    for name, features in featuresByName.items():
        shpOutPath = outputDir / name / "Inputs" / "REL" / name
        shpConv.writeShapefile(shpOutPath, fields, fieldNames, features, srs)
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
            featureName = featureProps.get("name", f"unnamed feature {i+1}").strip()

            message = f"Feature '{featureName}' in group '{groupName}' is isolated from all other features - consider assigning it to a different group"
            log.error(message)
            raise ValueError(message)


def clipDEMByReleaseGroup(dirList, inputDEM, outputDir, cfg):
    """Clip the DEM to include all features in each release group. Returns an error if any feature in a group is isolated.

    Parameters
    ----------
    dirList : list
        List of dictionaries containing dirName, and geometries list
    inputDEM : pathlib.Path
        Path to input DEM file
    outputDir : pathlib.Path
        Path to output directory where clipped DEMs will be saved
    cfg : configparser object
        Configuration settings containing:
            - GENERAL.bufferSize : float
                Size of buffer to add around release areas   Configuration settings

    Returns
    -------
    groupExtents : dict
        Dictionary with dirName as key and (xMin, xMax, yMin, yMax) as value.
        The extents are reduced by one pixel on each side to ensure DEM extents
        are larger than clip extents of other input.
    """
    # Read input DEM
    demData = rasterUtils.readRaster(inputDEM)
    header = demData["header"]
    raster = demData["rasterData"]
    cellSize = header["cellsize"]
    xOrigin = header["xllcenter"]
    yOrigin = header["yllcenter"]
    nRows = header["nrows"]
    nCols = header["ncols"]

    # Process each group
    groupExtents = {}
    for entry in dirList:
        dirName = entry["dirName"]
        geometries = entry["geometries"]
        properties = entry["properties"]

        if not geometries:
            message = f"No geometries found for {dirName}"
            log.error(message)
            raise ValueError(message)

        # Check if any features in the group are isolated
        bufferSize = cfg["GENERAL"].getfloat("bufferSize")
        checkFeatureIsolation(geometries, properties, bufferSize, dirName)

        # Get extent of all geometries in group
        bounds = [geom.bounds for geom in geometries]
        xMins, yMins, xMaxs, yMaxs = zip(*bounds)

        # Calculate extent with buffer
        bufferSize = cfg["GENERAL"].getfloat("bufferSize")
        xMin = min(xMins) - bufferSize
        xMax = max(xMaxs) + bufferSize
        yMin = min(yMins) - bufferSize
        yMax = max(yMaxs) + bufferSize
        groupExtents[dirName] = (xMin, xMax, yMin, yMax)  # Store extent for this group

        # Convert extent to grid indices (using top-left origin)
        colStart = max(0, int((xMin - xOrigin) / cellSize))
        colEnd = min(nCols, int((xMax - xOrigin) / cellSize) + 1)

        # Convert y-coordinates to row indices (flipped for bottom-left origin)
        rowStart = max(0, int((yOrigin + nRows * cellSize - yMax) / cellSize))
        rowEnd = min(nRows, int((yOrigin + nRows * cellSize - yMin) / cellSize) + 1)

        # Ensure valid row indices
        if rowEnd <= rowStart:
            log.warning(f"Invalid row indices calculated for {dirName}: start={rowStart}, end={rowEnd}")
            continue

        # Flip row indices for bottom-left origin
        rowStart, rowEnd = nRows - rowEnd, nRows - rowStart

        # Clip the DEM data
        clippedData = raster[rowStart:rowEnd, colStart:colEnd]

        # Create header for clipped DEM
        clippedHeader = header.copy()
        clippedHeader["ncols"] = colEnd - colStart
        clippedHeader["nrows"] = rowEnd - rowStart
        clippedHeader["xllcenter"] = xOrigin + colStart * cellSize
        clippedHeader["yllcenter"] = yOrigin + rowStart * cellSize

        # Update transformation matrix for clipped DEM
        clippedHeader["transform"] = rasterUtils.transformFromASCHeader(clippedHeader)

        # Write clipped DEM
        outputDEM = outputDir / dirName / "Inputs" / f"{dirName}_DEM"
        rasterUtils.writeResultToRaster(clippedHeader, clippedData, outputDEM, flip=True)
        log.debug(f"Clipped DEM saved to: {outputDEM}")

        # Store final DEM extents (reduced by one pixel on each side to ensure DEM > clip extents of other input)
        xMinDEM = clippedHeader["xllcenter"] + (cellSize * 0.5)
        yMinDEM = clippedHeader["yllcenter"] + (cellSize * 0.5)
        xMaxDEM = clippedHeader["xllcenter"] + (clippedHeader["ncols"] * cellSize) - (cellSize * 0.5)
        yMaxDEM = clippedHeader["yllcenter"] + (clippedHeader["nrows"] * cellSize) - (cellSize * 0.5)
        groupExtents[dirName] = (xMinDEM, xMaxDEM, yMinDEM, yMaxDEM)

    return groupExtents


def clipAndMoveOptionalInput(allSimInputFiles, outputDir, groupExtents):
    """Clip and move ENT and RES files based on group DEM extent.

    Parameters
    ----------
    allSimInputFiles: dict
        With all input information for a com1DFA sim
    outputDir : pathlib.Path
        Path to output directory where clipped files will be saved
    groupExtents : dict
        Dictionary with dirName as key and (xMin, xMax, yMin, yMax) as value,
        containing the DEM clipping extents for each group

    Returns
    -------
    groupFeatures : dict
        Dictionary containing clipped features for each group and type
        {dirName: {'ENT': [...], 'RES': [...]}}
    """
    groupFeatures = {}
    # Process ENT and RES directories
    for dirType in ["ENT", "RES"]:

        if dirType == "ENT":
            if allSimInputFiles["entFile"]:
                shpFile = allSimInputFiles["entFile"]
            else:
                log.info("No entrainment file found")
                continue

        if dirType == "RES":
            if allSimInputFiles["resFile"]:
                shpFile = allSimInputFiles["resFile"]
            else:
                log.info("No resistance file found")
                continue

        # Read shapefile
        fields, fieldNames, properties, geometries, srs = shpConv.readShapefile(shpFile)

        # Process each output directory that has extents
        for entry in outputDir.iterdir():
            if not entry.is_dir() or entry.name not in groupExtents:
                continue

            # Get extent
            xMin, xMax, yMin, yMax = groupExtents[entry.name]
            scenarioBbox = box(xMin, yMin, xMax, yMax)

            # Initialize group in dictionary if not exists
            if entry.name not in groupFeatures:
                groupFeatures[entry.name] = {"ENT": [], "RES": []}

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
            targetDir = entry / "Inputs" / dirType
            fU.makeADir(targetDir)
            outputPath = targetDir / f"{entry.name}_{dirType}.shp"
            shpConv.writeShapefile(outputPath, fields, fieldNames, clippedFeatures, srs)
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
    scenarios: dict
        Dictionary mapping scenario names to lists of shape records
    """
    scenarios = {}
    for shapeRecord in shapefile.Reader(str(inputShp)).iterShapeRecords():
        properties = {k.lower(): v for k, v in zip(fieldNames, shapeRecord.record)}
        scenarioValues = properties.get("scenario", "").split(",")
        for scenario in scenarioValues:
            # Check if scenario value is empty and set flag
            if scenario.strip() == "":
                scenario = "NULL"
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
    filteredFields = [field for field in fields if field[0].lower() != "scenario"]
    filteredFieldNames = [name for name in fieldNames if name.lower() != "scenario"]

    # Write the shapefile
    shpConv.writeShapefile(outputShp, filteredFields, filteredFieldNames, shapeFeatures, srs)


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
        inputShp = pathlib.Path(outputDir) / folder["dirName"] / "Inputs" / "REL" / folder["dirName"]
        fields, fieldNames, properties, geometries, srs = shpConv.readShapefile(inputShp)
        totalInputFiles += 1

        # Get the scenario attribute values
        if "scenario" in map(str.lower, fieldNames):
            # Group records by scenario
            scenarios = getScenarioGroups(inputShp, fieldNames)

            # Write a shapefile for each scenario
            for scenario, records in scenarios.items():
                if all(scenario == "NULL" for scenario in scenarios):
                    outputShp = (
                        pathlib.Path(outputDir)
                        / folder["dirName"]
                        / "Inputs"
                        / "REL"
                        / f"{folder['dirName']}_REL"
                    )
                elif scenario == "NULL":
                    outputShp = (
                        pathlib.Path(outputDir)
                        / folder["dirName"]
                        / "Inputs"
                        / "REL"
                        / f"{folder['dirName']}_NULL"
                    )
                else:
                    outputShp = (
                        pathlib.Path(outputDir)
                        / folder["dirName"]
                        / "Inputs"
                        / "REL"
                        / f"{folder['dirName']}_{scenario}"
                    )

                writeScenarioShapefile(outputShp, records, fields, fieldNames, srs)
                totalScenarioFiles += 1

            # Delete the intermediate shapefile
            for ext in [".shp", ".shx", ".dbf", ".prj"]:
                if (inputShp.with_suffix(ext)).exists():
                    (inputShp.with_suffix(ext)).unlink()
        else:
            # If no scenario attribute exists, rename the file (necessary for further processing)
            outputShp = (
                pathlib.Path(outputDir) / folder["dirName"] / "Inputs" / "REL" / f"{folder['dirName']}_REL"
            )
            shapeFeatures = [
                (dict(zip(fieldNames, record.record)), record.shape)
                for record in shapefile.Reader(str(inputShp)).iterShapeRecords()
            ]
            shpConv.writeShapefile(outputShp, fields, fieldNames, shapeFeatures, srs)
            for ext in [".shp", ".shx", ".dbf", ".prj"]:
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
    reportPath = outputDir / "splitInputs_scenarioReport.txt"

    with open(reportPath, "w") as f:
        f.write("SCENARIO REPORT\n")
        f.write("==============\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Process each group and their scenarios
        for group in sorted(dirListGrouped, key=lambda x: x["dirName"].lower()):
            dirName = group["dirName"]
            f.write(f"Group: {dirName}\n")
            f.write("-" * (len(dirName) + 7) + "\n\n")

            relPath = pathlib.Path(outputDir) / dirName / "Inputs" / "REL"
            scenarioFiles = sorted(relPath.glob(f"{dirName}_*.shp"), key=lambda x: x.stem.split("_")[-1])

            if not scenarioFiles:
                f.write("No scenarios found\n\n")
                continue

            # Write release areas for each scenario
            for scenFile in scenarioFiles:
                fields, fieldNames, properties, geometries, _ = shpConv.readShapefile(scenFile)
                scenName = scenFile.stem.split("_")[-1]

                f.write(f"Scenario: {scenName}\n")
                f.write(f"No. of release areas: {len(geometries)}\n")

                if "name" in map(str.lower, fieldNames):  # Handle case sensitivity
                    nameIdx = [i for i, name in enumerate(fieldNames) if name.lower() == "name"][0]
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
                totalEnt = sum(len(shpConv.readShapefile(ef)[3]) for ef in entFiles)
                if totalEnt > 0:
                    f.write(f"No. of entrainment areas: {totalEnt}\n")

            resFiles = list(resPath.glob(f"{dirName}_*.shp"))
            if resFiles:
                totalRes = sum(len(shpConv.readShapefile(rf)[3]) for rf in resFiles)
                if totalRes > 0:
                    f.write(f"No. of resistance areas: {totalRes}\n")

            f.write("\n")

    log.info(f"Scenario report written to: {reportPath}")


def createReportPlots(dirListGrouped, inputDEM, outputDir, groupExtents, groupFeatures):
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
    basicPath = createReportPlot(dirListGrouped, inputDEM, outputDir, groupExtents, groupFeatures, "basic")
    log.info(f"Visual report (basic) written to: {basicPath}")

    # Create optional features report
    optionalPath = createReportPlot(
        dirListGrouped, inputDEM, outputDir, groupExtents, groupFeatures, "optional"
    )
    log.info(f"Visual report (optional) written to: {optionalPath}")
