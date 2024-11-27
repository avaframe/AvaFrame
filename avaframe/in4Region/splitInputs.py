"""Helper functions for splitting avalanche input data into individual folders"""

import logging
import shapefile  # pyshp
from shapely.geometry import shape, Point, box
from avaframe.in2Trans import ascUtils

from avaframe.in3Utils.initializeProject import initializeFolderStruct

# setup local logger
log = logging.getLogger(__name__)

def createFolderList(inputShp):
    """Create a list of folders to be created based on the input shapefile."""
    folderList = []
    unnamedCount = 1

    with shapefile.Reader(inputShp) as src:
        fields = src.fields[1:]  # Skip deletion flag
        fieldNames = [field[0] for field in fields]

        for shapeRecord in src.iterShapeRecords():  # Use iterShapeRecords to access both geometry and attributes
            properties = dict(zip(fieldNames, shapeRecord.record))
            folderName = properties.get('name', '').strip() or f"unnamedRelScenario{str(unnamedCount).zfill(5)}"
            if not properties.get('name', '').strip():
                unnamedCount += 1

            geometry = shape(shapeRecord.shape.__geo_interface__)  # Access geometry
            folderList.append({
                'folderName': folderName,
                'properties': properties,
                'geometry': geometry,
            })

    log.info(f"Created folder list with '{len(folderList)}' entries.")
    return folderList

def createFoldersForReleaseAreas(folderList, outputDir):
    """Create folders for each release area."""
    for entry in folderList:
        folderName = entry['folderName']
        avaDirPath = outputDir / folderName

        initializeFolderStruct(str(avaDirPath), removeExisting=True)
        log.debug(f"Created folder structure for '{folderName}'.")


def splitAndMoveReleaseAreas(folderList, inputShp, outputDir):
    """Split release areas into individual shapefiles and move them to their respective folders."""
    with shapefile.Reader(inputShp) as src:
        fields = src.fields[1:]  # Skip deletion flag
        fieldNames = [field[0] for field in fields]

        for entry in folderList:
            folderName = entry['folderName']
            geometry = entry['geometry']
            avaDirPath = outputDir / folderName
            relDir = avaDirPath / 'Inputs' / 'REL'
            relDir.mkdir(parents=True, exist_ok=True)
            featureOutPath = relDir / f"{folderName}.shp"

            # Write shapefile for the feature
            with shapefile.Writer(str(featureOutPath)) as dst:
                dst.fields = fields  # Copy fields
                dst.record(*entry['properties'].values())  # Copy data
                dst.shape(geometry)  # Copy geometry

            log.debug(f"Moved release area '{folderName}' to '{featureOutPath}'.")


def splitDEMByCenterpointAndMove(folderList, inputDEM, outputDir, cfg):
    """Clip the DEM around each release area's centerpoint and move it to respective folders."""
    # Read DEM data
    demData = ascUtils.readRaster(inputDEM)
    header = demData['header']
    raster = demData['rasterData']

    cellsize = header["cellsize"]
    xOrigin = header["xllcenter"]
    yOrigin = header["yllcenter"]

    log.debug(f"Input DEM Header: {header}")

    for entry in folderList:
        folderName = entry['folderName']
        geometry = entry['geometry']
        center = geometry.centroid

        # Calculate bounding box indices for clipping
        bufferSize = float(cfg['DEM splitting']['bufferSize'])
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

        # Validate grid indices (happens e.g. if the bounding box falls completely outside the DEM)
        if colStart >= colEnd or rowStart >= rowEnd:
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