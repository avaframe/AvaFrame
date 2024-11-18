"""Helper functions for splitting avalanche input data into individual folders"""

import logging
import shapefile  # pyshp
from shapely.geometry import shape, Point, box
from avaframe.in2Trans import ascUtils

from avaframe.in3Utils.initializeProject import initializeFolderStruct

# setup local logger
log = logging.getLogger(__name__)

def createFolderList(input_shapefile):
    """Create a list of folders to be created based on the input shapefile."""
    folder_list = []
    unnamed_count = 1

    with shapefile.Reader(input_shapefile) as src:
        fields = src.fields[1:]  # Skip deletion flag
        field_names = [field[0] for field in fields]

        for shape_record in src.iterShapeRecords():  # Use iterShapeRecords to access both geometry and attributes
            properties = dict(zip(field_names, shape_record.record))
            folder_name = properties.get('name', '').strip() or f"unnamedRelScenario{unnamed_count}"
            if not properties.get('name', '').strip():
                unnamed_count += 1

            geometry = shape(shape_record.shape.__geo_interface__)  # Access geometry properly
            folder_list.append({
                'folder_name': folder_name,
                'properties': properties,
                'geometry': geometry,
            })

    log.info(f"Created folder list with '{len(folder_list)}' entries.")
    return folder_list

def createFoldersForReleaseAreas(folder_list, output_dir):
    """Create folders for each release area."""
    for entry in folder_list:
        folder_name = entry['folder_name']
        ava_folder_path = output_dir / folder_name

        initializeFolderStruct(str(ava_folder_path), removeExisting=True)
        #log.info(f"Created folder structure for '{folder_name}'.")


def splitAndMoveReleaseAreas(folder_list, input_shapefile, output_dir):
    """Split release areas into individual shapefiles and move them to their respective folders."""
    with shapefile.Reader(input_shapefile) as src:
        fields = src.fields[1:]  # Skip deletion flag
        field_names = [field[0] for field in fields]

        for entry in folder_list:
            folder_name = entry['folder_name']
            geometry = entry['geometry']
            ava_folder_path = output_dir / folder_name
            rel_dir = ava_folder_path / 'Inputs' / 'REL'
            rel_dir.mkdir(parents=True, exist_ok=True)
            feature_output_path = rel_dir / f"{folder_name}.shp"

            # Write shapefile for the feature
            with shapefile.Writer(str(feature_output_path)) as dst:
                dst.fields = fields  # Copy fields
                dst.record(*entry['properties'].values())  # Copy data
                dst.shape(geometry)  # Copy geometry

            log.info(f"Moved release area '{folder_name}' to '{feature_output_path}'.")


def splitDEMByCenterpointAndMove(folder_list, input_dem, output_dir):
    """Clip the DEM around each release area's centerpoint and move it to respective folders."""
    # Read DEM data
    dem_data = ascUtils.readRaster(input_dem)
    header = dem_data['header']
    raster = dem_data['rasterData']

    cellsize = header["cellsize"]
    x_origin = header["xllcenter"]
    y_origin = header["yllcenter"]

    log.info(f"DEM Header: {header}")

    for entry in folder_list:
        folder_name = entry['folder_name']
        geometry = entry['geometry']
        center = geometry.centroid

        # Calculate bounding box indices for clipping
        buffer_size = 2500  # buffer in each direction
        xmin, ymin, xmax, ymax = (
            center.x - buffer_size,
            center.y - buffer_size,
            center.x + buffer_size,
            center.y + buffer_size,
        )

        #log.info(f"Bounding box for '{folder_name}': xmin={xmin}, ymin={ymin}, xmax={xmax}, ymax={ymax}")

        # Convert bounding box to grid indices
        col_start = max(0, int((xmin - x_origin) / cellsize))
        col_end = min(header["ncols"], int((xmax - x_origin) / cellsize))
        # Calculate row indices
        row_start = max(0, int((y_origin + header["nrows"] * cellsize - ymax) / cellsize))
        row_end = min(header["nrows"], int((y_origin + header["nrows"] * cellsize - ymin) / cellsize))

        #log.info(f"Computed grid indices before adjustment for '{folder_name}': "
        #         f"col_start={col_start}, col_end={col_end}, row_start={row_start}, row_end={row_end}")

        # Flip rows for bottom-left origin
        row_start, row_end = header["nrows"] - row_end, header["nrows"] - row_start

        # Validate grid indices
        if col_start >= col_end or row_start >= row_end:
            log.warning(f"Invalid clipping bounds for '{folder_name}'. Skipping.")
            continue

        # Clip DEM section
        dem_clipped = raster[row_start:row_end, col_start:col_end]

        log.info(f"Clipped DEM for '{folder_name}' with '{len(dem_clipped)}' rows and '{len(dem_clipped[0])}' columns")

        # Prepare header for clipped DEM
        clipped_header = {
            "ncols": col_end - col_start,
            "nrows": row_end - row_start,
            "xllcenter": x_origin + col_start * cellsize,
            "yllcenter": y_origin + row_start * cellsize,
            "cellsize": cellsize,
            "nodata_value": header["nodata_value"],
        }

        # Define output path
        dem_output_path = output_dir / folder_name / "Inputs" / f"{folder_name}_DEM.asc"

        # Save clipped DEM
        try:
            ascUtils.writeResultToAsc(clipped_header, dem_clipped, dem_output_path, flip=True)
            log.info(f"Clipped DEM saved to '{dem_output_path}'.")
        except Exception as e:
            log.error(f"Failed to save clipped DEM for '{folder_name}': {e}")