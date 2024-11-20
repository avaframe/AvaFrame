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
            folder_name = properties.get('name', '').strip() or f"unnamed_rel_scenario_{unnamed_count}"
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

        initializeFolderStruct(str(ava_folder_path), removeExisting=False)
        log.info(f"Created folder structure for '{folder_name}'.")


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

    for entry in folder_list:
        folder_name = entry['folder_name']
        geometry = entry['geometry']
        center = geometry.centroid

        # Calculate bounding box indices for the clipping
        buffer_size = 2500  # 5x5 km -> 2.5 km buffer in each direction #ToDO: move to config
        xmin, ymin, xmax, ymax = (
            center.x - buffer_size,
            center.y - buffer_size,
            center.x + buffer_size,
            center.y + buffer_size,
        )

        # Convert bounding box to grid indices
        col_start = int((xmin - x_origin) / cellsize)
        col_end = int((xmax - x_origin) / cellsize)
        row_start = int((y_origin - ymax) / cellsize)
        row_end = int((y_origin - ymin) / cellsize)

        # Clip the raster data
        clipped_raster = raster[row_start:row_end, col_start:col_end]

        # Update header for clipped DEM
        clipped_header = header.copy()
        clipped_header["ncols"] = clipped_raster.shape[1]
        clipped_header["nrows"] = clipped_raster.shape[0]
        clipped_header["xllcenter"] = xmin + cellsize / 2
        clipped_header["yllcenter"] = ymin + cellsize / 2

        # Save clipped DEM
        ava_folder_path = output_dir / folder_name
        dem_output_path = ava_folder_path / 'Inputs' / 'DEM_clipped.asc'
        ascUtils.writeResultToAsc(clipped_header, clipped_raster, dem_output_path)

        log.info(f"Clipped and moved DEM for '{folder_name}' to '{dem_output_path}'.")