# -*- coding: utf-8 -*-

"""
    Functions to handle raster files.
"""

import rasterio
import sys
import logging

# create local logger
log = logging.getLogger(__name__)


def read_header(input_file):
    """
    Reads the header of the raster file
    raster file should be readable by rasterio (e.g. .tif, .asc)

    Parameters
    -----------
    input_file: str
        path to raster file

    Returns
    -----------
    header: dict
        header of raster file in style of ASCII-Rasters
    """

    raster = rasterio.open(input_file)
    if raster is None:
        print("Unable to open {}".format(input_file))
        sys.exit(1)

    header = {}
    header["ncols"] = raster.width
    header["nrows"] = raster.height
    header["xllcorner"] = (raster.transform * (0, 0))[0]
    header["yllcorner"] = (raster.transform * (0, raster.height))[1]
    header["cellsize"] = raster.transform[0]
    header["noDataValue"] = raster.nodata
    return header


def read_raster(input_file):
    """
    Reads in a raster file

    Parameters
    -----------
    input_file: str
        path to raster file

    Returns
    -----------
    my_array: np.array
        numpy array with values read in from the raster file
    header: dict
        header of raster file in style of ASCII-Rasters
    """

    header = read_header(input_file)
    raster = rasterio.open(input_file)
    my_array = raster.read(1)

    return my_array, header


def output_raster(referenceFile, file_out, raster):
    """
    Saves raster

    Parameters
    -----------
    referenceFile: str
        path to raster file to reference on, mostly DEM
    file_out: str
        path for the outputfile, possible extends are .asc or .tif
    raster: np.array
        raster (array) that is saved
    """

    raster_trans = rasterio.open(referenceFile)
    try:
        crs = rasterio.crs.CRS.from_dict(raster_trans.crs.data)
    except:
        # crs = rasterio.crs.CRS.from_epsg(4326)
        crs = None

    _success = True

    if file_out.suffix == ".asc":
        _driver = "AAIGrid"
    elif file_out.suffix == ".tif":
        _driver = "GTiff"

    try:
        with rasterio.open(
            file_out,
            "w",
            driver=_driver,
            height=raster.shape[0],
            width=raster.shape[1],
            count=1,
            dtype=raster.dtype,
            crs=crs,
            transform=raster_trans.transform,
            nodata=-9999,
        ) as new_dataset:
            new_dataset.write(raster, 1)
    except:
        _success = False
        log.error("could not write {} to {}".format(raster, file_out))

    try:
        if _success is True:
            log.info("wrote file: {}".format(file_out))
        else:
            log.info("failed to write file: {}".format(file_out))
    except:
        pass
