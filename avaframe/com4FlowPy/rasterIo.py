# -*- coding: utf-8 -*-

"""
    Functions to handle raster files.
"""

import rasterio
import sys
import logging

# create local logger
log = logging.getLogger(__name__)


def readHeader(inputFile):
    """
    Reads the header of the raster file
    raster file should be readable by rasterio (e.g. .tif, .asc)

    Parameters
    -----------
    inputFile: str
        path to raster file

    Returns
    -----------
    header: dict
        header of raster file in style of ASCII-Rasters
    """

    raster = rasterio.open(inputFile)
    if raster is None:
        print("Unable to open {}".format(inputFile))
        sys.exit(1)

    header = {}
    header["ncols"] = raster.width
    header["nrows"] = raster.height
    header["xllcorner"] = (raster.transform * (0, 0))[0]
    header["yllcorner"] = (raster.transform * (0, raster.height))[1]
    header["cellsize"] = raster.transform[0]
    header["noDataValue"] = raster.nodata
    return header


def readRaster(inputFile):
    """
    Reads in a raster file

    Parameters
    -----------
    inputFile: str
        path to raster file

    Returns
    -----------
    readArray: np.array
        numpy array with values read in from the raster file
    header: dict
        header of raster file in style of ASCII-Rasters
    """

    header = readHeader(inputFile)
    raster = rasterio.open(inputFile)
    readArray = raster.read(1)

    return readArray, header


def outputRaster(referenceFile, fileOut, raster):
    """
    Saves raster

    Parameters
    -----------
    referenceFile: str
        path to raster file to reference on, mostly DEM
    fileOut: str
        path for the outputfile, possible extends are .asc or .tif
    raster: np.array
        raster (array) that is saved
    """

    rasterTrans = rasterio.open(referenceFile)
    try:
        crs = rasterio.crs.CRS.from_dict(rasterTrans.crs.data)
    except:
        # crs = rasterio.crs.CRS.from_epsg(4326)
        crs = None

    _success = True

    if fileOut.suffix == ".asc":
        _driver = "AAIGrid"
    elif fileOut.suffix == ".tif":
        _driver = "GTiff"

    try:
        with rasterio.open(
            fileOut,
            "w",
            driver=_driver,
            height=raster.shape[0],
            width=raster.shape[1],
            count=1,
            dtype=raster.dtype,
            crs=crs,
            transform=rasterTrans.transform,
            nodata=-9999,
        ) as newDataset:
            newDataset.write(raster, 1)
    except:
        _success = False
        log.error("could not write {} to {}".format(raster, fileOut))

    try:
        if _success is True:
            log.info("wrote file: {}".format(fileOut))
        else:
            log.info("failed to write file: {}".format(fileOut))
    except:
        pass
