"""
    ASCII file reader and handler

"""

import logging
import rasterio
import numpy as np

# create local logger
log = logging.getLogger(__name__)


def readRaster(fname, noDataToNan=True):
    """Read raster file in .asc or .tif format

    Parameters
    -----------

    fname: pathlib object
        path to ascii/tif file
    noDataToNan: bool
        if True convert nodata_values to nan and set nodata_value to nan

    Returns
    --------
    data: dict
        -headerInfo: class
            information that is stored in header (ncols, nrows, xllcenter, yllcenter, nodata_value)
        -rasterdata : 2D numpy array
                2D numpy array of ascii matrix
    """

    log.debug("Reading raster file : %s", fname)

    raster = rasterio.open(fname)
    rasterData = raster.read(1).astype(np.float64)
    header = getHeaderFromRaster(raster)
    raster.close()

    # FSO: for testing whether we got all cases of llcenter vs llcornen. Header2 contains old setup
    header2 = readRasterHeader(fname)
    if (header["yllcenter"] - header2["yllcenter"]) > 0.1:
        message = "Header llcenter mismatch between direct file reading and rasterio"
        log.error(message)
        raise ValueError(message)

    data = {}
    data["header"] = header
    if noDataToNan:
        rasterData[rasterData == header["nodata_value"]] = np.nan
        data["header"]["nodata_value"] = np.nan
    data["rasterData"] = np.flipud(rasterData)

    return data


def getHeaderFromRaster(raster):
    """convert rasterio raster info to header info
    Parameters
    ----------
    raster: rasterio raster
        read by rasterio

    Returns
    -------
    header: dict
        header info
    """
    header = {}
    header["ncols"] = raster.width
    header["nrows"] = raster.height
    header["cellsize"] = raster.transform[0]
    # TODO maker sure xllcenter is correct here or if this needs transformation as in readRasterHeader
    header["xllcenter"] = (raster.transform * (0, 0))[0] + header["cellsize"] / 2.0
    header["yllcenter"] = (raster.transform * (0, raster.height))[1] + header["cellsize"] / 2.0
    header["nodata_value"] = raster.nodata
    header["crs"] = raster.crs
    header["driver"] = raster.driver
    header["transform"] = raster.transform

    return header


def transformFromASCHeader(header):
    """convert header info to raster transform info

    Parameters
    ----------
    header: dict
        header info

    Returns
    -------
    transfrom: dict
        rasterio transform info
    """
    # rasterio requires west, north
    # rasterio.transform.from_origin(west, north, xsize, ysize)
    xllCenter = header["xllcenter"]
    yllCenter = header["yllcenter"]
    cellSize = header["cellsize"]
    nRows = header["nrows"]

    transform = rasterio.transform.from_origin(xllCenter,
                                               yllCenter + nRows * cellSize,
                                               cellSize,
                                               cellSize)
    # crs = rasterio.crs.CRS.from_epsg(31287)

    return transform  # , crs

def readRasterHeader(fname):
    """return a class with information from an ascii file header

    Parameters
    -----------

    fname: str or pathlib object
        path to ascii file

    Returns
    --------
    headerInfo: class
        information that is stored in header (ncols, nrows, xllcenter, yllcenter, nodata_value)
    """

    # # read header
    # headerRows = 6  # six rows for header information
    # headerInfo = (
    #     {}
    # )  # store header information including ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value
    # rowCount = 1
    # with open(str(fname), "rt") as fileH:
    #     for line in fileH:
    #         if rowCount <= headerRows:
    #             line = line.split()
    #             headerInfo[line[0].lower()] = float(line[1])
    #         else:
    #             break
    #         rowCount = rowCount + 1
    #
    # if ("xllcenter" not in headerInfo and "xllcorner" in headerInfo) and (
    #     "yllcenter" not in headerInfo and "yllcorner" in headerInfo
    # ):
    #     headerInfo["xllcenter"] = headerInfo["xllcorner"] + headerInfo["cellsize"] / 2
    #     headerInfo["yllcenter"] = headerInfo["yllcorner"] + headerInfo["cellsize"] / 2
    #     # remove xllcorner, yllcorner
    #     headerInfo.pop("xllcorner")
    #     headerInfo.pop("yllcorner")
    #
    # # convert ncols and nrows to int
    # headerInfo["ncols"] = int(headerInfo["ncols"])
    # headerInfo["nrows"] = int(headerInfo["nrows"])
    #
    # headerItems = [item.lower() for item in list(headerInfo.keys())]

    raster = rasterio.open(fname)
    header = getHeaderFromRaster(raster)
    raster.close()

    headerItems = [item.lower() for item in list(header.keys())]
    # TODO: reanable check
    # if sorted(headerItems) != sorted(
    #     ["cellsize", "nrows", "ncols", "xllcenter", "yllcenter", "nodata_value"]
    # ):
    #     message = "DEM header is not in correct format - needs to contain values for: cellsize, nrows, ncols, xllcenter(-corner), yllcenter(-corner), nodata_value"
    #     log.error(message)
    #     raise ValueError(message)

    # fileH.close()

    return header


def isEqualASCheader(headerA, headerB):
    """Test if two headers (A,B) are the same (except for noData Values)

    Parameters
    -----------
    headerA: class
    headerB: class

    Returns
    --------
    boolean: True if header A and B are equal (disregrad the noData field)
    """
    a = headerA
    b = headerB
    return (
        (a["ncols"] == b["ncols"])
        and (a["nrows"] == b["nrows"])
        and (a["xllcenter"] == b["xllcenter"])
        and (a["yllcenter"] == b["yllcenter"])
        and (a["cellsize"] == b["cellsize"])
    )

def writeResultToRaster(header, resultArray, outFileName, flip=False):
    """Write 2D array to a raster file with header and save to location of outFileName

    Parameters
    ----------
    header : class
        class with methods that give cellsize, nrows, ncols, xllcenter
        yllcenter, nodata_value
    resultArray : 2D numpy array
        2D numpy array of values that shall be written to file
    outFileName : str
        path incl. name of file to be written
    flip: boolean
        if True, flip the rows of the resultArray when writing

    Returns
    -------
    outFile: path
        to file being written
    """

    if header["driver"] == "AAIGrid":
        outFile = outFileName.parent / (outFileName.name + ".asc")
    elif header["driver"] == "GTiff":
        outFile = outFileName.parent / (outFileName.name + ".tif")

    # try:
    rasterOut = rasterio.open(
        outFile,
        "w",
        driver=header["driver"],
        crs=header["crs"],
        nodata=header["nodata_value"],
        transform=header["transform"],
        height=resultArray.shape[0],
        width=resultArray.shape[1],
        count=1,
        dtype=resultArray.dtype,
    )
    if flip:
        rasterOut.write(np.flipud(resultArray), 1)
    else:
        rasterOut.write(resultArray, 1)
    rasterOut.close()
    # except:
    #     log.error("could not write {} to {}".format(resultArray, outFileName))
    return outFile
