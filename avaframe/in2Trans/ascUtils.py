"""
    ASCII file reader and handler

"""

import logging

import numpy as np

# create local logger
log = logging.getLogger(__name__)


def readASCheader(fname):
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

    # read header
    log.warning("This function is deprecated. Look at rasterUtils for replacements.")
    headerRows = 6  # six rows for header information
    headerInfo = (
        {}
    )  # store header information including ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value
    rowCount = 1
    with open(str(fname), "rt") as fileH:
        for line in fileH:
            if rowCount <= headerRows:
                line = line.split()
                headerInfo[line[0].lower()] = float(line[1])
            else:
                break
            rowCount = rowCount + 1

    if ("xllcenter" not in headerInfo and "xllcorner" in headerInfo) and (
            "yllcenter" not in headerInfo and "yllcorner" in headerInfo
    ):
        headerInfo["xllcenter"] = headerInfo["xllcorner"] + headerInfo["cellsize"] / 2
        headerInfo["yllcenter"] = headerInfo["yllcorner"] + headerInfo["cellsize"] / 2
        # remove xllcorner, yllcorner
        headerInfo.pop("xllcorner")
        headerInfo.pop("yllcorner")

    # convert ncols and nrows to int
    headerInfo["ncols"] = int(headerInfo["ncols"])
    headerInfo["nrows"] = int(headerInfo["nrows"])

    headerItems = [item.lower() for item in list(headerInfo.keys())]
    if sorted(headerItems) != sorted(
            ["cellsize", "nrows", "ncols", "xllcenter", "yllcenter", "nodata_value"]
    ):
        message = "DEM header is not in correct format - needs to contain values for: cellsize, nrows, ncols, xllcenter(-corner), yllcenter(-corner), nodata_value"
        log.error(message)
        raise ValueError(message)

    fileH.close()

    return headerInfo


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
    log.warning("This function is deprecated. Look at rasterUtils for replacements.")
    a = headerA
    b = headerB
    return (
            (a["ncols"] == b["ncols"])
            and (a["nrows"] == b["nrows"])
            and (a["xllcenter"] == b["xllcenter"])
            and (a["yllcenter"] == b["yllcenter"])
            and (a["cellsize"] == b["cellsize"])
    )


def readASCdata2numpyArray(fName):
    """Read ascii matrix as numpy array

    Parameters
    -----------

    fname: str or pathlib object
        path to ascii file

    Returns
    --------
    -rasterdata : 2D numpy array
            2D numpy array of ascii matrix
    """
    log.warning("This function is deprecated. Look at rasterUtils for replacements.")
    infile = open(fName, "r")
    rasterdata = np.loadtxt(fName, skiprows=6)
    infile.close()
    return rasterdata


def readRaster(fname, noDataToNan=True):
    """Read raster file (.asc)

    Parameters
    -----------

    fname: str or pathlib object
        path to ascii file
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

    log.warning("This function is deprecated. Look at rasterUtils for replacements.")
    log.debug("Reading dem : %s", fname)
    header = readASCheader(fname)
    rasterdata = readASCdata2numpyArray(fname)

    data = {}
    data["header"] = header
    if noDataToNan:
        rasterdata[rasterdata == header["nodata_value"]] = np.nan
        data["header"]["nodata_value"] = np.nan
    data["rasterData"] = np.flipud(rasterdata)

    return data


def writeResultToAsc(header, resultArray, outFileName, flip=False):
    """Write 2D array to an ascii file with header and save to location of outFileName

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
    """

    # Open outfile
    log.warning("This function is deprecated. Look at rasterUtils for replacements.")
    with open(outFileName, "w") as outFile:
        # write the header and array values to file
        outFile.write("ncols %d\n" % header["ncols"])
        outFile.write("nrows %d\n" % header["nrows"])
        outFile.write("xllcenter %.2f\n" % header["xllcenter"])
        outFile.write("yllcenter %.2f\n" % header["yllcenter"])
        outFile.write("cellsize %.2f\n" % header["cellsize"])
        outFile.write("nodata_value %.2f\n" % header["nodata_value"])

        M = resultArray.shape[0]
        for m in range(M):
            if flip:
                line = np.array([resultArray[M - m - 1, :]])
            else:
                line = np.array([resultArray[m, :]])
            np.savetxt(outFile, line, fmt="%.16g")

        outFile.close()
