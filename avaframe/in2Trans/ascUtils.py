"""
    ASCII file readers and handlers

    This file is part of Avaframe.
"""
import os
from decimal import *
import numpy as np
import logging


# create local logger
log = logging.getLogger(__name__)

class cASCheader:
    def __init__(self):
        self.nrows = None
        self.ncols = None
        self.cellsize = None
        self.xllcorner = None
        self.xllcenter = None
        self.yllcorner = None
        self.yllcenter = None
        self.noDataValue = None

    def __str__(self):
        rStr = ''
        rStr += "nrows:\t%s\nncols:\t%s\ncellsize:\t%s\n" % (self.nrows, self.ncols, self.cellsize)
        rStr += "xllcorner:\t%s\nyllcorner:\t%s\nnoDataValue:\t%s\n" % (
            self.xllcorner, self.yllcorner, self.noDataValue)
        return rStr


def readASCheader(fname):
    """input: Filename (z.B.: 'xx/xx/file.asc')
    returns: ascHeader for filename (see class cASCheader)
    """
    dhmInfo = cASCheader()
    infile = open(fname, "r")
    i = 0
    ln = 0
    for aline in infile:
        item = aline.split()
        if ln == 0:
            dhmInfo.ncols = int(item[1])
        elif ln == 1:
            dhmInfo.nrows = int(item[1])
        elif ln == 2:
            if item[0] == "xllcorner" or item[0] == "XLLCORNER" or item[0] == "Xllcorner":
                dhmInfo.xllcorner = float(item[1])
            elif item[0] == "xllcenter" or item[0] == "XLLCENTER" or item[0] == "Xllcenter":
                dhmInfo.xllcenter = float(item[1])
        elif ln == 3:
            if item[0] == "yllcorner" or item[0] == "YLLCORNER" or item[0] == "Yllcorner":
                dhmInfo.yllcorner = float(item[1])
            elif item[0] == "yllcenter" or item[0] == "YLLCENTER" or item[0] == "Yllcenter":
                dhmInfo.yllcenter = float(item[1])
        elif ln == 4:
            dhmInfo.cellsize = float(item[1])
        elif ln == 5:
            dhmInfo.noDataValue = float(item[1])
        ln += 1

    if dhmInfo.xllcenter is None and dhmInfo.yllcenter is None:
        dhmInfo.xllcenter = dhmInfo.xllcorner + dhmInfo.cellsize / 2
        dhmInfo.yllcenter = dhmInfo.yllcorner + dhmInfo.cellsize / 2
    else:
        dhmInfo.xllcorner = dhmInfo.xllcenter - dhmInfo.cellsize / 2
        dhmInfo.yllcorner = dhmInfo.yllcenter - dhmInfo.cellsize / 2

    infile.close()
    return dhmInfo


def isEqualASCheader(headerA, headerB):
    # test if two headers (A,B) are the same (except for noData Values)
    a = cASCheader()
    b = cASCheader()
    a = headerA
    b = headerB
    return (a.ncols == b.ncols) and (a.nrows == b.nrows) and (a.xllcenter == b.xllcenter) and\
           (a.xllcorner == b.xllcorner) and (a.yllcenter == b.yllcenter) and \
           (a.yllcorner == b.yllcorner) and (a.cellsize == b.cellsize)


def readASCdata2numpyArray(fName, headerFile=None):

    if headerFile is not None:
        # open input File
        infile = open(fName, "r")

        # create numpy array with dimesions based on ncols and nrows
        # initialize values of array with respective noDataValue obtained
        # from the header file
        data = np.zeros((headerFile.nrows, headerFile.ncols))
        data = np.where(data == 0, float(headerFile.noDataValue), data)

        i = 0  # index in vertical direction (row-count)
        j = 0  # index in horizontal direction (column-count)
        for aline in infile:
            item = aline.split()
            if i > 5:
                for x in range(headerFile.ncols):
                    data[i - 6, j] = float(item[x])
                    j += 1
                j = 0
            i += 1

    else:
        header = readASCheader(fName)
        infile = open(fName, "r")

        # create numpy array with dimesions based on ncols and nrows
        # initialize values of array with respective noDataValue obtained
        # from the header file
        data = np.zeros((header.nrows, header.ncols))
        data = np.where(data == 0, float(header.noDataValue), data)

        i = 0  # index in vertical direction (row-count)
        j = 0  # index in horizontal direction (column-count)
        for aline in infile:
            item = aline.split()
            if i > 5:
                for x in range(header.ncols):
                    data[i-6, j] = float(item[x])
                    j += 1
                j = 0
            i += 1

    return (data)


def readRaster(fname):
    """ Read raster file (.asc)"""

    log.debug('Reading dem : %s', fname)
    header = readASCheader(fname)
    rasterdata = readASCdata2numpyArray(fname, header)
    rasterdata[rasterdata == header.noDataValue] = np.NaN
    dem = {}
    dem['header'] = header
    dem['rasterData'] = np.flipud(rasterdata)
    return dem


def writeResultToAsc(header, resultArray, outFileName, flip=False):
    """ Write 2D array to an ascii file with header and save to location of outFileName

        Parameters
        ----------
        header : class
            class with methods that give cellsize, nrows, ncols, xllcenter
            yllcenter, noDataValue
        resultArray : numpy.ndarray
            2D numpy array of values that shall be written to file
        outFileName : str
            path incl. name of file to be written
        """

    # Open outfile
    with open(outFileName, 'w') as outFile:

        # write the header and array values to file
        outFile.write("ncols %d\n" % header.ncols)
        outFile.write("nrows %d\n" % header.nrows)
        outFile.write("xllcenter %.2f\n" % header.xllcenter)
        outFile.write("yllcenter %.2f\n" % header.yllcenter)
        outFile.write("cellsize %.2f\n" % header.cellsize)
        outFile.write("nodata_value %.2f\n" % header.noDataValue)

        M = resultArray.shape[0]
        for m in range(M):
            if flip:
                line = np.array([resultArray[M-m-1,:]])
            else:
                line = np.array([resultArray[m,:]])
            np.savetxt(outFile, line, fmt='%f')

        outFile.close()
