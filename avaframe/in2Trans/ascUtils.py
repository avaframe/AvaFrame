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
        self.xllcenter = None
        self.yllcenter = None
        self.noDataValue = None

    def __str__(self):
        rStr = ''
        rStr += "nrows:\t%s\nncols:\t%s\ncellsize:\t%s\n" % (self.nrows, self.ncols, self.cellsize)
        rStr += "xllcenter:\t%s\nyllcenter:\t%s\nnoDataValue:\t%s\n" % (
            self.xllcenter, self.yllcenter, self.noDataValue)
        return rStr


def readASCheader(fname):
    """ return a class with information from an ascii file header

    Parameters
    -----------

    fname: str
        path to ascii file

    Returns
    --------
    headerInfo: class
        information that is stored in header (ncols, nrows, xllcenter, yllcenter, noDataValue)
    """


    headerInfo = cASCheader()
    infile = open(fname, "r")
    i = 0
    ln = 0
    for aline in infile:
        item = aline.split()
        if ln == 0:
            if item[0] == 'ncols' or item[0] == 'nCols':
                headerInfo.ncols = int(item[1])
            else:
                log.error('DEM header is not in correct format - check line 1: should be ncols')
                raise ValueError
        elif ln == 1:
            if item[0] == 'nrows' or item[0] == 'nRows':
                headerInfo.nrows = int(item[1])
            else:
                log.error('DEM header is not in correct format - check line 2: should be nrows')
                raise ValueError
        elif ln == 2:
            if item[0] == "xllcorner" or item[0] == "XLLCORNER" or item[0] == "Xllcorner":
                headerInfo.xllcorner = float(item[1])
            elif item[0] == "xllcenter" or item[0] == "XLLCENTER" or item[0] == "Xllcenter":
                headerInfo.xllcenter = float(item[1])
            else:
                log.error('DEM header is not in correct format - check line 3: should be xllcenter or xllcorner')
                raise ValueError
        elif ln == 3:
            if item[0] == "yllcorner" or item[0] == "YLLCORNER" or item[0] == "Yllcorner":
                headerInfo.yllcorner = float(item[1])
            elif item[0] == "yllcenter" or item[0] == "YLLCENTER" or item[0] == "Yllcenter":
                headerInfo.yllcenter = float(item[1])
            else:
                log.error('DEM header is not in correct format - check line 4: should be yllcenter or yllcorner')
                raise ValueError
        elif ln == 4:
            if item[0] == 'cellSize' or item[0] == 'cellsize':
                headerInfo.cellsize = float(item[1])
            else:
                log.error('DEM header is not in correct format - check line 5: should be cellsize')
                raise ValueError
        elif ln == 5:
            if 'nodata' in item[0] or 'NODATA' in item[0]:
                headerInfo.noDataValue = float(item[1])
            else:
                log.error('DEM header is not in correct format - check line 6: should be NODATA_value')
                raise ValueError
        ln += 1

    if headerInfo.xllcenter is None:
        headerInfo.xllcenter = headerInfo.xllcorner + headerInfo.cellsize / 2
        headerInfo.yllcenter = headerInfo.yllcorner + headerInfo.cellsize / 2
        # remove xllcorner, yllcorner
        delattr(headerInfo, 'xllcorner')
        delattr(headerInfo, 'yllcorner')

    infile.close()

    return headerInfo


def isEqualASCheader(headerA, headerB):
    # test if two headers (A,B) are the same (except for noData Values)
    a = cASCheader()
    b = cASCheader()
    a = headerA
    b = headerB
    return (a.ncols == b.ncols) and (a.nrows == b.nrows) and (a.xllcenter == b.xllcenter) and\
            (a.yllcenter == b.yllcenter) and (a.cellsize == b.cellsize)


def readASCdata2numpyArray(fName, headerFile=None):
    infile = open(fName, "r")
    rasterdata = np.loadtxt(fName, skiprows=6)
    infile.close()
    return (rasterdata)


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
            np.savetxt(outFile, line, fmt='%.16g')

        outFile.close()
