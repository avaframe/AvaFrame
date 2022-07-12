"""
    ASCII file reader and handler

"""

import numpy as np
import logging


# create local logger
log = logging.getLogger(__name__)


def readASCheader(fname):
    """ return a class with information from an ascii file header

    Parameters
    -----------

    fname: str or pathlib object
        path to ascii file

    Returns
    --------
    headerInfo: class
        information that is stored in header (ncols, nrows, xllcenter, yllcenter, noDataValue)
    """

    headerInfo = {}
    infile = open(fname, "r")
    i = 0
    ln = 0
    for aline in infile:
        item = aline.split()
        if ln == 0:
            if item[0].lower() == 'ncols':
                headerInfo['ncols'] = int(item[1])
            else:
                message = 'DEM header is not in correct format - check line 1: should be ncols'
                log.error(message)
                raise ValueError(message)
        elif ln == 1:
            if item[0].lower() == 'nrows':
                headerInfo['nrows'] = int(item[1])
            else:
                message = ('DEM header is not in correct format - check line 2: should be nrows')
                log.error(message)
                raise ValueError(message)
        elif ln == 2:
            if item[0].lower() == "xllcorner":
                headerInfo['xllcorner'] = float(item[1])
            elif item[0].lower() == "xllcenter":
                headerInfo['xllcenter'] = float(item[1])
            else:
                message = ('DEM header is not in correct format - check line 3: should be xllcenter or xllcorner')
                log.error(message)
                raise ValueError(message)
        elif ln == 3:
            if item[0].lower() == "yllcorner":
                headerInfo['yllcorner'] = float(item[1])
            elif item[0].lower() == "yllcenter":
                headerInfo['yllcenter'] = float(item[1])
            else:
                message = ('DEM header is not in correct format - check line 4: should be yllcenter or yllcorner')
                log.error(message)
                raise ValueError(message)
        elif ln == 4:
            if item[0].lower() == 'cellsize':
                headerInfo['cellsize'] = float(item[1])
            else:
                message = ('DEM header is not in correct format - check line 5: should be cellsize')
                log.error(message)
                raise ValueError(message)
        elif ln == 5:
            if 'nodata' in item[0].lower():
                headerInfo['noDataValue'] = float(item[1])
            else:
                message = ('DEM header is not in correct format - check line 6: should be noDataValue')
                log.error(message)
                raise ValueError(message)
        ln += 1

    if 'xllcenter' not in headerInfo:
        headerInfo['xllcenter'] = headerInfo['xllcorner'] + headerInfo['cellsize'] / 2
        headerInfo['yllcenter'] = headerInfo['yllcorner'] + headerInfo['cellsize'] / 2
        # remove xllcorner, yllcorner
        headerInfo.pop('xllcorner')
        headerInfo.pop('yllcorner')

    infile.close()

    return headerInfo


def isEqualASCheader(headerA, headerB):
    """ Test if two headers (A,B) are the same (except for noData Values)

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
    return (a['ncols'] == b['ncols']) and (a['nrows'] == b['nrows']) and (a['xllcenter'] == b['xllcenter']) and\
            (a['yllcenter'] == b['yllcenter']) and (a['cellsize'] == b['cellsize'])


def readASCdata2numpyArray(fName):
    """ Read ascii matrix as numpy array

    Parameters
    -----------

    fname: str or pathlib object
        path to ascii file

    Returns
    --------
    -rasterdata : 2D numpy array
            2D numpy array of ascii matrix
    """
    infile = open(fName, "r")
    rasterdata = np.loadtxt(fName, skiprows=6)
    infile.close()
    return (rasterdata)


def readRaster(fname, noDataToNan=True):
    """ Read raster file (.asc)

    Parameters
    -----------

    fname: str or pathlib object
        path to ascii file
    noDataToNan: bool
        if True convert noDataValues to nan and set noDataValue to nan

    Returns
    --------
    data: dict
        -headerInfo: class
            information that is stored in header (ncols, nrows, xllcenter, yllcenter, noDataValue)
        -rasterdata : 2D numpy array
                2D numpy array of ascii matrix
    """

    log.debug('Reading dem : %s', fname)
    header = readASCheader(fname)
    rasterdata = readASCdata2numpyArray(fname)

    data = {}
    data['header'] = header
    if noDataToNan:
        rasterdata[rasterdata == header['noDataValue']] = np.NaN
        data['header']['noDataValue'] = np.nan
    data['rasterData'] = np.flipud(rasterdata)

    return data


def writeResultToAsc(header, resultArray, outFileName, flip=False):
    """ Write 2D array to an ascii file with header and save to location of outFileName

        Parameters
        ----------
        header : class
            class with methods that give cellsize, nrows, ncols, xllcenter
            yllcenter, noDataValue
        resultArray : 2D numpy array
            2D numpy array of values that shall be written to file
        outFileName : str
            path incl. name of file to be written
        flip: boolean
            if True, flip the rows of the resultArray when writing
        """

    # Open outfile
    with open(outFileName, 'w') as outFile:

        # write the header and array values to file
        outFile.write("ncols %d\n" % header['ncols'])
        outFile.write("nrows %d\n" % header['nrows'])
        outFile.write("xllcenter %.2f\n" % header['xllcenter'])
        outFile.write("yllcenter %.2f\n" % header['yllcenter'])
        outFile.write("cellsize %.2f\n" % header['cellsize'])
        outFile.write("noDataValue %.2f\n" % header['noDataValue'])

        M = resultArray.shape[0]
        for m in range(M):
            if flip:
                line = np.array([resultArray[M-m-1,:]])
            else:
                line = np.array([resultArray[m,:]])
            np.savetxt(outFile, line, fmt='%.16g')

        outFile.close()
