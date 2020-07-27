#####################################################
# Input-Output functions  for GraviProMod
# Einlesen und Ausschreiben von Hohenmodellen, Start-
# zellen bzw. Ergebnisrastern, etc.
#####################################################
# from: 11.03.2013
# by: Huber A.
# last modified: 23.07.2020
# by: Tonnel M.
#####################################################
# Import numpypy
import math
import os
import sys
import copy
import numpy
from tkinter import filedialog  # was from tkFileDialog import * in python 2
from decimal import *
from numpy import *
#from topslide_classes_functions import *


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
    '''input: Filename (z.B.: 'xx/xx/file.asc')
    returns: ascHeader for filename (see class cASCheader)'''

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
        data = zeros((headerFile.nrows, headerFile.ncols))
        data = where(data == 0, float(headerFile.noDataValue), data)

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
        data = zeros((header.nrows, header.ncols))
        data = where(data == 0, float(header.noDataValue), data)

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


def writeResultToAsc(header, resultArray, outType=None):
    # print "choose output file location and name"

    if outType is not None:
        outfile = asksaveasfile(mode='w',
                                defaultextension='.asc',
                                title='Save Output File %s' % (outType),
                                filetypes=[('ASCII', '*.asc')],
                                initialfile='path.asc',
                                initialdir=os.getcwd())
    else:
        outfile = asksaveasfile(mode='w',
                                defaultextension='.asc',
                                title='Save Output File',
                                filetypes=[('ASCII', '*.asc')],
                                initialfile='path.asc',
                                initialdir=os.getcwd())

    outfile.write("ncols %d\n" % header.ncols)
    outfile.write("nrows %d\n" % header.nrows)
    outfile.write("xllcorner %.2f\n" % header.xllcorner)
    outfile.write("yllcorner %.2f\n" % header.yllcorner)
    outfile.write("cellsize %.2f\n" % header.cellsize)
    outfile.write("nodata_value %.2f\n" % header.noDataValue)

    for i in range(resultArray.shape[0]):
        for j in range(resultArray.shape[1]):
            outfile.write("%f " % resultArray[i][j])
        outfile.write("\n")

    outfile.close()

def readpathdata2numpyArray(fName):
    infile = open(fName, "r")
    i = 0
    for aline in infile:
        item = aline.split()
        test = True
        try:
            float(item[0])
        except:
            test = False
        if test:
            break
        i += 1
    if not test:
        raise ValueError('[ReadPD2NA] wrong type of path input file')

    data = numpy.transpose(numpy.loadtxt(fName, skiprows=i))

    return (data)
