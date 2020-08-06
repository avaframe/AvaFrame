#! /usr/bin/python

import numpy
import re

def esriOpen(fileName):
    """
        ersiOpen(fileName) -> EsriGridFile

        Opens the filename and creates a EsriGridFile datastructure
    """
    f = open(fileName, 'r')
    header = {}
    for mustfind in ['ncols', 'nrows', 'xll.*er', 'yll.*er', 'cellsize', 'NODATA_value']:
        try:
            header[mustfind] = float(re.search( r'%s *(-?\d+)'%(mustfind), f.readline(), re.M|re.I).group(1))
        except :
            raise RuntimeError('unable to find %s in header of %s'%(mustfind, fileName))
    data = EsriGridFile(header['ncols'], header['nrows'], header['xll.*er'], header['yll.*er'], header['cellsize'], header['NODATA_value'])
    row = 0
    for line in f:
        col = 0
        for s in line.rstrip().split(' '):
            data.data[row, col] = float(s)
            col = col+1
        row = row+1
    return data

class EsriGridCell(object):

    def __init__(self, parent, col, row):

        self.parent = parent
        self.col = col
        self.row = row

    def center(self):
        return self.parent.cellCenter(self.col, self.row)

    def value(self):
        return self.parent.data[self.row, self.col]

class EsriGridFile(object):

    def __init__(self, ncols, nrows, xllcenter, yllcenter, cellsize, NODATA_value):

        self.nrows = int(nrows)
        self.ncols = int(ncols)
        self.data = numpy.empty((nrows,ncols))
        self.xllcenter = xllcenter
        self.yllcenter = yllcenter
        self.cellsize = cellsize
        self.NODATA_value = NODATA_value

    def left(self):

        return self.xllcenter

    def right(self):

        return self.xllcenter+self.cellsize*self.ncols

    def top(self):

        return self.yllcenter

    def bottom(self):

        return self.yllcenter+self.cellsize*self.nrows

    def valueAtXYCell(self, X, Y):
        return self.data[self.nrows-Y, X]

    def cellCenter(self, col, row):
        x = self.xllcenter + col * self.cellsize
        y = self.yllcenter + (self.nrows - row) * self.cellsize
        return numpy.array([x, y])

    def __iter__(self):
        self.__counter__ = 0
        return self

    def next(self):
        if self.__counter__ + 2 > self.ncols*self.nrows:
            raise StopIteration
        else:
            self.__counter__ += 1
            return EsriGridCell(self, (self.__counter__-1)%self.ncols, (self.__counter__-1)/self.ncols)


if __name__ == '__main__':
    pressureData = esriOpen('./radartestdata/pressure/000100.txt')

    print(pressureData.xllcenter, pressureData.yllcenter)
    print(pressureData.cellCenter(0, pressureData.nrows))
    print(pressureData.cellCenter(10, 10))
    print(pressureData.cellCenter(10, 11))
