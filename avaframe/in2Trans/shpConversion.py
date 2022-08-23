"""
    Conversion functions to read/ write Shape files
"""


import pathlib
import shapefile
import copy
import numpy as np
import logging

# create local logger
log = logging.getLogger(__name__)


def SHP2Array(infile, defname=None):
    """ Read shapefile and convert it to a python dictionary

    The dictionary contains the name of the paths in the shape file, the np array with
    the coordinates of the feature points (all stacked in the same array)
    and information about the starting index and length of each feature

    Parameters
    ----------
    infile: str
        path to shape file
    defname: str
        name to give to the feature in the shape file

    Returns
    -------
    SHPdata : dict
        sks :
            projection information
        Name : str
            list of feature names
        x : 1D numpy array
            np array of the x coord of the points in the features
        y : 1D numpy array
            np array of the y coord of the points in the features
        z : 1D numpy array
            np array of zeros and same size as the x an y coordinates array
        Start : 1D numpy array
            np array with the starting index of each feature in the coordinates
            arrays (as many indexes as features)
        Lengt : 1D numpy array
            np array with the length of each feature in the coordinates
            arrays (as many indexes as features)
        thickness (optional) : 1D numpy array
            np array with the (release or entrainment) thickness of each feature (as many values as features)
        id : list
            list of oid as string for each feature
        ci95: list
            list of 95% confidence interval of thickness value


    """
    #  Input shapefile
    sf = shapefile.Reader(str(infile))
    infile = pathlib.Path(infile)
    # set defaults for variables
    layername = None
    thickness = None
    rho = None
    sks = None
    iso = None
    id = None
    ci95 = None

    # get coordinate system
    sks = getSHPProjection(infile)

    # Start reading the shapefile
    records = sf.shapeRecords()
    shps = sf.shapes()

    SHPdata = {}
    SHPdata['sks'] = sks
    Name = []
    thicknessList = []
    idList = []
    ci95List = []
    Length = np.empty((0))
    Start = np.empty((0))
    Coordx = np.empty((0))
    Coordy = np.empty((0))
    Coordz = np.empty((0))
    start = 0

    for n, item in enumerate(shps):
        pts = item.points
        zs = [0.0] * len(pts)

        # check if records are available and extract
        if records:
            # loop through fields
            for (name, typ, size, deci), value in zip(sf.fields[1:], records[n].record):
                # get entity name
                name = name.lower()
                if (name == 'name'):
                    layername = str(value)
                if (name == 'thickness') or (name == 'd0'):
                    thickness = value
                if (name == 'ci95'):
                    ci95 = value
                if (name == 'rho'):
                    rho = value
                if (name == 'sks'):
                    sks = value
                if (name == 'iso'):
                    iso = value
            # if name is still empty go through file again and take Layer instead
            if ((type(layername) is bytes) or (layername is None)):
                for (name, typ, size, deci), value in zip(sf.fields[1:], records[n].record):
                    if (name == 'Layer'):
                        layername = value

        # if layer still not defined, use generic
        if layername is None:
            layername = defname

        Name.append(layername)
        log.debug('SHPConv: Found layer %s', layername)
        thicknessList.append(str(thickness))
        ci95List.append(str(ci95))

        Start = np.append(Start, start)
        length = len(pts)
        Length = np.append(Length, length)
        start += length

        for (pt, z) in zip(pts, zs):
            Coordx = np.append(Coordx, pt[0])
            Coordy = np.append(Coordy, pt[1])
            Coordz = np.append(Coordz, z)

    # get unique ID of features in shapefile
    for rec in sf.records():
        id = rec.oid
        idList.append(str(id))

    SHPdata['Name'] = Name
    SHPdata['thickness'] = thicknessList
    SHPdata['Start'] = Start
    SHPdata['Length'] = Length
    SHPdata['x'] = Coordx
    SHPdata['y'] = Coordy
    SHPdata['z'] = Coordz
    SHPdata['id'] = idList
    SHPdata['ci95'] = ci95List

    sf.close()

    return SHPdata


def getSHPProjection(infile):
    """ Fetch projection from shp file

    Parameters
    ----------
    infile: str
        path to shape file

    Returns
    -------
    sks: str
        projection string (if available, None if not)
    """
    # get coordinate system
    prjfile = infile.with_suffix('.prj')
    if prjfile.is_file():
        prjf = open(prjfile, 'r')
        sks = prjf.readline()
        prjf.close()
        return sks
    else:
        message = ('No projection layer for shp file %s' % infile)
        log.warning(message)
        return None


def readThickness(infile, defname=None):
    """ Read shapefile and fetch info on features' ids and thickness values

    Parameters
    ----------
    infile: str
        path to shape file
    defname: str
        name to give to the feature in the shape file

    Returns
    -------
    thickness: list
        list of strings with the (release or entrainment) thickness of each feature (as many values as features)
    id : list
        list of strings for oid of each feature in shp file
    ci95: list
        list of all ci95 values if provided

    """
    #  Input shapefile
    sf = shapefile.Reader(str(infile))
    infile = pathlib.Path(infile)

    thickness = None
    id = None
    ci95 = None

    # Start reading the shapefile
    records = sf.shapeRecords()
    shps = sf.shapes()

    thicknessList = []
    idList = []
    ci95List = []

    for n, item in enumerate(shps):
        pts = item.points
        zs = [0.0] * len(pts)

        # check if records are available and extract
        if records:
            # loop through fields
            for (name, typ, size, deci), value in zip(sf.fields[1:], records[n].record):
                # get entity name
                name = name.lower()
                if (name == 'thickness') or (name == 'd0'):
                    thickness = value
                if (name == 'ci95'):
                    ci95 = value

        thicknessList.append(str(thickness))
        ci95List.append(str(ci95))

    # get unique ID of features in shapefile
    for rec in sf.records():
        id = rec.oid
        idList.append(str(id))

    sf.close()

    return thicknessList, idList, ci95List


def readLine(fname, defname, dem):
    """ Read line from  .shp
    Use SHP2Array to read the shape file.
    Check if the lines are laying inside the dem extend

    Parameters
    ----------
    fname: str
        path to shape file
    defname: str
        name to give to the line in the shape file
    dem: dict
        dem dictionary

    Returns
    -------
    Line : dict
        Line['Name'] : list of lines names
        Line['Coord'] : np array of the coords of points in lines
        Line['Start'] : list of starting index of each line in Coord
        Line['Length'] : list of length of each line in Coord
    """

    log.debug('Reading avalanche path : %s ', str(fname))
    header = dem['header']
    rasterDEM = dem['rasterData']
    Line = SHP2Array(fname, defname)
    coordx = Line['x']
    coordy = Line['y']
    for i in range(len(coordx)):
        Lx = (coordx[i] - header['xllcenter']) / header['cellsize']
        Ly = (coordy[i] - header['yllcenter']) / header['cellsize']
        if ((Ly < 0) or (Ly > header['nrows']-1) or (Lx < 0) or (Lx > header['ncols']-1)):
            raise ValueError('The avalanche path exceeds dem extent. Try with another path')
        elif np.isnan(rasterDEM[int(np.floor(Ly)), int(np.floor(Lx))]):
            raise ValueError('Nan Value encountered. Try with another path')
    return Line


def readPoints(fname, dem):
    """ Read points from  .shp
    Use SHP2Array to read the shape file.
    Check if the points are laying inside the dem extend

    Parameters
    ----------
    fname: str
        path to shape file
    defname: str
        name to give to the points in the shape file
    dem: dict
        dem dictionary

    Returns
    -------
    Line : dict
        Line['Name'] : list of points names
        Line['Coord'] : np array of the coords of points in points
        Line['Start'] : list of starting index of each point in Coord
        Line['Length'] : list of length of each point in Coord
    """

    log.debug('Reading split point : %s ', str(fname))
    header = dem['header']
    rasterDEM = dem['rasterData']
    defname = 'SHP'
    Points = SHP2Array(fname, defname)
    Pointx = Points['x']
    Pointy = Points['y']
    for i in range(len(Pointx)):
        Lx = (Pointx[i] - header['xllcenter']) / header['cellsize']
        Ly = (Pointy[i] - header['yllcenter']) / header['cellsize']
        if (Ly < 0 or Ly > header['nrows']-1 or Lx < 0 or Lx > header['ncols']-1):
            raise ValueError('The split point is not on the dem. Try with another split point')
        elif np.isnan(rasterDEM[int(np.floor(Ly)), int(np.floor(Lx))]):
            raise ValueError('Nan Value encountered. Try with another split point')
    return Points


def removeFeature(featureIn, nFeature2Remove):
    """ Remove feature number nFeature2Remove from featureIn

    Parameters
    ----------
    featureIn: dict
        shape file dicionary (structure produced by SHP2Array, readLine or readPoint)
    nFeature2Remove: int
        index of feature to remove from featureIn

    Returns
    -------
    featureOut : dict
        shape file dicionary without feature nFeature2Remove
    """
    StartRel = featureIn['Start']
    LengthRel = featureIn['Length']
    thickness = featureIn['thickness']
    featureOut = copy.deepcopy(featureIn)
    start = StartRel[nFeature2Remove]
    end = start + LengthRel[nFeature2Remove]
    # remove feature
    featureOut['x'] = np.delete(featureIn['x'], np.arange(int(start), int(end)))
    featureOut['y'] = np.delete(featureIn['y'], np.arange(int(start), int(end)))
    if 'z' in featureIn.keys():
        featureOut['z'] = np.delete(featureIn['z'], np.arange(int(start), int(end)))

    del featureOut['Name'][nFeature2Remove]
    StartRel = featureOut['Start']
    StartRel[nFeature2Remove:] = StartRel[nFeature2Remove:] - LengthRel[nFeature2Remove]
    featureOut['Start'] = np.delete(StartRel, nFeature2Remove)
    featureOut['Length'] = np.delete(LengthRel, nFeature2Remove)
    featureOut['thickness'] = np.delete(thickness, nFeature2Remove)

    return featureOut


def extractFeature(featureIn, nFeature2Extract):
    """ Extract feature nFeature2Extract from featureIn

    Parameters
    ----------
    featureIn: dict
        shape file dicionary (structure produced by SHP2Array, readLine or readPoint)
    nFeature2Extract: int
        index of feature to extract from featureIn

    Returns
    -------
    featureOut : dict
        shape file dicionary with feature nFeature2Extract
    """
    NameRel = featureIn['Name']
    StartRel = featureIn['Start']
    LengthRel = featureIn['Length']
    thickness = featureIn['thickness']
    featureOut = copy.deepcopy(featureIn)
    # extract feature
    featureOut['Name'] = [NameRel[nFeature2Extract]]
    featureOut['Start'] = np.array([0])
    featureOut['Length'] = np.array([LengthRel[nFeature2Extract]])
    featureOut['thickness'] = np.array([thickness[nFeature2Extract]])
    start = StartRel[nFeature2Extract]
    end = start + LengthRel[nFeature2Extract]
    featureOut['x'] = featureIn['x'][int(start):int(end)]
    featureOut['y'] = featureIn['y'][int(start):int(end)]
    if 'z' in featureIn.keys():
        featureOut['z'] = featureIn['z'][int(start):int(end)]

    return featureOut


def writeLine2SHPfile(lineDict, lineName, fileName):
    """ write a line to shapefile

    Parameters
    ----------
    lineDict: dict
        line dict
    lineName: str
        line name
    fileName: str or pathlib path
        path where the line will be saved
    Returns
    -------
    fileName : str
        path where the line has been saved
    """
    fileName = str(fileName)
    line = np.zeros((np.size(lineDict['x']), 2))
    line[:, 0] = lineDict['x']
    line[:, 1] = lineDict['y']
    w = shapefile.Writer(fileName)
    w.field('name', 'C')
    w.line([line])
    w.record(lineName)
    w.close()
    return fileName


def writePoint2SHPfile(pointDict, pointName, fileName):
    """ write a point to shapefile

    Parameters
    ----------
    pointDict: dict
        point dict
    pointName: str
        point name
    fileName: str or pathlib path
        path where the point will be saved
    Returns
    -------
    fileName : str
        path where the point has been saved
    """
    fileName = str(fileName)
    w = shapefile.Writer(fileName)
    w.field('name', 'C')
    w.point(pointDict['x'], pointDict['y'])
    w.record(pointName)
    w.close()
    return fileName
