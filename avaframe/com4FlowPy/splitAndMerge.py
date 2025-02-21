#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Functions to handle the raster tiles.
     Tiling is intended to manage processing of large(r) computational domains.
"""

import logging
import pickle
import gc
import numpy as np
import avaframe.in2Trans.rasterUtils as IOf

# create local logger
log = logging.getLogger(__name__)


def tileRaster(fNameIn, fNameOut, dirName, xDim, yDim, U, isInit=False):
    """
    divides a raster into tiles and saves the tiles

    Parameters
    -----------
    fNameIn : str
        path to raster that is tiled
    fNameOut: str
        name of saved raster file
    dirName: str
        path to folder, where tiled raster is saved (temp - folder)
    xDim: int
        size of one tile in x dimension (number of raster columns)
    yDim: int
        size of one tile in y dimension (number of raster rows)
    U: int
        size of tile overlapping (number of raster cells)
    isInit: bool
        if isInit is True, edges are assigned to -9999 (default: False)
    """

    # if not os.path.exists(dirName):
    #    os.makedirs(dirName)

    # largeRaster, largeHeader = iof.f_readASC(fNameIn, dType='float')
    largeData = IOf.readRaster(fNameIn, noDataToNan=False)
    largeRaster = largeData["rasterData"]
    # einlesen des Rasters und der Header

    i, j, imax, jmax = 0, 0, 0, 0
    sX, sY, eX, eY = 0, 0, 0, 0
    # starte mit den tiles in der NW-Ecke sX,eX = cols; sY,eY = rows;
    # cs=largeHeader['cellsize']
    # xllc = largeHeader['xllcorner']
    # yllc = largeHeader['yllcorner']

    nrows, ncols = largeRaster.shape[0], largeRaster.shape[1]
    pickle.dump((nrows, ncols), open(dirName / "extentLarge", "wb"))

    # print ncols, nrows

    I, J, IMAX, JMAX = 0, 0, 0, 0

    while eY < nrows:
        eY = sY + yDim
        while eX < ncols:
            eX = sX + xDim

            # rangeRowsCols = ((sY,eY),(sX,eX))
            # pickle.dump(rangeRowsCols, open("%s/ext_%i_%i"%(dirName,i,j),"wb"))

            # headerTile = {}
            # headerTile['ncols'] = eX-sX
            # headerTile['nrows'] = eY-sY
            # headerTile['xllcorner'] = xllc + sX*cs
            # headerTile['yllcorner'] = yllc + nrows*cs - eY*cs
            # headerTile['cellsize'] = cs
            # headerTile['noDataValue'] = largeHeader['noDataValue']

            # pickle.dump( headerTile, open( "temp/header%d_%d.p"%(i,j), "wb" ) )
            # np.save("%s/%s_%i_%i"%(dirName,fNameOut, i, j), largeRaster[sY:eY,sX:eX])
            # pickle.dump(, open( "temp/header_large.p"%(fNameOut, i,j), "wb" ) )
            # log.info("saved %s - TileNr.: %i_%i"%(fNameOut,i,j))

            sX = eX - 2 * U
            JMAX = max(J, JMAX)
            J += 1
        sX, J, eX = 0, 0, 0
        sY = eY - 2 * U
        IMAX = max(I, IMAX)
        I += 1

    sX, sY, eX, eY = 0, 0, 0, 0

    if isInit is False:
        while eY < nrows:
            eY = sY + yDim
            while eX < ncols:
                eX = sX + xDim
                rangeRowsCols = ((sY, eY), (sX, eX))
                pickle.dump(rangeRowsCols, open(dirName / ("ext_%i_%i" % (i, j)), "wb"))

                # headerTile = {}
                # headerTile['ncols'] = eX-sX
                # headerTile['nrows'] = eY-sY
                # headerTile['xllcorner'] = xllc + sX*cs
                # headerTile['yllcorner'] = yllc + nrows*cs - eY*cs
                # headerTile['cellsize'] = cs
                # headerTile['noDataValue'] = largeHeader['noDataValue']

                # pickle.dump( headerTile,
                # open( "temp/header%d_%d.p"%(i,j), "wb" ) )
                np.save(dirName / ("%s_%i_%i" % (fNameOut, i, j)), largeRaster[sY:eY, sX:eX])
                # pickle.dump(,
                # open( "temp/header_large.p"%(fNameOut, i,j), "wb" ) )
                log.info("saved %s - TileNr.: %i_%i", fNameOut, i, j)

                sX = eX - 2 * U
                jmax = max(j, jmax)
                j += 1
            sX, j, eX = 0, 0, 0
            sY = eY - 2 * U
            imax = max(i, imax)
            i += 1
    else:
        while eY < nrows:
            eY = sY + yDim
            while eX < ncols:
                eX = sX + xDim

                rangeRowsCols = ((sY, eY), (sX, eX))
                pickle.dump(rangeRowsCols, open(dirName / ("ext_%i_%i" % (i, j)), "wb"))

                # headerTile = {}
                # headerTile['ncols'] = eX-sX
                # headerTile['nrows'] = eY-sY
                # headerTile['xllcorner'] = xllc + sX*cs
                # headerTile['yllcorner'] = yllc + nrows*cs - eY*cs
                # headerTile['cellsize'] = cs
                # headerTile['noDataValue'] = largeHeader['noDataValue']
                # pickle.dump(headerTile,
                # open( "temp/hd_%s%_d_%d.p"%(fNameOut, i,j), "wb" ) )

                initRas = largeRaster[sY:eY, sX:eX].copy()
                # shapeX = np.shape(initRas)[1]
                # shapeY = np.shape(initRas)[0]
                if j != JMAX:
                    initRas[:, -U:] = -9999  # Rand im Osten
                if i != 0:
                    initRas[0:U, :] = -9999  # Rand im Norden
                if j != 0:
                    initRas[:, 0:U] = -9999  # Rand im Westen
                if i != IMAX:
                    initRas[-U:, :] = -9999  # Rand im Sueden

                # log.info("%i_%i"%(shapeX-U, shapeX))
                np.save(dirName / ("%s_%i_%i" % (fNameOut, i, j)), initRas)
                del initRas
                # pickle.dump(,
                # open( "temp/header_large.p"%(fNameOut, i,j), "wb" ) )
                log.info("saved %s - TileNr.: %i_%i", fNameOut, i, j)

                sX = eX - 2 * U
                jmax = max(j, jmax)
                j += 1
            sX, j, eX = 0, 0, 0
            sY = eY - 2 * U
            imax = max(i, imax)
            i += 1

    pickle.dump((imax, jmax), open(dirName / "nTiles", "wb"))
    log.info("finished tiling %s: nTiles=%s" % (fNameOut, (imax + 1) * (jmax + 1)))
    log.info("----------------------------")
    # del largeRaster, largeHeader
    del largeRaster
    gc.collect()
    # return largeRaster


def mergeRaster(inDirPath, fName, method='max'):
    """
    Merges the results for each tile to one array using the
    method provided through the function parameters

    Parameters
    ----------
    inDirPath: str
        Path to the temporary files, that are results for each tile
    fName : str
        file name of the parameter which should be merged from tile-results
    method : str
        method, how the tiles should be merged (default: max)
        method 'min' calculates the minimum of input raster tiles,
        if the minimum is < 0, then 0 is used
        method 'sum' calculates the sum of the raster tiles

    Returns
    -------
    mergedRas : numpy array
        merged raster
    """

    # os.chdir(inDirPath)

    extL = pickle.load(open(inDirPath / "extentLarge", "rb"))
    # print extL
    nTiles = pickle.load(open(inDirPath / "nTiles", "rb"))

    mergedRas = np.zeros((extL[0], extL[1]))
    # create Raster with original size
    if method != 'sum':
        mergedRas[:, :] = np.nan

    for i in range(nTiles[0] + 1):
        for j in range(nTiles[1] + 1):
            smallRas = np.load(inDirPath / ("%s_%i_%i.npy" % (fName, i, j)))
            # print smallRas
            pos = pickle.load(open(inDirPath / ("ext_%i_%i" % (i, j)), "rb"))
            # print pos

            if method == 'max':
                mergedRas[pos[0][0]: pos[0][1], pos[1][0]: pos[1][1]] = np.fmax(
                    mergedRas[pos[0][0]: pos[0][1], pos[1][0]: pos[1][1]], smallRas
                )
            elif method == 'min':
                mergedRas[pos[0][0]:pos[0][1], pos[1][0]:pos[1][1]] =\
                    np.where((mergedRas[pos[0][0]:pos[0][1], pos[1][0]:pos[1][1]] >= 0) & (smallRas >= 0),
                    np.fmin(mergedRas[pos[0][0]:pos[0][1], pos[1][0]:pos[1][1]], smallRas),
                    np.fmax(mergedRas[pos[0][0]:pos[0][1], pos[1][0]:pos[1][1]], smallRas))
            if method == 'sum':
                mergedRas[pos[0][0]: pos[0][1], pos[1][0]: pos[1][1]] = np.add(
                    mergedRas[pos[0][0]: pos[0][1], pos[1][0]: pos[1][1]], smallRas
                )
            del smallRas
            log.info("appended result %s_%i_%i", fName, i, j)
    return mergedRas
