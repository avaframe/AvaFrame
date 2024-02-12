"""
    ASCII file reader and handler

"""

import numpy as np
import logging
import pathlib

# local imports
import avaframe.com1DFA.DFAtools as DFAtls
from avaframe.in3Utils import geoTrans as gT
import avaframe.in2Trans.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)


def convertDepthToThickness(depthDict, demDict):
    """convert depthField to thicknessField, using a DEM to compute the slope angle required for transformation
    also writes field to a new directory called transformed where depthFile is located
    Parameters
    -----------
    demDict: dict
        dictionary with dem header and rasterData (numpy nd array of z values)
    depthDict: dict
        dictionary with depthField header and rasterData (numpy nd array of depth values)

    Returns
    --------
    thicknessField: dict
        dictionary with header and thickness field numpy array as rasterData key

    """
    # get normal vector of the grid mesh
    demDict = gT.getNormalMesh(demDict)
    _, _, NzNormed = DFAtls.normalize(demDict["Nx"], demDict["Ny"], demDict["Nz"])

    # if resType field cellSize/extent is different to DEM reproject raster on a grid of shape DEM
    depthRasterNew, demData = gT.resizeData(depthDict, demDict)
    demDict["header"]["nodata_value"] = depthDict["header"]["nodata_value"]

    # multiply depth with cos(slopeAngle)
    thickness = depthRasterNew * NzNormed
    slopeAngleField = np.rad2deg(np.arccos(NzNormed))

    # create thickness dict
    thicknessDict = {'header': demDict['header'], 'rasterData': thickness}

    return thicknessDict, depthRasterNew, slopeAngleField


def fetchPointValuesFromField(dataDF, xyPoints, resType, interpMethod='bilinear'):
    """ derive field values at xyPoints using a interpMethod (options: nearest and bilinear)

        Parameters
        -----------
        dataDF: pandas dataFrame
            dataframe with info on simulations (including result file paths,configuration)
        xyPoints: dict
            dictionary with keys x, y and point of interest coordinates
        resType: str
            name of result type (pfv, pft, ...)
        interpMethod: str
            interpolation method to derive values of field at xyPoints

        Returns
        ---------
        dataDF: pandas dataFrame
            updated pandas dataFrame with new column pointValues_resType providing value of resType at xyPoint

    """

    for index, row in dataDF.iterrows():

        # read field
        field = IOf.readRaster(row[resType])

        value, _ = gT.projectOnRaster(field, xyPoints, interp=interpMethod, inData="rasterData", outData="value")
        dataDF.loc[index, ('pointValues_%s' % resType)] = value['value'][0]

    return dataDF


