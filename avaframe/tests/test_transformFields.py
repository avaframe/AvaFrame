"""
    Pytest for transforming fields
"""

#  Load modules
import numpy as np
import logging
import pytest
import configparser
import pathlib


import avaframe.in2Trans.transformFields as tF


def test_convertDepthThickness():
    """ test transforming fields """

    demHeader = {}
    demHeader["xllcenter"] = 1.0
    demHeader["yllcenter"] = 2.0
    demHeader["cellsize"] = 1.0
    demHeader["nodata_value"] = -9999
    demHeader["nrows"] = 4
    demHeader["ncols"] = 7

    # Compute coordinate grid
    xv = np.arange(1.0, 8., 1)
    yv = np.arange(2, 6, 1)
    x, y = np.meshgrid(xv, yv)

    # Set surface elevation from slope and max. elevation
    zv = 40 - np.tan(np.radians(34.)) * x
    demDict = {"header": demHeader, "rasterData": zv}

    print('zv', zv)
    print('x', xv)
    print('y', yv)

    inputDict = {'header': demHeader, 'rasterData': np.ones((4,7))}

    outputDict, inputRasterNew, slopeAngleField = tF.convertDepthThickness(inputDict, demDict, typeOfInput='depth')

    depthValue = (1. * np.cos(np.deg2rad(34.)))
    print('depth', depthValue)
    assert np.isclose(outputDict['rasterData'][2,2], depthValue)

    outputDict, inputRasterNew, slopeAngleField = tF.convertDepthThickness(inputDict, demDict, typeOfInput='thickness')

    thicknessValue = (1. / np.cos(np.deg2rad(34.)))
    print('depth', thicknessValue)
    assert np.isclose(outputDict['rasterData'][2, 2], thicknessValue)

    with pytest.raises(AssertionError) as e:
        assert tF.convertDepthThickness(inputDict, demDict, typeOfInput='test')
    assert "Type for thickness/depth conversion is" in str(e.value)





