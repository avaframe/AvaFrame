"""Tests for module in2Trans rasterUtils"""

import avaframe.in2Trans.rasterUtils as IOf
import pathlib
import pytest


def test_readASCheader(capfd):
    """Simple test for module readRasterHeader"""
    dirname = pathlib.Path(__file__).parents[0]
    DGMSource = dirname / ".." / "data" / "avaSlide" / "Inputs" / "slideTopo.asc"
    header = IOf.readRasterHeader(DGMSource)
    assert (header["ncols"] == 419) and (header["nrows"] == 201) and (header["cellsize"] == 5)


def test_isEqualASCheader(capfd):
    """Simple test for module readRasterHeader"""
    dirname = pathlib.Path(__file__).parents[0]
    DGMSource = dirname / ".." / "data" / "avaSlide" / "Inputs" / "slideTopo.asc"
    headerA = IOf.readRasterHeader(DGMSource)
    headerB = IOf.readRasterHeader(DGMSource)
    equal = IOf.isEqualASCheader(headerA, headerB)
    assert equal


def test_getRasterFileTypeFromHeader():
    """Test getRasterFileTypeFromHeader function with different driver types"""
    dirname = pathlib.Path(__file__).parents[0]

    # Test with AAIGrid (ASCII) file
    ascSource = dirname / ".." / "data" / "avaSlide" / "Inputs" / "slideTopo.asc"
    ascHeader = IOf.readRasterHeader(ascSource)
    fileType = IOf.getRasterFileTypeFromHeader(ascHeader)
    assert fileType == ".asc", f"Expected '.asc' but got '{fileType}'"
    assert ascHeader["driver"] == "AAIGrid", "Expected AAIGrid driver for .asc file"

    # Test with GTiff file
    tifSource = dirname / ".." / "data" / "avaAlr" / "Inputs" / "avaAlr.tif"
    tifHeader = IOf.readRasterHeader(tifSource)
    fileType = IOf.getRasterFileTypeFromHeader(tifHeader)
    assert fileType == ".tif", f"Expected '.tif' but got '{fileType}'"
    assert tifHeader["driver"] == "GTiff", "Expected GTiff driver for .tif file"

    # Test with unsupported driver - should raise AssertionError
    unsupportedHeader = {"driver": "HFA"}
    with pytest.raises(AssertionError, match="Unsupported driver for DEM"):
        IOf.getRasterFileTypeFromHeader(unsupportedHeader)
