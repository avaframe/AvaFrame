"""Tests for module com2AB"""
import avaframe.in2Trans.rasterUtils as IOf
import pathlib
import pytest


def test_readASCheader(capfd):
    '''Simple test for module readRasterHeader'''
    dirname = pathlib.Path(__file__).parents[0]
    DGMSource = dirname / '..'/'data'/'avaSlide'/'Inputs'/'slideTopo.asc'
    header = IOf.readRasterHeader(DGMSource)
    assert((header['ncols'] == 419) and (header['nrows'] == 201) and
           (header['cellsize'] == 5))

    with pytest.raises(ValueError) as e:
        DGMSource = dirname /'data'/'testReadAsc.asc'
        assert IOf.readRasterHeader(DGMSource)
    assert str(e.value) == ("DEM header is not in correct format - needs to contain values for: cellsize, nrows, ncols, xllcenter(-corner), yllcenter(-corner), nodata_value")


def test_isEqualASCheader(capfd):
    '''Simple test for module readRasterHeader'''
    dirname = pathlib.Path(__file__).parents[0]
    DGMSource = dirname / '..'/'data'/'avaSlide'/'Inputs'/'slideTopo.asc'
    headerA = IOf.readRasterHeader(DGMSource)
    headerB = IOf.readRasterHeader(DGMSource)
    equal = IOf.isEqualASCheader(headerA, headerB)
    assert equal


