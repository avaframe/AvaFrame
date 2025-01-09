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
