"""Tests for module com2AB"""
import avaframe.in2Trans.ascUtils as IOf
import pathlib
import pytest


def test_readASCheader(capfd):
    '''Simple test for module readASCheader'''
    dirname = pathlib.Path(__file__).parents[0]
    DGMSource = dirname / '..'/'data'/'avaSlide'/'Inputs'/'slideTopo.asc'
    header = IOf.readASCheader(DGMSource)
    assert((header['ncols'] == 419) and (header['nrows'] == 201) and
           (header['cellsize'] == 5))

    with pytest.raises(Exception) as e:
        DGMSource = dirname /'data'/'testReadAsc.asc'
        assert IOf.readASCheader(DGMSource)
    assert str(e.value) == 'DEM header is not in correct format - check line 4: should be yllcenter or yllcorner'


def test_isEqualASCheader(capfd):
    '''Simple test for module readASCheader'''
    dirname = pathlib.Path(__file__).parents[0]
    DGMSource = dirname / '..'/'data'/'avaSlide'/'Inputs'/'slideTopo.asc'
    headerA = IOf.readASCheader(DGMSource)
    headerB = IOf.readASCheader(DGMSource)
    equal = IOf.isEqualASCheader(headerA, headerB)
    assert equal


def test_readASCdata2numpyArray(capfd):
    '''Simple test for module readASCheader'''
    dirname = pathlib.Path(__file__).parents[0]
    DGMSource = dirname / '..'/'data'/'avaSlide'/'Inputs'/'slideTopo.asc'
    data = IOf.readASCdata2numpyArray(DGMSource)

    assert((data[0][0] == 1752.60) and (data[2][1] == 1749.10)
           and (data[0][3] == 1742.10))
