"""Tests for module com2AB"""
import numpy as np
import pathlib
import pytest

# Local imports
import avaframe.in2Trans.ascUtils as ascUtils
import avaframe.in2Trans.shpConversion as shpConv


def test_SHP2Array(capfd):
    '''Simple test for function SHP2Array'''
    dirname = pathlib.Path(__file__).parents[0]
    shpFileName = dirname / 'data' / 'testShpConv' / 'testLine.shp'
    SHPdata = shpConv.SHP2Array(shpFileName, defname=None)
    Name = SHPdata['Name']
    Start = SHPdata['Start']
    Length = SHPdata['Length']
    Coordx = SHPdata['x']
    Coordy = SHPdata['y']
    Coordz = SHPdata['z']

    # check that we read the proper values
    # check lines name
    atol = 1e-10
    assert Name == ['line1', 'line2']

    # check start index lines
    Sol = np.array([0, 3])
    testRes = np.allclose(Start, Sol, atol=atol)
    assert testRes

    # check length lines
    Sol = np.array([3, 3])
    testRes = np.allclose(Length, Sol, atol=atol)
    assert testRes

    # check line x coord
    Sol = np.array([15.7390148, 28.11881745, 110.48527428, 32.27956687,
                    55.29299685, 173.27281924])
    testRes = np.allclose(Coordx, Sol, atol=atol)
    assert testRes

    # check line y coord
    Sol = np.array([44.12575919, 51.67232464, 65.17502248, 29.69676696,
                    84.07879946, 69.05704811])
    testRes = np.allclose(Coordy, Sol, atol=atol)
    assert testRes

    # check line z coord
    Sol = np.array([0., 0., 0., 0., 0., 0.])
    testRes = np.allclose(Coordz, Sol, atol=atol)
    assert testRes


def test_extractFeature(capfd):
    '''Simple test for function extractFeature'''
    dirname = pathlib.Path(__file__).parents[0]
    shpFileName = dirname / 'data' / 'testShpConv' / 'testLine.shp'
    SHPdata = shpConv.SHP2Array(shpFileName, defname=None)

    atol = 1e-10
    # extract feature 0
    featureOut = shpConv.extractFeature(SHPdata, 0)

    assert featureOut['Name'] == ['line1']

    # check start index lines
    Sol = np.array([0])
    testRes = np.allclose(featureOut['Start'], Sol, atol=atol)
    assert testRes

    # check length lines
    Sol = np.array([3])
    testRes = np.allclose(featureOut['Length'], Sol, atol=atol)
    assert testRes

    # check line x coord
    Sol = np.array([15.7390148, 28.11881745, 110.48527428])
    testRes = np.allclose(featureOut['x'], Sol, atol=atol)
    assert testRes

    # check line y coord
    Sol = np.array([44.12575919, 51.67232464, 65.17502248])
    testRes = np.allclose(featureOut['y'], Sol, atol=atol)
    assert testRes

    # check line z coord
    Sol = np.array([0., 0., 0.])
    testRes = np.allclose(featureOut['z'], Sol, atol=atol)
    assert testRes

    featureOut = shpConv.extractFeature(SHPdata, 1)

    # extract feature 1
    assert featureOut['Name'] == ['line2']

    # check start index lines
    Sol = np.array([0])
    testRes = np.allclose(featureOut['Start'], Sol, atol=atol)
    assert testRes

    # check length lines
    Sol = np.array([3])
    testRes = np.allclose(featureOut['Length'], Sol, atol=atol)
    assert testRes

    # check line x coord
    Sol = np.array([32.27956687, 55.29299685, 173.27281924])
    testRes = np.allclose(featureOut['x'], Sol, atol=atol)
    assert testRes

    # check line y coord
    Sol = np.array([29.69676696, 84.07879946, 69.05704811])
    testRes = np.allclose(featureOut['y'], Sol, atol=atol)
    assert testRes

    # check line z coord
    Sol = np.array([0., 0., 0.])
    testRes = np.allclose(featureOut['z'], Sol, atol=atol)
    assert testRes


def test_removeFeature(capfd):
    '''Simple test for function removeFeature'''
    dirname = pathlib.Path(__file__).parents[0]
    shpFileName = dirname / 'data' / 'testShpConv' / 'testLine.shp'
    SHPdata = shpConv.SHP2Array(shpFileName, defname=None)
    atol = 1e-10
    # extract feature 0
    featureOut = shpConv.removeFeature(SHPdata, 1)

    assert featureOut['Name'] == ['line1']

    # check start index lines
    Sol = np.array([0])
    testRes = np.allclose(featureOut['Start'], Sol, atol=atol)
    assert testRes

    # check length lines
    Sol = np.array([3])
    testRes = np.allclose(featureOut['Length'], Sol, atol=atol)
    assert testRes

    # check line x coord
    Sol = np.array([15.7390148, 28.11881745, 110.48527428])
    testRes = np.allclose(featureOut['x'], Sol, atol=atol)
    assert testRes

    # check line y coord
    Sol = np.array([44.12575919, 51.67232464, 65.17502248])
    testRes = np.allclose(featureOut['y'], Sol, atol=atol)
    assert testRes

    # check line z coord
    Sol = np.array([0., 0., 0.])
    testRes = np.allclose(featureOut['z'], Sol, atol=atol)
    assert testRes

    featureOut = shpConv.removeFeature(SHPdata, 0)
    print(featureOut)
    # extract feature 1
    assert featureOut['Name'] == ['line2']

    # check start index lines
    Sol = np.array([0])
    testRes = np.allclose(featureOut['Start'], Sol, atol=atol)
    assert testRes

    # check length lines
    Sol = np.array([3])
    testRes = np.allclose(featureOut['Length'], Sol, atol=atol)
    assert testRes

    # check line x coord
    Sol = np.array([32.27956687, 55.29299685, 173.27281924])
    testRes = np.allclose(featureOut['x'], Sol, atol=atol)
    assert testRes

    # check line y coord
    Sol = np.array([29.69676696, 84.07879946, 69.05704811])
    testRes = np.allclose(featureOut['y'], Sol, atol=atol)
    assert testRes

    # check line z coord
    Sol = np.array([0., 0., 0.])
    testRes = np.allclose(featureOut['z'], Sol, atol=atol)
    assert testRes


def test_readLine(capfd):
    '''Simple test for function readLine'''
    dirname = pathlib.Path(__file__).parents[0]
    demFileName = dirname / 'data' / 'testShpConv' / 'testShpConv.asc'
    dem = ascUtils.readRaster(demFileName)
    shpFileName = dirname / 'data' / 'testShpConv' / 'testLine.shp'

    # do we react properly when the input line exceeds the dem?
    with pytest.raises(ValueError) as e:
        assert shpConv.readLine(shpFileName, '', dem)
    assert str(e.value) == "Nan Value encountered. Try with another path"

    # do we react properly when the input line exceeds the dem?
    shpFileName = dirname / 'data' / 'testShpConv' / 'testLineOut.shp'
    with pytest.raises(ValueError) as e:
        assert shpConv.readLine(shpFileName, '', dem)
    assert str(e.value) == "The avalanche path exceeds dem extent. Try with another path"

    shpFileName = dirname / 'data' / 'testShpConv' / 'testLineGood.shp'
    Line = shpConv.readLine(shpFileName, '', dem)

    # check lines name
    atol = 1e-10
    assert Line['Name'] == ['goodLine']

    # check start index lines
    Sol = np.array([0])
    testRes = np.allclose(Line['Start'], Sol, atol=atol)
    assert testRes

    # check length lines
    Sol = np.array([3])
    testRes = np.allclose(Line['Length'], Sol, atol=atol)
    assert testRes

    # check line x coord
    Sol = np.array([19.34206385, 35.20773381, 83.14231115])
    testRes = np.allclose(Line['x'], Sol, atol=atol)
    assert testRes

    # check line y coord
    Sol = np.array([83.06609712, 72.43272257, 71.42002023])
    testRes = np.allclose(Line['y'], Sol, atol=atol)
    assert testRes

    # check line z coord
    Sol = np.array([0., 0., 0.])
    testRes = np.allclose(Line['z'], Sol, atol=atol)
    assert testRes


def test_readPoints(capfd):
    '''Simple test for function readPoints'''
    dirname = pathlib.Path(__file__).parents[0]
    demFileName = dirname / 'data' / 'testShpConv' / 'testShpConv.asc'
    dem = ascUtils.readRaster(demFileName)

    # do we react properly when the input point exceeds the dem?
    shpFileName = dirname / 'data' / 'testShpConv' / 'testLine.shp'
    with pytest.raises(ValueError) as e:
        assert shpConv.readPoints(shpFileName, dem)
    assert str(e.value) == 'Nan Value encountered. Try with another split point'

    # do we react properly when the input point exceeds the dem?
    shpFileName = dirname / 'data' / 'testShpConv' / 'testLineOut.shp'
    with pytest.raises(ValueError) as e:
        assert shpConv.readPoints(shpFileName, dem)
    assert str(e.value) == 'The split point is not on the dem. Try with another split point'


def test_writeLine2SHPfile(tmp_path):
    '''Simple test for function writeLine2SHPfile'''
    temp = pathlib.Path(tmp_path, 'testSaveLine')
    linePath = temp / 'line2SHPfile.shp'
    lineDict = {'x': np.array([0, 1]), 'y': np.array([0, 1])}
    fileName = shpConv.writeLine2SHPfile(lineDict, 'myLine', linePath)
    lines = list(temp.glob(('*.shp')))
    print(lines)
    print(fileName)
    print(temp)
    assert str(lines[0]) == str(linePath)


def test_writePoint2SHPfile(tmp_path):
    '''Simple test for function writePoint2SHPfile'''
    temp = pathlib.Path(tmp_path, 'testSaveLine')
    pointpath = temp / 'line2SHPfile.shp'
    pointDict = {'x': np.array([0]), 'y': np.array([1])}
    fileName = shpConv.writePoint2SHPfile(pointDict, 'myPoint', pointpath)
    points = list(temp.glob(('*.shp')))
    print(points)
    print(fileName)
    print(temp)
    assert str(points[0]) == str(pointpath)
