"""Tests for plotUtils"""

import avaframe.out3Plot.plotUtils as pU
import pathlib
import datetime
import pytest
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm as cmapCameri
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.rasterUtils as IOf
from avaframe.in2Trans.rasterUtils import transformFromASCHeader
import rasterio


def test_constrainPlotsToData():
    """test constraining plots to data"""

    # setup required input
    inputData = np.zeros((10, 12))
    inputData[4:7, 3:9] = 4.0
    cellSize = 100

    # call function to be tested
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(inputData, cellSize, extentOption=False)

    #    print('rows', rowsMin, rowsMax)
    #    print('cols', colsMin, colsMax)

    assert rowsMin == 3
    assert rowsMax == 7
    assert colsMin == 2
    assert colsMax == 9

    # call function to be tested
    rowsMinPlot, rowsMaxPlot, colsMinPlot, colsMaxPlot, dataConstrained = pU.constrainPlotsToData(
        inputData, cellSize, extentOption=True, constrainedData=True
    )

    #    print('rows', rowsMinPlot, rowsMaxPlot)
    #    print('cols', colsMinPlot, colsMaxPlot)
    #    print('dataConstrained', dataConstrained)
    dataConstrainedTest = np.zeros((5, 8))
    dataConstrainedTest[1:4, 1:7] = 4.0

    assert rowsMinPlot == 300
    assert rowsMaxPlot == 700
    assert colsMinPlot == 200
    assert colsMaxPlot == 900
    assert np.array_equal(dataConstrained, dataConstrainedTest)

    # call function to be tested
    inputData = np.zeros((10, 12))
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(inputData, cellSize, extentOption=False)

    assert rowsMin == 0
    assert rowsMax == inputData.shape[0] - 1
    assert colsMin == 0
    assert colsMax == inputData.shape[1] - 1


def test_putAvaNameOnPlot(tmp_path):
    """test put avaName on a plot"""

    # setup required input
    avaDir = pathlib.Path(tmp_path, "avaTest")
    avaDir2 = pathlib.Path(tmp_path, "avaTest2")
    avaDirL = [avaDir, avaDir2]
    plt.figure()
    fig, ax = plt.subplots()

    # call function to be tested
    infoText = pU.putAvaNameOnPlot(ax, avaDir)
    infoText2 = pU.putAvaNameOnPlot(ax, avaDirL)
    #    print('info', infoText)
    #    print('info2', infoText2)

    assert infoText == (datetime.datetime.now().strftime("%d.%m.%y") + "; " + avaDir.stem)
    assert infoText2 == (datetime.datetime.now().strftime("%d.%m.%y") + ";" + "avaTest;avaTest2")

    # call function to be tested
    infoText = pU.putAvaNameOnPlot(ax, avaDir, date=False)
    infoText2 = pU.putAvaNameOnPlot(ax, avaDirL, date=False)
    #    print('info', infoText)
    #    print('info2', infoText2)

    assert infoText == (avaDir.stem)
    assert infoText2 == ("avaTest;avaTest2")


def test_makeColorMap():
    """test generating a customized colormap"""

    # setup required input
    colormapDict = cmapCameri.lapaz
    levMin = 1.0
    levMax = 11.0

    # call function to be tested
    cmap, colorsNew, levelsNew, norm = pU.makeColorMap(colormapDict, levMin, levMax, continuous=False)

    assert colorsNew == None
    assert levelsNew == None
    assert cmap == cmapCameri.lapaz

    # call function to be tested
    levMax = 400.0
    colormapDict3 = {
        "cmap": cmapCameri.hawaii.reversed(),
        "colors": ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"],
        "levels": [1.0, 10.0, 25.0, 50.0],
    }
    cmap2, colorsNew2, levelsNew2, norm2 = pU.makeColorMap(colormapDict3, levMin, levMax, continuous=False)
    assert colorsNew2 == ["#B0F4FA", "#75C165", "#A96C00", "#8B0069", "#8B0069"]
    assert levelsNew2 == [1.0, 10.0, 25.0, 50.0, 400]
    assert cmap2.colors == ["#B0F4FA", "#75C165", "#A96C00", "#8B0069", "#8B0069"]

    # call function to be tested
    colormapDict3 = {
        "cmap": cmapCameri.hawaii.reversed(),
        "colors": ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"],
        "levels": [1.0, 10.0, 25.0, 50.0, 100.0],
    }
    with pytest.raises(AssertionError) as e:
        assert pU.makeColorMap(colormapDict3, levMin, levMax, continuous=False)
    assert str(e.value) == "Number of levels is not allowed to exceed number of colors"

    # call function to be tested
    levMax = 400.0
    colormapDict4 = {"cmap": cmapCameri.hawaii.reversed(), "levels": [1.0, 10.0, 25.0, 50.0]}
    cmap4, colorsNew4, levelsNew4, norm4 = pU.makeColorMap(colormapDict4, levMin, levMax, continuous=False)

    assert np.array_equal(colorsNew4[0, :], np.asarray([0.703779, 0.948977, 0.993775, 1.0]))
    assert np.allclose(colorsNew4[1, :], np.asarray([0.638645, 0.946543, 0.974739, 1.0]), atol=1.0e-6)
    assert np.allclose(colorsNew4[2, :], np.asarray([0.529715, 0.938416, 0.932588, 1.0]), atol=1.0e-6)
    assert np.allclose(colorsNew4[3, :], np.asarray([0.398074, 0.911369, 0.826627, 1.0]), atol=1.0e-6)
    assert levelsNew4 == [1.0, 10.0, 25.0, 50.0, 400]

    # call function to be tested
    levMax = 400.0
    colormapDict5 = {
        "cmap": cmapCameri.hawaii.reversed(),
        "colors": ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"],
    }
    cmap5, colorsNew5, levelsNew5, norm5 = pU.makeColorMap(colormapDict5, levMin, levMax, continuous=False)
    assert levelsNew5 == [1.0, 100.75, 200.5, 300.25, 400]

    # call function to be tested
    levMax = 400.0
    colormapDict5 = {"cmap": cmapCameri.hawaii.reversed()}
    cmap5, colorsNew5, levelsNew5, norm5 = pU.makeColorMap(colormapDict5, levMin, levMax, continuous=False)
    assert levelsNew5 == list(np.linspace(1, 400, 6))

    # call function to be tested
    levMax = 400.0
    colormapDict6 = {
        "colors": ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"],
        "levels": [1.0, 10.0, 25.0, 50.0],
    }
    cmap6, colorsNew6, levelsNew6, norm6 = pU.makeColorMap(colormapDict6, levMin, levMax, continuous=True)
    assert colorsNew6 == ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"]
    assert levelsNew6 == None

    # call function to be tested
    colormapDict6 = {"levels": [1.0, 10.0, 25.0, 50.0]}
    with pytest.raises(FileNotFoundError) as e:
        assert pU.makeColorMap(colormapDict6, levMin, levMax, continuous=False)
    assert str(e.value) == "A `colors` list or a `cmap` is required to create the colormap"

    # call function to be tested
    colormapDict6 = {"levels": [1.0, 10.0, 25.0, 50.0]}
    with pytest.raises(FileNotFoundError) as e:
        assert pU.makeColorMap(colormapDict6, levMin, levMax, continuous=True)
    assert str(e.value) == "You need a `colors` list or a `cmap` to be able to create the colormap"


def test_getColors4Scatter(tmp_path):
    """test put avaName on a plot"""

    # no parameter variation, feed a none to the function
    values = None
    nSamples = np.size(values)
    unitSC = ""
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(
        values, nSamples, unitSC
    )
    #    print('colorSC', colorSC)
    #    print('displayColorBar', displayColorBar)

    assert np.allclose(colorSC, np.array([0.5, 0.5, 0.5, 0.5, 0.5]))
    assert displayColorBar is False

    # parameter variation, feed stings
    values = ["a", "b", "c", "c", "b", "b"]
    nSamples = np.size(values)
    unitSC = ""
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(
        values, nSamples, unitSC
    )
    #    print('itemsList', itemsList)
    #    print('colorSC', colorSC)
    #    print('ticksSC', ticksSC)
    #    print('normSC', normSC)
    #    print('displayColorBar', displayColorBar)

    assert np.allclose(colorSC, np.array([1.0, 2.0, 3.0, 3.0, 2.0, 2.0]))
    assert itemsList == ["a", "b", "c"]
    assert np.allclose(ticksSC, np.array([1.0, 2.0, 3.0]))
    assert displayColorBar is True

    # parameter variation, feed floats
    values = np.array([1, 2, 3, -1, 0, 5])
    nSamples = np.size(values)
    unitSC = "m"
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(
        values, nSamples, unitSC
    )
    #    print('itemsList', itemsList)
    #    print('colorSC', colorSC)
    #    print('ticksSC', ticksSC)
    #    print('normSC', normSC)
    #    print('displayColorBar', displayColorBar)

    assert np.allclose(colorSC, np.array([1, 2, 3, -1, 0, 5]))
    assert itemsList == ""
    assert ticksSC is None
    assert displayColorBar is True


def test_addHillShadeContours(tmp_path):

    # set up required inputs
    testDir = pathlib.Path(tmp_path, "testsPlots")
    fU.makeADir(testDir)

    c = np.array([1, 2, 3, 4])
    field1 = np.tile(c, (4, 1))
    for row in range(field1.shape[0]):
        field1[row, :] = field1[row, :] * (row + 1)
    header = {"cellsize": 1, "xllcenter": 0, "yllcenter": 0, "nrows": 4, "ncols": 4, "nodata_value": -9999}
    header["driver"] = "AAIGrid"
    transform = transformFromASCHeader(header)
    header["transform"] = transform
    header["crs"] = rasterio.crs.CRS()

    outFileName = testDir / "testDEM.asc"
    IOf.writeResultToRaster(header, field1, outFileName, flip=False)

    field2 = np.zeros((4, 4)) + 10.0
    outFileName2 = testDir / "testData.asc"
    IOf.writeResultToRaster(header, field2, outFileName2, flip=False)

    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(field2, header["cellsize"])
    xllcenter = header["xllcenter"]
    yllcenter = header["yllcenter"]
    rowsMinPlot = rowsMin * header["cellsize"] + yllcenter
    rowsMaxPlot = (rowsMax) * header["cellsize"] + yllcenter
    colsMinPlot = colsMin * header["cellsize"] + xllcenter
    colsMaxPlot = (colsMax) * header["cellsize"] + xllcenter
    extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]
    #    print('EXT', extent)
    fig, ax = plt.subplots()
    ls, CS = pU.addHillShadeContours(
        ax, field1, header["cellsize"], extent, colors=["gray"], onlyContours=False
    )

    assert CS != None

    # call function to be tested
    xmin = header["xllcenter"]
    xmax = header["xllcenter"] + header["cellsize"] * (header["ncols"])
    ymin = header["yllcenter"]
    ymax = header["yllcenter"] + header["cellsize"] * (header["nrows"])
    extent2 = [xmin, xmax, ymin, ymax]
    #    print('extent2 is', extent2)
    with pytest.raises(AssertionError) as e:
        assert pU.addHillShadeContours(
            ax, field1, header["cellsize"], extent2, colors=["gray"], onlyContours=False
        )
    assert "Extent of dem data for hillshade and provided extent" in str(e.value)


def test_checkMeshgridInputs():
    """test if inputs for meshgrid are correct shapewise"""

    # setup required inputs
    data = np.zeros((8, 10))
    cellSize = 5
    rowsMin, rowsMax, colsMin, colsMax, dataConstrained = pU.constrainPlotsToData(
        data, cellSize, extentOption=False, constrainedData=True
    )
    header = {"cellsize": 5, "xllcenter": 0, "yllcenter": 0, "nrows": 8, "ncols": 10, "nodata_value": -9999}
    rowsMinPlot = rowsMin * cellSize + header["yllcenter"]
    rowsMaxPlot = (rowsMax) * cellSize + header["yllcenter"]
    colsMinPlot = colsMin * cellSize + header["xllcenter"]
    colsMaxPlot = (colsMax) * cellSize + header["xllcenter"]
    # create meshgrid for contour plot also constrained to where there is data
    xx = np.linspace(colsMinPlot, colsMaxPlot, data.shape[1])
    yy = np.linspace(rowsMinPlot, rowsMaxPlot, data.shape[0])

    pU.checkMeshgridInputs(xx, yy, cellSize, data, plotName="")

    x2 = np.linspace(colsMinPlot, colsMaxPlot + 1, data.shape[1])
    y2 = np.linspace(rowsMinPlot, rowsMaxPlot + 1, data.shape[0])
    with pytest.raises(AssertionError) as e:
        assert pU.checkMeshgridInputs(x2, y2, cellSize, data, plotName="")
    assert (
            "Meshgrid coordinates created have a different spacing between points than cell centers of the corresponding data "
            in str(e.value)
    )
