"""Tests for plotUtils"""
import avaframe.out3Plot.plotUtils as pU
import pathlib
import datetime
import pytest
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm as cmapCameri


def test_constrainPlotsToData():
    """ test constraining plots to data """

    # setup required input
    inputData = np.zeros((10,12))
    inputData[4:7,3:9] = 4.
    cellSize = 100

    # call function to be tested
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(inputData, cellSize, extentOption=False)

    print('rows', rowsMin, rowsMax)
    print('cols', colsMin, colsMax)

    assert rowsMin == 3
    assert rowsMax == 7
    assert colsMin == 2
    assert colsMax == 9

    # call function to be tested
    rowsMinPlot, rowsMaxPlot, colsMinPlot, colsMaxPlot, dataConstrained = pU.constrainPlotsToData(
                                                                          inputData, cellSize,
                                                                          extentOption=True,
                                                                          constrainedData=True)

    print('rows', rowsMinPlot, rowsMaxPlot)
    print('cols', colsMinPlot, colsMaxPlot)
    print('dataConstrained', dataConstrained)
    dataConstrainedTest = np.zeros((5, 8))
    dataConstrainedTest[1:4, 1:7] = 4.

    assert rowsMinPlot == 300
    assert rowsMaxPlot == 800
    assert colsMinPlot == 200
    assert colsMaxPlot == 1000
    assert np.array_equal(dataConstrained, dataConstrainedTest)

    # call function to be tested
    inputData = np.zeros((10,12))
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(inputData, cellSize, extentOption=False)

    assert rowsMin == 0
    assert rowsMax == inputData.shape[0]
    assert colsMin == 0
    assert colsMax == inputData.shape[1]


def test_putAvaNameOnPlot(tmp_path):
    """ test put avaName on a plot """

    # setup required input
    avaDir = pathlib.Path(tmp_path, 'avaTest')
    avaDir2 = pathlib.Path(tmp_path, 'avaTest2')
    avaDirL = [avaDir, avaDir2]
    plt.figure()
    fig, ax = plt.subplots()

    # call function to be tested
    infoText = pU.putAvaNameOnPlot(ax, avaDir)
    infoText2 = pU.putAvaNameOnPlot(ax, avaDirL)
    print('info', infoText)
    print('info2', infoText2)

    assert infoText == (datetime.datetime.now().strftime("%d.%m.%y")  + '; ' + avaDir.stem)
    assert infoText2 == (datetime.datetime.now().strftime("%d.%m.%y") + ';' + 'avaTest;avaTest2')


    # call function to be tested
    infoText = pU.putAvaNameOnPlot(ax, avaDir, date=False)
    infoText2 = pU.putAvaNameOnPlot(ax, avaDirL, date=False)
    print('info', infoText)
    print('info2', infoText2)

    assert infoText == (avaDir.stem)
    assert infoText2 == ('avaTest;avaTest2')    


def test_makeColorMap():
    """ test generating a customized colormap """

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
    colormapDict3 = {'cmap': cmapCameri.hawaii.reversed(), 'colors': ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"],
                    'levels': [1.0, 10.0, 25.0, 50.0]}
    cmap2, colorsNew2, levelsNew2, norm2 = pU.makeColorMap(colormapDict3, levMin, levMax, continuous=False)
    assert colorsNew2 == ["#B0F4FA", "#75C165", "#A96C00", "#8B0069", "#8B0069"]
    assert levelsNew2 == [1.0, 10.0, 25.0, 50.0, 400]
    assert cmap2.colors == ["#B0F4FA", "#75C165", "#A96C00", "#8B0069", "#8B0069"]

    # call function to be tested
    colormapDict3 = {'cmap': cmapCameri.hawaii.reversed(), 'colors': ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"],
                    'levels': [1.0, 10.0, 25.0, 50.0, 100.]}
    with pytest.raises(AssertionError) as e:
        assert pU.makeColorMap(colormapDict3, levMin, levMax, continuous=False)
    assert str(e.value) == "Number of levels is not allowed to exceed number of colors"

    # call function to be tested
    levMax = 400.0
    colormapDict4 = {'cmap': cmapCameri.hawaii.reversed(), 'levels': [1.0, 10.0, 25.0, 50.0]}
    cmap4, colorsNew4, levelsNew4, norm4 = pU.makeColorMap(colormapDict4, levMin, levMax, continuous=False)

    assert np.array_equal(colorsNew4[0,:], np.asarray([0.703779, 0.948977, 0.993775, 1.]))
    assert np.allclose(colorsNew4[1,:], np.asarray([0.638645, 0.946543, 0.974739, 1.]), atol=1.e-6)
    assert np.allclose(colorsNew4[2,:], np.asarray([0.529715, 0.938416, 0.932588, 1.]), atol=1.e-6)
    assert np.allclose(colorsNew4[3,:], np.asarray([0.398074, 0.911369, 0.826627, 1.]), atol=1.e-6)
    assert levelsNew4 == [1.0, 10.0, 25.0, 50.0, 400]

    # call function to be tested
    levMax = 400.0
    colormapDict5 = {'cmap': cmapCameri.hawaii.reversed(), 'colors': ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"]}
    cmap5, colorsNew5, levelsNew5, norm5 = pU.makeColorMap(colormapDict5, levMin, levMax, continuous=False)
    assert levelsNew5 == [1.0, 100.75, 200.5, 300.25, 400]


    # call function to be tested
    levMax = 400.0
    colormapDict5 = {'cmap': cmapCameri.hawaii.reversed()}
    cmap5, colorsNew5, levelsNew5, norm5 = pU.makeColorMap(colormapDict5, levMin, levMax, continuous=False)
    assert levelsNew5 == list(np.linspace(1, 400, 6))

    # call function to be tested
    levMax = 400.0
    colormapDict6 = {'colors': ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"],
                    'levels': [1.0, 10.0, 25.0, 50.0]}
    cmap6, colorsNew6, levelsNew6, norm6 = pU.makeColorMap(colormapDict6, levMin, levMax, continuous=True)
    assert colorsNew6 == ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"]
    assert levelsNew6 == None

    # call function to be tested
    colormapDict6 = {'levels': [1.0, 10.0, 25.0, 50.0]}
    with pytest.raises(FileNotFoundError) as e:
        assert pU.makeColorMap(colormapDict6, levMin, levMax, continuous=False)
    assert str(e.value) == "A `colors` list or a `cmap` is required to create the colormap"

    # call function to be tested
    colormapDict6 = {'levels': [1.0, 10.0, 25.0, 50.0]}
    with pytest.raises(FileNotFoundError) as e:
        assert pU.makeColorMap(colormapDict6, levMin, levMax, continuous=True)
    assert str(e.value) == 'You need a `colors` list or a `cmap` to be able to create the colormap'


def test_getColors4Scatter(tmp_path):
    """ test put avaName on a plot """

    # no parameter variation, feed a none to the function
    values = None
    nSamples = np.size(values)
    unitSC = ''
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(values, nSamples, unitSC)
    print('colorSC', colorSC)
    print('displayColorBar', displayColorBar)

    assert np.allclose(colorSC, np.array([0.5, 0.5, 0.5, 0.5, 0.5]))
    assert displayColorBar is False

    # parameter variation, feed stings
    values = ['a', 'b', 'c', 'c', 'b', 'b']
    nSamples = np.size(values)
    unitSC = ''
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(values, nSamples, unitSC)
    print('itemsList', itemsList)
    print('colorSC', colorSC)
    print('ticksSC', ticksSC)
    print('normSC', normSC)
    print('displayColorBar', displayColorBar)

    assert np.allclose(colorSC, np.array([1., 2., 3., 3., 2., 2.]))
    assert itemsList == ['a', 'b', 'c']
    assert np.allclose(ticksSC, np.array([1., 2., 3.]))
    assert displayColorBar is True

    # parameter variation, feed floats
    values = np.array([1, 2, 3, -1, 0, 5])
    nSamples = np.size(values)
    unitSC = 'm'
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(values, nSamples, unitSC)
    print('itemsList', itemsList)
    print('colorSC', colorSC)
    print('ticksSC', ticksSC)
    print('normSC', normSC)
    print('displayColorBar', displayColorBar)

    assert np.allclose(colorSC, np.array([1,  2,  3, -1,  0,  5]))
    assert itemsList == ''
    assert ticksSC is None
    assert displayColorBar is True
