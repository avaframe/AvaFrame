"""
Functions to plot 2D simulation results: plot of all peak files at once

"""

import logging
import numpy as np
import matplotlib

matplotlib.use("agg")
from matplotlib import pyplot as plt
import pathlib

import avaframe.out3Plot.plotUtils as pU
import avaframe.in1Data.getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)


def plotAllPeakFields(avaDir, cfgFLAGS, modName, demData=""):
    """Plot all peak fields and return dictionary with paths to plots
    with DEM in background

    Parameters
    ----------
    avaDir : str
        path to avalanche directoy
    cfgFLAGS : str
        general configuration, required to define if plots saved to reports directoy
    modName : str
        name of module that has been used to produce data to be plotted
    demData: dictionary
        optional - if not the dem in the avaDir/Inputs folder has been used but a different one

    Returns
    -------
    plotDict : dict
        dictionary with info on plots, like path to plot
    """

    # Load all infos on simulations
    avaDir = pathlib.Path(avaDir)
    inputDir = avaDir / "Outputs" / modName / "peakFiles"
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    if demData == "":
        demFile = gI.getDEMPath(avaDir)
        demData = IOf.readRaster(demFile, noDataToNan=True)
        demDataField = demData["rasterData"]
    else:
        # check if noDataValue is found and if replace with nans for plotting
        demDataField = np.where(
            demData["rasterData"] == demData["header"]["noDataValue"], np.nan, demData["rasterData"]
        )
    demField = demDataField

    # Output directory
    if cfgFLAGS.getboolean("ReportDir"):
        outDir = avaDir / "Outputs" / modName / "reports"
        fU.makeADir(outDir)
    else:
        outDir = avaDir / "Outputs" / "out1Peak"
        fU.makeADir(outDir)

    # Initialise plot dictionary with simulation names
    plotDict = {}
    for sName in peakFilesDF["simName"]:
        plotDict[sName] = {}

    # Loop through peakFiles and generate plot
    for m in range(len(peakFilesDF["names"])):

        # Load names and paths of peakFiles
        name = peakFilesDF["names"][m]
        fileName = peakFilesDF["files"][m]
        avaName = peakFilesDF["avaName"][m]
        resType = peakFilesDF["resType"][m]
        log.debug("now plot %s:" % (fileName))

        plotName = outDir / ("%s.%s" % (name, pU.outputFormat))

        # only produce a plot if it does not already exists
        # make sure to remove the Outputs folder if you want to regenerate the plot
        # this enables to append simulations to an already existing output without regenerating all plots
        if not plotName.is_file():

            # Figure  shows the result parameter data
            fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))

            # fetch cellSize of peak field
            cellSize = peakFilesDF["cellSize"][m]

            # add peak field data now
            ax, rowsMinPlot, colsMinPlot = addConstrainedDataField(fileName, resType, demField, ax, cellSize)

            # add title, labels and ava Info
            title = str("%s" % name)
            ax.set_title(title)
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            pU.putAvaNameOnPlot(ax, avaDir)

            # save and or show figure
            plotPath = pU.saveAndOrPlot({"pathResult": outDir}, plotName.stem, fig)
            plotDict[peakFilesDF["simName"][m]].update({peakFilesDF["resType"][m]: plotPath})

    return plotDict


def addConstrainedDataField(fileName, resType, demField, ax, cellSize, alpha=1.0, setLimits=False, oneColor=''):
    """ find fileName data, constrain data and demField to where there is data,
        create colormap, define extent, add hillshade contours, add to axes
        and add colorbar

        Parameters
        -----------
        fileName: pathlib path
            path to data
        resType: str
            name of result variable type
        demField: numpy ndarray
            array of dem data
        ax: matplotlib axes object
            axes where to add data plot to
        cellSize: float
            cellSize of data
        alpha: float
            from 0 transparent to 1 opaque for plot of constrained data
        setLimits: bool
            if True set limits of constrained data to plot
        oneColor: str
            optional to add a color for a single color for field
            
        Return
        --------
        ax: matplotlib axes object
            axes updated
    """

    # Load data
    raster = IOf.readRaster(fileName, noDataToNan=True)
    data = raster["rasterData"]

    # constrain data to where there is data
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(data, cellSize)
    dataConstrained = data[rowsMin : rowsMax + 1, colsMin : colsMax + 1]
    demConstrained = demField[rowsMin : rowsMax + 1, colsMin : colsMax + 1]

    data = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)
    unit = pU.cfgPlotUtils["unit%s" % resType]

    # Set extent of peak file
    ny = data.shape[0]
    nx = data.shape[1]
    Ly = ny * cellSize
    Lx = nx * cellSize

    # choose colormap
    cmap, col, ticks, norm = pU.makeColorMap(
        pU.colorMaps[resType], np.amin(data), np.amax(data), continuous=pU.contCmap
        )
    cmap.set_bad(alpha=0)
    # uncomment this to set the under value for discrete cmap transparent
    # cmap.set_under(alpha=0)
    # add origin of data to extent and cellssize to rows and columns
    xllcenter = raster["header"]["xllcenter"]
    yllcenter = raster["header"]["yllcenter"]
    rowsMinPlot = rowsMin * cellSize + yllcenter
    rowsMaxPlot = (rowsMax + 1) * cellSize + yllcenter
    colsMinPlot = colsMin * cellSize + xllcenter
    colsMaxPlot = (colsMax + 1) * cellSize + xllcenter
    extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]

    # add DEM hillshade with contour lines
    ls, CS = pU.addHillShadeContours(ax, demConstrained, cellSize, extent)

    # add peak field data
    if oneColor != '':
        dataOneColor = np.where(data > 0.0, np.amax(data)*0.25,np.nan)
        im1 = ax.imshow(dataOneColor, cmap=oneColor, norm=norm, extent=extent, origin="lower", aspect="equal", zorder=4,
                        alpha=alpha)
    else:
        im1 = ax.imshow(data, cmap=cmap, norm=norm, extent=extent, origin="lower", aspect="equal", zorder=4,
            alpha=alpha)
        pU.addColorBar(im1, ax, ticks, unit)

    if setLimits:
        ax.set_xlim([colsMinPlot, colsMaxPlot])
        ax.set_ylim([rowsMinPlot, rowsMaxPlot])

    return ax, rowsMinPlot, colsMinPlot


def plotAllFields(avaDir, inputDir, outDir, unit="", constrainData=True):
    """Plot all fields within given directory and save to outDir

    Parameters
    ----------
    avaDir : str
        path to avalanche directoy
    inputDir : str
        path to input directoy
    outDir : str
        path to directoy where plots shall be saved to
    unit: str
        unit of result type

    """

    # Load all infos on simulations
    inputDir = pathlib.Path(inputDir)
    outDir = pathlib.Path(outDir)
    if outDir.is_dir() is False:
        # create out dir if not already existing
        outDir.mkdir()
    peakFiles = list(inputDir.glob("*.asc"))

    # Loop through peakFiles and generate plot
    for filename in peakFiles:

        # Load data
        raster = IOf.readRaster(filename)
        data = raster["rasterData"]
        data = np.ma.masked_where(data == 0.0, data)
        name = filename.stem

        # get header info for file writing
        header = raster["header"]
        cellSize = header["cellsize"]

        # Set extent of peak file
        ny = data.shape[0]
        nx = data.shape[1]
        Ly = ny * cellSize
        Lx = nx * cellSize

        # Figure  shows the result parameter data
        fig = plt.figure(figsize=(pU.figW, pU.figH))
        fig, ax = plt.subplots()
        # choose colormap
        cmap, _, ticks, norm = pU.makeColorMap(
            pU.cmapNN, np.amin(data), np.amax(data), continuous=pU.contCmap
        )
        cmap.set_bad("w")

        if constrainData:
            rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(data, cellSize)
            dataConstrained = data[rowsMin : rowsMax + 1, colsMin : colsMax + 1]
            data = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)
            rowsMinPlot = rowsMin * cellSize
            rowsMaxPlot = (rowsMax + 1) * cellSize
            colsMinPlot = colsMin * cellSize
            colsMaxPlot = (colsMax + 1) * cellSize
            im1 = ax.imshow(
                data,
                cmap=cmap,
                norm=norm,
                extent=[colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot],
                origin="lower",
                aspect=nx / ny,
            )
        else:
            im1 = ax.imshow(
                data, cmap=cmap, norm=norm, extent=[0, Lx, 0, Ly], origin="lower", aspect=nx / ny
            )

        pU.addColorBar(im1, ax, ticks, unit)

        title = str("%s" % name)
        ax.set_title(title)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")

        plotName = "%s" % (name)

        pU.putAvaNameOnPlot(ax, avaDir)

        # save and or show figure
        plotPath = pU.saveAndOrPlot({"pathResult": outDir}, plotName, fig)
