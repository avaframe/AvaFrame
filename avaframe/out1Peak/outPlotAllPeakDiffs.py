"""
Functions to plot 2D simulation results and area indicators

"""

import logging
import numpy as np
from matplotlib import pyplot as plt
import pathlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import avaframe.out3Plot.plotUtils as pU
import geopandas as gpd
import copy
from cmcrameri import cm as cmapCrameri

# create local logger
log = logging.getLogger(__name__)


def computeAreaDiff(
    refRaster,
    simRaster,
    thresholdValueReference,
    thresholdValueSimulation,
    dem,
    cropToArea=None,
):
    """compute difference between two rasters based on a thresholdValue

    Parameters
    ------------
    refRaster: numpy ndarray
        reference raster
    simRaster: numpy ndarray
        simulation result raster
    thresholdValueReference: float
        threshold value to mask reference raster
    thresholdValueSimlation: float
        threshold value to mask simulation raster
    dem: dict
        dictionary with info on DEM used for sims, here areaRaster is required to compute actual areas of indicators
    cropToArea: None or numpy ndarray
        None or raster that is used to define where comparison is computed (0 -no, 1-yes)

    Returns
    ---------
    refMask: numpy ndarray
        array with 1 indicating reference data found exceeding threshold and 0 not found
    compRaster: numpy ndarray
        array with 1 indicating simulation data found exceeding threshold and 0 not found
    indicatorDict: dict
        dictionary with info on true positive, false positive, false negative
                    arrays with 1 indicating true positive, false positive, false negative areas and 0 all other entries
                    number of cells where tp, fp, np
                    total area for tp, fp, fn
    """

    refMask = copy.deepcopy(refRaster)
    # set all nans to 0 values
    refMask = np.where(np.isnan(refMask), 0, refMask)
    # set to 0 when values <= threshold, and else if threshold exceeded
    refMask = np.where(refMask <= thresholdValueReference, 0, refMask)
    refMask = np.where(refMask > thresholdValueReference, 1, refMask)

    # set all nans to 0 values
    compRasterMask = copy.deepcopy(simRaster)
    # set to 0 when values <= threshold, and else if threshold exceeded
    compRasterMask = np.where(np.isnan(compRasterMask), 0, compRasterMask)
    compRasterMask = np.where(compRasterMask <= thresholdValueSimulation, 0, compRasterMask)
    compRasterMask = np.where(compRasterMask > thresholdValueSimulation, 1, compRasterMask)

    # if array to crop data is provided, perform analysis only for this area
    if isinstance(cropToArea, np.ndarray):
        refMask = np.where(cropToArea == 0, 0, refMask)
        compRasterMask = np.where(cropToArea == 0, 0, compRasterMask)

    # create true positive, false positive, false negative areas
    tp = np.where((compRasterMask == 1) & (refMask == 1), 1, 0)
    fp = np.where((compRasterMask == 1) & (refMask == 0), 1, 0)
    fn = np.where((compRasterMask == 0) & (refMask == 1), 1, 0)

    indicatorDict = {
        "truePositive": {
            "array": tp,
            "nCells": np.nansum(tp),
            "areaSum": np.nansum(tp * dem["areaRaster"]),
        },
        "falsePositive": {
            "array": fp,
            "nCells": np.nansum(fp),
            "areaSum": np.nansum(fp * dem["areaRaster"]),
        },
        "falseNegative": {
            "array": fn,
            "nCells": np.nansum(fn),
            "areaSum": np.nansum(fn * dem["areaRaster"]),
        },
    }

    return refMask, compRasterMask, indicatorDict


def plotAreaDiff(
    refRaster,
    refMask,
    simRaster,
    compRasterMask,
    resType,
    rasterHeader,
    thresholdValueSimulation,
    outDir,
    indicatorDict,
    simName,
    cropFile=None,
):
    """plot comparison of reference and simulation and area indicators

    Parameters
    -----------
    refRaster, simRaster: numpy ndarray
        reference raster, simulation raster
    refMask, compRasterMask: numpy ndarray
        array with 1 indicating reference (simulation) data found exceeding threshold and 0 not found
    resType: str
        result type of simulation
    rasterHeader: dict
        info on extent
    thresholdValueSimulation: float
        threshold value to mask simulation raster
    outDir: pathlib path
        path to save plot
    indicatorDict: dict
        dictionary with info on true positive, false positive, false negative
                    arrays with 1 indicating true positive, false positive, false negative areas and 0 all other entries
                    number of cells where tp, fp, np
                    total area for tp, fp, fn
    simName: str
        name of simulation
    cropFile: None or pathlib path
        None or path to cropfile with polygon used to perfrom analysis

    """

    data = compRasterMask - refMask
    # mask where both don't show values
    data = np.ma.masked_where((refMask == 0) & (compRasterMask == 0), data)

    # create figure
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(pU.figW * 3, pU.figH))
    # fetch extent of domain
    extentCellCenters, extentCellCorners = pU.createExtentMinMax(
        refRaster, rasterHeader, originLLCenter=True
    )
    cmapTF, _, ticks, normTF = pU.makeColorMap(cmapCrameri.hawaii, -0.5, 0.5, continuous=pU.contCmap)
    cmapTF.set_under(color="b")
    cmapTF.set_over(color="r")
    im1 = ax[0].imshow(
        data,
        cmap=cmapTF,
        norm=normTF,
        extent=extentCellCorners,
        origin="lower",
        aspect="equal",
        zorder=4,
    )
    # add area indicator info
    ax[0].text(
        0.01,
        0.85,
        (
            "false negative:\n %.0f cells, %.2f m2"
            % (
                indicatorDict["falseNegative"]["nCells"],
                indicatorDict["falseNegative"]["areaSum"],
            )
        ),
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax[0].transAxes,
        color="royalblue",
        alpha=0.75,
        zorder=5,
    )
    ax[0].text(
        0.01,
        0.95,
        (
            "false positive:\n %.0f cells, %.2f m2"
            % (
                indicatorDict["falsePositive"]["nCells"],
                indicatorDict["falsePositive"]["areaSum"],
            )
        ),
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax[0].transAxes,
        color="darkred",
        alpha=0.75,
        zorder=5,
    )
    ax[0].text(
        0.01,
        0.75,
        (
            "true positive:\n %.0f cells, %2.f m2"
            % (
                indicatorDict["truePositive"]["nCells"],
                indicatorDict["truePositive"]["areaSum"],
            )
        ),
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax[0].transAxes,
        color="darkgreen",
        alpha=0.75,
        zorder=5,
    )
    title = str(
        "Affected area based on %s %.2f %s"
        % (resType, thresholdValueSimulation, pU.cfgPlotUtils["unit" + resType])
    )

    if isinstance(cropFile, pathlib.Path):
        focus = gpd.read_file(cropFile)
        focus.plot(
            ax=ax[0],
            zorder=12,
            edgecolor="darkgrey",
            linewidth=2,
            facecolor="none",
            alpha=1.0,
        )

    # choose colormap for simulation result
    cmap1, col1, ticks1, norm1 = pU.makeColorMap(
        pU.colorMaps[resType],
        np.nanmin(refRaster),
        np.nanmax(refRaster),
        continuous=True,
    )
    cmap1.set_bad(alpha=0)
    im2 = ax[1].imshow(
        np.ma.masked_where(simRaster == 0.0, simRaster),
        cmap=cmap1,
        extent=extentCellCorners,
        origin="lower",
        aspect="equal",
    )
    # plot reference raster
    ax[2].imshow(refRaster, extent=extentCellCorners, origin="lower", aspect="equal")

    # add colorbar
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im2, cax=cax)
    cbar.set_label(resType)

    # add labels
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[0].set_title(title)
    ax[1].set_xlabel("x [m]")
    ax[1].set_ylabel("y [m]")
    ax[1].set_title("Simulation %s" % resType)
    ax[2].set_xlabel("x [m]")
    ax[2].set_ylabel("y [m]")
    ax[2].set_title("Reference mask")

    pU.saveAndOrPlot({"pathResult": outDir}, "areaIndicatorAnalysis_%s" % simName, fig)
