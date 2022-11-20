"""
    Plot settings for output figures

    This file is part of Avaframe.
"""

import os
import warnings
import seaborn as sns
import copy
import numpy as np
import matplotlib
import datetime
import pathlib
from matplotlib.image import NonUniformImage
from matplotlib import pyplot as plt
from matplotlib import colors as mplCol
import logging
from cmcrameri import cm as cmapCrameri
from matplotlib.colors import LightSource

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.out3Plot import plotUtils
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in1Data.getInput as gI
import avaframe.in3Utils.geoTrans as gT


# create local logger
log = logging.getLogger(__name__)

# read main configuration to get acces to the plotting flags
cfgMain = cfgUtils.getGeneralConfig()
cfgFlags = cfgMain['FLAGS']

# Load all input Parameters from config file
# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(plotUtils)
cfgPlotUtils = cfg['UNITS']
cfg = cfg['MAIN']

# define seaborn style and color maps
sns.set(font_scale=1)
sns.set_style("ticks", {'axes.linewidth': 1, 'axes.edgecolor': 'black',
                        'font.family': [cfg['fontFamily']]})

# define figure dimensions
figW = float(cfg['figW'])
figH = float(cfg['figH'])
# define lines and marker properties
markers = cfg['markerStyle']
ms = float(cfg['markerSize'])
matplotlib.rcParams['lines.linewidth'] = float(cfg['lineWidth'])
matplotlib.rcParams['lines.markersize'] = ms
# font size
fs = float(cfg['fontSize'])
matplotlib.rcParams['figure.titlesize'] = cfg['titleSize']
matplotlib.rcParams['figure.dpi'] = float(cfg['figResolution'])
matplotlib.rcParams['figure.autolayout'] = True
ls = ['-', '--', '-.']
matplotlib.rcParams['axes.titlesize'] = cfg['axesTitleSize']
matplotlib.rcParams['axes.labelsize'] = cfg['labelSize']
matplotlib.rcParams['axes.linewidth'] = 1.0
matplotlib.rcParams['axes.edgecolor'] = 'lightgrey'
matplotlib.rcParams['axes.labelcolor'] = 'grey'
matplotlib.rcParams['xtick.color'] = 'grey'
matplotlib.rcParams['xtick.major.width'] = 1.0
matplotlib.rcParams['xtick.major.size'] = 3
matplotlib.rcParams['xtick.labelsize'] = cfg['tickLabelSize']
matplotlib.rcParams['ytick.color'] = 'grey'
matplotlib.rcParams['ytick.major.width'] = 1.0
matplotlib.rcParams['ytick.major.size'] = 3
matplotlib.rcParams['ytick.labelsize'] = cfg['tickLabelSize']


# set output extension {png, ps, pdf, svg}
outputFormat = cfg['savefigFormat']
matplotlib.rcParams['savefig.format'] = outputFormat
matplotlib.rcParams['savefig.bbox'] = 'tight'

matplotlib.rcParams['legend.edgecolor'] = 'None'
matplotlib.rcParams['legend.fontsize'] = cfg['fontSize']
matplotlib.rcParams['font.family'] = cfg['fontFamily']
matplotlib.rcParams['text.usetex'] = cfg.getboolean('usetex')
matplotlib.rc('text.latex', preamble=r'\usepackage{cmbright}')

matplotlib.rcParams['grid.color'] = 'whitesmoke'
matplotlib.rcParams['grid.linestyle'] = ':'
matplotlib.rcParams['grid.linewidth'] = 0.3

# for hillshade
azimuthDegree = cfg.getfloat('azimuthDegree')
elevationDegree = cfg.getfloat('elevationDegree')
vertExag = cfg.getfloat('vertExag')
hillshadeContLevs = cfg.getint('hillshadeContLevs')

# define settings for colormaps creation
discreteLevels = cfg.getint('discreteLevels')
############################
# Color maps
############################
# hell white/green to dark blue
cmapGB = copy.copy(sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True))
cmapGB.set_bad(color='k')

cmapReds = copy.copy(matplotlib.cm.Reds)
cmapReds.set_bad(color='k')

cmapBlues = copy.copy(matplotlib.cm.Blues)
cmapBlues.set_bad(color='k')

cmapGreys = copy.copy(matplotlib.cm.get_cmap("Greys"))

cmapPlasma = copy.copy(matplotlib.cm.plasma)
cmapPlasma.set_bad(color='k')

cmapViridis = copy.copy(matplotlib.cm.viridis)
cmapViridis.set_bad(color='k')

# divergent color map
cmapdiv = copy.copy(cmapCrameri.broc)

# custom colomaps
# cmap based on avaframe logo colors
colorAvaframe = ['#0EF8EA', '#12E4E6', '#28D0DF', '#3CBCD5', '#4AA8C9',
                 '#5595B9', '#5C82A8', '#5F6F95', '#5E5E81', '#5A4D6C', '#523E58', '#483045']
cmapAvaframe = mplCol.ListedColormap(colorAvaframe)
cmapAvaframe.set_bad(color='k')
# add a continuous version
cmapAvaframeCont = mplCol.LinearSegmentedColormap.from_list('cmapAvaframeCont', colorAvaframe, N=256)


# for the choice of the colormaps, check https://www.fabiocrameri.ch/colourmaps/
# and http://hclwizard.org:3000/hclwizard/
# multi sequential colormap for pressure

levP = list(fU.splitIniValueToArraySteps(cfgPlotUtils['pressureColorLevels']))
# Hawaii color map
colorsP = ["#B0F4FA", "#75C165", "#A96C00", "#8B0069"]
cmapP = copy.copy(cmapCrameri.hawaii.reversed())

# multi sequential colormap for flow thickness
levT = list(fU.splitIniValueToArraySteps(cfgPlotUtils['thicknessColorLevels']))
# Lajolla color map
colorsT = ["#FCFFC9", "#EBCE7B", "#DE9529", "#BE5A32", "#7F2B3F", "#1D0B14"]
cmapT = copy.copy(cmapCrameri.lajolla)

# multi sequential colormap for speed
levS = list(fU.splitIniValueToArraySteps(cfgPlotUtils['speedColorLevels']))
# Batflow color map
colorsS = ['#FFCEF4', '#FFA7A8', '#C19A1B', '#578B21', '#007054', '#004960',
           '#201158']
cmapS = copy.copy(cmapCrameri.batlow.reversed())

# multi sequential colormap for Travel Angle
levTA = list(fU.splitIniValueToArraySteps(cfgPlotUtils['travelAngleColorLevels']))
# Batflow color map
colorsTA = ['#FFCEF4', '#FFA7A8', '#C19A1B', '#578B21', '#007054', '#004960',
            '#201158']
cmapTA = copy.copy(cmapCrameri.lapaz)

# colormap used if no resType provided
cmapNN = copy.copy(cmapCrameri.imola.reversed())

# colormap used for range time diagram - continuous
cmapRangeTime = copy.copy(cmapCrameri.batlowW_r)

# colormap used for radar field of view plot
cmapRadarFOV = copy.copy(cmapCrameri.davos)

# multi sequential colormap for kinetic energy
levE = list(fU.splitIniValueToArraySteps(cfgPlotUtils['energyColorLevels']))
# Tokyo color map
colorsE = ['#aed6a2', '#949e8c', '#8b6f7f', '#673762', '#1a0e34']
# colormap used for peak kinetic energy
cmapE = copy.copy(cmapCrameri.tokyo.reversed())

# colormap for probabilities
levProb = list(fU.splitIniValueToArraySteps(cfgPlotUtils['probaColorLevels']))
# lapaz color map
colorsProb = ['#FEF1F1', '#B2AB96', '#5B8BA3', '#2D5393', '#1A0C64']
cmapProbmap = copy.copy(cmapCrameri.lapaz.reversed())

###############################################
# Set colormaps to use
###############################################
# ploting with a descrete (contCmap = continuousCmap = False) or continuous colormap (contCmap = True)?
# if continuous, only the cmap argument in the cmapTictionnary maters
# replace it with the wanted colormap
contCmap = cfg.getboolean('contCmap')
# for pressure
cmapPres = {}
cmapPres['cmap'] = cmapP
cmapPres['colors'] = colorsP
cmapPres['levels'] = levP


cmapThickness = {}
cmapThickness['cmap'] = cmapT
cmapThickness['colors'] = colorsT
cmapThickness['levels'] = levT

cmapSpeed = {}
cmapSpeed['cmap'] = cmapS
cmapSpeed['colors'] = colorsS
cmapSpeed['levels'] = levS

cmapTravelAngle = {}
cmapTravelAngle['cmap'] = cmapTA
cmapTravelAngle['colors'] = colorsTA
cmapTravelAngle['levels'] = levTA

cmapProb = {}
cmapProb['cmap'] = cmapProbmap
cmapProb['colors'] = colorsProb
cmapProb['levels'] = levProb

cmapEnergy = {}
cmapEnergy['cmap'] = cmapE
cmapEnergy['colors'] = colorsE
cmapEnergy['levels'] = levE

colorMaps = {'ppr': cmapPres, 'pfv': cmapSpeed, 'pft': cmapThickness, 'P': cmapPres, 'FV': cmapSpeed,
             'FM': cmapThickness, 'Vx': cmapSpeed, 'Vy': cmapSpeed, 'Vz': cmapSpeed, 'FTV': cmapSpeed,
             'FT': cmapThickness, 'prob': cmapProb, 'pta': cmapTravelAngle, 'TA': cmapTravelAngle,
             'pke': cmapEnergy}

cmapDEM = cmapGreys

cmapDEM2 = {}
cmapDEM2['cmap'] = cmapGreys
cmapDEM2['colors'] = colorsS
cmapDEM2['levels'] = None

cmapAimec = cmapAvaframe

cmapVar = {}
# cmapVar['cmap'] = cmapAvaframeCont
cmapVar['colors'] = colorAvaframe
cmapVar['levels'] = None

cmapGreys1 = {}
cmapGreys1['cmap'] = cmapGreys


def makeColorMap(colormapDict, levMin, levMax, continuous=False):
    """Get colormap for plot

    get the colormap, norm, levels... for ploting results


    1) colormapDict = dict (details see Parameters) - depending on the continuous flag and what is
    given in the colormapDict dictionary, the cmap is created and the levelsNew and
    colorsNew are created.
    2) colormapDict = matplotlib cmap - the continuous flag is ignored and a continuous cmap is
    returned with levelsNew set to None.


    Parameters
    ----------
        colormapDict: dict or cmap
            if colormap= dict :
                cmap: matplotlib colormap
                    continuous colormap
                colors: list
                    list of hex or rgb or dec colors
                levels: list
                    levels corresponding to the colors (same number of levels as colors)
        levMin : float
            minimum value of the colormap
        levMax: float
            maximum value of the colormap
        continuous: boolean
            True for continuous cmaps

    Returns
    -------
        cmap: matplotlib colormap
            new color map
        colorsNew: list
            new list of hex or rgb or dec colors
        levelsNew: list
            new levels corresponding to the colors (one more level than colors)
        norm: matplotlib norm
            norm associated to the levels and the colors (includes the vmin and vmax for
            continuous colormaps)
    """
    if type(colormapDict) in [matplotlib.colors.LinearSegmentedColormap,
                              matplotlib.colors.ListedColormap]:
        cmap = colormapDict
        colorsNew = None
        norm = mplCol.Normalize(vmin=levMin, vmax=levMax, clip=False)
        levelsNew = None
    elif continuous:
        # make a continuous color map
        # check if cmap is provided
        if 'cmap' in colormapDict.keys():
            cmap = colormapDict['cmap']
            colorsNew = None
        # check if list of colors is provided
        elif 'colors' in colormapDict.keys():
            colorsNew = colormapDict['colors']
            cmap = mplCol.LinearSegmentedColormap.from_list('myCmap', colorsNew, N=256)
        # Houston ve have a problem
        else:
            message = 'You need a `colors` list or a `cmap` to be able to create the colormap'
            log.error(message)
            raise FileNotFoundError(message)

        norm = mplCol.Normalize(vmin=levMin, vmax=levMax, clip=False)
        levelsNew = None
    else:
        # make a discrete color map
        if 'levels' in colormapDict:
            levels = colormapDict['levels']
        else:
            if 'colors' in colormapDict:
                levels = list(np.linspace(levMin, levMax,
                                          len(colormapDict['colors'])+1))
                log.warning('No `levels` list is provided to generate a discrete colormap, \
                            creating %d levels ranging from %.2f to %.2f' %
                            (len(colormapDict['colors']), levMin, levMax))
            else:
                levels = list(np.linspace(levMin, levMax, discreteLevels))
                log.warning('No `levels` list is provided to generate a discrete colormap, \
                            creating %d levels ranging from %.2f to %.2f' %
                            (discreteLevels, levMin, levMax))
        try:
            indEnd = np.where(np.asarray(levels) >= levMax)[0][0]
        except IndexError:
            indEnd = len(levels)
        levelsNew = levels[:indEnd]
        levelsNew.append(levMax)

        # check if list of colors is provided
        if 'colors' in colormapDict.keys():
            if indEnd > len(colormapDict['colors']):
                message = 'Number of levels is not allowed to exceed number of colors'
                log.error(message)
                raise AssertionError(message)
            colors = colormapDict['colors']
            colorsNew = colors[:indEnd]
            colorsNew.append(colors[indEnd-1])
            cmap = mplCol.ListedColormap(colorsNew)
            cmap.set_bad(color='w')
            cmap.set_under(color='w')
        # check if a cmap is provided
        elif 'cmap' in colormapDict.keys():
            cmap = colormapDict['cmap']
            colorsNew = cmap(np.asarray(levelsNew[:-1])/levelsNew[-1])
        # Houston ve have a problem
        else:
            message = 'A `colors` list or a `cmap` is required to create the colormap'
            log.error(message)
            raise FileNotFoundError(message)

        if len(levelsNew) == 1:
            levelsNew = [0] + levelsNew
        norm = mplCol.BoundaryNorm(levelsNew, cmap.N)
    return cmap, colorsNew, levelsNew, norm


###################################
# shortcut plot functions
###################################
def NonUnifIm(ax, x, y, z, xlab, ylab, **kwargs):
    im = NonUniformImage(ax, **kwargs)
    im.set_data(x, y, z)
    # im.set_clim(vmin=vmin, vmax=vmax)
    ref = ax.images.append(im)
    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min(), y.max()])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ref, im


def saveAndOrPlot(pathDict, outFileName, fig):
    """
    Receive a plot handle and config and check whether to save and or plot
    closes it afterwards.
    If saved, the plot will be saved in pathDict['pathResult']/outFileName.extension and the path will be returned

    Parameters
    ----------
        pathDict : dict
            with field "pathResult" (path to output folder)
        outFileName: str
            output file name
        fig: matplotlib figure

    """

    if cfgFlags.getboolean('showPlot'):
        plt.show()
    else:
        plt.ioff()

    outPath = None
    if cfgFlags.getboolean('savePlot'):
        outFileName = outFileName.replace(".", "_")
        outDir = pathlib.Path(pathDict['pathResult'])
        outPath = outDir / (outFileName + '.' + outputFormat)
        outName = outDir / outFileName
        if not outDir.is_dir():
            fU.makeADir(outDir)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.savefig(outName)
            log.info('saved to : %s ' % outName)

    plt.close(fig)

    return outPath


def constrainPlotsToData(inputData, cellSize, extentOption=False, constrainedData=False, buffer=''):
    """ constrain inut raster dataset to where there is data plus buffer zone

        Parameters
        -----------
        inputData: numpy array
            raster data to be plotted
        cellSize: float
            cellsize of raster data
        extentOption: bool
            if True rows and columns limits converted to acutal extent in meters
        buffer: float
            buffer for constraining data in meters - optional if not provided read from ini file

        Returns
        --------
        rowsMin, rowsMax: int or float
            indices where there is data in x direction (or min/max extent in meters)
        colsMin, colsMax: int or float
            indices where there is data in y direction (or min/max extent in meters)
        dataConstrained: numpy array
            constrained array where there is data
        """

    # check if buffer is given as input or needs to be read from ini file
    if buffer != '':
        plotBuffer = int(buffer / cellSize)
    else:
        plotBuffer = int(cfg.getfloat('plotBuffer') / cellSize)

    ind = np.where(inputData > 0)
    if len(ind[0]) > 0:
        rowsMin = max(np.amin(ind[0])-plotBuffer, 0)
        rowsMax = min(np.amax(ind[0])+plotBuffer, inputData.shape[0]-1)
        colsMin = max(np.amin(ind[1])-plotBuffer, 0)
        colsMax = min(np.amax(ind[1])+plotBuffer, inputData.shape[1]-1)
    else:
        rowsMin = 0
        rowsMax = inputData.shape[0]
        colsMin = 0
        colsMax = inputData.shape[1]

    if extentOption:
        rowsMinPlot = rowsMin*cellSize
        rowsMaxPlot = (rowsMax+1)*cellSize
        colsMinPlot = colsMin*cellSize
        colsMaxPlot = (colsMax+1)*cellSize
        if constrainedData:
            dataConstrained = inputData[rowsMin:rowsMax+1, colsMin:colsMax+1]
            return rowsMinPlot, rowsMaxPlot, colsMinPlot, colsMaxPlot, dataConstrained
        else:
            return rowsMinPlot, rowsMaxPlot, colsMinPlot, colsMaxPlot
    else:
        if constrainedData:
            dataConstrained = inputData[rowsMin:rowsMax+1, colsMin:colsMax+1]
            return rowsMin, rowsMax, colsMin, colsMax, dataConstrained
        else:
            return rowsMin, rowsMax, colsMin, colsMax


def addColorBar(im, ax2, ticks, myUnit, title='', extend='neither', pad=0.05, tickLabelsList='', location='right'):
    '''
    Adds, styles and labels a colorbar to the given image and axes
    '''


    cbar = ax2.figure.colorbar(im, ax=ax2, ticks=ticks, extend=extend, pad=pad, shrink=0.9, location=location)
    cbar.outline.set_visible(False)
    # make sure the cbar title does not overlap with the cbar itself
    if extend in ['both', 'max']:
        pad = 25
    else:
        pad = 10
    if myUnit != None:
        cbar.ax.set_title('[' + myUnit + ']', pad=pad)
    if title != '':
        cbar.set_label(title)
    if len(tickLabelsList) > 0:
        cbar.ax.set_yticklabels(tickLabelsList)
    return cbar


def putAvaNameOnPlot(ax, avaDir, date=True):
    '''
    Puts the date and avalanche name (or a list of ava names) in the lower left corner of the given
    matplotlib axes, if date=False only avalanche name is put
    '''

    # if avaDir is just a single avaDir or a list of avaDirs
    if isinstance(avaDir, str) or isinstance(avaDir, pathlib.Path):
        avaName = pathlib.PurePath(avaDir).name
        if date:
            infoText = datetime.datetime.now().strftime("%d.%m.%y") + '; ' + str(avaName)
        else:
            infoText = str(avaName)
    else:
        if date:
            infoText = datetime.datetime.now().strftime("%d.%m.%y") + ';'
        else:
            infoText = ''
        for ava in avaDir[:-1]:
            avaName = pathlib.PurePath(ava).name
            infoText = infoText + str(avaName) + ';'
        infoText = infoText + str(pathlib.PurePath(avaDir[-1]).name)

    ax.annotate(infoText, xy=(0.01, 0.01), xycoords='axes fraction',
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.5))

    return infoText


def putInfoBox(ax, infoText, location='lowerRight', color='black', hAlignment='right', alphaF=0.5):
    '''
    Puts the infoBox in the lower right or upper left or upper right or lowerLeft corner of the given
    matplotlib axes

    Parameters
    ----------
    ax: matplotlib axes
    infoText: str
        text of infobox
    location: str
        default: lowerRight, options: upperRight, upperLeft
    color: str
        color of font
    hAlignment: str
        horizontal alignment
    alphaF: float
        alpha value of text

    '''

    # if avaDir is just a single avaDir or a list of avaDirs
    if isinstance(infoText, str) is False:
        log.warning('Info text to be added to plot is not a string - ignored')

    if location == 'upperRight':
        xy = (0.97, 0.99)
    elif location == 'upperLeft':
        xy = (0.01, 0.97)
    elif location == 'lowerLeft':
        xy = (0.01, 0.01)
    else:
        xy = (0.99, 0.01)

    ax.annotate(infoText, xy=xy, xycoords='axes fraction', fontsize=8, horizontalalignment=hAlignment,
        verticalalignment='top', color=color, alpha=alphaF,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.5))


def constrainToMinElevation(avaDir, data, cfg, cellSize, extentOption=False):
    """ constrain data array to bufferzone around min elevation of dem where there is data in data array

        Parameters
        -----------
        avaDir: pathlib path of str
            path to avalanche directory
        data: numpy array
            data array of equal shape as dem data
        cfg: configparser object
            configuration settings for buffer zone
        extentOption: bool
            if True in meters if False in rows and cols

        Returns
        --------
        dataCut : numpy array
            data constrained to a bufferzone
        xOrigin: float
            origin of x axis
        yOrigin: float
            origin of y axis
    """

    # load dem to identify runout area according to min elevation where peak result != 0
    dem = gI.readDEM(avaDir)

    # mask dem to where there is data in result file
    demCut = np.where(data > 0, dem['rasterData'], np.nan)

    # identify min elevation and cut data to buffer zone around min elevation
    indMin = np.where(demCut == np.nanmin(demCut))
    nrowsMin = indMin[0][0]
    ncolsMin = indMin[1][0]
    rangePlot = int(cfg.getfloat('zoomBuffer') / cellSize)
    dataCut = data[nrowsMin-rangePlot:nrowsMin+rangePlot, ncolsMin-rangePlot:ncolsMin+rangePlot]

    # to determine the extent for plotting
    if extentOption:
        yOrigin = (nrowsMin-rangePlot) * cellSize
        xOrigin = (ncolsMin-rangePlot) * cellSize
    else:
        yOrigin = nrowsMin-rangePlot
        xOrigin = ncolsMin-rangePlot

    return dataCut, xOrigin, yOrigin


def getColorbarTicksForStrings(varVal):
    """ if values for colorbar are strings convert to numbers and provide labels

        Parameters
        -----------
        varVal: list
            list of strings

        Returns
        --------
        itemsList: list
            list of unique parameter values (strings)
        ticksList: numpy array
            array of unique parameter values (floats)
        varValV: numpy array
            array with parameter values (floats)
    """

    countItems = 0
    itemsList = []
    itemDict = {}
    for index, item in enumerate(varVal):
        if item not in itemsList:
            countItems = countItems + 1
            itemDict[item] = countItems
            itemsList.append(item)

    varValV = np.array([])
    for item in varVal:
        varValV = np.append(varValV, itemDict[item])
    ticksList = np.linspace(1, countItems, countItems)

    return itemsList, ticksList, varValV


def getColors4Scatter(values, nSamples, unitSC):
    """ provide cMap, colors, ticks... for a scatter plot



        Parameters
        -----------
        values: list
            list of strings
        nSamples: int
            number of samples in the scatter
        unitSC: str
            unit for the colorbar

        Returns
        --------
        cmapSC: matplotlib colormap
            color map
        colorSC: numpy array
            values for coding the color of each sample
        ticksSC: list
            ticks levels for the colorbar
        normSC: matplotlib norm
            norm associated to the levels and the colors
        unitSC: str
            unit for the values displayed in the colorbar
        itemsList: list
            list of unique parameter values (strings). In case the color variation is done for strig values
        displayColorBar: boolean
            True if the colorBar should be displayed (if a color variation is applied)
    """
    itemsList = ''
    if values is None:
        displayColorBar = False
        colorSC = 0.5 * np.ones(nSamples)
        cmapSC, _, ticksSC, normSC = makeColorMap(cmapVar, None, None, continuous=True)
    else:
        typeCP = type(values[0])
        if typeCP == str:
            itemsList, ticksSC, colorSC = getColorbarTicksForStrings(values)
            cmapSC, _, _, normSC = makeColorMap(cmapVar, np.amin(colorSC), np.amax(colorSC), continuous=True)
            displayColorBar = True
            unitSC = '-'
        else:
            colorSC = values
            cmapSC, _, ticksSC, normSC = makeColorMap(cmapVar, np.nanmin(colorSC), np.nanmax(colorSC), continuous=True)
            displayColorBar = True
    return cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar


def addHillShadeContours(ax, data, cellSize, extent, colors=['gray'], onlyContours=False):
    """ add hillshade and contours for given DEM data

        Parameters
        -----------
        ax: matplotlib axes
            axes of plot
        data: numpy array
            dem data
        cellSize: float
            cell size of data
        extent: list
            extent [x0, x1, y0, y1] x0, y0 lower left corner and extent for imshow plot
        colors: list
            optional, colors for elevation contour lines
        onlyContours: bool
            if True add only contour lines but no hillshade
    """

    if onlyContours:
        ls = None
    else:
        # create lightSource
        ls = LightSource(azdeg=azimuthDegree, altdeg=elevationDegree)

        # add hillshade to axes
        im0 = ax.imshow(ls.hillshade(data, vert_exag=vertExag, dx=data.shape[1], dy=data.shape[0]), cmap='gray',
                        extent=extent, origin='lower', aspect='equal', zorder=1)

    # create x,y coors for data array
    nrows, ncols = data.shape
    X, Y = gT.makeCoordinateGrid(extent[0], extent[2], cellSize, ncols, nrows)

    # add contour lines
    CS = ax.contour(X, Y, data, colors=colors, levels=hillshadeContLevs, alpha=1.,
                    linewidths=0.5, zorder=2)

    # add labels
    ax.clabel(CS, CS.levels, inline=True, fontsize=8, zorder=3)

    # add info box with indication of label meaning
    putInfoBox(ax, '- elevation [m]', location='upperLeft', color='gray',
               hAlignment='left', alphaF=1.0)

    return ls, CS
