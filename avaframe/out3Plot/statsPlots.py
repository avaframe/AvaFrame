"""
    plot statistics of simulations

"""

# load python modules
import os
import numpy as np
import logging
import matplotlib.pyplot as plt
import numpy.ma as ma
import seaborn as sns
import pandas as pd
import pathlib

# local imports
import avaframe.out3Plot.plotUtils as pU
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in1Data.getInput as gI
from avaframe.in3Utils import cfgHandling


# create local logger
log = logging.getLogger(__name__)


def plotValuesScatter(peakValues, resType1, resType2, cfg, avalancheDir, statsMeasure='max', flagShow=False):
    """ Produce scatter plot of statistical measures (eg. max values) (resType1 und resType2),
        for one set of simulations or multiple

        Parameters
        -----------
        peakDictList: dict
            peakValues dictionary that contain max values of peak parameters and parameter variation info
        resType1: str
            result parameter 1, 'ppr', 'pft', 'pfv'
        resType2: str
            result parameter 1, 'ppr', 'pft', 'pfv'
        cfg: dict
            configuration, for now contains output location and varPar: parameter that is varied
            to perfom a set of simulations
        statsMeasure: str
            statistical measure for plotting, options: max, mean, min, std
        flagShow: bool
            if True show plot
        """

    varPar = cfg['varPar']
    # extract values from dictionaries
    varVal = []
    values1 = []
    values2 = []
    for key in peakValues:
        values1.append(peakValues[key][resType1][statsMeasure])
        values2.append(peakValues[key][resType2][statsMeasure])
        varVal.append(peakValues[key]['varPar'])

    log.info('Number of simulations is: %d' % (len(varVal)))

    # Get name and units for resTypes and parameter variation to annotate plots
    name1 = pU.cfgPlotUtils['name%s' % resType1]
    name2 = pU.cfgPlotUtils['name%s' % resType2]
    unit1 = pU.cfgPlotUtils['unit%s' % resType1]
    unit2 = pU.cfgPlotUtils['unit%s' % resType2]
    nameVar = cfg['varParName']
    unitVar = cfg['varParUnit']

    # if the varied parameter used for colorcoding is a string
    if isinstance(varVal[0], str):
        itemsList, ticksList, varValV = pU.getColorbarTicksForStrings(varVal)
    else:
        varValV = np.array(varVal)
        itemsList = ''

    # load variation colormap
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapVar, np.amin(varValV), np.amax(varValV), continuous=True)
    if isinstance(varVal[0], str):
        ticks = ticksList

    fig, ax = plt.subplots()
    plt.title('%s vs. %s' % (name1, name2))
    # set range and steps of colormap
    cc = varValV
    sc = ax.scatter(values1, values2, c=cc, cmap=cmap)
    ax.set_xlabel('%s %s [%s]' % (statsMeasure, name1, unit1))
    ax.set_ylabel('%s %s [%s]' % (statsMeasure, name2, unit2))
    pU.addColorBar(sc, ax, ticks, unitVar, nameVar, tickLabelsList=itemsList)
    pU.putAvaNameOnPlot(ax, avalancheDir)

    # shpw plot
    if flagShow:
        plt.show()

    # save fig
    outDir = cfg['outDir']
    fig.savefig(os.path.join(outDir, 'Scatter_%s_vs_%s_dist%s_%s.png' %
                (resType1, resType2, cfg['distType'], cfg['scenario'])))
    plt.close(fig)


def plotValuesScatterHist(peakValues, resType1, resType2, cfg, avalancheDir,
                          statsMeasure='max', flagShow=False, flagHue=False):
    """ Produce scatter and marginal kde plot of max values, for one set of simulations or multiple

        Parameters
        -----------
        peakValues: dict
            peakValues dictionary that contain max values of peak parameters and parameter variation info
        resType1: str
            result parameter 1, 'ppr', 'pft', 'pfv'
        resType2: str
            result parameter 1, 'ppr', 'pft', 'pfv'
        cfg: dict
            configuration, for now contains output location and varPar: parameter that is varied
            to perfom a set of simulations
        statsMeasure: str
            statistical measure for plotting, options: max, mean, min, std
        flagShow: bool
            if True show plot

    """

    varPar = cfg['varPar']
    # extract values from dictionaries
    varVal = []
    values1 = []
    values2 = []
    scenario = []

    for key in peakValues:
        values1.append(peakValues[key][resType1][statsMeasure])
        values2.append(peakValues[key][resType2][statsMeasure])
        varVal.append(peakValues[key]['varPar'])
        if 'scenario' in peakValues[key] and flagHue:
            scenario.append(peakValues[key]['scenario'])

    log.info('Number of simulations is: %d' % (len(varVal)))

    # Get name and units for resTypes and parameter variation to annotate plots
    name1 = pU.cfgPlotUtils['name%s' % resType1]
    name2 = pU.cfgPlotUtils['name%s' % resType2]
    unit1 = pU.cfgPlotUtils['unit%s' % resType1]
    unit2 = pU.cfgPlotUtils['unit%s' % resType2]
    nameVar = cfg['varParName']
    unitVar = cfg['varParUnit']
    varValV = np.array(varVal)
    title1 = statsMeasure + ' ' + name1+' [' + unit1 + ']'
    title2 = statsMeasure + ' ' + name2+' [' + unit2 + ']'

    if flagHue:
        # create pandas data frame reqiured for seaborn jointplot
        maxVals = {name1: values1, name2: values2, 'scenario': scenario}
        maxData = pd.DataFrame(maxVals)
        fig1 = sns.jointplot(data=maxData, x=name1, y=name2, hue="scenario")
    else:
        fig1 = sns.jointplot(x=values1, y=values2)

    # add title and text box
    fig1.ax_joint.set_xlabel(title1)
    fig1.ax_joint.set_ylabel(title2)
    pU.putAvaNameOnPlot(fig1.ax_joint, avalancheDir)

    # save fig
    outDir = cfg['outDir']
    plt.savefig(os.path.join(outDir, 'Scatterkde_%s_vs_%s_dist%s_%s.png' %
                (resType1, resType2, cfg['distType'], cfg['scenario'])))

    # shpw plot
    if flagShow:
        plt.show()

    plt.close()


def plotHistCDFDiff(dataDiffPlot, ax1, ax2, insert='True', title=['', '']):
    """ Produce histogram and CDF plot of the raster difference of two simulations

    Parameters
    -----------
    dataDiffPlot: 2D numpy array
        raster of the difference of the two simulations
    ax1: axes
        axes for the histogram plot
    ax2: axes
        axes for the CDF plot
    insert: boolean
        true if the plots are in inserted axes (size of the lables is then smaller)
    title: list
        if not inserts, title for the plots

    """
    # Difference between datasets
    diffMax = np.nanmax(dataDiffPlot)
    diffMin = np.nanmin(dataDiffPlot)

    sortedDiffPlot = np.sort(np.abs(dataDiffPlot))
    nSample = np.size(sortedDiffPlot)
    hist = np.array(range(nSample))/float(nSample)
    ticks = []
    centile1 = float(pU.cfg['centile1'])
    centile2 = float(pU.cfg['centile2'])
    for val in [centile1, centile2]:
        ind = np.searchsorted(hist, val)
        ind = min(ind, np.size(hist)-1)
        ax1.plot(sortedDiffPlot, hist)
        ax1.hlines(hist[ind], 0, sortedDiffPlot[ind], linestyles='--', linewidths=0.5)
        ax1.vlines(sortedDiffPlot[ind], 0, hist[ind], linestyles='--', linewidths=0.5)
        ticks.append(sortedDiffPlot[ind])

    ax2.set_xlim([-sortedDiffPlot[ind], sortedDiffPlot[ind]])
    width = diffMax - diffMin
    stepsInterval = int(pU.cfg['steps2Centile2'])
    stepWidth = 2*sortedDiffPlot[ind]/stepsInterval    # stepsInterval bins in the [-2sigma,+2sigma] interval
    bins = int(width/stepWidth)

    # reduce bins to a sensible size
    if bins > 1000:
        bins = 1000

    ax2.hist(dataDiffPlot, bins=bins, density=True, histtype="stepfilled")
    ax2.get_yaxis().set_ticks([])
    if insert:
        ax2.tick_params(axis='x', which='major', labelsize=8, rotation=45)
        ax2.set_title(title[0])

    ticks.append(np.floor(np.nanmax(np.abs(dataDiffPlot))))
    ax1.set_xticks(ticks)
    if insert:
        ax1.tick_params(axis='both', which='major', labelsize=8, rotation=45)
        ax1.set_title(title[1])

    return ticks


def plotProbMap(avaDir, inDir, cfgFull, demPlot=False):
    """ plot probability maps including contour lines

        Parameters
        ----------
        avaDir: str
            path to avalanche directory
        inDir: str
            path to datasets that shall be plotted
        cfgFull: configParser object
            configuration settings for probAna
            keys used: name, cmapType, levels, unit
        demPlot: bool
            if True plot dem in background with contourlines for elevation that is found in avaDir/Inputs

    """

    # configuration settings
    cfg = cfgFull['PLOT']

    if demPlot:
        demFile = gI.getDEMPath(avaDir)
        demData = IOf.readRaster(demFile, noDataToNan=True)
        demField = demData['rasterData']

    # fetch probabiltiy map datasets in inDir
    dataFiles = list(inDir.glob('*.asc'))
    if dataFiles == []:
        message = 'No probability map dataset found in: %s' % (inDir)
        log.error(message)
        raise FileNotFoundError(message)
        log.info('Probability maps found for generating plots: %d' % len(dataFiles))

    # load levels and define colors
    levels = fU.splitIniValueToArraySteps(cfg['levels'])
    multLabel = False
    if len(levels) > 2:
        multLabel = True
        cmapType = pU.cmapGreys1
        colorBackGround = 'seashell'
    else:
        cmapType = pU.colorMaps[cfg['cmapType']]
        colorsP = ['black', 'white']
        colorBackGround = 'gainsboro'
        if len(levels) == 1 and levels[0] > 0.5:
            colorsP = ['white']
    unit = cfg['unit']

    # set colors for contourlines according to levels
    labels = []
    for lev in levels:
        labels.append('%.0f %%' % (lev*100))

    # loop over all datasets and create plots
    for data in dataFiles:

        raster = IOf.readRaster(data, noDataToNan=True)
        dataPlot = raster['rasterData']
        header = IOf.readASCheader(data)
        cellSize = header['cellsize']

        # Set dimensions of plots
        ny = dataPlot.shape[0]
        nx = dataPlot.shape[1]
        Ly = ny*cellSize
        Lx = nx*cellSize

        # constrain plot to where there is data
        rowsMin, rowsMax, colsMin, colsMax, dataConstrained = pU.constrainPlotsToData(dataPlot, cellSize,
            extentOption=False, constrainedData=True)
        dataPlot = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)

        # create figure and add title
        fig = plt.figure(figsize=(pU.figW*2, pU.figH))
        fullTitle = cfg['name'] + ' based on %s $>$ %s %s' % (cfgFull['GENERAL']['peakVar'], cfgFull['GENERAL']['peakLim'], cfgFull['GENERAL']['unit'])
        suptitle = fig.suptitle(fullTitle, fontsize=14, color='0.5')
        ax1 = fig.add_subplot(121)

        # set extent in meters using cellSize
        rowsMinPlot = rowsMin*cellSize
        rowsMaxPlot = (rowsMax+1)*cellSize
        colsMinPlot = colsMin*cellSize
        colsMaxPlot = (colsMax+1)*cellSize
        extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]

        # for now continuous color map is desired
        cmap, _, ticks, norm = pU.makeColorMap(cmapType, np.nanmin(dataPlot), np.nanmax(dataPlot),
            continuous=True)

        if demPlot and demData['header']['cellsize'] == cellSize:
            # also constrain DEM to data constrained
            demConstrained = demField[rowsMin:rowsMax+1, colsMin:colsMax+1]
            # add DEM hillshade with contour lines
            ls, CS = pU.addHillShadeContours(ax1, demConstrained, cellSize, extent)
            dataPlot = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)
            cmap.set_bad(alpha=0)
        else:
            cmap.set_bad(colorBackGround)

        # add data plot
        im1 = ax1.imshow(dataPlot, cmap=cmap, extent=extent, origin='lower', aspect=nx/ny, norm=norm,
            zorder=3)

        # create meshgrid for contour plot also constrained to where there is data
        xx = np.arange(colsMinPlot, colsMaxPlot, cellSize)
        yy = np.arange(rowsMinPlot, rowsMaxPlot, cellSize)
        X, Y = np.meshgrid(xx, yy)

        # add contourlines for levels
        if multLabel:
            CS = ax1.contour(X, Y, dataPlot, levels=levels, cmap=pU.cmapT.reversed(), linewidths=1, zorder=4)
        else:
            CS = ax1.contour(X, Y, dataPlot, levels=levels, colors=colorsP, linewidths=1, zorder=4)
        for i in range(len(labels)):
            CS.collections[i].set_label(labels[i])

        pU.addColorBar(im1, ax1, ticks, unit)
        title = str('%s' % cfg['name'])
        ax1.set_xlabel('x [m]')
        ax1.set_ylabel('y [m]')

        # add zoom plot of runout area
        ax2 = fig.add_subplot(122)

        # determine zoom in runout area
        dataCut, xOrigin, yOrigin = pU.constrainToMinElevation(avaDir, raster['rasterData'], cfg, cellSize,
            extentOption=True)

        # constrain to where there is data
        rowsMinPlot, rowsMaxPlot, colsMinPlot, colsMaxPlot, dataCutConstrained = pU.constrainPlotsToData(dataCut,
            cellSize, extentOption=True, constrainedData=True, buffer=cfg.getfloat('constrainBuffer'))
        dataCutConstrained = np.ma.masked_where(dataCutConstrained == 0.0, dataCutConstrained)

        # set extent of zoom Plot
        x0 = xOrigin + colsMinPlot
        x1 = xOrigin + colsMaxPlot
        y0 = yOrigin + rowsMinPlot
        y1 = yOrigin + rowsMaxPlot

        # add plot
        im2 = ax2.imshow(dataCutConstrained, cmap=cmap, extent=[x0, x1, y0, y1],
                         origin='lower', aspect=nx/ny, norm=norm)

        # create meshgrid for contour plot also constrained to where there is data
        xx2 = np.arange(x0, x1, cellSize)
        yy2 = np.arange(y0, y1, cellSize)
        X2, Y2 = np.meshgrid(xx2, yy2)

        # add contourlines for levels
        if multLabel:
            CS2 = ax2.contour(X2, Y2, dataCutConstrained, levels=levels, cmap=pU.cmapT.reversed(), linewidths=1)
        else:
            CS2 = ax2.contour(X2, Y2, dataCutConstrained, levels=levels, colors=colorsP, linewidths=1)

        # Get the handles for the legend elements
        handles, _ = CS2.legend_elements()

        ax2.set_xlabel('x [m]')
        ax2.set_ylabel('y [m]')

        plt.legend(handles, labels, facecolor = 'black', framealpha = 0.04)
        pU.addColorBar(im2, ax2, ticks, unit)

        outDir = inDir / 'plots'
        fU.makeADir(outDir)
        avaName = pathlib.PurePath(avaDir).name
        outFile = outDir / ('%s_probMap_lim%s.%s' % (avaName, cfgFull['GENERAL']['peakLim'], pU.outputFormat))
        fig.savefig(outFile)
        plt.close(fig)


def resultHistPlot(cfg, dataDF, xName='', scenario='', stat='count', parametersDict=''):
    """ create a histogram of values and optional colorcode using scenario name
        and option to filter simulations using parametersDict

        Parameters
        -----------
        cfg: configparser object
            configuration info here used outDir
        dataDF: dataFrame
            dataFrame with info on simulation results and configuration one line per simulation
            (e.g. aimec resAnalysisDF)
        xName: str
            column name for x axis
        scenario: str
            column name used to colorcode values
        stat: str
            statistical measure to show (percent, probability, density, count, frequency), default count
        parametersDict: dict
            optional - dictionary filter criteria, parameter name and list of values


        Returns
        --------
        plotPath: pathlib path
            path to figure

    """

    # filter DF with parametersDict
    if parametersDict != '':
        simNameList = cfgHandling.filterSims(cfg['avalancheDir'], parametersDict, specDir='', simDF=dataDF)
        dataDF = dataDF[dataDF['simName'].isin(simNameList)]

    # initialize figure
    fig, ax = plt.subplots()

    # create histogram
    bars = sns.histplot(data=dataDF, x=xName, hue="scenario", stat=stat, ax=ax)

    # create second y axis for ecdf
    ax2 = ax.twinx()
    cdf = sns.ecdfplot(data=dataDF, x=xName, hue='scenario', stat='count', ax=ax2)

    ax.set_ylabel('histogram ' + stat)
    ax2.set_ylabel('ecdf count')

    outFileName = '%s_' % xName + 'histogram'
    plotPath = pU.saveAndOrPlot({'pathResult': cfg['outDir']}, outFileName, fig)

    return plotPath


def plotDistFromDF(cfg, dataDF, name1, name2, scenario='', parametersDict='', type=''):
    """ create a dist plot from dataframe for name1 on x axis and name2 on y axis, optionally
        colorcoded with scenario name and filtered with parametersDict

        Parameters
        -----------
        cfg: configparser object
            configuration settings here outDir
        dataDF: dataframe
            dataframe with one line per simulation and info on model parameters and results
        name1: str
            column name of dataDF to use for plot x axis
        name2: str
            column name of dataDF to use for plot y axis
        scenario: str
            optional name of column used to colorcode points in plots for type=scatter
            or kde for type=dist
        parametersDict: dict
            optional - dictionary filter criteria
        type: str
            optional - type of plot dist or scatter

        Returns
        --------
        plotPath: pathlib path
            path to figure

    """

    # filter DF using parametersDict
    if parametersDict != '':
        simNameList = cfgHandling.filterSims(cfg['avalancheDir'], parametersDict, specDir='', simDF=dataDF)
        dataDF = dataDF[dataDF['simName'].isin(simNameList)]

    # # create figure
    if scenario !='':
        if type == 'scatter':
            fig, ax = plt.subplots()
            ax = sns.scatterplot(data=dataDF[dataDF['simName'].isin(simNameList)], x=name1, y=name2, hue=scenario)
        else:
            dist = sns.jointplot(data=dataDF[dataDF['simName'].isin(simNameList)], x=name1, y=name2, hue=scenario)
            ax = dist.ax_joint
            fig = dist.fig
    else:
        dist = sns.jointplot(data=dataDF[dataDF['simName'].isin(simNameList)], x=name1, y=name2)
        ax = dist.ax_joint
        fig = dist.fig

    # put ava name on plot and save figure
    pU.putAvaNameOnPlot(ax, cfg['avalancheDir'])
    outFileName = '%s_vs_%s_' % (name1, name2) + 'distplot'
    plotPath = pU.saveAndOrPlot({'pathResult': cfg['outDir']}, outFileName, fig)

    return plotPath
