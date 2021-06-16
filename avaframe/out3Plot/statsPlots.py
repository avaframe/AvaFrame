"""
    plot statistics of simulations

"""

# load python modules
import os
import numpy as np
import logging
import matplotlib.pyplot as plt
import avaframe.out3Plot.plotUtils as pU
import avaframe.out3Plot.makePalette as makePalette
import seaborn as sns
import pandas as pd

# create local logger
log = logging.getLogger(__name__)


def plotValuesScatter(peakValues, resType1, resType2, varPar, cfg, avalancheDir, flagShow=False):
    """ Produce scatter plot of max values (resType1 und resType2), for one set of simulations or multiple

        Parameters
        -----------
        peakDictList: dict
            peakValues dictionary that contain max values of peak parameters and parameter variation info
        resType1: str
            result parameter 1, 'ppr', 'pfd', 'pfv'
        resType2: str
            result parameter 1, 'ppr', 'pfd', 'pfv'
        varPar: str
            parameter that is varied to perfom a set of simulations
        cfg: dict
            configuration, for now contains output location
        flagShow: bool
            if True show plot
        """

    # extract values from dictionaries
    varVal = []
    values1 = []
    values2 = []
    for key in peakValues:
        values1.append(peakValues[key][resType1])
        values2.append(peakValues[key][resType2])
        varVal.append(peakValues[key]['varPar'])

    log.info('Number of simulations is: %d' % (len(varVal)))

    # Get name and units for resTypes and parameter variation to annotate plots
    name1 = pU.cfgPlotUtils['name%s' % resType1]
    name2 = pU.cfgPlotUtils['name%s' % resType2]
    unit1 = pU.cfgPlotUtils['unit%s' % resType1]
    unit2 = pU.cfgPlotUtils['unit%s' % resType2]
    nameVar = pU.cfgPlotUtils['name%s' % varPar.lower()]
    unitVar = pU.cfgPlotUtils['unit%s' % varPar.lower()]
    varValV = np.array(varVal)

    # load variation colormap
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        pU.cmapVar, np.amin(varValV), np.amax(varValV), continuous=True)

    fig, ax = plt.subplots()
    plt.title('%s vs. %s' % (name1, name2))
    # set range and steps of colormap
    cc = varValV
    sc = ax.scatter(values1, values2, c=cc, cmap=cmap)
    ax.set_xlabel('%s [%s]' % (name1, unit1))
    ax.set_ylabel('%s [%s]' % (name2, unit2))
    pU.addColorBar(sc, ax, ticks, unitVar, nameVar)
    pU.putAvaNameOnPlot(ax, avalancheDir)

    # shpw plot
    if flagShow:
        plt.show()

    # save fig
    outDir = cfg['outDir']
    fig.savefig(os.path.join(outDir, 'Scatter_%s_vs_%s_dist%s_%s.png' % (resType1, resType2, cfg['distType'], cfg['scenario'])))
    plt.close('all')


def plotValuesScatterHist(peakValues, resType1, resType2, varPar, cfg, avalancheDir, flagShow=False, flagHue=False):
    """ Produce scatter and marginal kde plot of max values, for one set of simulations or multiple

        Parameters
        -----------
        peakDictList: dict
            peakValues dictionary that contain max values of peak parameters and parameter variation info
        resType1: str
            result parameter 1, 'ppr', 'pfd', 'pfv'
        resType2: str
            result parameter 1, 'ppr', 'pfd', 'pfv'
        varPar: str
            parameter that is varied to perfom a set of simulations
        cfg: dict
            configuration, for now contains output location
        flagShow: bool
            if True show plot

    """

    # extract values from dictionaries
    varVal = []
    values1 = []
    values2 = []
    scenario = []

    for key in peakValues:
        values1.append(peakValues[key][resType1])
        values2.append(peakValues[key][resType2])
        varVal.append(peakValues[key]['varPar'])
        if 'scenario' in peakValues[key] and flagHue:
            scenario.append(peakValues[key]['scenario'])

    log.info('Number of simulations is: %d' % (len(varVal)))

    # Get name and units for resTypes and parameter variation to annotate plots
    name1 = pU.cfgPlotUtils['name%s' % resType1]
    name2 = pU.cfgPlotUtils['name%s' % resType2]
    unit1 = pU.cfgPlotUtils['unit%s' % resType1]
    unit2 = pU.cfgPlotUtils['unit%s' % resType2]
    nameVar = pU.cfgPlotUtils['name%s' % varPar.lower()]
    unitVar = pU.cfgPlotUtils['unit%s' % varPar.lower()]
    varValV = np.array(varVal)
    title1 = name1+' [' + unit1 + ']'
    title2 = name2+' [' + unit2 + ']'

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
    plt.savefig(os.path.join(outDir, 'Scatterkde_%s_vs_%s_dist%s_%s.png' % (resType1, resType2, cfg['distType'], cfg['scenario'])))

    # shpw plot
    if flagShow:
        plt.show()

    plt.close('all')


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
    for val in [0.95, 0.99]:
        ind = np.searchsorted(hist, val)
        ind = min(ind, np.size(hist)-1)
        ax1.plot(sortedDiffPlot, hist)
        ax1.hlines(hist[ind], 0, sortedDiffPlot[ind], linestyles='--', linewidths=0.5)
        ax1.vlines(sortedDiffPlot[ind], 0, hist[ind], linestyles='--', linewidths=0.5)
        ticks.append(sortedDiffPlot[ind])

    ax2.set_xlim([-sortedDiffPlot[ind], sortedDiffPlot[ind]])
    width = diffMax - diffMin
    stepWidth = 2*sortedDiffPlot[ind]/50     # 50 bins in the [-2sigma,+2sigma] interval
    bins = int(width/stepWidth)
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
