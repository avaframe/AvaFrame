"""
    plot statistics of simulations

"""

# load python modules
import os
import numpy as np
import matplotlib.pyplot as plt
from avaframe.out3Plot.plotUtils import *
import avaframe.out3Plot.makePalette as makePalette
import seaborn as sns
import pandas as pd


def plotValuesScatter(peakDictList, resType1, resType2, varPar, cfg, avalancheDir, flagShow=False):
    """ Produce scatter plot of max values (resType1 und resType2), for one set of simulations or multiple

        Parameters
        -----------
        peakDictList: list
            list of peakValues dictionaries that contain max values of peak parameters and parameter variation info
        resType1: str
            result parameter 1, 'ppr', 'pfd', 'pv'
        resType2: str
            result parameter 1, 'ppr', 'pfd', 'pv'
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
    if len(peakDictList) > 1:
        for peakValues in peakDictList:
            for key in peakValues:
                values1.append(peakValues[key][resType1])
                values2.append(peakValues[key][resType2])
                varVal.append(peakValues[key]['varPar'])
    else:
        peakValues = peakDictList[0]
        for key in peakValues:
            values1.append(peakValues[key][resType1])
            values2.append(peakValues[key][resType2])
            varVal.append(peakValues[key]['varPar'])

    log.info('Number of simulations is: %d' % (len(varVal)))

    # Get name and units for resTypes and parameter variation to annotate plots
    name1 = cfgPlotUtils['name%s' % resType1]
    name2 = cfgPlotUtils['name%s' % resType2]
    unit1 = cfgPlotUtils['unit%s' % resType1]
    unit2 = cfgPlotUtils['unit%s' % resType2]
    nameVar = cfgPlotUtils['name%s' % varPar]
    unitVar = cfgPlotUtils['unit%s' % varPar]
    varValV = np.array(varVal)

    # load variation colormap
    cmap, _, _, norm, ticks = makePalette.makeColorMap(
        cmapVar, np.amin(varValV), np.amax(varValV), continuous=True)

    fig, ax = plt.subplots()
    plt.title('%s vs. %s' % (name1, name2))
    # set range and steps of colormap
    cc = varValV
    sc = ax.scatter(values1, values2, c=cc, cmap=cmap)
    ax.set_xlabel('%s [%s]' % (name1, unit1))
    ax.set_ylabel('%s [%s]' % (name2, unit2))
    addColorBar(sc, ax, ticks, unitVar, nameVar)
    putAvaNameOnPlot(ax, avalancheDir)

    # shpw plot
    if flagShow:
        plt.show()

    # save fig
    outDir = cfg['outDir']
    fig.savefig(os.path.join(outDir, 'Scatter_%s_vs_%s_dist%s_%s.png' % (resType1, resType2, cfg['distType'], cfg['scenario'])))
    plt.close('all')


def plotValuesScatterHist(peakDictList, resType1, resType2, varPar, cfg, avalancheDir, flagShow=False, flagHue=False):
    """ Produce scatter and marginal kde plot of max values, for one set of simulations or multiple

        Parameters
        -----------
        peakDictList: list
            list of peakValues dictionaries that contain max values of peak parameters and parameter variation info
        resType1: str
            result parameter 1, 'ppr', 'pfd', 'pv'
        resType2: str
            result parameter 1, 'ppr', 'pfd', 'pv'
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
    if len(peakDictList) > 1:
        for peakValues in peakDictList:
            for key in peakValues:
                values1.append(peakValues[key][resType1])
                values2.append(peakValues[key][resType2])
                varVal.append(peakValues[key]['varPar'])
                if flagHue:
                    scenario.append(peakValues[key]['scenario'])

    else:
        peakValues = peakDictList[0]
        for key in peakValues:
            if 'scenario' in peakValues:
                flagHue = True
            values1.append(peakValues[key][resType1])
            values2.append(peakValues[key][resType2])
            varVal.append(peakValues[key]['varPar'])
            if flagHue:
                scenario.append(peakValues[key]['scenario'])

    log.info('Number of simulations is: %d' % (len(varVal)))

    # Get name and units for resTypes and parameter variation to annotate plots
    name1 = cfgPlotUtils['name%s' % resType1]
    name2 = cfgPlotUtils['name%s' % resType2]
    unit1 = cfgPlotUtils['unit%s' % resType1]
    unit2 = cfgPlotUtils['unit%s' % resType2]
    nameVar = cfgPlotUtils['name%s' % varPar]
    unitVar = cfgPlotUtils['unit%s' % varPar]
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
    putAvaNameOnPlot(fig1.ax_joint, avalancheDir)

    # save fig
    outDir = cfg['outDir']
    plt.savefig(os.path.join(outDir, 'Scatterkde_%s_vs_%s_dist%s_%s.png' % (resType1, resType2, cfg['distType'], cfg['scenario'])))

    # shpw plot
    if flagShow:
        plt.show()

    plt.close('all')
