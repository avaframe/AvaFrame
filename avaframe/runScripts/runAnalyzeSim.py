"""
    Compare the com1DFA kernel results to the similarity solution
    This script computes the similarity solution for a gliding avalanche on
    a inclined plane according to similarity solution from :
    Hutter, K., Siegel, M., Savage, S.B. et al.
    Two-dimensional spreading of a granular avalanche down an inclined plane
    Part I. theory. Acta Mechanica 100, 37–68 (1993).
    https://doi.org/10.1007/BF01176861
    and compares it to the DFA kernel com1DFA
"""

import pathlib
import pandas as pd

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in1Data import getInput as gI
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.outAna1Plots as outAna1Plots
import avaframe.ana4Stats.anaCompare as anaCompare


import numpy as np
import matplotlib.pyplot as plt

# local imports
import avaframe.out3Plot.plotUtils as pU


def plotAimecRes(simDF, outDirTest, xField, yFieldArray, coloredBy, sizedBy, logScaleX=False, logScaleY=False, fit=False):
    """plot error between all com1DFA sol and analytic sol
    function of whatever you want

    Parameters
    -----------
    simDF: dataFrame
        the simulation data with the postprocessing results
    outDirTest: str or pathlib
        output directory
    cfgSimi: configparser
        the cfg
    xField: str
        column of the simDF to use for the x axis
    yFieldArray: list
        list of max 2 column of the simDF to use for the y axis
    coloredBy: str
        column of the simDF to use for the colors
    sizedBy: str
        column of the simDF to use for the marker size
    logScale: boolean
        If you want a loglog scale
    """
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapAvaframeCont, min(simDF[coloredBy]), max(simDF[coloredBy]), continuous=pU.contCmap)
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(2*pU.figW, 1.5*pU.figH))
    # fig1, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(4*pU.figW, 2*pU.figH))
    # get the sizing function
    sizeList = simDF[sizedBy].unique()
    lenSize = len(sizeList)
    minSize = np.nanmin(sizeList)
    maxSize = np.nanmax(sizeList)
    if lenSize > 1:
        sizeList = (simDF[sizedBy].to_numpy() - minSize) / (maxSize - minSize) * 70 + 10
    else:
        sizeList = np.array([100])
    # make the scatter plot
    scatter1 = ax1.scatter(simDF[xField], simDF[yFieldArray[0]], c=simDF[coloredBy], s=sizeList, cmap=cmap, norm=norm,
                          marker=pU.markers[0], alpha=1)#, edgecolors='k')
    # scatter2 = ax2.scatter(simDF[xField], simDF[yFieldArray[1]], c=simDF[coloredBy], s=sizeList, cmap=cmap, norm=norm,
    #                        marker=pU.markers[0], alpha=1)#, edgecolors='k')

    if logScaleX:
        ax1.set_xscale('log')
        # ax2.set_xscale('log')
    if logScaleY:
        ax1.set_yscale('log')
        # ax2.set_yscale('log')

    ax1.set_title(yFieldArray[0])
    # ax2.set_title(yFieldArray[1])
    ax1.set_xlabel(xField)
    # ax2.set_xlabel(xField)
    ax1.set_ylabel(yFieldArray[0])
    # ax2.set_ylabel(yFieldArray[1])
    if len(simDF[sizedBy].unique())<=10:
        lenColor = None
    legend1 = ax1.legend(*scatter1.legend_elements(num=lenColor), loc="upper center", title=coloredBy)
    ax1.add_artist(legend1)
    # legend2 = ax2.legend(*scatter2.legend_elements(num=lenColor), loc="upper center", title=coloredBy)
    # ax2.add_artist(legend2)

    # produce a legend with a cross section of sizes from the scatter
    if lenSize<=10:
        lenSize = None
    kw = dict(prop="sizes", color=scatter1.cmap(0.7),
          func=lambda s: (s-10)*(maxSize - minSize)/70 + minSize)
    legend3 = ax1.legend(*scatter1.legend_elements(num=lenSize, **kw), loc="upper right", title=sizedBy)
    # kw = dict(prop="sizes", color=scatter2.cmap(0.7),
          # func=lambda s: (s-10)*(maxSize - minSize)/70 + minSize)
    # ax2.legend(*scatter2.legend_elements(num=lenSize, **kw), loc="upper right", title=sizedBy)
    ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    ax1.grid(color='grey', which='minor', linestyle='--', linewidth=0.25, alpha=0.5)
    # ax2.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    # ax2.grid(color='grey', which='minor', linestyle='--', linewidth=0.25, alpha=0.5)
    pU.saveAndOrPlot({'pathResult': outDirTest / 'pics'}, 'aimecResPlot', fig1)

    return fig1, ax1




# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# Load general configuration
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('Make posprocessing plots')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for similarity solution test
cfg = cfgUtils.getModuleConfig(anaCompare)

# create output directory for test result plots
outDirTest = pathlib.Path(avalancheDir, 'Outputs', 'outPlots')
fU.makeADir(outDirTest)


simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# Define release thickness distribution
pathToAimecRes = pathlib.Path(avalancheDir, 'Outputs', 'ana3AIMEC', 'com1DFA')
configFiles = list(pathToAimecRes.glob('*.csv'))
with open(configFiles[1], 'rb') as file:
    resultDF = pd.read_csv(file, index_col=0, keep_default_na=False)

simDF = simDF.reset_index().merge(resultDF, on='simName').set_index('index')

# filter on viscosity parameter
# simDF = simDF[simDF['subgridMixingFactor'].isin([10])]
# filter on time stepping parameter
simDF = simDF[simDF['cMax'].isin([0.01])]
# filter on nPPK0
# simDF = simDF[simDF['aPPK']==-2]
simDF = simDF[simDF['nPPK0']==15]
simDF = simDF[simDF['distReproj']==0]
simDF = simDF[simDF['massPerParticleDeterminationMethod']=='MPPKR']

# now do some plotting
# compare the simulations to the reference
# outAna1Plots.plotErrorRef(simDF, outDirTest, cfg['SIMISOL'], 'subgridMixingFactor', ['hErrorL2', 'vhErrorL2'],
#                           'deltaTh', 'dt', logScale=False)

# make convergence plot
plotAimecRes(simDF, outDirTest, 'nPart', ['maxpfvCrossMax', 'sRunout'],
                          'aPPK', 'nPPK0', logScaleX=True, logScaleY=False, fit=False)

outAna1Plots.plotTimeCPULog(simDF, outDirTest, cfg['MAIN'], 'nPart', 'aPPK', 'nPPK0')
