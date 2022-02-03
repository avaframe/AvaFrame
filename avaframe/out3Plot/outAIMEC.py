"""
    Plotting and saving AIMEC results

    This file is part of Avaframe.
"""

import os
import logging
import math
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib
from matplotlib import cm


# Local imports
import avaframe.out3Plot.plotUtils as pU
from avaframe.out3Plot import statsPlots as sPlot
from avaframe.in3Utils import cfgUtils

# create local logger
log = logging.getLogger(__name__)


def visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict):
    """
    Plot and save the domain transformation figure
    The left subplot shows the reference result raster with the outline of the
    new domain. The second one shows this same data in the (s,l) coordinate
    system define by the outline in the first plot.

    Parameters
    ----------
    rasterTransfo: dict
        domain transformation information
    inputData : dict
        inputData dictionary:
            slRaster: numpy array with (s,l) raster
            xyRaster: numpy array with (x,y) raster
            headerXY: header corresponding to xyRaster
    cfgSetup : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini
        'resType'
    pathDict : dict
        dictionary with path to data to analyze

    """
    ####################################
    # Get input data
    resType = cfgSetup['resType']
    unit = pU.cfgPlotUtils['unit' + resType]
    # read paths
    projectName = pathDict['projectName']
    # read rasterdata
    slRaster = inputData['slRaster']
    xyRaster = inputData['xyRaster']
    headerXY = inputData['xyHeader']
    cellsize = headerXY['cellsize']
    n = headerXY['nrows']
    m = headerXY['ncols']
    xllc = headerXY['xllcenter']
    yllc = headerXY['yllcenter']
    # read avaPath with scale
    xPath = rasterTransfo['x']
    yPath = rasterTransfo['y']
    # read domain boundarries with scale
    cellSizeSL = rasterTransfo['cellSizeSL']
    DBXl = rasterTransfo['DBXl']*cellSizeSL
    DBXr = rasterTransfo['DBXr']*cellSizeSL
    DBYl = rasterTransfo['DBYl']*cellSizeSL
    DBYr = rasterTransfo['DBYr']*cellSizeSL

    ############################################
    # prepare for plot
    x = np.arange(m)*cellsize + xllc
    y = np.arange(n)*cellsize + yllc
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    xx = rasterTransfo['x'][indStartOfRunout]
    yy = rasterTransfo['y'][indStartOfRunout]

    l = rasterTransfo['l']
    s = rasterTransfo['s']

    maskedArray = np.ma.masked_where(xyRaster == 0, xyRaster)
    maskedArraySL = np.ma.masked_where(slRaster == 0, slRaster)

    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin(
        maskedArray), np.nanmax(maskedArray), continuous=pU.contCmap)
    cmap.set_under(color='w')

    ############################################
    # Figure: Raster transformation
    fig = plt.figure(figsize=(pU.figW*2, pU.figH))

    ax1 = plt.subplot(121)

    ref0, im = pU.NonUnifIm(ax1, x, y, maskedArray, 'x [m]', 'y [m]',
                            extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, norm=norm)
    plt.plot(xx, yy, 'k+', label='start of run-out area point : %.1f °' %
             rasterTransfo['startOfRunoutAreaAngle'])
    plt.plot(xPath, yPath, 'k--', label='flow path')
    plt.plot(DBXl, DBYl, 'k-', label='domain')
    plt.plot(DBXr, DBYr, 'k-')
    plt.plot([DBXl, DBXr], [DBYl, DBYr], 'k-')

    ax1.set_title('XY Domain')
    ax1.legend(loc=4)
    pU.putAvaNameOnPlot(ax1, pathDict['projectName'])

    ax2 = plt.subplot(122)

    ref0, im = pU.NonUnifIm(ax2, l, s, maskedArraySL, 'l [m]', 's [m]',
                            extent=[l.min(), l.max(), s.min(), s.max()],
                            cmap=cmap, norm=norm)
    ax2.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                label='start of run-out area point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])

    ax2.set_title('sl Domain' + '\n' + 'Black = out of raster')
    ax2.legend(loc=4)
    pU.addColorBar(im, ax2, ticks, unit)

    outFileName = '_'.join([projectName, 'DomainTransformation'])
    pU.saveAndOrPlot(pathDict, outFileName, fig)


def visuRunoutComp(rasterTransfo, resAnalysisDF, cfgSetup, pathDict):
    """
    Plot and save the Peak Fields distribution (max mean per cross section)
    after coordinate transformation

    Parameters
    ----------
    rasterTransfo: dict
        domain transformation information
    resAnalysis: dict
        results from Aimec analysis (for ppr, pfd and pfv):
            PPRCrossMax: numpy array with max peak field along path for each file to analyse
            PPRCrossMean: numpy array with mean peak field along path for each file to analyse
    cfgSetup : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini
        'resType' and 'thresholdValue'
    pathDict : dict
        dictionary with path to data to analyze
    """
    ####################################
    # Get input data
    resType = cfgSetup['resType']
    thresholdValue = cfgSetup['thresholdValue']
    # read paths
    projectName = pathDict['projectName']
    refSimName = pathDict['refSimulation']
    simName = pathDict['compSimulation']
    # read data
    s = rasterTransfo['s']
    l = rasterTransfo['l']

    ############################################
    # prepare for plot
    title = ['Pressure ', 'Flow Depth ', 'Flow Velocity ']
    unit = ['$P(s)$ [kPa]', '$fd(s)$ [m]', '$v(s) [m.s^{-1}]$']
    peakList = ['ppr', 'pfd', 'pfv']

    ############################################
    # Figure: Pressure depth speed

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(pU.figW*3, pU.figH))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.3)

    for ax, peak, titleVal, unitVal in zip(axes.flatten(), peakList, title, unit):
        ax.plot(resAnalysisDF.loc[refSimName, peak + 'CrossMax'], s, '--k', label='Max Reference')
        ax.plot(resAnalysisDF.loc[refSimName, peak + 'CrossMean'], s, '-k', label='Mean Reference')
        ax.plot(resAnalysisDF.loc[simName, peak + 'CrossMax'], s, '--b', label='Max Simulation')
        ax.plot(resAnalysisDF.loc[simName, peak + 'CrossMean'], s, '-b', label='Mean Simulation')

        ax.set_title(titleVal + 'distribution along path')
        ax.legend(loc='best')
        ax.set_ylabel('s [m]')
        ax.set_ylim([s.min(), s.max()])
        ax.set_xlim(auto=True)
        ax.set_xlabel(unitVal)
    pU.putAvaNameOnPlot(ax, projectName)

    outFileName = '_'.join([projectName, resType,
                            'thresholdValue', str(thresholdValue).replace('.', 'p'), 'slComparison'])
    outFilePath = pU.saveAndOrPlot(pathDict, outFileName, fig)

    return outFilePath


def visuRunoutStat(rasterTransfo, inputsDF, resAnalysisDF, newRasters, cfgSetup, pathDict):
    """
    Plot and save the Peak field  distribution after coord transfo
    used when more then 2 simulations are compared

    Parameters
    ----------
    rasterTransfo: dict
        domain transformation information
    inputsDF: dataFrame
        aimec inputs DF
    resAnalysis: dataFrame
        results from Aimec analysis:
            numpy array with max peak field (of the 'resType') along path for each file to analyse
            runout for 'resType'
    newRasters: dict
        dictionary with new (s, l) raster for the 'resType'
    cfgSetup : configparser
        configparser with ana3AIMEC settings defined in ana3AIMECCfg.ini
        'percentile', 'resType' and 'thresholdValue'
    pathDict : dict
        dictionary with path to data to analyze
    """
    ####################################
    # Get input data
    varParList = cfgSetup['varParList'].split('|')
    paraVar = varParList[0]
    percentile = cfgSetup.getfloat('percentile')
    resType = cfgSetup['resType']
    thresholdValue = cfgSetup['thresholdValue']
    unit = pU.cfgPlotUtils['unit' + resType]
    name = pU.cfgPlotUtils['name' + resType]
    # read paths
    projectName = pathDict['projectName']
    # read data
    s = rasterTransfo['s']
    l = rasterTransfo['l']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    rasterdataPres = newRasters['newRefRaster' + resType.upper()]
    runout = resAnalysisDF['sRunout'].to_numpy()
    pprCrossMax = np.stack(resAnalysisDF[resType.lower() + 'CrossMax'].to_numpy())
    ############################################
    # prepare for plot
    pMean = np.mean(pprCrossMax, axis=0)
    pMedian = np.median(pprCrossMax, axis=0)
    pPercentile = np.percentile(pprCrossMax, [percentile/2, 50, 100-percentile/2], axis=0)
    maskedArray = np.ma.masked_where(rasterdataPres <= float(thresholdValue), rasterdataPres)

    # get plots limits
    indYMin = max(0, indStartOfRunout-5)
    yMin = s[indYMin]
    yMax = max(runout) + 25
    indXMin = max(0, np.min(np.nonzero(np.any(maskedArray[indStartOfRunout:, :] > 0, axis=0))[0])-5)
    xMin = l[indXMin]
    indXMax = min(np.max(np.nonzero(np.any(maskedArray[indStartOfRunout:, :] > 0, axis=0))[0])+5, len(l)-1)
    xMax = l[indXMax]

    # get colormap for raster plot
    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin(
        maskedArray[indYMin:, indXMin:indXMax]), np.nanmax(maskedArray[indYMin:, indXMin:indXMax]), continuous=pU.contCmap)
    cmap.set_bad('w', 1.)

    # Get colors for scatter
    unitSC = cfgSetup['unit']
    nSamples = np.size(runout)
    if 'colorParameter' in pathDict:
        if pathDict['colorParameter'] is False:
            values = None
        else:
            values = inputsDF[varParList[0]]
    else:
        values = None
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(values, nSamples, unitSC)

    ############################################
    # Figure: Analysis runout
    fig = plt.figure(figsize=(pU.figW*2, pU.figH))
    ax1 = plt.subplot(121)

    ax1.axhline(y=np.max(runout), color='k', linestyle='-.', label='runout max %.0f m' % np.max(runout))
    ax1.axhline(y=np.average(runout), color='k', linestyle='-', label='runout mean %.0f m' % np.mean(runout))
    ax1.axhline(y=np.min(runout), color='k', linestyle=':', label='runout min %.0f m' % np.min(runout))

    ax1.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                label='start of run-out area point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])
    ref5, im = pU.NonUnifIm(ax1, l, s, maskedArray, 'l [m]', 's [m]',
                            extent=[xMin, xMax, yMin, yMax],
                            cmap=cmap, norm=norm)
    sc = ax1.scatter(resAnalysisDF['lRunout'], resAnalysisDF['sRunout'],
                     c=colorSC, cmap=cmapSC, norm=normSC, marker=pU.markers[0], label='runout points')

    if displayColorBar:
        pU.addColorBar(sc, ax1, ticksSC, unitSC, title=paraVar, pad=0.08, tickLabelsList=itemsList)
    ax1.set_xlim([xMin, xMax])
    ax1.set_ylim([yMin, yMax])
    ax1.set_title('Peak Pressure 2D plot for the reference')
    ax1.legend(loc=4)
    pU.putAvaNameOnPlot(ax1, projectName)

    pU.addColorBar(im, ax1, ticks, unit)

    ax2 = plt.subplot(122)

    ax2.fill_betweenx(s, pPercentile[2], pPercentile[0], facecolor=[.8, .8, .8], alpha=0.5,
                      label=('[%.2f, %.2f]%% interval' % (percentile/2, 100-percentile/2)))
    matplotlib.patches.Patch(alpha=0.5, color=[.8, .8, .8])
    ax2.plot(pMedian, s, color='r', label='median')
    ax2.plot(pMean, s, color='b', label='mean')

    ax2.set_title('%s distribution along the path between runs' % name)
    ax2.legend(loc='upper right')
    ax2.set_ylabel('s [m]')
    ax2.set_ylim([s.min(), s.max()])
    ax2.set_xlim(auto=True)
    ax2.set_xlabel('$P_{max}(s)$ [%s]' % unit)

    outFileName = '_'.join([projectName, resType, 'thresholdValue', str(thresholdValue).replace('.', 'p'),
                           'slComparisonStat'])

    pU.saveAndOrPlot(pathDict, outFileName, fig)


def visuMass(resAnalysisDF, pathDict, simName, refSimulationName, timeMass):
    """
    Plot and save the results from mass analysis

    Parameters
    ----------
    resAnalysis: dict
        mass results from Aimec analysis:
            entMassFlowArray: entrained mass array corresponding to the time array for each simulation to analyse
            totalMassArray: entrained mass array corresponding to the time array for each simulation to analyse
            entMass: final entrained for each simulation
            finalMass: final mass for each simulation
            time: time array
    pathDict : dict
        dictionary with path to data to analyze
    """
    ####################################
    # Get input data
    # read paths
    projectName = pathDict['projectName']
    # read data
    entMassRef = resAnalysisDF.loc[refSimulationName, 'entMass']
    finalMassRef = resAnalysisDF.loc[refSimulationName, 'finalMass']
    entMass = resAnalysisDF.loc[simName, 'entMass']
    finalMass = resAnalysisDF.loc[simName, 'finalMass']

    ############################################
    # prepare for plot
    Title = ['Entrained Mass Flow', 'Total Mass']
    Unit = ['Entrained Mass Flow [$kg.s{^-1}$]', 'Total Mass [kg]']
    fieldList = ['entMassFlowArray', 'totalMassArray']
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(pU.figW*2, pU.figH))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.3)

    for ax, field, title, unit in zip(axes.flatten(), fieldList, Title, Unit):
        refArray = resAnalysisDF.loc[refSimulationName, field]
        simArray = resAnalysisDF.loc[simName, field]
        ax.plot(timeMass, refArray, '-k', label='Reference : %s ' % refSimulationName)
        ax.plot(timeMass, simArray, '-b', label='Simulation : %s ' % simName)

        ax.set_title(title + ' function of time')
        ax.legend(loc=1)
        ax.set_xlabel('t [s]')
        ax.set_ylabel(unit)

    ax2 = axes.flatten()[1].twinx()
    ax2.spines['right'].set_color('r')
    ax2.tick_params(axis='y', colors='r')
    ax2.plot(timeMass, (simArray-refArray) / refArray*100, 'r', label='total mass')

    if np.any(entMass):
        axes.flatten()[1].text(timeMass[-1]/4, (np.nanmin(refArray) + np.nanmax(refArray))/2,
                               'Entrained Mass Difference : %.2f kg \n Relative to total mass : %.2f %% ' %
                               ((entMassRef-entMass), (entMassRef-entMass)/finalMassRef*100),
                               bbox=dict(boxstyle="square", ec='white', fc='white'),
                               horizontalalignment='left', verticalalignment='bottom')

    ax2.set_ylabel('Entrained Mass Difference relative to total mass[%]', color='r')

    outFileName = '_'.join([projectName, 'massAnalysis', str(simName)])
    pU.putAvaNameOnPlot(ax2, pathDict['projectName'])
    outFilePath = pU.saveAndOrPlot(pathDict, outFileName, fig)

    return outFilePath


def visuSimple(cfgSetup, rasterTransfo, resAnalysisDF, newRasters, pathDict):
    """
    Plot and save the Peak Pressure Peak Flow depth and Peak speed
    fields after coord transfo

    Parameters
    ----------
    rasterTransfo: dict
        domain transformation information
    resAnalysis: dict
        results from Aimec analysis:
            numpy array with the 'runout' for each simulation and the 'thresholdValue'
    newRasters: dict
        dictionary with new (s, l) raster for ppr, pfd and pfd
    pathDict : dict
        dictionary with path to data to analyze
    """
    ####################################
    # Get input data
    # read paths
    projectName = pathDict['projectName']
    refSimulationName = pathDict['refSimulation']
    # read data
    plim = cfgSetup['thresholdValue']
    s = rasterTransfo['s']
    l = rasterTransfo['l']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    rasterdataPres = newRasters['newRefRasterPPR']
    rasterdataDepth = newRasters['newRefRasterPFD']
    rasterdataSpeed = newRasters['newRefRasterPFV']
    runout = resAnalysisDF.loc[refSimulationName, 'sRunout']

    ############################################
    # prepare for plot
    Cmap = [pU.cmapPres, pU.cmapDepth, pU.cmapSpeed]
    Title = ['Peak Pressure', 'Peak Flow Depth', 'Peak Speed']
    Unit = [pU.cfgPlotUtils['unitppr'], pU.cfgPlotUtils['unitpfd'], pU.cfgPlotUtils['unitpfv']]
    Data = np.array(([None] * 3))
    Data[0] = rasterdataPres
    Data[1] = rasterdataDepth
    Data[2] = rasterdataSpeed

    ############################################
    # Figure: Pressure depth speed

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(pU.figW*3, pU.figH))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.3)

    for ax, cmap, data, title, unit in zip(axes.flatten(), Cmap, Data, Title, Unit):
        maskedArray = np.ma.masked_where(data == 0, data)
        cmap, _, ticks, norm = pU.makeColorMap(cmap, np.nanmin(maskedArray), np.nanmax(maskedArray),
                                               continuous=pU.contCmap)
        cmap.set_bad('w', 1.)
        ax.axhline(y=runout, color='k', linestyle='-', label='runout')

        ax.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                   label='Start or run-out point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])
        ref3, im = pU.NonUnifIm(ax, l, s, maskedArray, 'l [m]', 's [m]',
                                extent=[l.min(), l.max(), s.min(), s.max()],
                                cmap=cmap, norm=norm)

        ax.set_title(title)
        ax.legend(loc=4)
        pU.addColorBar(im, ax, ticks, unit)
        pU.putAvaNameOnPlot(ax, pathDict['projectName'])

    outFileName = '_'.join([projectName, 'plim', str((plim)).replace('.', 'p'), 'referenceFields'])

    pU.saveAndOrPlot(pathDict, outFileName, fig)


def visuComparison(rasterTransfo, inputs, pathDict):
    """
    Plot and save the comparison between current simulation and Reference
    in the run-out area

    Parameters
    ----------
    rasterTransfo: dict
        domain transformation information
    inputs: dict
        input data for plot:
            'refData' and 'compData' arrays
            'refRasterMask' and 'compRasterMask' arrays
            'thresholdArray'
            'i' id of the simulation
            'resType' result type analyzed
            'diffLim'
    pathDict : dict
        dictionary with path to data to analyze
    """
    ####################################
    # Get input data
    # read paths
    projectName = pathDict['projectName']
    # read data
    s = rasterTransfo['s']
    l = rasterTransfo['l']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    refData = inputs['refData']
    compData = inputs['compData']
    refRasterMask = inputs['refRasterMask']
    compRasterMask = inputs['compRasterMask']
    simName = inputs['simName']
    refSimulation = pathDict['refSimulation']
    resType = inputs['resType']
    unit = pU.cfgPlotUtils['unit' + resType]
    name = pU.cfgPlotUtils['name' + resType]
    thresholdArray = inputs['thresholdArray']
    thresholdValue = thresholdArray[-1]

    yLimRef = s[np.max(np.nonzero(np.any(refData > 0, axis=1))[0])] + 20
    yLim = s[max(np.max(np.nonzero(np.any(refData > 0, axis=1))[0]),
              np.max(np.nonzero(np.any(compData > 0, axis=1))[0]))] + 20
    ############################################
    # Figure: Raster comparison (mask for the pThreshold given in the ini file)
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))
    ax1 = plt.subplot2grid((1, 2), (0, 0))

    # get color map
    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin(
        (refData)), np.nanmax((refData)), continuous=pU.contCmap)
    cmap.set_bad(color='w')
    refDataPlot = np.ma.masked_where(refData == 0.0, refData)
    ref0, im = pU.NonUnifIm(ax1, l, s, refDataPlot, 'l [m]', 's [m]',
                            extent=[l.min(), l.max(), 0, yLimRef],
                            cmap=cmap, norm=norm)
    ax1.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                label='start of run-out area point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])
    ax1.set_title('Reference %s in the RunOut area' % name
                  + '\n' + '%s threshold: %.1f %s' % (name, thresholdValue, unit))
    pU.addColorBar(im, ax1, ticks, unit)

    ax1.set_ylim([0, yLimRef])
    ax1.legend(loc='lower right')
    pU.putAvaNameOnPlot(ax1, projectName)

    ax2 = plt.subplot2grid((1, 2), (0, 1))
    cmap.set_under(color='b')
    cmap.set_over(color='r')
    cmap.set_bad(alpha=0)
    data = compRasterMask-refRasterMask
    data = np.ma.masked_where(data == 0.0, data)
    if compRasterMask[indStartOfRunout:, :].any() or refRasterMask[indStartOfRunout:, :].any():
        ref1, im1 = pU.NonUnifIm(ax2, l, s, data, 'l [m]', 's [m]',
                                 extent=[l.min(), l.max(), s[indStartOfRunout], yLim], cmap=cmap)
        im1.set_clim(vmin=-0.5, vmax=0.5)
        ax2.set_ylim([s[indStartOfRunout], yLim])
    else:
        ax2.text(.5, .5, 'No data in the run out area!', fontsize=18, color='red',
                 bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')
    if pathDict['compType'][0] == 'comModules':
        namePrint = 'refMod:' + pathDict['compType'][1] + '_' + 'compMod:' + pathDict['compType'][2]
        pU.putAvaNameOnPlot(ax2, namePrint)
    else:
        namePrint = 'ref:' + str(refSimulation) + '_' + 'sim:' + str(simName)
        pU.putAvaNameOnPlot(ax2, namePrint)

    ax2.set_title('Difference %s current - reference in runout area' % resType + '\n' + 'Blue = FN, Red = FP')

    outFileName = '_'.join([projectName, 'thresholdValue', str(
                  thresholdValue).replace('.', 'p'),  'sim', str(simName), 'AreaComparisonToReference'])
    pU.saveAndOrPlot(pathDict, outFileName, fig)

    ############################################
    # Figure: Raster comparison
    # , constrained_layout=True)
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))
    ax1 = plt.subplot2grid((3, 3), (0, 0), rowspan=3)
    # get color map
    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin(
        (refData)), np.nanmax((refData)), continuous=pU.contCmap)
    cmap.set_bad(color='w')
    refDataPlot = np.ma.masked_where(refData == 0.0, refData)
    ref0, im = pU.NonUnifIm(ax1, l, s, refDataPlot, 'l [m]', 's [m]',
                            extent=[l.min(), l.max(), 0, yLimRef],
                            cmap=cmap, norm=norm)
    ax1.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                label='start of run-out area point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])
    ax1.set_title('Reference %s' % name)
    pU.addColorBar(im, ax1, ticks, unit)
    ax1.set_ylim([0, yLimRef])
    ax1.legend(loc='lower right')
    pU.putAvaNameOnPlot(ax1, projectName)

    compData = compData[indStartOfRunout:, :]
    refData = refData[indStartOfRunout:, :]
    dataDiff = compData - refData
    dataDiff = np.where((refData == 0) & (compData == 0), np.nan, dataDiff)
    dataDiffPlot = np.where((refData < thresholdArray[-1]) & (compData < thresholdArray[-1]), np.nan, dataDiff)
    dataDiffPlot = dataDiffPlot[~np.isnan(dataDiffPlot)]

    if dataDiffPlot.size:
        # only add the second axis if one of the two avalanches reached the run out area
        indDiff = np.abs(dataDiffPlot) > 0

        if indDiff.any():
            # only plot hist and CDF if there is a difference in the data
            ax2 = plt.subplot2grid((3, 3), (0, 1), rowspan=2, colspan=2)
        else:
            ax2 = plt.subplot2grid((3, 3), (0, 1), rowspan=3, colspan=3)

        cmap = pU.cmapdiv
        cmap.set_bad(color='w')
        elev_max = inputs['diffLim']
        ref0, im3 = pU.NonUnifIm(ax2, l, s[indStartOfRunout:], (dataDiff), 'l [m]', 's [m]',
                                 extent=[l.min(), l.max(), s[indStartOfRunout], yLim], cmap=cmap)
        im3.set_clim(vmin=-elev_max, vmax=elev_max)

        # print contour lines only if the thre threshold is reached
        L, S = np.meshgrid(l, s[indStartOfRunout:])
        colorsP = pU.colorMaps['pfd']['colors'][1:]
        if (np.where(refData > thresholdArray[-1], True, False)).any():
            contourRef = ax2.contour(L, S, refData, levels=thresholdArray[:-1], linewidths=2, colors=colorsP)
            # generate corresponding labels
            labels = [str(level) for level in thresholdArray[:-1]]
            labels = labels[0:len(contourRef.collections)]
            # add legend associated to the contour plot
            handles, _ = contourRef.legend_elements()
            legend2 = ax2.legend(title=name + '\ncontour lines\n[' + unit + ']', handles=handles, labels=labels)
            plt.setp(legend2.get_title(), multialignment='center')
        else:
            log.warning('Reference %s did not reach the run out area!' % refSimulation)
            ax2.text(0, (s[indStartOfRunout] + yLim)/2, 'Reference %s did not reach the run out area!' % refSimulation,
                     fontsize=24, color='red',
                     bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')
        if (np.where(compData > thresholdArray[-1], True, False)).any():
            contourComp = ax2.contour(
                L, S, compData, levels=thresholdArray[:-1], linewidths=2, colors=colorsP, linestyles='dashed')
        else:
            log.warning('Simulation %d did not reach the run out area!' % simName)
            ax2.text(0, (s[indStartOfRunout] + yLim)/2, 'Simulation %d did not reach the run out area!' % simName,
                     fontsize=24, color='red',
                     bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')

        if pathDict['compType'][0] == 'comModules':
            namePrint = 'refMod:' + pathDict['compType'][1] + '_' + 'compMod:' + pathDict['compType'][2]
            pU.putAvaNameOnPlot(ax2, namePrint)
        else:
            namePrint = 'ref:' + str(refSimulation) + '_' + 'sim:' + str(simName)
            pU.putAvaNameOnPlot(ax2, namePrint)

        if indDiff.any():
            # only plot hist and CDF if there is a difference in the data
            ax3 = plt.subplot2grid((3, 3), (2, 1))
            ax4 = plt.subplot2grid((3, 3), (2, 2))
            # there is data to compare in the run out area
            centiles = sPlot.plotHistCDFDiff(dataDiffPlot, ax4, ax3, insert='False', title=['%s diff histogram' % name,
                                             '%s diff CDF (95%% and 99%% centiles)' % name])

        ax2.set_ylim([s[indStartOfRunout], yLim])
        pU.addColorBar(im3, ax2, None, unit, title=name, extend='both')
    else:
        # if no avalanche reached the run out area print a warning on the second plot
        ax2 = plt.subplot2grid((3, 3), (0, 1), rowspan=3, colspan=3)
        log.warning('No data in the run out area!')
        ax2.text(.5, .5, 'No data in the run out area!', fontsize=24, color='red',
                 bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')

    if pathDict['compType'][0] == 'comModules':
        ax2.set_title('%s difference and contour lines' % name + '\n' + 'refMod = full, compMod = dashed line')
    else:
        ax2.set_title('%s difference and contour lines' % name + '\n' + 'ref = full, sim = dashed line')

    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    outFileName = '_'.join([projectName, 'plim', str(thresholdValue).replace(
                  '.', 'p'),  'sim', str(simName), 'ContourComparisonToReference'])
    outFilePath = pU.saveAndOrPlot(pathDict, outFileName, fig)

    return outFilePath


def resultWrite(pathDict, cfg, rasterTransfo, resAnalysisDF):
    """
    This function writes the main Aimec results to a file (outputFile)
    in pathDict
    """
    ####################################
    # Get input data
    projectName = pathDict['projectName']
    pathResult = pathDict['pathResult']
    pathName = pathDict['pathName']
    refSimName = pathDict['refSimulation']
    demName = os.path.basename(pathDict['demSource'])
    # dataName = [os.path.basename(name) for name in pathDict['ppr']]
    cfgSetup = cfg['AIMECSETUP']
    cfgFlags = cfg['FLAGS']
    flagMass = cfgFlags.getboolean('flagMass')
    domainWidth = cfgSetup['domainWidth']
    thresholdValue = cfgSetup['thresholdValue']
    resType = cfgSetup['resType']
    unit = pU.cfgPlotUtils['unit' + resType]
    name = pU.cfgPlotUtils['name' + resType]

    startOfRunoutAreaAngle = rasterTransfo['startOfRunoutAreaAngle']
    areaSum = resAnalysisDF.loc[refSimName, 'TP'] + resAnalysisDF.loc[refSimName, 'FN']
    if areaSum == 0:
        log.warning('Reference did not reach the run-out area. Not normalizing area indicators')
    else:
        resAnalysisDF['TP'] = resAnalysisDF['TP'] / areaSum
        resAnalysisDF['FN'] = resAnalysisDF['FN'] / areaSum
        resAnalysisDF['FP'] = resAnalysisDF['FP'] / areaSum
        resAnalysisDF['TN'] = resAnalysisDF['TN'] / areaSum
    # do some statistics
    forStats = ['xRunout', 'yRunout', 'sRunout', 'elevRel', 'deltaH', 'maxpprCrossMax', 'maxpfdCrossMax',
                'maxpfvCrossMax', 'TP', 'FN', 'FP', 'TN']
    if flagMass:
        forStats = forStats + ['relMass', 'entMass', 'finalMass', 'relativMassDiff', 'growthIndex', 'growthGrad']
    forStats = set(forStats) & set(resAnalysisDF.columns)
    # compute some statistics
    resAnalysisStatsDF = resAnalysisDF[forStats].describe(percentiles=None)

    header = ''.join(['projectName: ', projectName, '\n',
                      'path: ', pathName, '\n',
                      'dhm: ', demName, '\n',
                      'domain_width: ', str(domainWidth), ' m\n',
                      'run out based on ', name, ' results\n',
                      name, ' limit: ', str(thresholdValue), unit, '\n',
                      'start of runout area Angle (SROA angle): ', str(round(startOfRunoutAreaAngle, 2)), ' °\n'])

    outFileName = '_'.join(['Results', projectName, str(resType), 'lim', str(thresholdValue), 'w',
                           str(domainWidth)]) + '.csv'
    outname = os.path.join(pathResult, outFileName)
    outFileNameStats = '_'.join(['Results', projectName, str(resType), 'lim', str(thresholdValue), 'w',
                                str(domainWidth)]) + 'stats.csv'
    outnameStats = os.path.join(pathResult, outFileNameStats)

    # check if folder exists / create
    if not os.path.exists(os.path.dirname(outname)):
        os.makedirs(os.path.dirname(outname))
    resAnalysisDF.to_csv(outname)

    if not os.path.exists(os.path.dirname(outnameStats)):
        os.makedirs(os.path.dirname(outnameStats))
    resAnalysisStatsDF.to_csv(outnameStats)

    log.info('File written: %s' % outname)
    log.info('File written: %s' % outnameStats)


def resultVisu(cfgSetup, inputsDF, pathDict, cfgFlags, rasterTransfo, resAnalysisDF):
    """
    Visualize results in a nice way
    """
    ####################################
    # Get input data
    resType = cfgSetup['resType']
    varParList = cfgSetup['varParList'].split('|')
    unit = cfgSetup['unit']
    paraVar = cfgSetup['varParList'].split('|')[0]
    name = pU.cfgPlotUtils['name' + resType]
    nSim = len(resAnalysisDF.index)
    refSimName = pathDict['refSimulation']
    maxMaxDPPRRef = resAnalysisDF.loc[refSimName, 'max' + resType.lower() + 'CrossMax']
    maxMaxDPPR = resAnalysisDF['max' + resType.lower() + 'CrossMax'].to_numpy()

    thresholdValue = cfgSetup['thresholdValue']

    flag = float(cfgFlags['typeFlag'])

    zPath = rasterTransfo['z']
    sPath = rasterTransfo['s']

    runout = resAnalysisDF['sRunout'].to_numpy()
    runoutRef = resAnalysisDF.loc[refSimName, 'sRunout']
    areaSum = resAnalysisDF.loc[refSimName, 'TP'] + resAnalysisDF.loc[refSimName, 'FN']
    if areaSum == 0:
        log.warning('Reference did not reach the run-out area. Not normalizing area indicators')
        areaSum = 1
    rTP = resAnalysisDF['TP'].to_numpy() / areaSum
    rFP = resAnalysisDF['FP'].to_numpy() / areaSum
    rTPRef = resAnalysisDF.loc[refSimName, 'TP'] / areaSum
    rFPRef = resAnalysisDF.loc[refSimName, 'FP'] / areaSum

    title = 'Visualizing max ' + name + ' data'
    tipo = 'relMax' + resType + '_thresholdValue' + str(thresholdValue).replace('.', 'p')
    data = maxMaxDPPR / maxMaxDPPRRef
    dataRef = maxMaxDPPRRef / maxMaxDPPRRef
    yaxis_label = 'relative max ' + name + ' [-]'

    log.info(title)

    # If more than 100 files are provided, add a density plot
    plotDensity = 0
    if (nSim > cfgFlags.getfloat('nDensityPlot')):
        plotDensity = 1
    # Get colors for scatter
    unitSC = cfgSetup['unit']
    nSamples = np.size(runout)
    if 'colorParameter' in pathDict:
        if pathDict['colorParameter'] is False:
            values = None
        else:
            values = inputsDF[varParList[0]]
    else:
        values = None
    cmapSC, colorSC, ticksSC, normSC, unitSC, itemsList, displayColorBar = pU.getColors4Scatter(values, nSamples, unitSC)
    #######################################
    # Final result diagram - z_profile+data

    fig = plt.figure(figsize=(pU.figW*2, pU.figH))

    # show flow path
    ax1 = fig.add_subplot(111)
    ax1.set_title(title)
    ax1.set_ylabel(yaxis_label, color=pU.cmapVar['colors'][-1])
    ax1.spines['left'].set_color(pU.cmapVar['colors'][-1])
    ax1.tick_params(axis='y', colors=pU.cmapVar['colors'][-1])
    ax1.set_xlabel(''.join(['s [m] - runout with ', str(thresholdValue), ' kPa threshold']))
    pU.putAvaNameOnPlot(ax1, pathDict['projectName'])
    if plotDensity:  # estimate 2D histogram --> create pcolormesh
        nbins = 10
        H, xedges, yedges = np.histogram2d(runout, data, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
        dataDensity = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.Blues)
        cbar = plt.colorbar(dataDensity, orientation='horizontal')
        cbar.ax.set_ylabel('Counts')

    ax2 = ax1.twinx()
    ax2.set_ylabel('z [m]')
    ax2.spines['left'].set_color(pU.cmapVar['colors'][-1])
    ax2.tick_params(axis='y', colors='k')
    ax2.plot(sPath, zPath, color='k', label='path', linestyle='--')
    plt.xlim([0, max(sPath) + 50])
    plt.ylim([math.floor(min(zPath)/10)*10, math.ceil(max(zPath)/10)*10])

    if not plotDensity:
        sc = ax1.scatter(resAnalysisDF['sRunout'], data, c=colorSC, cmap=cmapSC, norm=normSC, marker=pU.markers[0],
                         label='runout points')

        if displayColorBar:
            pU.addColorBar(sc, ax2, ticksSC, unitSC, title=paraVar, pad=0.08, tickLabelsList=itemsList)

        ax1.plot(runoutRef, dataRef, color='g', label='Reference', marker='+', markersize=2*pU.ms, linestyle='None')
        ax1.legend(loc=4)

    ax1.grid('on')

    outFileName = '_'.join([pathDict['projectName'], tipo])

    pU.saveAndOrPlot(pathDict, outFileName, fig)

    ############################################
    # Final result diagram - roc-plots

    fig = plt.figure(figsize=(pU.figW, pU.figH))
    ax1 = fig.add_subplot(111)
    ax1.set_title('Normalized difference compared to reference')
    ax1.set_ylabel('True positive rate')
    ax1.set_xlabel('False positive rate')
    pU.putAvaNameOnPlot(ax1, pathDict['projectName'])
    if plotDensity:  # estimate 2D histogram --> create pcolormesh
        nbins = 100
        H, xedges, yedges = np.histogram2d(rFP, rTP, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
        dataDensity = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.cool)
        cbar = plt.colorbar(dataDensity, orientation='horizontal')
        cbar.ax.set_ylabel('hit rate density')
    else:
        sc = ax1.scatter(rFP, rTP, c=colorSC, cmap=cmapSC, norm=normSC, marker=pU.markers[0], label='runout points')

        if displayColorBar:
            pU.addColorBar(sc, ax1, ticksSC, unitSC, title=paraVar, pad=0.08, tickLabelsList=itemsList)

        ax1.plot(rFPRef, rTPRef, color='g', label='Reference', marker='+', markersize=2*pU.ms, linestyle='None')
        ax1.legend(loc=4)

    plt.xlim([-0.03, max(1, max(rFP)+0.03)])
    plt.ylim([-0.03, 1.03])
    plt.grid('on')

    outFileName = '_'.join([pathDict['projectName'], resType,
                            'thresholdValue', str(thresholdValue).replace('.', 'p'), 'ROC'])

    pU.saveAndOrPlot(pathDict, outFileName, fig)

    return
