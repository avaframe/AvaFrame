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

# create local logger
log = logging.getLogger(__name__)


def visuTransfo(rasterTransfo, inputData, cfgSetup, pathDict, cfgFlags):
    """
    Plot and save the domain transformation figure
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
    # read avaPath with scale
    xPath = rasterTransfo['x']
    yPath = rasterTransfo['y']
    # read domain boundarries with scale
    xllc = rasterTransfo['xllc']
    yllc = rasterTransfo['yllc']
    cellSize = rasterTransfo['cellSize']
    DBXl = rasterTransfo['DBXl']*cellSize+xllc
    DBXr = rasterTransfo['DBXr']*cellSize+xllc
    DBYl = rasterTransfo['DBYl']*cellSize+yllc
    DBYr = rasterTransfo['DBYr']*cellSize+yllc

    ############################################
    # prepare for plot
    n, m = np.shape(xyRaster)
    x = np.arange(m)*cellSize+xllc
    y = np.arange(n)*cellSize+yllc
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    xx = rasterTransfo['x'][indStartOfRunout]
    yy = rasterTransfo['y'][indStartOfRunout]

    l = rasterTransfo['l']
    s = rasterTransfo['s']

    maskedArray = np.ma.masked_where(xyRaster == 0, xyRaster)
    maskedArraySL = np.ma.masked_where(slRaster == 0, slRaster)

    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin(maskedArray), np.nanmax(maskedArray), continuous=pU.contCmap)
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
    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)


def visuRunoutComp(rasterTransfo, resAnalysis, newRasters, cfgSetup, pathDict, cfgFlags):
    """
    Plot and save the Peak Fields distribution (max mean per cross section)
    after coordinate transformation
    """
    ####################################
    # Get input data
    resType = cfgSetup['resType']
    thresholdValue = cfgSetup['thresholdValue']
    # read paths
    projectName = pathDict['projectName']
    # read data
    s = rasterTransfo['s']
    l = rasterTransfo['l']
    PPRCrossMax = resAnalysis['PPRCrossMax']
    PPRCrossMean = resAnalysis['PPRCrossMean']
    PFDCrossMax = resAnalysis['PFDCrossMax']
    PFDCrossMean = resAnalysis['PFDCrossMean']
    PFVCrossMax = resAnalysis['PFVCrossMax']
    PFVCrossMean = resAnalysis['PFVCrossMean']

    ############################################
    # prepare for plot
    title = ['Pressure ', 'Flow Depth ', 'Flow Velocity ']
    unit = ['$P(s)$ [kPa]', '$fd(s)$ [m]', '$v(s) [m.s^{-1}]$']
    dataMax = np.array(([None] * 3))
    dataMax[0] = PPRCrossMax
    dataMax[1] = PFDCrossMax
    dataMax[2] = PFVCrossMax

    dataMean = np.array(([None] * 3))
    dataMean[0] = PPRCrossMean
    dataMean[1] = PFDCrossMean
    dataMean[2] = PFVCrossMean

    ############################################
    # Figure: Pressure depth speed

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(pU.figW*3, pU.figH))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.3)

    for ax, maxVal, meanVal, titleVal, unitVal in zip(axes.flatten(), dataMax, dataMean, title, unit):
        ax.plot(maxVal[0, :], s, '--k', label='Max Reference')
        ax.plot(meanVal[0, :], s, '-k', label='Mean Reference')
        ax.plot(maxVal[1, :], s, '--b', label='Max Simulation')
        ax.plot(meanVal[1, :], s, '-b', label='Mean Simulation')

        ax.set_title(titleVal + 'distribution along path')
        ax.legend(loc='best' )
        ax.set_ylabel('s [m]')
        ax.set_ylim([s.min(), s.max()])
        ax.set_xlim(auto=True)
        ax.set_xlabel(unitVal)
    pU.putAvaNameOnPlot(ax, projectName)

    outFileName = '_'.join([projectName, resType,
                            'thresholdValue', str(thresholdValue).replace('.', 'p'), 'slComparison'])
    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)

    outFilePath = os.path.join(pathDict['pathResult'], 'pics', outFileName + '.png')

    return outFilePath


def visuRunoutStat(rasterTransfo, resAnalysis, newRasters, cfgSetup, pathDict, cfgFlags):
    """
    Plot and save the Peak field  distribution after coord transfo
    used when more then 2 simulations are compared
    """
    ####################################
    # Get input data
    percentile = cfgSetup.getfloat('percentile')
    resType = cfgSetup['resType']
    thresholdValue = cfgSetup['thresholdValue']
    unit = pU.cfgPlotUtils['unit' + resType]
    name = pU.cfgPlotUtils['name' + resType]
    # read paths
    projectName = pathDict['projectName']
    nRef = pathDict['referenceFile']
    # read data
    s = rasterTransfo['s']
    l = rasterTransfo['l']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    dataPressure = newRasters['newRaster' + resType.upper()]
    rasterdataPres = dataPressure[nRef]
    runout = resAnalysis['runout'][0]
    PPRCrossMax = resAnalysis[resType.upper() + 'CrossMax']

    ############################################
    # prepare for plot
    pMean = np.mean(PPRCrossMax, axis=0)
    pMedian = np.median(PPRCrossMax, axis=0)
    pPercentile = np.percentile(PPRCrossMax, [percentile/2, 50, 100-percentile/2], axis=0)

    maskedArray = np.ma.masked_where(rasterdataPres == 0, rasterdataPres)

    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin(maskedArray), np.nanmax(maskedArray), continuous=pU.contCmap)
    cmap.set_bad('w', 1.)

    ############################################
    # Figure: Analysis runout
    fig = plt.figure(figsize=(pU.figW*2, pU.figH))
    ax1 = plt.subplot(121)

    ax1.axhline(y=np.max(runout), color='k', linestyle='-.', label='runout max')
    ax1.axhline(y=np.average(runout), color='k', linestyle='-', label='runout mean')
    ax1.axhline(y=np.min(runout), color='k', linestyle=':', label='runout min')

    ax1.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                label='start of run-out area point : %.1f °' % resAnalysis['startOfRunoutAreaAngle'])
    ref5, im = pU.NonUnifIm(ax1, l, s, maskedArray, 'l [m]', 's [m]',
                         extent=[l.min(), l.max(), s.min(), s.max()],
                         cmap=cmap, norm=norm)

    ax1.set_title('Peak Pressure 2D plot for the reference')
    ax1.legend(loc=4)
    pU.putAvaNameOnPlot(ax1, projectName)

    pU.addColorBar(im, ax1, ticks, unit)

    ax2 = plt.subplot(122)

    ax2.fill_betweenx(s, pPercentile[2], pPercentile[0],
                      facecolor=[.8, .8, .8], alpha=0.5, label=('[%.2f, %.2f]%% interval' % (percentile/2, 100-percentile/2)))
    matplotlib.patches.Patch(alpha=0.5, color=[.8, .8, .8])
    ax2.plot(pMedian, s, color='r', label='median')
    ax2.plot(pMean, s, color='b', label='mean')

    ax2.set_title('%s distribution along the path between runs' % name)
    ax2.legend(loc=4)
    ax2.set_ylabel('s [m]')
    ax2.set_ylim([s.min(), s.max()])
    ax2.set_xlim(auto=True)
    ax2.set_xlabel('$P_{max}(s)$ [%s]' % unit)

    outFileName = '_'.join([projectName, resType,
                            'thresholdValue', str(thresholdValue).replace('.', 'p'), 'slComparisonStat'])

    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)


def visuMass(resAnalysis, pathDict, cfgFlags):
    """
    Plot and save the results from mass analysis
    """
    ####################################
    # Get input data
    # read paths
    projectName = pathDict['projectName']
    # read data
    entMassFlowArray = resAnalysis['entMassFlowArray']
    totalMassArray = resAnalysis['totalMassArray']
    entMass = resAnalysis['entMass']
    finalMass = resAnalysis['finalMass']
    time = resAnalysis['time']

    ############################################
    # prepare for plot
    Title = ['Entrained Mass Flow', 'Total Mass']
    Unit = ['Entrained Mass Flow [$kg.s{^-1}$]', 'Total Mass [kg]']
    DataMass = np.array(([None] * 2))
    DataMass[0] = entMassFlowArray
    DataMass[1] = totalMassArray
    ############################################
    # Figure: Pressure depth speed

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(pU.figW*2, pU.figH))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.3)

    for ax, dataMass, title, unit in zip(axes.flatten(), DataMass, Title, Unit):
        ax.plot(time, dataMass[0, :], '-k', label='Reference')
        ax.plot(time, dataMass[1, :], '-b', label='Simulation')

        ax.set_title(title + ' function of time')
        ax.legend(loc=4)
        ax.set_xlabel('t [s]')
        ax.set_ylabel(unit)

    ax2 = axes.flatten()[1].twinx()
    # ax2.set_ylabel('z [m]')
    ax2.spines['right'].set_color('r')
    ax2.tick_params(axis='y', colors='r')
    ax2.plot(time, (dataMass[1, :]-dataMass[0, :])/dataMass[0, :]*100, 'r', label='total mass')

    if np.any(entMass):
        axes.flatten()[1].text(time[-1]/4, (np.nanmin(dataMass[0, :])+np.nanmax(dataMass[0, :]))/2, 'Entrained Mass Difference : %.2f kg \n Relative to entrained mass : %.2f %% \n Relative to total mass : %.2f %% ' %
             ((entMass[0]-entMass[1]), (entMass[0]-entMass[1])/entMass[0]*100, (entMass[0]-entMass[1])/finalMass[0]*100),
             bbox=dict(boxstyle="square", ec='white', fc='white'),
             horizontalalignment='left', verticalalignment='bottom')

    ax2.set_ylabel('Entrained Mass Difference relative to total mass[%]', color='r')

    outFileName = '_'.join([projectName, 'massAnalysis'])
    pU.putAvaNameOnPlot(ax2, pathDict['projectName'])
    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)

    outFilePath = os.path.join(pathDict['pathResult'], 'pics', outFileName + '.png')

    return outFilePath


def visuSimple(rasterTransfo, resAnalysis, newRasters, pathDict, cfgFlags):
    """
    Plot and save the Peak Pressure Peak Flow depth and Peak speed
    fields after coord transfo
    """
    ####################################
    # Get input data
    # read paths
    projectName = pathDict['projectName']
    nRef = pathDict['referenceFile']
    # read data
    plim = resAnalysis['thresholdValue']
    s = rasterTransfo['s']
    l = rasterTransfo['l']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    dataPressure = newRasters['newRasterPPR']
    rasterdataPres = dataPressure[nRef]
    dataDepth = newRasters['newRasterPFD']
    rasterdataDepth = dataDepth[nRef]
    dataSpeed = newRasters['newRasterPFV']
    rasterdataSpeed = dataSpeed[nRef]
    runout = resAnalysis['runout'][0]

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
        cmap, _, ticks, norm = pU.makeColorMap(cmap, np.nanmin(maskedArray), np.nanmax(maskedArray), continuous=pU.contCmap)
        cmap.set_bad('w', 1.)
        ax.axhline(y=runout[0], color='k', linestyle='-', label='runout')

        ax.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                   label='Start or run-out point : %.1f °' % resAnalysis['startOfRunoutAreaAngle'])
        ref3, im = pU.NonUnifIm(ax, l, s, maskedArray, 'l [m]', 's [m]',
                             extent=[l.min(), l.max(), s.min(), s.max()],
                             cmap=cmap, norm=norm)

        ax.set_title(title)
        ax.legend(loc=4)
        pU.addColorBar(im, ax, ticks, unit)
        pU.putAvaNameOnPlot(ax, pathDict['projectName'])

    outFileName = '_'.join([projectName, 'plim',
                            str((plim)).replace('.', 'p'), 'referenceFields'])

    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)


def visuComparison(rasterTransfo, inputs, pathDict, cfgFlags):
    """
    Plot and save the comparison between current simulation and Reference
    in the run-out area
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
    newRasterMask = inputs['newRasterMask']
    i = inputs['i']
    resType = inputs['resType']
    unit = pU.cfgPlotUtils['unit' + resType]
    name = pU.cfgPlotUtils['name' + resType]
    thresholdArray = inputs['thresholdArray']
    thresholdValue = thresholdArray[-1]
    contCmap = pathDict['contCmap']

    ############################################
    # Figure: Raster comparison (mask for the pThreshold given in the ini file)
    fig = plt.figure(figsize=(pU.figW*2, pU.figH))
    ax1 = plt.subplot2grid((1,2), (0,0))

    # get color map
    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin((refData)), np.nanmax((refData)), continuous=contCmap)
    cmap.set_bad(color='w')
    refDataPlot = np.ma.masked_where(refData == 0.0, refData)
    ref0, im = pU.NonUnifIm(ax1, l, s, refDataPlot, 'l [m]', 's [m]',
                         extent=[l.min(), l.max(), s.min(), s.max()],
                         cmap=cmap, norm=norm)
    ax1.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                label='start of run-out area point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])
    ax1.set_title('Reference %s in the RunOut area' % name +
                  '\n' + '%s threshold: %.1f %s' % (name, thresholdValue, unit))
    pU.addColorBar(im, ax1, ticks, unit)

    yLim = s[max(np.max(np.nonzero(np.any(refData > 0, axis=1))[0]), np.max(np.nonzero(np.any(compData > 0, axis=1))[0]))]
    ax1.set_ylim([0, yLim])
    ax1.legend(loc='lower right')
    pU.putAvaNameOnPlot(ax1, projectName)

    ax2 = plt.subplot2grid((1,2), (0,1))
    colorsList = [[0, 0, 1], [1, 1, 1], [1, 0, 0]]
    cmap = matplotlib.colors.ListedColormap(colorsList)
    cmap.set_under(color='b')
    cmap.set_over(color='r')
    cmap.set_bad(alpha=0)
    data = newRasterMask-refRasterMask
    data = np.ma.masked_where(data == 0.0, data)
    if newRasterMask[indStartOfRunout:, :].any() or refRasterMask[indStartOfRunout:, :].any():
        ref1, im1 = pU.NonUnifIm(ax2, l, s, data, 'l [m]', 's [m]',
                         extent=[l.min(), l.max(), s[indStartOfRunout], yLim], cmap=cmap)
        im1.set_clim(vmin=-0.5, vmax=0.5)
        ax2.set_ylim([s[indStartOfRunout], yLim])
    else:
        ax2.text(.5,.5, 'No data in the run out area!', fontsize=18, color='red',
        bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')
    if pathDict['compType'][0] == 'comModules':
        namePrint = 'refMod:' + pathDict['compType'][1] +'_' + 'compMod:' +pathDict['compType'][2]
        pU.putAvaNameOnPlot(ax2, namePrint)

    ax2.set_title('Difference %s current - reference in runout area' % resType + '\n' + 'Blue = FN, Red = FP')

    outFileName = '_'.join([projectName, 'thresholdValue', str(thresholdValue).replace('.', 'p'),  'sim', str(i), 'AreaComparisonToReference'])
    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)

    ############################################
    # Figure: Raster comparison
    fig = plt.figure(figsize=(pU.figW*3, pU.figH*2))#, constrained_layout=True)
    ax1 = plt.subplot2grid((3,3), (0,0), rowspan=3)
    # get color map
    cmap, _, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.nanmin((refData)), np.nanmax((refData)), continuous=pU.contCmap)
    cmap.set_bad(color='w')
    refDataPlot = np.ma.masked_where(refData == 0.0, refData)
    ref0, im = pU.NonUnifIm(ax1, l, s, refDataPlot, 'l [m]', 's [m]',
                         extent=[l.min(), l.max(), 0, yLim],
                         cmap=cmap, norm=norm)
    ax1.axhline(y=s[indStartOfRunout], color='k', linestyle='--',
                label='start of run-out area point : %.1f °' % rasterTransfo['startOfRunoutAreaAngle'])
    ax1.set_title('Reference %s' % name)
    pU.addColorBar(im, ax1, ticks, unit)
    ax1.set_ylim([0, yLim])
    ax1.legend(loc='lower right')
    pU.putAvaNameOnPlot(ax1, projectName)

    compData = compData[indStartOfRunout:, :]
    refData = refData[indStartOfRunout:, :]
    dataDiff = compData - refData
    dataDiff = np.where((refData==0) & (compData==0), np.nan, dataDiff)
    dataDiffPlot = np.where((refData<thresholdArray[-1]) & (compData<thresholdArray[-1]), np.nan, dataDiff)
    dataDiffPlot = dataDiffPlot[np.isnan(dataDiffPlot) == False]

    if dataDiffPlot.size:
        # only add the second axis if one of the two avalanches reached the run out area
        indDiff = np.abs(dataDiffPlot) > 0

        if indDiff.any():
            # only plot hist and CDF if there is a difference in the data
            ax2 = plt.subplot2grid((3,3), (0,1), rowspan=2, colspan=2)
        else:
            ax2 = plt.subplot2grid((3,3), (0,1), rowspan=3, colspan=3)

        cmap = pU.cmapdiv
        cmap.set_bad(color='w')
        elev_max = inputs['diffLim']
        ref0, im3 = pU.NonUnifIm(ax2, l, s[indStartOfRunout:], (dataDiff), 'l [m]', 's [m]',
                             extent=[l.min(), l.max(), s[indStartOfRunout], yLim],
                             cmap=cmap)
        im3.set_clim(vmin=-elev_max, vmax=elev_max)

        # print contour lines only if the thre threshold is reached
        L, S = np.meshgrid(l, s[indStartOfRunout:])
        colorsP = pU.colorMaps['pfd']['colors'][1:]
        if (np.where(refData > thresholdArray[-1], True, False)).any():
            contourRef = ax2.contour(L, S, refData, levels=thresholdArray[:-1], linewidths=2, colors=colorsP)
            labels = [str(level) + unit for level in thresholdArray[:-1]]
            for j in range(len(contourRef.collections)):
                contourRef.collections[j].set_label(labels[j])
        else:
            log.warning('Reference did not reach the run out area!')
            ax2.text(0, (s[indStartOfRunout] + yLim)/2, 'Reference did not reach the run out area!', fontsize=24, color='red',
            bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')
        if (np.where(compData > thresholdArray[-1], True, False)).any():
            contourComp = ax2.contour(L, S, compData, levels=thresholdArray[:-1], linewidths=2, colors=colorsP, linestyles= 'dashed')
        else:
            log.warning('Simulation did not reach the run out area!')
            ax2.text(0, (s[indStartOfRunout] + yLim)/2, 'Simulation did not reach the run out area!', fontsize=24, color='red',
            bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')

        if pathDict['compType'][0] == 'comModules':
            namePrint = 'refMod:' + pathDict['compType'][1] +'_' + 'compMod:' +pathDict['compType'][2]
            pU.putAvaNameOnPlot(ax2, namePrint)

        if indDiff.any():
            # only plot hist and CDF if there is a difference in the data
            ax3 = plt.subplot2grid((3,3), (2, 1))
            ax4 = plt.subplot2grid((3,3), (2, 2))
            # there is data to compare in the run out area
            centiles = sPlot.plotHistCDFDiff(dataDiffPlot, ax4, ax3, insert='False',
                                                 title=['%s diff histogram' % name, '%s diff CDF (95%% and 99%% centiles)' % name])

        ax2.set_ylim([s[indStartOfRunout], yLim])
        ax2.legend(loc='lower right')
        pU.addColorBar(im3, ax2, None, unit, title=name, extend='both')
    else:
        # if no avalanche reached the run out area print a warning on the second plot
        ax2 = plt.subplot2grid((3,3), (0,1), rowspan=3, colspan=3)
        log.warning('No data in the run out area!')
        ax2.text(.5,.5, 'No data in the run out area!', fontsize=24, color='red',
        bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=1'), ha='center', va='center')

    ax2.set_title(
        '%s difference and contour lines' % name + '\n' + 'refMod = full, compMod = dashed line')

    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    outFileName = '_'.join([projectName, 'plim', str(thresholdValue).replace('.', 'p'),  'sim', str(i), 'ContourComparisonToReference'])
    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)

    outFilePath = os.path.join(pathDict['pathResult'], 'pics', outFileName + '.png')

    return outFilePath


def resultWrite(pathDict, cfgSetup, flagMass, rasterTransfo, resAnalysis):
    """
    This function writes the main Aimec results to a file (outputFile)
    in pathDict
    """
    ####################################
    # Get input data
    projectName = pathDict['projectName']
    pathResult = pathDict['pathResult']
    pathName = pathDict['pathName']
    nRef = pathDict['referenceFile']
    demName = os.path.basename(pathDict['demSource'])
    dataName = [os.path.basename(name) for name in pathDict['ppr']]
    domainWidth = cfgSetup['domainWidth']
    thresholdValue = cfgSetup['thresholdValue']
    resType = cfgSetup['resType']
    unit = pU.cfgPlotUtils['unit' + resType]
    name = pU.cfgPlotUtils['name' + resType]

    startOfRunoutAreaAngle = resAnalysis['startOfRunoutAreaAngle']
    indStartOfRunout = rasterTransfo['indStartOfRunout']
    s = rasterTransfo['s']
    sStart = s[indStartOfRunout]
    runoutFromMid = resAnalysis['runout'][0] - sStart
    runout = resAnalysis['runout']
    MMPPR = resAnalysis['MMPPR']
    MMPFD = resAnalysis['MMPFD']
    MMPFV = resAnalysis['MMPFV']
    deltaH = resAnalysis['deltaH']
    elevRel = resAnalysis['elevRel']
    if flagMass:
        relMass = resAnalysis['relMass']
        entMass = resAnalysis['entMass']
        relativMassDiff = resAnalysis['relativMassDiff']
        finalMass = resAnalysis['finalMass']
        GI = resAnalysis['growthIndex']
        GR = resAnalysis['growthGrad']
    TP = resAnalysis['TP']
    FN = resAnalysis['FN']
    FP = resAnalysis['FP']
    TN = resAnalysis['TN']
    areaSum = TP[nRef] + FN[nRef]
    if areaSum == 0:
        log.warning('Reference did not reach the run-out area. Not normalizing area indicators')
        areaSum = np.ones(np.shape(TP))

    ############################################
    # prepare for writing
    legend = ['fileNr', 'Xrunout', 'Yrunout', 'Lrunout', 'runoutFromSROA', 'elevRel', 'deltaH',
              'MMPPR', 'MMPFD', 'MMPFV', 'TP ', 'FN ', 'FP ', 'TN']
    resfile = [runout[1], runout[2], runout[0], runoutFromMid, elevRel, deltaH, MMPPR,
               MMPFD, MMPFV, TP/areaSum, FN/areaSum, FP/areaSum, TN/areaSum]

    if flagMass:
        legend = legend + ['relMass', 'entMass', 'finalMass', 'rMassDif', 'GI', 'GR']
        resfile = resfile + [relMass, entMass, finalMass, relativMassDiff, GI, GR]

    header = ''.join(['projectName: ', projectName, '\n',
                      'path: ', pathName, '\n',
                      'dhm: ', demName, '\n',
                      'domain_width: ', str(domainWidth), ' m\n',
                      'run out based on ', name, ' results\n',
                      name, ' limit: ', str(thresholdValue), unit, '\n',
                      'start of runout area Angle (SROA angle): ', str(round(startOfRunoutAreaAngle, 2)), ' °\n'])

    outFileName = '_'.join(['Results', projectName, '', '', 'plim',
                            str(thresholdValue), 'w', str(domainWidth)]) + '.txt'
    outname = os.path.join(pathResult, outFileName)

    # check if folder exists / create
    if not os.path.exists(os.path.dirname(outname)):
        os.makedirs(os.path.dirname(outname))

    ############################################
    # write results to file
    output = resfile
    log.info('write output file: %s' % outname)
    fid = open(outname, 'w')
    fid.write(header)
    # write table legend
    for j in range(len(legend)):
        fid.write('{:<15s}'.format(legend[j]))
    fid.write('\n')
    # write table values
    for i in range(np.shape(output)[1]):
        fid.write('{:<15d}'.format(i))
        for j in range(np.shape(output)[0]):
            try:
                fid.write('{:<15.3f}'.format(output[j][i]))
            except:
                fid.write('{:<15}'.format('NaN'))
        fid.write('\n')
    # compute and write min, max, mean, std to file
    fid.write('Range of the output parameters\n')
    for j in range(len(legend)):
        fid.write('{:<15s}'.format(legend[j]))
    fid.write('\n')

    names = ['Minimum', 'Maximum', 'Mean', 'STD']
    operators = ['nanmin', 'nanmax', 'nanmean', 'nanstd']
    for name, operator in zip(names, operators):
        fid.write('{:<15s}'.format(name))
        for j in range(np.shape(output)[0]):
            fid.write('{:<15.3f}'.format(getattr(np, operator)(output[j][:], axis=0)))
        fid.write('\n')

    for i in range(np.shape(output)[1]):
        tmp = os.path.basename(dataName[i])
        name = os.path.splitext(tmp)[0]
        fid.write('file number: %d = %s \n' %(i, name))
    fid.close()

    log.info('File written: %s' % outname)


def resultVisu(cfgSetup, pathDict, cfgFlags, rasterTransfo, resAnalysis):
    """
    Visualize results in a nice way
    """
    ####################################
    # Get input data
    resType = cfgSetup['resType']
    unit = cfgSetup['unit']
    paraVar = cfgSetup['varParList'].split('|')[0]
    name = pU.cfgPlotUtils['name' + resType]
    fnames = pathDict[resType]
    maxMaxDPPR = resAnalysis['MM' + resType.upper()]

    thresholdValue = cfgSetup['thresholdValue']

    nRef = pathDict['referenceFile']

    flag = float(cfgFlags['typeFlag'])

    zPath = rasterTransfo['z']
    sPath = rasterTransfo['s']

    runout = resAnalysis['runout'][0]
    areaSum = resAnalysis['TP'][nRef] + resAnalysis['FN'][nRef]
    if areaSum == 0:
        log.warning('Reference did not reach the run-out area. Not normalizing area indicators')
        areaSum = 1
    rTP = resAnalysis['TP'] / areaSum
    rFP = resAnalysis['FP'] / areaSum

    # prepare for plot
    if flag == 2:
        title = 'Visualizing EGU growth index data'
        tipo = 'growthInd'
        GI = resAnalysis['growthIndex']
        data = GI
        yaxis_label = 'growth index [GI]'

    elif flag == 3:
        title = 'Visualizing max ' + name + ' data'
        tipo = 'relMax' + resType + '_thresholdValue' + str(thresholdValue).replace('.', 'p')
        data = maxMaxDPPR / maxMaxDPPR[nRef]
        yaxis_label = 'relative max ' + name + ' [-]'

    else:
        log.error('Wrong flag')
        return None

    log.info(title)

    # If more than 100 files are provided, add a density plot
    plotDensity = 0
    if (len(fnames) > 100):
        plotDensity = 1

    if 'colorParameter' in pathDict:
        if pathDict['colorParameter'] == []:
            nSamples = np.size(runout)
            colors = np.zeros(nSamples)
            cmap, _, ticks, norm = pU.makeColorMap(pU.cmapVar, None, None, continuous=True)
            displayColorBar = False
            dataFrame = False
        else:
            typeCP = type(pathDict['colorParameter'][0])
            if typeCP == str:
                keys = list(set(pathDict['colorParameter']))
                nKeys = len(keys)
                cmap = pU.cmapAimec(np.linspace(0, 1, nKeys, dtype=float))
                df = pd.DataFrame(dict(runout=runout, data=data, rFP=rFP, rTP=rTP, colorParameter = pathDict['colorParameter']))
                dataFrame = True
                displayColorBar = False
            else:
                colors = pathDict['colorParameter']
                cmap, _, ticks, norm = pU.makeColorMap(pU.cmapVar, np.nanmin(pathDict['colorParameter']), np.nanmax(pathDict['colorParameter']), continuous=True)
                displayColorBar = True
                dataFrame = False
    else:
        displayColorBar = False
        dataFrame = False
        nSamples = np.size(runout)
        colors = np.zeros(nSamples)
        cmap, _, ticks, norm = pU.makeColorMap(pU.cmapVar, None, None, continuous=True)
    #######################################
    # Final result diagram - z_profile+data

    fig = plt.figure(figsize=(pU.figW*2, pU.figH))

    # show flow path
    ax1 = fig.add_subplot(111)
    ax1.set_title(title)
    ax1.set_ylabel(yaxis_label, color=pU.cmapVar['colors'][-1])
    ax1.spines['left'].set_color(pU.cmapVar['colors'][-1])
    ax1.tick_params(axis='y', colors=pU.cmapVar['colors'][-1])
    ax1.set_xlabel(''.join(['s [m] - runout with ', str(thresholdValue),
                            ' kPa threshold']))
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
        if dataFrame:
            sns.scatterplot('runout', 'data', marker=pU.markers, data=df, hue='colorParameter', palette=cmap, ax=ax1)
        else:
            sc = ax1.scatter(runout, data, marker=pU.markers, c=colors, cmap=cmap)
            if displayColorBar:
                pU.addColorBar(sc, ax2, ticks, unit, title=paraVar, pad=0.08)
        ax1.plot(runout[nRef], data[nRef], color='g', label='Reference', marker='+',
                         markersize=2*pU.ms, linestyle='None')
        ax1.legend(loc=4)

    ax1.grid('on')

    outFileName = '_'.join([pathDict['projectName'], tipo])

    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)

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
        if dataFrame:
            sns.scatterplot('rFP', 'rTP', marker=pU.markers, data=df, hue='colorParameter', palette=cmap, ax=ax1)
        else:
            sc = ax1.scatter(rFP, rTP, marker=pU.markers, c=colors, cmap=cmap)
            if displayColorBar:
                pU.addColorBar(sc, ax1, ticks, unit, title=paraVar)
        ax1.plot(rFP[nRef], rTP[nRef], color='g', label='Reference', marker='+',
                         markersize=2*pU.ms, linestyle='None')
        ax1.legend(loc=4)

    plt.xlim([-0.03, max(1, max(rFP)+0.03)])
    plt.ylim([-0.03, 1.03])
    plt.grid('on')

    outFileName = '_'.join([pathDict['projectName'], resType,
                            'thresholdValue', str(thresholdValue).replace('.', 'p'), 'ROC'])

    pU.saveAndOrPlot(pathDict, cfgFlags, outFileName, fig)

    return
