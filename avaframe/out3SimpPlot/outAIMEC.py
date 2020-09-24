"""
    Plotting and saving AIMEC results

    This file is part of Avaframe.
"""

import sys
import os
import logging
import math
import numpy as np
import scipy as sp
import copy
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.image import NonUniformImage


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in3Utils.ascUtils as IOf
from avaframe.out3SimpPlot.plotSettings import *

# create local logger
log = logging.getLogger(__name__)


def visuTransfo(rasterTransfo, inputData, cfgPath, cfgFlags):
    """
    Plot and save the domain transformation figure
    """
    # read paths
    pathResult = cfgPath['pathResult']
    projectName = cfgPath['dirName']
    # read rasterdata
    sourceData = inputData['sourceData']
    header = sourceData['header']
    xllc = rasterTransfo['xllc']
    yllc = rasterTransfo['yllc']
    cellsize = rasterTransfo['cellsize']
    rasterdata = sourceData['rasterData']
    # read avaPath with scale
    Avapath = inputData['Avapath']
    xPath = Avapath['x']
    yPath = Avapath['y']
    # read domain boundarries with scale
    DBXl = rasterTransfo['DBXl']*cellsize+xllc
    DBXr = rasterTransfo['DBXr']*cellsize+xllc
    DBYl = rasterTransfo['DBYl']*cellsize+yllc
    DBYr = rasterTransfo['DBYr']*cellsize+yllc

    fig = plt.figure(figsize=(figW*2, figH))

#    for figure: referenz-simulation bei pLim=1
    ax1 = plt.subplot(121)
    indRunoutPoint = rasterTransfo['indRunoutPoint']
    xx = rasterTransfo['x'][indRunoutPoint]
    yy = rasterTransfo['y'][indRunoutPoint]
    newRasterdata = rasterdata
    maskedArray = newRasterdata
    # get color map
    ticks = ticksPres
    cmap, _, _, norm = makeColorMap(colorsPres, levPres, 0.0, np.nanmax(maskedArray))
    cmap.set_under(color='w')

    n, m = np.shape(newRasterdata)
    x = np.arange(m)*cellsize+xllc
    y = np.arange(n)*cellsize+yllc
    im0 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap, norm=norm)
    # im.set_interpolation('bilinear')
    im0.set_clim(vmin=0.000000001)
    im0.set_data(x, y, maskedArray)
    ref1 = ax1.images.append(im0)
    # cbar = ax1.figure.colorbar(im0, ax=ax1, use_gridspec=True)
    plt.autoscale(False)
    ref0 = plt.plot(xx, yy, 'ro', label='Beta point : %.1f °' %
                    rasterTransfo['runoutAngle'])
    ref2 = plt.plot(xPath, yPath, 'b-', label='flow path')
    ref3 = plt.plot(DBXl, DBYl, 'g-', label='domain')
    ref3 = plt.plot(DBXr, DBYr, 'g-', label='domain')
    ref3 = plt.plot([DBXl, DBXr], [DBYl, DBYr], 'g-', label='domain')
    refs = [ref0[0], ref2[0], ref3[0]]

    labels = ['Beta point : %.1f °' % rasterTransfo['runoutAngle'], 'flow path', 'domain']
    ax1.set_title('XY Domain')
    ax1.legend(refs, labels, loc=0)
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')

    ax2 = plt.subplot(122)
    ax2.set_title('sl Domain' + '\n' + 'Black = out of raster')
    isosurf = copy.deepcopy(inputData['avalData'])
    lcoord = rasterTransfo['l']
    scoord = rasterTransfo['s']
    ref1 = ax2.axhline(y=scoord[indRunoutPoint], color='r',
                       label='Beta point : %.1f °' % rasterTransfo['runoutAngle'])
    maskedArray = isosurf  # np.ma.array(isosurf,mask=np.isnan(isosurf))
    im = NonUniformImage(ax2, extent=[lcoord.min(), lcoord.max(),
                                      scoord.min(), scoord.max()], cmap=cmap, norm=norm)
    # im.set_interpolation('bilinear')
    im.set_clim(vmin=0.000000001)
    im.set_data(lcoord, scoord, maskedArray)
    ref0 = ax2.images.append(im)
    cbar = ax2.figure.colorbar(im, ax=ax2, ticks=ticks)
    cbar.ax.set_ylabel('peak pressure [kPa]')
    ax2.set_xlim([lcoord.min(), lcoord.max()])
    ax2.set_ylim([scoord.min(), scoord.max()])
    ax2.set_xlabel('l [m]')
    ax2.set_ylabel('s [m]')
    ax2.legend()

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()
    if cfgFlags.getboolean('savePlot'):
        outFileName = projectName + '_domTransfo'
        outname = os.path.join(pathResult, 'pics', outFileName)
        if not os.path.exists(os.path.dirname(outname)):
            os.makedirs(os.path.dirname(outname))
        fig.savefig(outname)

    plt.close(fig)


def visuRunout(rasterTransfo, resAnalysis, pLim, newRasters, cfgPath, cfgFlags):
    """
    Plot and save the Peak Pressure  distribution after coord transfo
    """
    # read paths
    pathResult = cfgPath['pathResult']
    projectName = cfgPath['dirName']
    # read data
    scoord = rasterTransfo['s']
    lcoord = rasterTransfo['l']
    indRunoutPoint = rasterTransfo['indRunoutPoint']
    sBeta = scoord[indRunoutPoint]
    rasterArea = rasterTransfo['rasterArea']
    dataPressure = newRasters['newRasterPressure']
    rasterdataPres = dataPressure[0]  # ana3AIMEC.makeRasterAverage(dataPressure)
    runout = resAnalysis['runout'][0] + sBeta
    runoutMean = resAnalysis['runoutMean'][0] + sBeta
    pCrossAll = resAnalysis['pCrossAll']

    # prepare for plot
    pMean = np.mean(pCrossAll, axis=0)
    pMedian = np.median(pCrossAll, axis=0)
    pPercentile = sp.percentile(pCrossAll, [2.5, 50, 97.5], axis=0)

    fig = plt.figure(figsize=(figW*2, figH))
    ax1 = plt.subplot(121)
    ax1.set_title('Peak Pressure 2D plot for the reference')
    ref1 = ax1.axhline(y=scoord[indRunoutPoint], color='k',
                       label='Beta point : %.1f °' % resAnalysis['runoutAngle'])
    ref1 = ax1.axhline(y=np.max(runout), color='r', label='runout max')
    ref2 = ax1.axhline(y=np.average(runout), color='y', label='runout mean')
    ref3 = ax1.axhline(y=np.min(runout), color='g', label='runout min')
    # ref3 = ax1.plot(np.zeros(np.shape(scoord)), scoord,'.r', linewidth=0.1)
    isosurf = copy.deepcopy(rasterdataPres)
    xx, yy = np.meshgrid(lcoord, scoord)
    maskedArray = np.ma.masked_where(isosurf == 0, isosurf)
    # get color map
    ticks = ticksPres
    cmap, _, _, norm = makeColorMap(colorsPres, levPres, 0.0, np.nanmax(maskedArray))
    cmap.set_bad('w', 1.)
    im = NonUniformImage(ax1, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=cmap, norm=norm)
    # im.set_interpolation('bilinear')
    im.set_clim(vmin=maskedArray.min(), vmax=maskedArray.max())
    im.set_data(lcoord, scoord, maskedArray)
    ref0 = ax1.images.append(im)
    cbar = ax1.figure.colorbar(im, ax=ax1, ticks=ticks)
    cbar.ax.set_ylabel('peak pressure [kPa]')
    plt.autoscale(False)
    ax1.set_xlim([xx.min(), xx.max()])
    ax1.set_ylim([yy.min(), yy.max()])
    ax1.set_xlabel('l [m]')
    ax1.set_ylabel('s [m]')
    ax1.legend(loc=0)

    ax2 = plt.subplot(122)
    ax2.set_title('Peak Pressure distribution along the path between runs')
    ax2.fill_betweenx(scoord, pPercentile[2], pPercentile[0],
                      facecolor=[.8, .8, .8], alpha=0.5, label='quantiles')
    ref1 = mpatches.Patch(alpha=0.5, color=[.8, .8, .8])
    ref2 = ax2.plot(pMedian, scoord, color='r', label='median')
    ref3 = ax2.plot(pMean, scoord, color='b', label='mean')
    # ref3 = mlines.Line2D([], [], color='b', linewidth=2)
    ax2.set_ylabel('l [m]')
    ax2.set_ylim([yy.min(), yy.max()])
    ax2.set_xlim(auto=True)
    ax2.set_xlabel('$P_{max}(s)$ [kPa]')
    ax2.legend(loc=0)

    if cfgFlags.getboolean('savePlot'):
        outFileName = projectName + '_dptr' + str(int(pLim)) + '_slComparison'
        outname = os.path.join(pathResult, 'pics', outFileName)
        if not os.path.exists(os.path.dirname(outname)):
            os.makedirs(os.path.dirname(outname))
        fig.savefig(outname)

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()

    plt.close(fig)


def visuSimple(rasterTransfo, resAnalysis, newRasters, cfgPath, cfgFlags):
    """
    Plot and save the Peak Pressure Peak Flow depth and Peak speed
    fields after coord transfo
    """
    # read paths
    pathResult = cfgPath['pathResult']
    projectName = cfgPath['dirName']
    # read data
    scoord = rasterTransfo['s']
    lcoord = rasterTransfo['l']
    indRunoutPoint = rasterTransfo['indRunoutPoint']
    sBeta = scoord[indRunoutPoint]
    rasterArea = rasterTransfo['rasterArea']
    dataPressure = newRasters['newRasterPressure']
    rasterdataPres = dataPressure[0]  # ana3AIMEC.makeRasterAverage(dataPressure)
    dataDepth = newRasters['newRasterDepth']
    rasterdataDepth = dataDepth[0]
    dataSpeed = newRasters['newRasterSpeed']
    rasterdataSpeed = dataSpeed[0]
    runout = resAnalysis['runout'][0] + sBeta
    runoutMean = resAnalysis['runoutMean'][0] + sBeta

    fig = plt.figure(figsize=(figW*3, figH))
    ax1 = plt.subplot(131)
    ax1.set_title('Peak Pressure')
    ref1 = ax1.axhline(y=scoord[indRunoutPoint], color='k',
                       label='Beta point : %.1f °' % resAnalysis['runoutAngle'])
    ref2 = ax1.axhline(y=runout[0], color='b', label='runout')
    # ref3 = ax1.plot(np.zeros(np.shape(scoord)), scoord,'.r', linewidth=0.1)
    isosurf = copy.deepcopy(rasterdataPres)
    xx, yy = np.meshgrid(lcoord, scoord)
    maskedArray = np.ma.masked_where(isosurf == 0, isosurf)
    # get color map
    ticks = ticksPres
    cmap, _, _, norm = makeColorMap(colorsPres, levPres, 0.0, np.nanmax(maskedArray))
    cmap.set_bad('w', 1.)
    im = NonUniformImage(ax1, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=cmap, norm=norm)
    # im.set_interpolation('bilinear')
    im.set_clim(vmin=maskedArray.min(), vmax=maskedArray.max())
    im.set_data(lcoord, scoord, maskedArray)
    ref0 = ax1.images.append(im)
    cbar = ax1.figure.colorbar(im, ax=ax1, ticks=ticks)
    cbar.ax.set_ylabel('peak pressure [kPa]')
    plt.autoscale(False)
    ax1.set_xlim([xx.min(), xx.max()])
    ax1.set_ylim([yy.min(), yy.max()])
    ax1.set_xlabel('l [m]')
    ax1.set_ylabel('s [m]')
    ax1.legend(loc=0)

    ax2 = plt.subplot(132)
    ax2.set_title('Peak Flow Depth')
    ref1 = ax2.axhline(y=scoord[indRunoutPoint], color='k',
                       label='Beta point : %.1f °' % resAnalysis['runoutAngle'])
    ref2 = ax2.axhline(y=runout[0], color='b', label='runout')
    # ref3 = ax1.plot(np.zeros(np.shape(scoord)), scoord,'.r', linewidth=0.1)
    isosurf = copy.deepcopy(rasterdataDepth)
    xx, yy = np.meshgrid(lcoord, scoord)
    maskedArray = np.ma.masked_where(isosurf == 0, isosurf)
    cmap = cmapDepth
    cmap.set_bad('w', 1.)
    im = NonUniformImage(ax2, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    im.set_data(lcoord, scoord, maskedArray)
    ref0 = ax2.images.append(im)
    cbar = ax2.figure.colorbar(im, ax=ax2, use_gridspec=True)
    cbar.ax.set_ylabel('peak flow depth $[m]$')
    plt.autoscale(False)
    ax2.set_xlim([xx.min(), xx.max()])
    ax2.set_ylim([yy.min(), yy.max()])
    ax2.set_xlabel('l [m]')
    ax2.set_ylabel('s [m]')
    ax2.legend(loc=0)

    ax3 = plt.subplot(133)
    ax3.set_title('Peak Speed')
    ref1 = ax3.axhline(y=scoord[indRunoutPoint], color='k',
                       label='Beta point : %.1f °' % resAnalysis['runoutAngle'])
    ref2 = ax3.axhline(y=runout[0], color='b', label='runout')
    # ref3 = ax1.plot(np.zeros(np.shape(scoord)), scoord,'.r', linewidth=0.1)
    isosurf = copy.deepcopy(rasterdataSpeed)
    xx, yy = np.meshgrid(lcoord, scoord)
    maskedArray = np.ma.masked_where(isosurf == 0, isosurf)
    cmap = cmapSpeed
    cmap.set_bad('w', 1.)
    im = NonUniformImage(ax3, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    im.set_data(lcoord, scoord, maskedArray)
    ref0 = ax3.images.append(im)
    cbar = ax3.figure.colorbar(im, ax=ax3, use_gridspec=True)
    cbar.ax.set_ylabel('peak speed [m/s]')
    plt.autoscale(False)
    ax3.set_xlim([xx.min(), xx.max()])
    ax3.set_ylim([yy.min(), yy.max()])
    ax3.set_xlabel('l [m]')
    ax3.set_ylabel('s [m]')
    ax3.legend(loc=0)

    if cfgFlags.getboolean('savePlot'):
        outFileName = projectName + '_referenceFields'
        outname = os.path.join(pathResult, 'pics', outFileName)
        if not os.path.exists(os.path.dirname(outname)):
            os.makedirs(os.path.dirname(outname))
        fig.savefig(outname)

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()

    plt.close(fig)


def resultWrite(cfgPath, cfgSetup, resAnalysis):
    """
    This function writes the main Aimec results to a file (outputFile)
    in cfgPath
    """

    projectName = cfgPath['projectName']
    pathResult = cfgPath['pathResult']
    pathName = cfgPath['pathName']
    demName = os.path.basename(cfgPath['demSource'])
    dataName = [os.path.basename(name) for name in cfgPath['pressurefileList']]
    domainWidth = cfgSetup['domainWidth']
    pressureLimit = cfgSetup['pressureLimit']

    runoutAngle = resAnalysis['runoutAngle']

    runout = resAnalysis['runout']
    AMPP = resAnalysis['AMPP']
    MMPP = resAnalysis['MMPP']
    AMD = resAnalysis['AMD']
    MMD = resAnalysis['MMD']
    AMS = resAnalysis['AMS']
    MMS = resAnalysis['MMS']
    deltaH = resAnalysis['deltaH']
    elevRel = resAnalysis['elevRel']
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
    areaSum = TP + FN

    legend = ['fileNr', 'Xrunout', 'Yrunout', 'Lrunout', 'elevRel', 'deltaH', 'AMPP',
              'MMPP', 'AMD', 'MMD', 'AMS', 'MMS', 'relMass', 'entMass', 'finalMass', 'rMassDif', 'GI', 'GR', 'TP ', 'FN ', 'FP ', 'TN']
    resfile = [runout[1], runout[2], runout[0],
               elevRel, deltaH, AMPP, MMPP, AMD, MMD, AMS, MMS, relMass, entMass, finalMass, relativMassDiff, GI, GR, TP/areaSum, FN/areaSum, FP/areaSum, TN/areaSum]

    header = ''.join(['projectName: ',  projectName, '\n',
                      'path: ', pathName, '\n',
                      'dhm: ', demName, '\n',
                      'domain_width: ', str(domainWidth), ' m\n',
                      'pressure_limit: ', str(pressureLimit), ' kPa\n',
                      'start of runout area Angle: ', str(round(runoutAngle, 2)), ' °\n'])

    outFileName = 'Results_pl' + str(pressureLimit) + '_w' + str(domainWidth) + '.txt'
    outname = os.path.join(pathResult, outFileName)


#    chekc if folder exists / create
    if not os.path.exists(os.path.dirname(outname)):
        os.makedirs(os.path.dirname(outname))

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
        tmp = os.path.basename(dataName[i])
        name = os.path.splitext(tmp)[0]
        fid.write('{:<15s}'.format(name))
        for j in range(np.shape(output)[0]):
            try:
                fid.write('{:<15.3f}'.format(output[j][i]))
            except:
                fid.write('{:<15}'.format('NaN'))
        fid.write('\n')
    fid.close()

    log.info('File written: %s' % outname)


def resultVisu(cfgPath, rasterTransfo, resAnalysis, plim):
    """
    Visualize results in a nice way
    Jan-Thomas Fischer BFW 2010-2012
    AK BFW 2014-2015
    """

    fnames = cfgPath['pressurefileList']
    rasterSource = cfgPath['demSource']
    ProfileLayer = cfgPath['profileLayer']
    outpath = cfgPath['pathResult']
    DefaultName = cfgPath['projectName']

    xPath = rasterTransfo['x']
    yPath = rasterTransfo['y']
    zPath = rasterTransfo['z']
    sPath = rasterTransfo['s']

    indRunoutPoint = rasterTransfo['indRunoutPoint']
    sBeta = sPath[indRunoutPoint]
    runout = resAnalysis['runout'][0] + sBeta
    meanMaxDPP = resAnalysis['AMPP']
    maxMaxDPP = resAnalysis['MMPP']
    GI = resAnalysis['growthIndex']

    # includes flag for y axis -
    # 1 = mean pressure data
    # 2 = groth index
    # 3 = max pressure data
    flag = 3
    if (len(fnames) > 100):
        plotDensity = 1
    else:
        plotDensity = 0

    if flag == 1:
        title = 'Visualizing mean peak pressure data'
        log.info(title)
        tipo = 'rapp'
        data = meanMaxDPP / meanMaxDPP[0]
        yaxis_label = 'rAPP [-]'
        ytick_increment = 0.25
        ymax = 3
    elif flag == 2:
        title = 'Visualizing EGU growth index data'
        log.info(title)
        tipo = 'GI'
        data = GI
        yaxis_label = 'growth index [GI]'
        ytick_increment = 2
    elif flag == 3:
        title = 'Visualizing max peak pressure data'
        log.info(title)
        tipo = 'rmpp'
        data = maxMaxDPP / maxMaxDPP[0]
        yaxis_label = 'rMPP [-]'
        ytick_increment = 0.1
        ymax = max(data[:])+(max(data[:])-min(data[:]))*0.1
        ymin = min(data[:])-(max(data[:])-min(data[:]))*0.1
    else:
        log.error('Wrong flag')
        return None

    xlimProfAxis = max(sPath) + 50
    # color = cm.get_cmap('autumn', len(runout) + 3)
    color = cmapAimec(np.linspace(1, 0, len(runout) + 3, dtype=float))
    # Final result diagram - z_profile+data
    fig = plt.figure(figsize=(figW*2, figH))
    # fig.suptitle(title)
    mk = 0
#    show flow path
    ax1 = fig.add_subplot(111)
    ax1.set_title(title)
    ax1.set_ylabel(yaxis_label, color=color[1])
    ax1.spines['left'].set_color(color[1])
    ax1.tick_params(axis='y', colors=color[1])
    ax1.set_xlabel(''.join(['s [m] - runout with ', str(plim),
                            ' kPa threshold']), color='black')
    if plotDensity:  # estimate 2D histogram --> create pcolormesh
        nbins = 100
        H, xedges, yedges = np.histogram2d(runout, data, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
        dataDensity = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.Blues)
#        dataDensity = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.cool)
        cbar = plt.colorbar(dataDensity, orientation='horizontal')
        cbar.ax.set_ylabel('Counts')
    ax2 = ax1.twinx()
    ax2.set_ylabel('z [m]', color='k')
    ax2.spines['right'].set_color('k')
    ax2.spines['left'].set_color(color[1])
    ax2.tick_params(axis='y', colors='k')
    ax2.plot(sPath, zPath, color='k', label='path', linestyle='--')
    plt.xlim([0, xlimProfAxis])
    plt.ylim([math.floor(min(zPath)/10)*10, math.ceil(max(zPath)/10)*10])
    if not plotDensity:
        for k in range(len(runout)):
            topoName = cfgPath['projectName']
            pfarbe = color[k+1]  # (float(k), len(runout), colorflag)
            if k == 0:
                ax1.plot(runout[k], data[k], marker='+', linestyle='None',
                         markersize=2*ms, color='g', label='Reference')
            elif k == 1:
                ax1.plot(runout[k], data[k], marker=markers, label='sims',
                         color=pfarbe, linestyle='None')
            else:
                ax1.plot(runout[k], data[k], marker=markers, color=pfarbe,
                         linestyle='None')
            mk = mk+1
            if mk == len(markers):
                mk = 1
        ax1.legend(loc='lower left')
    plt.grid('on')
    outFileName = ''.join([cfgPath['dirName'], '_dptr',
                           str(int(plim)), '_', tipo])
    outname = os.path.join(cfgPath['pathResult'], 'pics', outFileName)

    if not os.path.exists(os.path.dirname(outname)):
        os.makedirs(os.path.dirname(outname))
    fig.savefig(outname)

    plt.close(fig)

    # Final result diagram - roc-plots
    rTP = resAnalysis['TP'] / (resAnalysis['TP'][0] + resAnalysis['FN'][0])
    rFP = resAnalysis['FP'] / (resAnalysis['TP'][0] + resAnalysis['FN'][0])

    fig = plt.figure(figsize=(figW, figH))
    # fig.suptitle('Normalized difference compared to reference')
    mk = 0
    ax1 = fig.add_subplot(111)
    ax1.set_title('Normalized difference compared to reference')
    ax1.set_ylabel('True positive rate')
    ax1.set_xlabel('False positive rate')
    if plotDensity:  # estimate 2D histogram --> create pcolormesh
        nbins = 100
        H, xedges, yedges = np.histogram2d(rFP, rTP, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
#        dataDensity = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.Blues)
        dataDensity = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.cool)
        cbar = plt.colorbar(dataDensity, orientation='horizontal')
        cbar.ax.set_ylabel('hit rate density')
    if not plotDensity:
        for k in range(len(rTP)):
            topoName = cfgPath['projectName']
            pfarbe = color[k+1]  # colorvar(float(k), len(rTP), colorflag)
            if k == 0:
                ax1.plot(rFP[k], rTP[k], color='g', label='Reference', marker='+',
                         markersize=2*ms, linestyle='None')
            elif k == 1:
                ax1.plot(rFP[k], rTP[k], marker=markers, label='sims',
                         color=pfarbe, linestyle='None')
            else:
                ax1.plot(rFP[k], rTP[k], marker=markers,
                         color=pfarbe, linestyle='None')
            mk = mk+1
            if mk == len(markers):
                mk = 0
        ax1.legend(loc='lower left')
    plt.xlim([-0.01, max(1, max(rFP))])
    plt.ylim([0, 1.01])
    plt.grid('on')

    outFileName = ''.join([cfgPath['dirName'], '_dptr', str(int(plim)), '_ROC'])
    outname = os.path.join(cfgPath['pathResult'], 'pics', outFileName)

    if not os.path.exists(os.path.dirname(outname)):
        os.makedirs(os.path.dirname(outname))
    fig.savefig(outname)

    plt.close(fig)

    return
