# -*- coding: utf-8 -*-
"""
Created on Mon Wed 12 2022

@author: dicko

Tools to extract information on the avalanche simulations run in the Output files

modified by AvaFrame
"""

# Python imports
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import configparser
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import math
import seaborn as sns
import logging
from cmcrameri import cm
import matplotlib.patheffects as pe

# Local imports
import avaframe.out3Plot.plotUtils as pU
import avaframe.in1Data.getInput as gI
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana5Utils.distanceTimeAnalysis as dtAna
import avaframe.in1Data.getInput as gI
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.outDistanceTimeAnalysis as dtAnaPlots
from avaframe.Tools import PostProcessingTools
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.geoTrans as gT
from avaframe.Tools import PlotTools
import avaframe.out1Peak.outPlotAllPeak as oP

log = logging.getLogger(__name__)


#%% Function to plot the peak flow quantities and the velocity thalweg envelope
def plotPeakVelVelThalwegEnvelope(avalancheDir, simIndex, simDF, rasterTransfo, dictVelAltThalweg,
    resTypePlots, modName, demData):
    """ plot peak flow fields and velocity thalweg envelope

        Parameters
        ------------
        avalancheDir: pathlib path or str
            path to avalanche directory
        simIndex: str
            index of current simulation
        rasterTransfo: dict
            info on domain transformation from xy to sl thalweg
        dictVelAltThalweg: dict
            info on velocity and altitude of particles along thalweg s
        resTypePlots: list
            list of result types that shall be plotted, one plot for each resType
        modName: str
            name of com module used to perform sims
        demData: dict
            dict of dem from sim
    """

    # convert avalancheDir to pathlib path
    avaDir = pathlib.Path(avalancheDir)
    # fetch name of simulation
    simName = simDF['simName'].loc[simIndex]

    # create a plot for result variables in resTypePlots
    for resType in resTypePlots:

        # initialize figure
        fig, ax = plt.subplots(1, 2, figsize=(pU.figW+10, pU.figH+3))

        # ax[0]
        # add peak field
        ax[0], rowsMinPlot, colsMinPlot = addPeakFieldConstrained(avaDir, modName, simName, resType, demData, ax[0], alpha=1.0)
        ax[0].plot(rasterTransfo['avaPath']['x'], rasterTransfo['avaPath']['y'], '-y', zorder=20, linewidth=1.0, path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()], label='thalweg')
        ax[0].legend(loc='upper right')

        # ax[1]
        # plot the thalweg diagram
        # First Y axis
        cmap = cm.vik  # colormap
        ax[1].plot(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['maxVelocity'], color=cmap(
            0.01), label='Max and min velocity')  # max
        ax[1].plot(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['minVelocity'], color=cmap(0.01))  # min
        ax[1].plot(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['meanVelocity'], linestyle='dashed',
                 color=cmap(0.7), markersize=0.5, label='Mean velocity', linewidth=0.5)  # mean
        ax[1].plot(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['medianVelocity'], linestyle='dotted',
                 color=cmap(1.0), markersize=0.7, label='Median velocity', linewidth=0.7)  # median
        # filling the space between the max and min
        ax[1].fill_between(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['minVelocity'],
                         dictVelAltThalweg['maxVelocity'], color=cmap(0.2), alpha=0.2, zorder=1)

        # y label
        ax[1].set_ylabel("Velocity [m/s]", fontsize=22)
        ax[1].tick_params(axis='both', labelsize=15)

        # Second Y axis
        ax2 = ax[1].twinx()
        # Plotting the altitude envelope - dashed lines
        ax2.plot(dictVelAltThalweg['sXYThalweg'],  dictVelAltThalweg['maxZ'], color='red',
                 linestyle='dashed', linewidth=0.4,label='Max and min altitude')
        ax2.plot(dictVelAltThalweg['sXYThalweg'],  dictVelAltThalweg['minZ'], color='red',
                 linestyle='dashed', linewidth=0.4)
        ax2.fill_between(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['minZ'], dictVelAltThalweg['maxZ'],hatch='|||||||', facecolor='black')
        # X and Y labels
        ax[1].set_xlabel('$s_{xy}$[m]\n\n', fontsize=22)
        ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax2.set_ylabel('Altitude [m]', fontsize=22)
        ax2.tick_params(axis='both', labelsize=15)

        # add legend
        fig.legend(loc='lower center', ncol=4, fancybox=True, shadow=True, fontsize=15)

        # save and or plot
        outDir = avaDir / 'Outputs' / modName / 'Plots'/'Peak flow quantities and Velocity altitude thalweg'
        fU.makeADir(outDir)
        name = 'PfVelAltThalweg_'+'_'+resType
        plotName = ('%s_%s' % (name, simIndex))
        plotPath = pU.saveAndOrPlot({"pathResult": outDir}, plotName, fig)


# %% Plot the peak flow velocity, thalweg time diagram and energy line
def plotPeakQuantThalTimeEnergyLine(avalancheDir, simIndex, simDF, rasterTransfo,
    dictRaster, modName, demSim, mtiInfo, cfgRangeTime):
    """ Create plot showing the resType peak field with thalweg,
        thalweg vs altitude with max peak field values along thalweg derived from peak fields
        and the tt-diagram

        Parameters
        -----------
        avalancheDir: str or pathlib path
            path to avalanche directory
        simIndex: str
            index of current sim in simDF
        simDF: dataFrame
            dataframe with one row per simulation, and all model config parameters
        rasterTransfo: dict
            dict with info on transformation from cartesian to thalweg coordinate system
        dictRaster: dict
            dict with info on peak fields
        modName: str
            name of com module used to perform sims
        demSim: dict
            dict with info on dem used for sims
    """

    # Load all infos from the peak files
    avaDir = pathlib.Path(avalancheDir)
    # fetch name of simulation
    simName = simDF['simName'].loc[simIndex]

    # initialize figure
    fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
    gs = fig.add_gridspec(2,2)
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[:, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    # ax1 Energy line
    # fetch energy line info for plot
    # TODO: shall we leave g definition here?
    g = 9.81 # gravitation constant
    EnergyLineRaster = rasterTransfo['z'] + pd.Series(dictRaster['maxPeakFlowVelocity']).values[0]**2 / (2*g) # energy line
    velocity = pd.Series(dictRaster['maxPeakFlowVelocity']).values[0]
    scat = ax1.scatter(rasterTransfo['s'], EnergyLineRaster, marker='s', cmap=pU.cmapRangeTime, s=8*pU.ms, c=velocity)
    ax1.plot(rasterTransfo['s'], rasterTransfo['z'], '-y.', zorder=20, linewidth=0.1, markersize=0.8,label='Thalweg altitude')

    # add colorbar
    cbar2 = ax1.figure.colorbar(scat, ax=ax1, use_gridspec=True)
    cbar2.ax.set_title('[m/s]', pad=10)
    cbar2.ax.set_ylabel('Max peak flow velocity')

    # draw the horizontal and vertical bars
    zLim = ax1.get_ylim()
    sLim = ax1.get_xlim()
    ax1.vlines(x=rasterTransfo['s'][0], ymin=rasterTransfo['z'][-1], ymax=rasterTransfo['z'][0],
               color='r', linestyle='--')
    ax1.hlines(y=rasterTransfo['z'][-1], xmin=0, xmax=rasterTransfo['s'][-1],
               color='r', linestyle='--')
    deltaz = rasterTransfo['z'][0] - rasterTransfo['z'][-1]
    deltas = rasterTransfo['s'][-1] - rasterTransfo['s'][0]
    alpha = np.arctan(deltaz/deltas)*(180/math.pi)
    # add textbox with angles, delta values
    textString = ('$\Delta z$=%s m\n$\Delta s_{xy}$=%s m\n' % (str(round(deltaz,1)), str(round(deltas,1)))) + r'$\alpha$=' + str(round(alpha,2)) + '°'
    ax1.text(0.98,0.9, textString, horizontalalignment='right',
        verticalalignment='top', fontsize=10, transform=ax1.transAxes, multialignment='left')
    X = [0,rasterTransfo['s'][-1]]
    Y = [rasterTransfo['z'][0],rasterTransfo['z'][-1]]
    ax1.plot(X,Y,color='black', linestyle='dashdot', linewidth=0.8)

    # Labels
    ax1.set_xlabel('$s_{xy}$ [m]', fontsize = 20)
    ax1.set_ylabel('z [m]', fontsize = 20)
    ax2.tick_params(axis='both', labelsize=15)
    ax1.set_xlim(sLim)
    ax1.set_ylim(zLim)
    ax1.legend(loc='upper right', fontsize = 10)
    ax1.set_title('Thalweg-Altitude')

    # ax2
    # add peak file plot
    ax2, rowsMinPlot, colsMinPlot = addPeakFieldConstrained(avaDir, modName, simName, 'pfv', demSim, ax2, alpha=1.0)
    # add the thalweg
    ax2.plot(rasterTransfo['avaPath']['x'], rasterTransfo['avaPath']['y'], '-y', zorder=20, linewidth=1.0, path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()], label='thalweg')

    # # add labels
    ax2.set_xlabel('x [m] \n\n', fontsize=20)
    ax2.set_ylabel('y [m] \n\n', fontsize=20)
    ax2.tick_params(axis='both', labelsize=15)

    # ax3 thalweg time diagram
    ax3, rangeTimeResType = dtAnaPlots.addRangeTimePlotToAxes(mtiInfo, cfgRangeTime, ax3)

    # save and or plot
    avaDir = pathlib.Path(avalancheDir)
    outDir = avaDir / 'Outputs' / modName / 'Plots'/ 'Energy line, Peak flow velocity and Thalweg time diagram'
    fU.makeADir(outDir)
    plotName = ('EnVelAltTT_%s' % (simName))
    plotPath = pU.saveAndOrPlot({"pathResult": outDir}, plotName, fig)
    log.info("Plot for %s successfully saved at %s" % (plotName, str(plotPath)))


# %% Plot the peak flow quantities with tracked particles and velocities
def plotPeakQuantTrackedPartVel(avalancheDir, simName, dictVelAltThalweg,
    dictVelAltThalwegPart, trackedPartProp, dictVelEnvelope, demSim, modName, rasterTransfo):
    """ Create plot showing particle properties over time and along avalanche thalweg
        in light blue envelope for all particles (filled between min and max values)
        in dark blue the values for tracked particles

        panel 1: map view of flow variable peak field
        panel 2: particle trajectoryLengthXYZ vs time
        panel 3: particle velocityMagnitude vs time
        panel 4: particle acceleration vs time
        panel 5: particle trajectoryLengthXYZ vs thalweg Sxy
        panel 6: particle velocity vs thalweg Sxy
        panel 7: particle acceleration vs thalweg Sxy

        Parameters
        -----------
        avalancheDir: pathlib path or str
            path to avalanche directory
        simName: str
            name of simulation
        dictVelAltThalweg: dict
            dict with velocity and altitude envelope info for all particles
        dictVelAltThalwegPart: dict
            dict with velocity and altitude envelope info for tracked particles
        trackedPartProp: dict
            dict with time series of tracked particle properties
        dictVelEnvelope: dict
            dict with velocity envelope info
        demSim: dict
            dict with sim dem info
        modName: str
            name of computational module that has been used to produce the sims
    """

    # create pathlib path
    avaDir = pathlib.Path(avalancheDir)

    # setup figure with subplots 7 panels
    fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
    gs = fig.add_gridspec(3,3)
    ax1 = fig.add_subplot(gs[0:2, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2, 1])
    ax5 = fig.add_subplot(gs[0, 2])
    ax6 = fig.add_subplot(gs[1, 2])
    ax7 = fig.add_subplot(gs[2, 2])
    ax8 = fig.add_subplot(gs[2, 0])

    # ax1
    # add peak file plot
    ax1, rowsMinPlot, colsMinPlot = addPeakFieldConstrained(avaDir, modName, simName, 'pfv', demSim, ax1, alpha=0.5)

    # add tracked particles locations over time
    xllcenter = demSim['header']['xllcenter']
    yllcenter = demSim['header']['yllcenter']
    cmap = cm.vik
    ax1.plot(trackedPartProp['x'][:,0]+xllcenter,trackedPartProp['y'][:,0]+yllcenter, zorder=1, linewidth=1.0, color=cmap(0.25), label='tracked particles')
    # ax1.plot(trackedPartProp['x']+xllcenter,trackedPartProp['y']+yllcenter, '-y', zorder=20, linewidth=1.0, path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()])
    ax1.plot(trackedPartProp['x']+xllcenter,trackedPartProp['y']+yllcenter, zorder=2, linewidth=1.0, color=cmap(0.25))

    # labels and ticks
    ax1.set_xlabel('x [m] \n\n', fontsize=22)
    ax1.set_ylabel('y [m]', fontsize=22)
    ax1.tick_params(axis='both', labelsize=13)
    ax1.legend(loc='upper right')

    # ax8
    ax8, rowsMinPlot, colsMinPlot = addPeakFieldConstrained(avaDir, modName, simName, 'pfv', demSim, ax8, alpha=0.5)
    ax8.plot(rasterTransfo['avaPath']['x'], rasterTransfo['avaPath']['y'], '-y', zorder=20, linewidth=1.0, path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()], label='thalweg')
    # labels and ticks
    ax8.set_xlabel('x [m] \n\n', fontsize=22)
    ax8.set_ylabel('y [m]', fontsize=22)
    ax8.tick_params(axis='both', labelsize=13)
    ax8.legend(loc='upper right')

    # ax2
    # # plot the trajectoryLengthXYZ vs time
    ax2.plot(trackedPartProp['t'], trackedPartProp['trajectoryLengthXYZ'], zorder=2, linewidth=1.0, color=cmap(0.25))
    ax2.plot(trackedPartProp['t'][:], trackedPartProp['trajectoryLengthXYZ'][:,0], zorder=1, linewidth=1.0, color=cmap(0.25), label='tracked particles')
    ax2.fill_between(dictVelEnvelope['Time'], dictVelEnvelope['SxyzMin'],
                    dictVelEnvelope['SxyzMax'], color=cmap(0.2), alpha=0.2, zorder=0)
    # labels and ticks
    ax2.set_xlabel('Time [s] \n\n', fontsize=15)
    ax2.set_ylabel('$trajectory_{XYZ}$ [m]', fontsize=15)
    ax2.tick_params(axis='both', labelsize=13)

    # ax3
    # plot velocity Magnitude of particles vs time
    cmap = cm.vik
    ax3.plot(trackedPartProp['t'], trackedPartProp['velocityMag'], zorder=2, linewidth=1.0, color=cmap(0.25), alpha=0.5)
    ax3.fill_between(dictVelEnvelope['Time'], dictVelEnvelope['Min'],
                    dictVelEnvelope['Max'], color=cmap(0.2), alpha=0.2, zorder=0)
    # labels and ticks
    ax3.set_xlabel('Time [s] \n\n', fontsize=15)
    ax3.set_ylabel('Velocity [m/s]', fontsize=15)
    ax3.tick_params(axis='both', labelsize=13)

    # ax4
    # plot acceleration of particles vs time
    ax4.fill_between(dictVelEnvelope['Time'], np.nanmin(dictVelEnvelope['Acc'], axis=1),
            np.nanmax(dictVelEnvelope['Acc'], axis=1), color=cmap(0.2), alpha=0.2, zorder=0, label='all particles')
    ax4.plot(trackedPartProp['t'][:], trackedPartProp['uAcc'][:,0], zorder=1, linewidth=1.0, color=cmap(0.25), label='tracked particles')
    ax4.plot(trackedPartProp['t'], trackedPartProp['uAcc'], zorder=1, linewidth=1.0, color=cmap(0.25))
    # labels and ticks
    ax4.set_xlabel('Time [s] \n\n', fontsize=15)
    ax4.set_ylabel('Acceleration [m/s²]', fontsize=15)
    ax4.tick_params(axis='both', labelsize=13)

    # add legend
    ax4.legend(loc='upper left', bbox_to_anchor=(0., -0.3), ncol=4, fancybox=True, shadow=True, fontsize=13)

    # ax5
    # plot travel length along thalweg
    ax5.fill_between(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['minSxyz'],
                     dictVelAltThalweg['maxSxyz'], color=cmap(0.2), alpha=0.2, zorder=0)
    ax5.plot(trackedPartProp['sAimec'], trackedPartProp['trajectoryLengthXYZ'], color=cmap(0.25), zorder=1)
    # labels and ticks
    ax5.set_xlabel('$S_{xy}$ (thalweg) [m]', fontsize=15)
    ax5.set_ylabel('$trajectory_{XYZ}$ [m]', fontsize=15)
    ax5.tick_params(axis='both', labelsize=13)

    # ax6
    # plot velocity along thalweg
    l2 = ax6.fill_between(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['minVelocity'],
                     dictVelAltThalweg['maxVelocity'], color=cmap(0.2), alpha=0.2, zorder=0, label='all particles')
    ax6.plot(trackedPartProp['sAimec'], trackedPartProp['velocityMag'], color=cmap(0.25), zorder=1)
    # labels and ticks
    ax6.set_xlabel('$S_{xy}$ [m]', fontsize=15)
    ax6.set_ylabel('Velocity [m/s]', fontsize=15)
    ax6.tick_params(axis='both', labelsize=13)

    # ax7
    # plot acceleration along the thalweg
    ax7.fill_between(dictVelAltThalweg['sXYThalweg'], dictVelAltThalweg['maxAcc'],
                     dictVelAltThalweg['minAcc'], color=cmap(0.2), alpha=0.2, zorder=0)
    ax7.plot(trackedPartProp['sAimec'], trackedPartProp['uAcc'], color=cmap(0.25), zorder=1)  # max
    # labels and ticks
    ax7.set_xlabel('$S_{xy}$ [m]', fontsize=15)
    ax7.set_ylabel('Acceleration [m/s²]', fontsize=15)
    ax7.tick_params(axis='both', labelsize=13)

    # save and or plot
    outDir = avaDir / 'Outputs' / modName / 'Plots'/'Peak flow quantities and Velocity of tracked particles'
    fU.makeADir(outDir)
    plotName = ('PfVelTrackedParticles_%s' % (simName))
    plotPath = pU.saveAndOrPlot({"pathResult": outDir}, plotName, fig)
    log.info("Plot for %s successfully saved at %s" % (plotName, str(plotPath)))


def addPeakFieldConstrained(avaDir, modName, simName, resType, demData, ax, alpha):
    """ use out1Peak functions to add plot of a peak field

        Parameters
        ------------
        avaDir: pathlib path
            path to avalanche dir
        modName: str
            name of computational module used to produce result peak fields
        simName: str
            name of simulation
        resType: str
            result variable name to look for peak field
        demData: dict
            dict with info on dem used to run sim of peak field data
        ax: matplotlib axes object
            axes where plot shall be added to
        alpha: float
            value for transparency from 0-1
    """

    # prepare dem data
    demField = np.where(demData["rasterData"] == demData["header"]["noDataValue"], np.nan,
        demData["rasterData"])

    # set input dir
    inputDir = avaDir / 'Outputs' / modName / 'peakFiles'

    # fetch all sims in inputDir with info on result files
    simResDF, resTypeList = fU.makeSimFromResDF(avaDir, modName, inputDir=inputDir)

    # find peakFile for resType
    peakFilePath = simResDF[resType].loc[simResDF['simName']==simName].values[0]

    # fetch cell size
    cellSize = simResDF['cellSize'].loc[simResDF['simName']==simName].values[0]

    # add constrained peak result field data plot to axes
    ax = oP.addConstrainedDataField(peakFilePath, resType, demField, ax, cellSize, alpha=alpha,
        setLimits=True)

    return ax
