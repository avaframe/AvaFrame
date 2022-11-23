# -*- coding: utf-8 -*-
"""
Created on Mon Wed 12 2022

@author: dicko

Tools to extract information on the avalanche simulations run in the Output files  

"""

# Python imports 
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
import pathlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib import pyplot as plt
import configparser
from matplotlib.ticker import FormatStrFormatter 
from scipy.interpolate import interp1d
import pandas as pd 
import math 

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
from avaframe.Tools import NodeTools
from avaframe.Tools import PlotTools



#%%  Function to plot the velocity altitude thalweg diagram 

def plotVelocityAltitudeThalweg(simu_number,dictVelAltThalweg,dictNodes,avaDict,avalancheDir,rasterTransfo,dictRaster,Title,Save,AvaNodes,TrackedPart,Raster,EnergyLine,modName):

    # figure 
    fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
    ax = fig.add_subplot(111)
    
    # First Y axis 
    cmap = cm.vik  # colormap
    ax.plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.01), label='Max and min velocity')  # max
    ax.plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number], color=cmap(0.01))  # min
    #ax.plot(sSortedDeduplicated[simu_number], simMeanVelThalweg[simu_number], linestyle='dashed',
    #         color=cmap(0.7), markersize=0.5, label='Mean velocity', linewidth=1.5)  # mean
    #ax.plot(sSortedDeduplicated[simu_number], simMedianVelThalweg[simu_number], linestyle='dotted',
    #         color=cmap(1.0), markersize=0.7, label='Median velocity', linewidth=1.5)  # median
    # filling the space between the max and min
    ax.fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                     dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1)
    # AvaNodes data
    if AvaNodes: 
        ax.plot(dictNodes['thalwegDomain']['sC07'], dictNodes['temporalDomain']['C07']['Velocity'], color='brown', label='AvaNode C07')
        ax.plot(dictNodes['thalwegDomain']['sC09'], dictNodes['temporalDomain']['C09']['Velocity'], color='green', label='AvaNode C09')
        ax.plot(dictNodes['thalwegDomain']['sC10'], dictNodes['temporalDomain']['C10']['Velocity'], color='orange', label='AvaNode C10')
         
    # raster data 
    if Raster:
        ax.plot(rasterTransfo['s'],pd.Series(dictRaster['maxPeakFlowVelocity']).values[0], color=cmap(0.3))
        ax.plot(rasterTransfo['s'],dictRaster['minPeakFlowVelocity'], color=cmap(0.3))
        ax.plot(rasterTransfo['s'],pd.Series(dictRaster['meanPeakFlowVelocity']).values[0], linestyle='dashed',
                 color=cmap(0.8), markersize=0.5, linewidth=1.5, label='Mean velocity of the raster file')
        ax.fill_between(rasterTransfo['s'],dictRaster['minPeakFlowVelocity'],pd.Series(dictRaster['maxPeakFlowVelocity']).values[0],hatch='|||||||', facecolor='black', label='velocity envelope of the pfv raster data') 
        #ax.set_xlim(0,rasterTransfo['s'][len(simMaxVelThalweg[simu_number])-1])  
        
    # y label 
    ax.set_ylabel("Velocity [m/s]", fontsize=30)
    ax.tick_params(axis='both', labelsize=18)
    
    # Second Y axis 
    ax2 = ax.twinx()
    # Plotting the altitude boxplot 
    medianprops = {}
    medianprops['linewidth'] = 0.7
    #ax2.boxplot(x=listZsorted[simu_number], positions=sSortedDeduplicated[simu_number],
    #           showfliers=False, medianprops=medianprops, widths=2.8) 
    
    #ax2.fill_between(sSortedDeduplicated[simu_number],simMinAltThalweg[simu_number],simMaxAltThalweg[simu_number],hatch='|||||||', facecolor='black')
    # Plotting the altitude envelope - dashed lines 
    #ax2.plot(sSortedDeduplicated[simu_number], simMaxAltThalweg[simu_number], color='red',
    #         linestyle='dashed', linewidth=0.5,label='Max and min altitude')
    #ax2.plot(sSortedDeduplicated[simu_number], simMinAltThalweg[simu_number], color='red',
    #         linestyle='dashed', linewidth=0.5)
    ax2.plot(rasterTransfo['s'], rasterTransfo['z'], color='black', linestyle='dashdot', linewidth=1.5,label='Thalweg')
    # raster 
    ax2.plot(rasterTransfo['s'], dictRaster['AltitudeCrossMax'], color='red', linestyle='dashdot', linewidth=0.5)
    ax2.plot(rasterTransfo['s'], dictRaster['AltitudeCrossMin'], color='red', linestyle='dashdot', linewidth=0.5,label='Altitude envelope from the raster file')
    # EnergyLines 
    index = len(dictVelAltThalweg['maxVelocity'][simu_number])
    EnergyLineParticles = rasterTransfo['z'][0:index] + np.array(dictVelAltThalweg['maxVelocity'][simu_number])**2/ (2*9.81) # energy line
    EnergyLineRaster = rasterTransfo['z'] + pd.Series(dictRaster['maxPeakFlowVelocity']).values[0]**2 / (2*9.81) # energy line
    ax2.plot(rasterTransfo['s'],EnergyLineRaster, color='red', linestyle='dashed', linewidth=0.8,label='Energy line from the raster files')
    ax2.plot(rasterTransfo['s'][0:index],EnergyLineParticles, color='brown', linestyle='dotted', linewidth=2.0,label='Energy line from the particle files')
    ax2.plot(rasterTransfo['s'], rasterTransfo['z'], color='black', linestyle='dashdot', linewidth=0.8,label='Altitude from the thalweg')

    
    # X and Y labels 
    ax.set_xlabel('$s_{xy}$ [m]\n', fontsize=30)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax2.set_ylabel('Altitude [m]', fontsize=30)
    ax2.tick_params(axis='both', labelsize=20)
    
    # legend     
    fig.legend(loc='lower center', ncol=7, fancybox=True, shadow=True, fontsize=14)   
    
    # saving the plot
    if Save:
        avaDir = pathlib.Path(avalancheDir)
        if TrackedPart==False:
            outDir = avaDir / 'Outputs' / modName / 'Plots'/'Velocity Altitude Thalweg' / 'Whole avalanche'
        else:
            outDir = avaDir / 'Outputs' / modName / 'Plots'/'Velocity Altitude Thalweg' / 'Tracked numerical particles'
        fU.makeADir(outDir)
        name = 'VelAltThalweg_'+avaDict[simu_number][0]['simName']
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName
        print("Plot for "+str(name)+" successfully saved at"+str(plotPath))
    
    
#%%  Function to plot the thalweg time diagram for each particle 

def plotThalwegTimeParticles(simu_number,dictVelEnvelope,dictVelAltThalweg,avaDict,avaDir,Save,Show,modName):
    fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
    ax = fig.add_subplot(111)
    sParticle = [None]*len(dictVelEnvelope['Time'][simu_number])

    for i in range (0,len(dictVelEnvelope['Time'][simu_number])):
        sParticle[i] = np.array(dictVelAltThalweg['sXYPart'][simu_number][i] - dictVelAltThalweg['sBetaPoint'][simu_number][i]) 
        
    plt.plot(dictVelEnvelope['Time'][simu_number],sParticle)
    ax.set_xlabel('Time[s]', fontsize=20)
    ax.set_ylabel('Projected distance to beta point ($s_{xy}$-sBetaPoint)[m]', fontsize=15)
    plt.gca().invert_yaxis()  
     
    # saving the plot
    if Save:
        outDir = avaDir / 'Outputs' / modName / 'Plots'/'Velocity Time' / 'Whole avalanche particles'
        fU.makeADir(outDir)
        name = 'VelAltThalweg_'+avaDict[simu_number][0]['simName']
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName
        print("Plot for "+str(name)+" successfully saved at"+str(plotPath))

    if Show:
        plt.show()


#%% Function to plot the boxplot of the difference between the mean simulation flow velocity and the average velocity of the three AvaNodes

def plotBoxplot(avalancheDir,Sim,number_ava,dictNodes,dictVelEnvelope,TrackedPart,Save,Show,modName):
    # Calculating the average velocity of the nodes 
    NumberNodes = [7, 9, 10]  # Numbers of the nodes you want to average
    meanNodeVelocity = NodeTools.averageNodesVelocity(NumberNodes,avalancheDir)
    
    # Calculating the difference between the mean avaframe velocity and the mean sensor velocity
    boxplot_data = [None]*number_ava
    for i in range(0, number_ava):
        # interpolating the node velocity
        f = interp1d(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], meanNodeVelocity)
        if max(dictVelEnvelope['Time'][i]) > max(dictNodes['temporalDomain']['C09']['TimeAfterBlasting']):
            index = dictVelEnvelope['Time'][i].index(
                next((j for j in dictVelEnvelope['Time'][i] if int(j) > max(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'])), None))
            vel_experimental_interpolated = f(dictVelEnvelope['Time'][i][0:index-1])
            boxplot_data[i] = dictVelEnvelope['Mean'][i][0:index-1] - vel_experimental_interpolated
        else:
            vel_experimental_interpolated = f(dictVelEnvelope['Time'][i])
            boxplot_data[i] = dictVelEnvelope['Mean'][i] - vel_experimental_interpolated
                
    # figure 
    fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
    ax = fig.add_subplot(111)
    # boxplot 
    plt.boxplot(boxplot_data, whis=200, zorder=2)
    # add a horizontal line at y=0 
    ax.axhline(y=0, color='black', linestyle='dashed', linewidth=0.8, zorder=0)
    # title and x labels  
    fig.suptitle('Velocity difference between the mean avaframe velocity and the Node velocity\nSeilbahnrinne', fontsize=22) 
    labels = [ PlotTools.labelFrictParamBoxplot(Sim, i) for i in range(0, number_ava) ]
    # tick labels 
    ax.set_xticklabels(labels, color='black', fontsize=15)
    ax.set_ylabel("Velocity difference [m/s]", color='black', fontsize=22)
    ax.tick_params(axis='both',colors='black', labelsize=15)
    
    if Show:
        plt.show()
    
    # saving the plot
    if Save:
        avaDir = pathlib.Path(avalancheDir)
        if TrackedPart==False:
            outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Comparison Plots'/'Boxplot'/'Temporal domain'/'Whole avalanche'
        else:
            outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Comparison Plots'/'Boxplot'/'Temporal domain'/'Tracked particles'
        fU.makeADir(outDir)
        name = 'Boxplot_'+'velocity_time_comparison'
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName
        print("Plot for "+str(name)+" successfully saved at"+str(plotPath))
        

#%% Function to plot the velocity time envelopes for all the simulations in the Outputs file 

def plotVelocityTimeEnvelope(number_ava,dictVelEnvelope,dictNodes,Sim,avaDir,AvaNodes,Show,modName): 

    # Choose your plot 
    Title = False # choose if you want a title or not 
    Save = True # choose if you want to save the plot or not 
    
    # Preparing the subplots regarding the number of avalanches found
    nrow,ncol,fig,ax = PlotTools.subplot(number_ava)
    
    # max simulation time, sometimes needed 
    #max_time = max(max(simTime[:]))
    
    # Colormap
    cmap = cm.vik
    rgba = cmap(0.2)
    
    # For 1D plots
    if nrow == 1 and ncol == 1:# add title, labels and ava Info
        ax.plot(dictVelEnvelope['Time'][0], dictVelEnvelope['Max'][0], color=cmap(0.01))  # max
        ax.plot(dictVelEnvelope['Time'][0], dictVelEnvelope['Min'][0], color=cmap(0.01), label='Maximum and minimum values')  # min
        ax.plot(dictVelEnvelope['Time'][0], dictVelEnvelope['Mean'][0],
                linestyle='dashed', color=cmap(0.7), markersize=1.8,  label='Mean')  # mean
        ax.plot(dictVelEnvelope['Time'][0], dictVelEnvelope['Median'][0],
                linestyle='dotted', color=cmap(0.7), markersize=1.8,  label='Median')  # median
        # filling the space between the max and min
        ax.fill_between(dictVelEnvelope['Time'][0], dictVelEnvelope['Min'][0],
                        dictVelEnvelope['Max'][0], color=cmap(0.2), alpha=0.2)
        # Experiment (AvaRange)
        if AvaNodes: 
            ax.plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']['VelocityAfterBlasting'],
                    color=cmap(0.80), label='AvaNode C07')  # AvaNode velocity
            ax.plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']['VelocityAfterBlasting'],
                    color=cmap(0.90), label='AvaNode C09')  # AvaNode velocity
            ax.plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']['VelocityAfterBlasting'],
                    color=cmap(0.95), label='AvaNode C10')  # AvaNode velocity
        # if AvaNodes: 
        #     ax.plot(velNodes['C07']['TimeAfterBlasting'], velNodes['C07']['VelocityAfterBlasting'], color='brown', label='AvaNode C07')
        #     ax.plot(velNodes['C09']['TimeAfterBlasting'], velNodes['C09']['VelocityAfterBlasting'], color='green', label='AvaNode C09')
        #     ax.plot(velNodes['C10']['TimeAfterBlasting'], velNodes['C10']['VelocityAfterBlasting'], color='orange', label='AvaNode C10')
        # title
        if Title:
            title = PlotTools.titleFrictParam(Sim, 0) 
            ax.set_title(title, fontsize=18)
        # ax.set_xlim(right=max_time)
    
    elif nrow == 1 or ncol == 1:
        simu_number = 0
        for i in range(0, max(nrow, ncol)):
            if i != max(nrow, ncol)-1:
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(0.01))  # max
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number], color=cmap(0.01))  # min
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Mean'][simu_number],
                           linestyle='dashed', color=cmap(0.7), markersize=1.8)  # mean
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Median'][simu_number],
                           linestyle='dotted', color=cmap(0.7), markersize=1.8)  # median
                # filling the space between the max and min
                ax[i].fill_between(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number],
                                   dictVelEnvelope['Max'][simu_number], color=cmap(0.2), alpha=0.2)
                # Experiment (AvaRange)
                if AvaNodes: 
                    ax[i].plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']
                               ['VelocityAfterBlasting'], color=cmap(0.80))  # AvaNode velocity
                    ax[i].plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']
                               ['VelocityAfterBlasting'], color=cmap(0.90))  # AvaNode velocity
                    ax[i].plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']
                               ['VelocityAfterBlasting'], color=cmap(0.95))  # AvaNode velocity
                # title
                if Title:
                    title = PlotTools.titleFrictParam(Sim, simu_number) 
                    ax[i].set_title(title, fontsize=18)  
                # ax[i].set_xlim(right=max_time)
                simu_number += 1
            else: 
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(0.01))  # max
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number], color=cmap(0.01), label='Maximum and minimum values')  # min
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Mean'][simu_number],
                           linestyle='dashed', color=cmap(0.7), markersize=1.8, label='Mean')  # mean
                ax[i].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Median'][simu_number],
                           linestyle='dotted', color=cmap(0.7), markersize=1.8, label='Median')  # median
                # filling the space between the max and min
                ax[i].fill_between(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number],
                                   dictVelEnvelope['Max'][simu_number], color=cmap(0.2), alpha=0.2)
                # Experiment (AvaRange)
                if AvaNodes: 
                    ax[i].plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']
                               ['VelocityAfterBlasting'], color=cmap(0.80), label='AvaNode C07')  # AvaNode velocity
                    ax[i].plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']
                              ['VelocityAfterBlasting'], color=cmap(0.90), label='AvaNode C09')  # AvaNode velocity
                    ax[i].plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']
                               ['VelocityAfterBlasting'], color=cmap(1.0), label='AvaNode C10')  # AvaNode velocity
                # title
                if Title:
                    title = PlotTools.titleFrictParam(Sim, simu_number) 
                    ax[i].set_title(title, fontsize=18)  
                #ax[i].set_xlim(right=120)
                simu_number += 1
                
    # For 2D plots
    else:
        simu_number = 0
        for i in range(0, nrow):
            for j in range(0, ncol):
                if i != nrow-1 or j!=ncol-1:
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(0.01))  # max
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number], color=cmap(0.01))  # min
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Mean'][simu_number], color=cmap(
                        0.7), linestyle='dashed', markersize=1.8)  # mean
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Median'][simu_number], linestyle='dotted', color=cmap(
                        1.0), markersize=1.8)  # median
                    # filling the space between the max and min
                    ax[i, j].fill_between(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(0.2), alpha=0.2)
                    # Experiment (AvaRange)
                    if AvaNodes: 
                        ax[i,j].plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']
                                     ['VelocityAfterBlasting'], color='brown')  # AvaNode velocity  #color=cmap(0.80)
                        ax[i,j].plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']
                                     ['VelocityAfterBlasting'], color='green')  # AvaNode velocity  #color=cmap(0.90)
                        ax[i,j].plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']
                                     ['VelocityAfterBlasting'], color='orange')  # AvaNode velocity   #color=cmap(1.0)
                    # title
                    title = PlotTools.titleFrictParam(Sim, simu_number) 
                    ax[i,j].set_title(title, fontsize=18)  
                    #ax[i,j].set_xlim(right=max_time)
                    simu_number += 1
                else: 
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(0.01))  # max
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number], color=cmap(0.01), label='Maximum and minimum velocity')  # min
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Mean'][simu_number], color=cmap(
                        0.7), linestyle='dashed', markersize=1.8, label='Mean velocity')  # mean
                    ax[i, j].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Median'][simu_number], linestyle='dotted', color=cmap(
                        1.0), markersize=1.8, label='Median velocity')  # median
                    # filling the space between the max and min
                    ax[i, j].fill_between(
                        dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1)
                    # Experiment (AvaRange)
                    if AvaNodes: 
                        ax[i,j].plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']
                                     ['VelocityAfterBlasting'], color='brown', label='AvaNode C07')  # AvaNode velocity  #color=cmap(0.80)
                        ax[i,j].plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']
                                     ['VelocityAfterBlasting'], color='green', label='AvaNode C09')  # AvaNode velocity  #color=cmap(0.90)
                        ax[i,j].plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']
                                     ['VelocityAfterBlasting'], color='orange', label='AvaNode C10')  # AvaNode velocity  #color=cmap(1.0)
                    # title
                    title = PlotTools.titleFrictParam(Sim, simu_number) 
                    ax[i,j].set_title(title, fontsize=18)  
                    #ax[i,j].set_xlim(right=max_time)
                    simu_number += 1
                    
    # set legend 
    fig.legend(loc='lower center', ncol=6, fancybox=True, shadow=True, fontsize=15)
    
    # Hide x labels and tick labels for top plots and y ticks for right plots.
    if nrow != 1 or ncol != 1:
        for axs in ax.flat:
            axs.set_ylabel("Velocity[m/s]", fontsize=22)
            axs.set_xlabel("Time[s]\n\n", fontsize=22)
            axs.tick_params(axis='both', labelsize=15)
            axs.label_outer()
    elif nrow==1 and ncol==1:
        ax.set_ylabel("Velocity[m/s]", fontsize=22)
        ax.set_xlabel("Time[s]\n\n", fontsize=22)
        ax.tick_params(axis='both', labelsize=15)
    else:
        for axs in ax.flat:
            axs.set_ylabel("Velocity[m/s]", fontsize=22)
            axs.set_xlabel("Time[s]\n\n", fontsize=22)
            axs.tick_params(axis='both', labelsize=15)
    
    # set spacing to subplots
    fig.tight_layout()
    
    # saving the plot
    if Save:
        outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Comparison Plots' / 'Velocity Envelopes' / 'Whole avalanche'
        fU.makeADir(outDir)
        name = 'Velocity_envelope'
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName
        print("Plot for "+str(name)+" successfully saved at"+str(plotPath))
    
    if Show:
        plt.show()


#%% Function to plot the range time diagram with peak flow velocity and the velocity time envelope 

def plotRangeTimePeakVelVelTimeEnvelope(Sim,avalancheDir,number_ava,dictVelEnvelope,dictNodes,Save,Show,AvaNodes,modName):
  
    avaDir = pathlib.Path(avalancheDir)
    # preparing simulation data
    simDF = cfgUtils.createConfigurationInfo(
        avalancheDir, standardCfg='', writeCSV=False, specDir='')
    
    for simu_number in range(0, number_ava):
    
        # Load all infos from the peak files
        inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
        peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)
        
        # simulation index
        sim = Sim.index[simu_number]
        simu_numb = np.where(peakFilesDF.simID == sim)
    
        for k in range(0, 3):
    
            fig, ax = plt.subplots(1, 3, figsize=(pU.figW+10, pU.figH+3))
    
            # ax[0]
            # preparing range time data for the simulation
            mtiInfo, cfgRangeTime = PlotTools.RangeTimeDiagram(
                avalancheDir, sim, simDF)
            # preparing range time data for the node
            NumberNodes = [10]  # Numbers of the nodes you want to plot
            DictNodRangeTime = NodeTools.findNodeRangeTime(
                NumberNodes, avalancheDir)
            # fetch required input info
            mti = mtiInfo['mti']
            rangeGates = mtiInfo['rangeGates']
            timeList = mtiInfo['timeList']
            rangeList = mtiInfo['rangeList']
            rangeTimeResType = cfgRangeTime['GENERAL']['rangeTimeResType']
            maxVel, rangeVel, timeVel = dtAna.approachVelocity(mtiInfo)
            # in case time steps are not ordered - the colormesh x and y need to be ordered
            timeIndex = np.argsort(np.array(timeList))
            timeListNew = np.array(timeList)[timeIndex]
            mti = mti[:, timeIndex]
            # fetch velocity legend style info
            width = cfgRangeTime['PLOTS'].getfloat('width')
            height = cfgRangeTime['PLOTS'].getfloat('height')
            lw = cfgRangeTime['PLOTS'].getfloat('lw')
            textsize = cfgRangeTime['PLOTS'].getfloat('textsize')
    
            # plotting
            pc = ax[0].pcolormesh(timeListNew, rangeGates,
                                  mti, cmap=pU.cmapRangeTime)
            ax[0].plot(timeList, rangeList, '.', color='black',
                        markersize=4, label='avalanche front')
            ax[0].set_xlabel('Time [s]')
            # add y label axis
            if mtiInfo['type'] == 'thalwegTime':
                sTypeCapital = mtiInfo['sType'][0].upper() + mtiInfo['sType'][1:]
                ax[0].set_ylabel('%s distance to %s [m]' %
                                  (sTypeCapital, mtiInfo['referencePointName']))
            else:
                ax[0].set_ylabel('Distance to %s [m]' %
                                  mtiInfo['referencePointName'])
    
            # add colorbar and infobox
            unit = pU.cfgPlotUtils['unit' + rangeTimeResType]
            if mtiInfo['type'] == 'thalwegTime' and cfgRangeTime['GENERAL']['maxOrMean'].lower() == 'max':
                avgType = 'max'
            else:
                avgType = 'avg.'
            cName = '%s ' % avgType + pU.cfgPlotUtils['name' + rangeTimeResType]
            pU.addColorBar(pc, ax[0], None, unit, title=cName)
            pU.putAvaNameOnPlot(ax[0], cfgRangeTime['GENERAL']['avalancheDir'])
            # add range time velocity legend
            dtAnaPlots.rangeTimeVelocityLegend(
                ax[0], maxVel, width, height, lw, textsize)
    
            # add max velocity location
            ax[0].plot(timeVel, rangeVel, 'r*', label='max velocity location')
    
            # add experimental AvaNode data
            ax[0].plot(DictNodRangeTime['C10']['Time'], DictNodRangeTime['C10']['RangeTime'],
                       color='red', linestyle='dashdot', linewidth=1, label='AvaNode location')  # velocity
    
            # add info on avalanche front in legend
            ax[0].legend(facecolor='grey', framealpha=0.2,
                          loc='upper left', fontsize=8)
    
            # ax[1]
            # Generate data for the simulation peakFile
            plotDict = PlotTools.PeakFields(
                avalancheDir, peakFilesDF, simu_numb[0][2+k], demData='')
            # choose colormap
            cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[plotDict['resType']], np.amin(
                plotDict['data']), np.amax(plotDict['data']), continuous=pU.contCmap)
            cmap.set_bad(alpha=0)
            # uncomment this to set the under value for discrete cmap transparent
            # cmap.set_under(alpha=0)
            xllcenter = plotDict['raster']['header']['xllcenter']
            yllcenter = plotDict['raster']['header']['yllcenter']
            rowsMinPlot = plotDict['rowsMin']*plotDict['cellSize'] + yllcenter
            rowsMaxPlot = (plotDict['rowsMax']+1)*plotDict['cellSize'] + yllcenter
            colsMinPlot = plotDict['colsMin']*plotDict['cellSize'] + xllcenter
            colsMaxPlot = (plotDict['colsMax']+1)*plotDict['cellSize'] + xllcenter
    
            extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]
    
            # add DEM hillshade with contour lines
            ls, CS = pU.addHillShadeContours(
                ax[1], plotDict['demConstrained'], plotDict['cellSize'], extent)
    
            # add peak field data
            im1 = ax[1].imshow(plotDict['data'], cmap=cmap, norm=norm,
                                extent=extent, origin='lower', aspect='equal', zorder=2)
            pU.addColorBar(im1, ax[1], ticks, plotDict['unit'])
    
            # add AvaNode data
            if AvaNodes:
                ax[1].plot(dictNodes['Coordinates']['e10'], dictNodes['Coordinates']['n10'], color='orange', label='AvaNode C10')
                ax[1].plot(dictNodes['Coordinates']['e09'], dictNodes['Coordinates']['n09'], color='green', label='AvaNode C09')
                ax[1].plot(dictNodes['Coordinates']['e07'], dictNodes['Coordinates']['n07'], color='brown', label='AvaNode C07')
    
            # # add title, labels and ava Info
            ax[1].set_xlabel('x [m] \n\n')
            ax[1].set_ylabel('y [m]')
    
            # ax[2]
            # Colormap
            cmap = cm.vik
            # Plotting the velocity envelope
            ax[2].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(
                0.01), label='Maximum and Minimum values')  # max
            ax[2].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number],
                        color=cmap(0.01))  # min
            ax[2].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Mean'][simu_number], linestyle='dashed',
                        color=cmap(0.7), markersize=1.8, label='Mean')  # mean
            ax[2].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Median'][simu_number], linestyle='dotted',
                        color=cmap(0.7), markersize=1.8, label='Median')  # median
            # filling the space between the max and min
            ax[2].fill_between(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number],
                                dictVelEnvelope['Max'][simu_number], color=cmap(0.2), alpha=0.2)
            # Experiment (AvaRange)
            ax[2].plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']['VelocityAfterBlasting'],
                      color='brown')  # AvaNode velocity
            ax[2].plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']['VelocityAfterBlasting'],
                       color='green')  # AvaNode velocity
            ax[2].plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']['VelocityAfterBlasting'],
                       color='orange')  # AvaNode velocity
            ax[2].set_ylabel("Velocity[m/s]")
            ax[2].set_xlabel("Time[s]\n\n")
    
            # Add title, legends and save the plot
            # title
            PlotTools.titleFrictParam(Sim, simu_number) 
            # add legend
            fig.legend(loc='lower center', ncol=3,
                        fancybox=True, shadow=True, fontsize=10)
            
            # show plot
            if Show:
                plt.show()
                
            # save plot
            if Save:
                avaDir = pathlib.Path(avalancheDir)
                outDir = avaDir / 'Outputs' / modName / 'Plots' /'Range time diagram, Peak flow quantities and Velocity envelope'
                fU.makeADir(outDir)
                plotName = outDir / ('%s.%s' % (plotDict['name'], pU.outputFormat))
                fig.savefig(plotName)
                plotPath = pathlib.Path.cwd() / plotName
                print("Plot for "+str(plotDict['name'])+" successfully saved at"+str(plotPath))


#%% Function to plot the peak flow quantities and the velocity thalweg envelope 

def plotPeakVelVelThalwegEnvelope(Sim,avalancheDir,number_ava,dictNodes,rasterTransfo,dictVelAltThalweg,Save,AvaNodes,modName):

    for simu_number in range(0, number_ava):

        # Load all infos from the peak files
        avaDir = pathlib.Path(avalancheDir)
        inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
        peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

        sim = Sim.index[simu_number]
        simu_numb = np.where(peakFilesDF.simID == sim)

        for k in range(0, 3):

            fig, ax = plt.subplots(1, 2, figsize=(pU.figW+10, pU.figH+3))

            # ax[0]
            # Generate data for the simulation peakFile
            plotDict = PlotTools.PeakFields(
                avalancheDir, peakFilesDF, simu_numb[0][3+k], demData='')
            # choose colormap
            cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[plotDict['resType']], np.amin(
                plotDict['data']), np.amax(plotDict['data']), continuous=pU.contCmap)
            cmap.set_bad(alpha=0)
            # uncomment this to set the under value for discrete cmap transparent
            # cmap.set_under(alpha=0)
            xllcenter = plotDict['raster']['header']['xllcenter']
            yllcenter = plotDict['raster']['header']['yllcenter']
            rowsMinPlot = plotDict['rowsMin']*plotDict['cellSize'] + yllcenter
            rowsMaxPlot = (plotDict['rowsMax']+1)*plotDict['cellSize'] + yllcenter
            colsMinPlot = plotDict['colsMin']*plotDict['cellSize'] + xllcenter
            colsMaxPlot = (plotDict['colsMax']+1)*plotDict['cellSize'] + xllcenter

            extent = [0, colsMaxPlot-colsMinPlot, 0, rowsMaxPlot-rowsMinPlot]

            # add DEM hillshade with contour lines
            ls, CS = pU.addHillShadeContours(
                ax[0], plotDict['demConstrained'], plotDict['cellSize'], extent)

            # add peak field data
            im1 = ax[0].imshow(plotDict['data'], cmap=cmap, norm=norm,
                               extent=extent, origin='lower', aspect='equal', zorder=2)
            pU.addColorBar(im1, ax[0], ticks, plotDict['unit'])

            # add AvaNode data
            if AvaNodes: 
                ax[0].plot(dictNodes['Coordinates']['e10']-colsMinPlot, dictNodes['Coordinates']['n10']-rowsMinPlot, color='orange')
                ax[0].plot(dictNodes['Coordinates']['e09']-colsMinPlot, dictNodes['Coordinates']['n09']-rowsMinPlot, color='green')
                ax[0].plot(dictNodes['Coordinates']['e07']-colsMinPlot, dictNodes['Coordinates']['n07']-rowsMinPlot, color='brown')
            
            # add center of mass path 
            ax[0].plot(rasterTransfo['x']-colsMinPlot, rasterTransfo['y']-rowsMinPlot, '-y.', zorder=20, linewidth=0.3, markersize=2.5, label='Avalanche thalweg')

            # labels and ticks 
            ax[0].set_xlabel('x [m] \n\n', fontsize=22)
            ax[0].set_ylabel('y [m]', fontsize=22) 
            ax[0].tick_params(axis='both', labelsize=13)
            
            # ax[1]
            # plot the thalweg diagram
            # First Y axis 
            cmap = cm.vik  # colormap
            ax[1].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(
                0.01), label='Max and min velocity')  # max
            ax[1].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number], color=cmap(0.01))  # min
            ax[1].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['meanVelocity'][simu_number], linestyle='dashed',
                     color=cmap(0.7), markersize=0.5, label='Mean velocity', linewidth=0.5)  # mean
            ax[1].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['medianVelocity'][simu_number], linestyle='dotted',
                     color=cmap(1.0), markersize=0.7, label='Median velocity', linewidth=0.7)  # median
            # filling the space between the max and min
            ax[1].fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                             dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1)
            # AvaNode
            if AvaNodes: 
                ax[1].plot(dictNodes['thalwegDomain']['sC07'],dictNodes['temporalDomain']['C07']['Velocity'], color='brown', label='AvaNode C07')
                ax[1].plot(dictNodes['thalwegDomain']['sC09'],dictNodes['temporalDomain']['C09']['Velocity'], color='green', label='AvaNode C09')
                ax[1].plot(dictNodes['thalwegDomain']['sC10'],dictNodes['temporalDomain']['C10']['Velocity'], color='orange', label='AvaNode C10')
            # y label 
            ax[1].set_ylabel("Velocity [m/s]", fontsize=22)
            ax[1].tick_params(axis='both', labelsize=15)
            
            # Second Y axis 
            ax2 = ax[1].twinx()
            # Plotting the altitude boxplot 
            #medianprops = {}
            #medianprops['linewidth'] = 0.7
            #ax2.boxplot(x=listZsorted[simu_number], positions=sSortedDeduplicated[simu_number],
            #           showfliers=False, medianprops=medianprops, widths=2.8) 
            # Plotting the altitude envelope - dashed lines 
            ax2.plot(dictVelAltThalweg['sXYThalweg'][simu_number],  dictVelAltThalweg['maxZ'][simu_number], color='red',
                     linestyle='dashed', linewidth=0.4,label='Max and min altitude')
            ax2.plot(dictVelAltThalweg['sXYThalweg'][simu_number],  dictVelAltThalweg['minZ'][simu_number], color='red',
                     linestyle='dashed', linewidth=0.4)
            ax2.fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minZ'][simu_number], dictVelAltThalweg['maxZ'][simu_number],hatch='|||||||', facecolor='black')
            # X and Y labels 
            ax[1].set_xlabel('$s_{xy}$[m]\n\n', fontsize=22)
            ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax2.set_ylabel('Altitude [m]', fontsize=22)
            ax2.tick_params(axis='both', labelsize=15)

            # add legend
            fig.legend(loc='lower center', ncol=4, fancybox=True, shadow=True, fontsize=15)
                
            # title 
            PlotTools.titleFrictParam(Sim, simu_number)
            
            # saving the plot
            if Save:
                outDir = avaDir / 'Outputs' / modName / 'Plots'/'Peak flow quantities and Velocity altitude thalweg'
                fU.makeADir(outDir)
                name = 'PfVelAltThalweg_'+'_'+peakFilesDF['resType'][7*simu_number+k+3]
                plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
                fig.savefig(plotName)
                plotPath = pathlib.Path.cwd() / plotName
                print("Plot for "+str(name)+" successfully saved at"+str(plotPath))


#%% Function to plot the peak flow quantities and the velocity time envelope 

def plotPeakVelVelEnvelope(Sim,avalancheDir,number_ava,dictNodes,dictVelEnvelope,Save,AvaNodes,Show,modName):

    for simu_number in range(0, number_ava):

        # Load all infos from the peak files
        avaDir = pathlib.Path(avalancheDir)
        inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
        peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

        sim = Sim.index[simu_number]
        simu_numb = np.where(peakFilesDF.simID == sim)

        for k in range(0, 3):

            fig, ax = plt.subplots(1, 2, figsize=(pU.figW+10, pU.figH+3))

            # ax[0]
            # Generate data for the simulation peakFile
            plotDict = PlotTools.PeakFields(
                avalancheDir, peakFilesDF, simu_numb[0][3+k], demData='')
            # choose colormap
            cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[plotDict['resType']], np.amin(
                plotDict['data']), np.amax(plotDict['data']), continuous=pU.contCmap)
            cmap.set_bad(alpha=0)
            # uncomment this to set the under value for discrete cmap transparent
            # cmap.set_under(alpha=0)
            xllcenter = plotDict['raster']['header']['xllcenter']
            yllcenter = plotDict['raster']['header']['yllcenter']
            rowsMinPlot = plotDict['rowsMin']*plotDict['cellSize'] + yllcenter
            rowsMaxPlot = (plotDict['rowsMax']+1)*plotDict['cellSize'] + yllcenter
            colsMinPlot = plotDict['colsMin']*plotDict['cellSize'] + xllcenter
            colsMaxPlot = (plotDict['colsMax']+1)*plotDict['cellSize'] + xllcenter

            #extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]
            extent = [0, colsMaxPlot-colsMinPlot, 0, rowsMaxPlot-rowsMinPlot]

            # add DEM hillshade with contour lines
            ls, CS = pU.addHillShadeContours(
                ax[0], plotDict['demConstrained'], plotDict['cellSize'], extent)

            # add peak field data
            im1 = ax[0].imshow(plotDict['data'], cmap=cmap, norm=norm,
                               extent=extent, origin='lower', aspect='equal', zorder=2)
            pU.addColorBar(im1, ax[0], ticks, plotDict['unit'])
            
            # add AvaNode data
            if AvaNodes: 
                ax[0].plot(dictNodes['Coordinates']['e10']-colsMinPlot, dictNodes['Coordinates']['n10']-rowsMinPlot, color='orange', label='AvaNode C10')
                ax[0].plot(dictNodes['Coordinates']['e09']-colsMinPlot, dictNodes['Coordinates']['n09']-rowsMinPlot, color='green', label='AvaNode C09')
                ax[0].plot(dictNodes['Coordinates']['e07']-colsMinPlot, dictNodes['Coordinates']['n07']-rowsMinPlot, color='brown', label='AvaNode C07')

            # # add title, labels and ava Info
            # #title = str('%s' % name)
            # #ax.set_title(title +'\n')
            ax[0].set_xlabel('x [m] \n\n', fontsize=22)
            ax[0].set_ylabel('y [m] \n\n', fontsize=22)
            ax[0].tick_params(axis='both', labelsize=15)
            # ax[1].set_ylabel('y [m]')

            # ax[1]
            # Colormap
            cmap = cm.vik
            # Plotting the velocity envelope
            ax[1].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Max'][simu_number], color=cmap(
                0.01))  # max
            ax[1].plot(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number],
                       color=cmap(0.01))  # min
            #ax[1].plot(simTime[simu_number], simMeanVel[simu_number], linestyle='dashed',
            #           color=cmap(0.7), markersize=1.8, label='Mean')  # mean
            #ax[1].plot(simTime[simu_number], simMedianVel[simu_number], linestyle='dotted',
            #           color=cmap(0.7), markersize=1.8, label='Median')  # median
            # filling the space between the max and min
            ax[1].fill_between(dictVelEnvelope['Time'][simu_number], dictVelEnvelope['Min'][simu_number],
                               dictVelEnvelope['Max'][simu_number], color=cmap(0.2), alpha=0.2, label='Velocity envelope')
            # Experiment (AvaRange)
            if AvaNodes: 
                ax[1].plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']['VelocityAfterBlasting'],
                           color='brown')  # AvaNode velocity
                ax[1].plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']['VelocityAfterBlasting'],
                           color='green')  # AvaNode velocity
                ax[1].plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']['VelocityAfterBlasting'],
                           color='orange')  # AvaNode velocity
            ax[1].set_ylabel("Velocity[m/s]", fontsize=22)
            ax[1].set_xlabel("Time[s]\n\n", fontsize=22)
            ax[1].tick_params(axis='both', labelsize=15)

            # title 
            PlotTools.titleFrictParam(Sim, simu_number)
                
            # add legend
            fig.legend(loc='lower center', ncol=4,
                       fancybox=True, shadow=True, fontsize=15)
            
            if Show:
                plt.show()
            else:
                plt.close() 
            
            # save data
            if Save: 
                avaDir = pathlib.Path(avalancheDir)
                outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Peak flow quantities and Velocity envelope'
                fU.makeADir(outDir)
                name = 'PfVel_'+peakFilesDF['resType'][7*simu_number+k+3]
                plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
                fig.savefig(plotName)
                plotPath = pathlib.Path.cwd() / plotName 
                print("Plot for "+str(str(name))+" successfully saved at"+str(plotPath))


# %% Plotting velocity-altitude-thalweg envelope for all the different friction parameters

def plotVelocityAltitudeThalwegAllSim(dictVelAltThalweg,Sim,dictNodes,avaDir,number_ava,TrackedPart,Title,Save,Show,AvaNodes,modName):

    # Preparing the subplots regarding the number of avalanches found
    nrow,ncol,fig,ax = PlotTools.subplot(number_ava)

    # Colormap
    cmap = cm.vik

    # For 1D plots
    if nrow == 1 and ncol == 1:# add title, labels and ava Info
        ax.plot(dictVelAltThalweg['sXYThalweg'][0], dictVelAltThalweg['maxVelocity'][0], color=cmap(0.01), label='Maximum and minimum velocity')  # max
        ax.plot(dictVelAltThalweg['sXYThalweg'][0], dictVelAltThalweg['minVelocity'][0], color=cmap(0.01))  # min
        ax.plot(dictVelAltThalweg['sXYThalweg'][0], dictVelAltThalweg['meanVelocity'][0], linestyle='dashed',
                 color=cmap(0.7), markersize=0.5, label='Mean velocity', linewidth=1.5)  # mean
        ax.plot(dictVelAltThalweg['sXYThalweg'][0], dictVelAltThalweg['medianVelocity'][0], linestyle='dotted',
                 color=cmap(1.0), markersize=0.7, label='Median velocity', linewidth=1.5)  # median
        # filling the space between the max and min
        ax.fill_between(dictVelAltThalweg['sXYThalweg'][0], dictVelAltThalweg['minVelocity'][0],
                         dictVelAltThalweg['maxVelocity'][0], color=cmap(0.2), alpha=0.2, zorder=1)
        # AvaNode
        if AvaNodes: 
            ax.plot(dictNodes['thalwegDomain']['sC07'],dictNodes['temporalDomain']['C07']['Velocity'], color='brown', label='AvaNode C07')
            ax.plot(dictNodes['thalwegDomain']['sC09'],dictNodes['temporalDomain']['C09']['Velocity'], color='green', label='AvaNode C09')
            ax.plot(dictNodes['thalwegDomain']['sC10'],dictNodes['temporalDomain']['C10']['Velocity'], color='orange', label='AvaNode C10')
        # title
        title = PlotTools.titleFrictParam(Sim, 0) 
        ax.set_title(title, fontsize=18)
        # ax.set_xlim(right=max_time)

    elif nrow == 1 or ncol == 1:
        simu_number = 0
        for i in range(0, max(nrow, ncol)):
            if i != max(nrow, ncol)-1:
                ax[i].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.01))  # max
                ax[i].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number], color=cmap(0.01))  # min
                #ax[i].plot(sSortedDeduplicated[]simu_number, simMeanVelThalweg[simu_number], linestyle='dashed',
                #         color=cmap(0.7), markersize=0.5, linewidth=1.5)  # mean
                #ax[i].plot(sSortedDeduplicated[simu_number], simMedianVelThalweg[simu_number], linestyle='dotted',
                #         color=cmap(1.0), markersize=0.7, linewidth=1.5)  # median
                # filling the space between the max and min
                ax[i].fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                                 dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1)
                # AvaNode
                if AvaNodes: 
                    ax[i].plot(dictNodes['thalwegDomain']['sC07'],dictNodes['temporalDomain']['C07']['Velocity'], color='brown')
                    ax[i].plot(dictNodes['thalwegDomain']['sC09'],dictNodes['temporalDomain']['C09']['Velocity'], color='green')
                    ax[i].plot(dictNodes['thalwegDomain']['sC10'],dictNodes['temporalDomain']['C10']['Velocity'], color='orange')
                # title
                title = PlotTools.titleFrictParam(Sim, simu_number) 
                ax[i].set_title(title, fontsize=18)  
                # ax[i].set_xlim(right=max_time)
                simu_number += 1
            else: 
                ax[i].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.01), label='Maximum and minimum velocity')  # max
                ax[i].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number], color=cmap(0.01))  # min
                #ax[i].plot(sSortedDeduplicated[simu_number], simMeanVelThalweg[simu_number], linestyle='dashed',
                #         color=cmap(0.7), markersize=0.5, label='Mean velocity', linewidth=1.5)  # mean
                #ax[i].plot(sSortedDeduplicated[simu_number], simMedianVelThalweg[simu_number], linestyle='dotted',
                #         color=cmap(1.0), markersize=0.7, label='Median velocity', linewidth=1.5)  # median
                # filling the space between the max and min
                ax[i].fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                                 dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1)
                # AvaNode
                if AvaNodes: 
                    ax[i].plot(dictNodes['thalwegDomain']['sC07'],dictNodes['temporalDomain']['C07']['Velocity'], color='brown', label='AvaNode C07')
                    ax[i].plot(dictNodes['thalwegDomain']['sC09'],dictNodes['temporalDomain']['C09']['Velocity'], color='green', label='AvaNode C09')
                    ax[i].plot(dictNodes['thalwegDomain']['sC10'],dictNodes['temporalDomain']['C10']['Velocity'], color='orange', label='AvaNode C10')
                # title
                title = PlotTools.titleFrictParam(Sim, simu_number) 
                ax[i].set_title(title, fontsize=18)  
                #ax[i].set_xlim(right=120)
                simu_number += 1
                
    # For 2D plots
    else:
        simu_number = 0
        for i in range(0, nrow):
            for j in range(0, ncol):
                if i != nrow-1 or j!=ncol-1:
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.01))  # max
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number], color=cmap(0.01))  # min
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['meanVelocity'][simu_number], linestyle='dashed',
                             color=cmap(0.7), markersize=0.5, linewidth=1.0)  # mean
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['medianVelocity'][simu_number], linestyle='dotted',
                             color=cmap(1.0), markersize=0.5, linewidth=1.0)  # median
                    # filling the space between the max and min
                    ax[i,j].fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                                     dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1)
                    # AvaNode
                    if AvaNodes: 
                        ax[i,j].plot(dictNodes['thalwegDomain']['sC07'],dictNodes['temporalDomain']['C07']['Velocity'], color='brown')
                        ax[i,j].plot(dictNodes['thalwegDomain']['sC09'],dictNodes['temporalDomain']['C09']['Velocity'], color='green')
                        ax[i,j].plot(dictNodes['thalwegDomain']['sC10'],dictNodes['temporalDomain']['C10']['Velocity'], color='orange')
                    # title
                    title = PlotTools.titleFrictParam(Sim, simu_number) 
                    ax[i,j].set_title(title, fontsize=18)  
                    #ax[i,j].set_xlim(right=max_time)
                    simu_number += 1
                else: 
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.01), label='Maximum and minimum velocity')  # max
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number], color=cmap(0.01))  # min
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['meanVelocity'][simu_number], linestyle='dashed',
                             color=cmap(0.7), markersize=0.5, label='Mean velocity', linewidth=1.0)  # mean
                    ax[i,j].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['medianVelocity'][simu_number], linestyle='dotted',
                             color=cmap(1.0), markersize=0.5, label='Median velocity', linewidth=1.0)  # median
                    # filling the space between the max and min
                    ax[i,j].fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                                     dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1)
                    # AvaNode
                    if AvaNodes: 
                        ax[i,j].plot(dictNodes['thalwegDomain']['sC07'],dictNodes['temporalDomain']['C07']['Velocity'], color='brown', label='AvaNode C07')
                        ax[i,j].plot(dictNodes['thalwegDomain']['sC09'],dictNodes['temporalDomain']['C09']['Velocity'], color='green', label='AvaNode C09')
                        ax[i,j].plot(dictNodes['thalwegDomain']['sC10'],dictNodes['temporalDomain']['C10']['Velocity'], color='orange', label='AvaNode C10')
                    # title
                    title = PlotTools.titleFrictParam(Sim, simu_number) 
                    ax[i,j].set_title(title, fontsize=18)  
                    #ax[i,j].set_xlim(right=max_time)
                    simu_number += 1

    # set title 
    if Title: 
        fig.suptitle('Velocity envelope for the whole avalanche')
    # set legend 
    fig.legend(loc='lower center', ncol=6, fancybox=True, shadow=True, fontsize=15)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    if nrow != 1 or ncol != 1:
        for axs in ax.flat:
            axs.set_ylabel("Velocity[m/s]", fontsize=22)
            axs.set_xlabel("$S_{xy}$[m]\n\n", fontsize=22)
            axs.tick_params(axis='both', labelsize=15)
            axs.label_outer()
    elif nrow==1 and ncol==1:
        ax.set_ylabel("Velocity[m/s]", fontsize=22)
        ax.set_xlabel("$S_{xy}$[m]\n\n", fontsize=22)
        ax.tick_params(axis='both', labelsize=15)
    else:
        for axs in ax.flat:
            axs.set_ylabel("Velocity[m/s]", fontsize=22)
            axs.set_xlabel("$S_{xy}$[m]\n\n", fontsize=22)
            axs.tick_params(axis='both', labelsize=15)

    # set spacing to subplots
    fig.tight_layout()

    # saving the plot
    if Save:
        if TrackedPart==False:
            outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Comparison Plots' / 'Velocity Thalweg Envelopes' / 'Whole avalanche'
        else:
            outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Comparison Plots' / 'Velocity Thalweg Envelopes' / 'Tracked numerical particles'
        fU.makeADir(outDir)
        name = 'Velocity_envelopes_all_simulations'
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName
        print("Plot for "+str(name)+" successfully saved at"+str(plotPath))

    if Show:
        plt.show()




# %% Function to generate the boxplot of the difference between the mean simulation flow velocity and the average velocity
# of the three AvaNodes along the thalweg 

def plotBoxplotThalweg(avalancheDir,number_ava,Sim,dictNodes,dictVelAltThalweg,TrackedPart,Save,modName):

    # Calculating the average velocity of the nodes 
    NumberNodes = [7, 9, 10]  # Numbers of the nodes you want to average
    meanNodeVelocity = NodeTools.averageNodesVelocity(NumberNodes,avalancheDir)
    
    # Calculating the difference between the mean avaframe velocity and the mean sensor velocity
    boxplot_data = [None]*number_ava
    for i in range(0, number_ava):
        # interpolating the node velocity
        f = interp1d(dictNodes['thalwegDomain']['sC09AfterBlasting'], meanNodeVelocity)
        if max(dictVelAltThalweg['sXYThalweg'][i]) > max(dictNodes['thalwegDomain']['sC09AfterBlasting']) or min(dictVelAltThalweg['sXYThalweg'][i]) < min(dictNodes['thalwegDomain']['sC09AfterBlasting']):
            indexMin = 0 
            indexMax = len(dictVelAltThalweg['sXYThalweg'][i])-1
            if min(dictVelAltThalweg['sXYThalweg'][i]) < min(dictNodes['thalwegDomain']['sC09AfterBlasting']):
                indexMin = dictVelAltThalweg['sXYThalweg'][i].index(
                    next((j for j in dictVelAltThalweg['sXYThalweg'][i] if int(j) > min(dictNodes['thalwegDomain']['sC09AfterBlasting'])), None))
            if max(dictVelAltThalweg['sXYThalweg'][i]) > max(dictNodes['thalwegDomain']['sC09AfterBlasting']):
                indexMax = dictVelAltThalweg['sXYThalweg'][i].index(
                    next((j for j in dictVelAltThalweg['sXYThalweg'][i] if int(j) > max(dictNodes['thalwegDomain']['sC09AfterBlasting'])), None))
            vel_experimental_interpolated = f(dictVelAltThalweg['sXYThalweg'][i][indexMin:indexMax-1])
            boxplot_data[i] = dictVelAltThalweg['meanVelocity'][i][indexMin:indexMax-1] - vel_experimental_interpolated
        
        else:
            vel_experimental_interpolated = f(dictVelAltThalweg['sXYThalweg'][i])
            boxplot_data[i] = dictVelAltThalweg['meanVelocity'][i] - vel_experimental_interpolated
                
    # figure 
    fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
    ax = fig.add_subplot(111)
    # boxplot 
    plt.boxplot(boxplot_data, whis=200, zorder=2)
    # add a horizontal line at y=0 
    ax.axhline(y=0, color='black', linestyle='dashed', linewidth=0.8, zorder=0)
    # title and x labels  
    fig.suptitle('Velocity difference between the mean avaframe velocity and the Node velocity along the thalweg') 
    labels = [ PlotTools.labelFrictParamBoxplot(Sim, i) for i in range(0, number_ava) ]
    # tick labels 
    ax.set_xticklabels(labels, color='black', fontsize=15)
    ax.set_ylabel("Velocity difference [m/s]", color='black', fontsize=22)
    ax.tick_params(axis='both',colors='black', labelsize=15)
    
    plt.show()
    
    # saving the plot
    if Save:
        avaDir = pathlib.Path(avalancheDir)
        if TrackedPart==False:
            outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Comparison Plots'/'Boxplot'/'Along the thalweg'/'Whole avalanche'
        else:
            outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Comparison Plots'/'Boxplot'/'Along the thalweg'/'Tracked particles'
        fU.makeADir(outDir)
        name = 'Boxplot_Thalweg'
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName
        print("Plot for "+str(name)+" successfully saved at"+str(plotPath))
    
    
# %% Plot the energy line 

def plotEnergyLine(avalancheDir,avaDict,simu_number,Sim,dictVelAltThalweg,rasterTransfo,dictRaster,Show,Save,modName):

    # figure 
    fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
    ax = fig.add_subplot(111)

    index = len(dictVelAltThalweg['maxVelocity'][simu_number])

    # First Y axis 
    g = 9.8 # gravitation constant 
    EnergyLineParticles = rasterTransfo['z'][0:index] + np.array(dictVelAltThalweg['maxVelocity'][simu_number])**2/ (2*g) # energy line
    EnergyLineRaster = rasterTransfo['z'] + pd.Series(dictRaster['maxPeakFlowVelocity']).values[0]**2 / (2*g) # energy line
    ax.plot(rasterTransfo['s'],EnergyLineRaster, color='red', linestyle='dashed', linewidth=0.8,label='Energy line from the raster files')
    ax.plot(rasterTransfo['s'][0:index],EnergyLineParticles, color='brown', linestyle='dotted', linewidth=2.0,label='Energy line from the particle files')
    ax.plot(rasterTransfo['s'], rasterTransfo['z'], color='black', linestyle='dashdot', linewidth=0.8,label='Altitude from the thalweg')
    # x and y label
    ax.set_xlabel("$S_{xy}$ [m]\n\n", fontsize=30)  
    ax.set_ylabel("Energy [kg.m²/s²]", fontsize=30)
    ax.tick_params(axis='both', labelsize=18)
    ax.set_xlim(0,rasterTransfo['s'][index-1])

    # legend     
    fig.legend(loc='lower center', ncol=3, fancybox=True, shadow=True, fontsize=14)   

    # title 
    title = PlotTools.titleFrictParam(Sim, simu_number)  
    ax.set_title(title, fontsize=18)
        
    # saving the plot
    if Save:
        avaDir = pathlib.Path(avalancheDir)
        outDir = avaDir / 'Outputs' / modName / 'Plots'/'Energy Line' / 'Whole avalanche'
        fU.makeADir(outDir)
        name = 'VelAltThalweg_'+'_'+avaDict[simu_number][0]['simName']
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName 
        print("Plot for "+str(name)+" successfully saved at"+str(plotPath))

    if Show:
        plt.show() 
    
    
# %% Plot the velocity altitude thalweg and thalweg time diagram
 
def plotVelAltThalTimeDiag(avalancheDir,number_ava,simu_number,avaDict,Sim,dictVelAltThalweg,rasterTransfo,dictNodes,Show,Save,AvaNodes,modName):

    # preparing simulation data for the thalweg time diagram 
    simDF = cfgUtils.createConfigurationInfo(
        avalancheDir, standardCfg='', writeCSV=False, specDir='')

    for simu_number in range(0, number_ava):

        # simulation index 
        sim = Sim.index[simu_number]

        fig, ax = plt.subplots(1, 2, figsize=(pU.figW+10, pU.figH+3))

        # ax[0]
        # First Y axis 
        cmap = cm.vik  # colormap
        ax[0].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.01))  # max
        ax[0].plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number], color=cmap(0.01))  # min
        #ax[0].plot(sSortedDeduplicated[simu_number], simMeanVelThalweg[simu_number], linestyle='dashed',
        #         color=cmap(0.7), markersize=0.5, label='Mean velocity', linewidth=1.5)  # mean
        #ax[0].plot(sSortedDeduplicated[simu_number], simMedianVelThalweg[simu_number], linestyle='dotted',
        #         color=cmap(1.0), markersize=0.7, label='Median velocity', linewidth=1.5)  # median
        # filling the space between the max and min
        ax[0].fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                         dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=1, label='Velocity envelope')
        # AvaNode
        if AvaNodes: 
            ax[0].plot(dictNodes['thalwegDomain']['sC07'],dictNodes['temporalDomain']['C07']['Velocity'], color='brown', label='AvaNode C07')
            ax[0].plot(dictNodes['thalwegDomain']['sC09'],dictNodes['temporalDomain']['C09']['Velocity'], color='green', label='AvaNode C09')
            ax[0].plot(dictNodes['thalwegDomain']['sC10'],dictNodes['temporalDomain']['C10']['Velocity'], color='orange', label='AvaNode C10')
                
                
        # y label 
        ax[0].set_ylabel("Velocity [m/s]", fontsize=20)
        ax[0].tick_params(axis='both', labelsize=15)
      
        # Second Y axis 
        ax2 = ax[0].twinx()
        # Plotting the altitude boxplot 
        medianprops = {}
        medianprops['linewidth'] = 0.7
        #ax2.boxplot(x=listZsorted[simu_number], positions=sSortedDeduplicated[simu_number],
        #           showfliers=False, medianprops=medianprops, widths=2.8) 
      
        ax2.fill_between(dictVelAltThalweg['sXYThalweg'][simu_number],dictVelAltThalweg['minZ'][simu_number],dictVelAltThalweg['maxZ'][simu_number],hatch='|||||||', facecolor='black',label='Altitude envelope')
        # Plotting the altitude envelope - dashed lines 
        ax2.plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['maxZ'][simu_number], color='red',
                 linestyle='dashed', linewidth=0.5)
        ax2.plot(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minZ'][simu_number], color='red',
                 linestyle='dashed', linewidth=0.5)
        ax2.plot(rasterTransfo['s'], rasterTransfo['z'], color='black', linestyle='dashdot', linewidth=1.5,label='Thalweg')

        # X and Y labels np.array(
        ax[0].set_xlabel('$s_{xy}$ [m]\n', fontsize=20)
        ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax2.set_ylabel('Altitude [m]', fontsize=20)
        ax2.tick_params(axis='both', labelsize=15)
        # add legend
        fig.legend(loc='lower center', ncol=6,
                   fancybox=True, shadow=True, fontsize=15)
        
        # ax[1]
        # preparing range time data for the simulation
        mtiInfo, cfgRangeTime = PlotTools.ThalwegTimeDiagram(avalancheDir, sim, simDF)
        
        # fetch required input info
        mti = mtiInfo['mti']
        rangeGates = mtiInfo['rangeGates']
        timeList = mtiInfo['timeList']
        rangeList = mtiInfo['rangeList']
        rangeTimeResType = cfgRangeTime['GENERAL']['rangeTimeResType']
        maxVel, rangeVel, timeVel = dtAna.approachVelocity(mtiInfo)
        # in case time steps are not ordered - the colormesh x and y need to be ordered
        timeIndex = np.argsort(np.array(timeList))
        timeListNew = np.array(timeList)[timeIndex]
        mti = mti[:, timeIndex]
        # fetch velocity legend style info
        width = cfgRangeTime['PLOTS'].getfloat('width')
        height = cfgRangeTime['PLOTS'].getfloat('height')
        lw = cfgRangeTime['PLOTS'].getfloat('lw')
        textsize = cfgRangeTime['PLOTS'].getfloat('textsize')

        # plotting
        pc = ax[1].pcolormesh(timeListNew, rangeGates,
                              mti, cmap=pU.cmapRangeTime)
        ax[1].plot(timeList, rangeList, '.', color='black',
                    markersize=4, label='avalanche front')
        ax[1].set_xlabel('Time [s]\n\n', fontsize=20)
        # add y label axis
        mtiInfo['type'] == 'thalwegTime'
        ax[1].set_ylabel('$s_{xy}$ [m]', fontsize=20)

        # add colorbar and infobox
        unit = pU.cfgPlotUtils['unit' + rangeTimeResType]
        if mtiInfo['type'] == 'thalwegTime' and cfgRangeTime['GENERAL']['maxOrMean'].lower() == 'max':
            avgType = 'max'
        else:
            avgType = 'avg.'
        cName = '%s ' % avgType + pU.cfgPlotUtils['name' + rangeTimeResType]
        pU.addColorBar(pc, ax[1], None, unit, title=cName)
        pU.putAvaNameOnPlot(ax[1], cfgRangeTime['GENERAL']['avalancheDir'])
        # add range time velocity legend
        dtAnaPlots.rangeTimeVelocityLegend(
            ax[1], maxVel, width, height, lw, textsize)

        # add max velocity location
        ax[1].plot(timeVel, rangeVel, 'r*', label='max velocity location')
        
        # add info on avalanche front in legend
        ax[1].legend(facecolor='grey', framealpha=0.2,
                      loc='lower right', fontsize=8)

        # add experimental AvaNode data    
        if AvaNodes: 
            ax[1].plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['thalwegDomain']['sC07AfterBlasting'], color='brown', label='AvaNode C07')  # AvaNode velocity
            ax[1].plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['thalwegDomain']['sC09AfterBlasting'], color='green', label='AvaNode C09')  # AvaNode velocity
            ax[1].plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['thalwegDomain']['sC10AfterBlasting'], color='orange', label='AvaNode C10')  # AvaNode velocity

        
        # if tt-diagram add beta point info
        #if mtiInfo['type'] == 'thalwegTime':
        # invert y axis as ava flow starts from minus distance to beta point
        ax[1].invert_yaxis()
        ax[1].axhline(y=0.0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax[1].tick_params(axis='both', labelsize=15)
       
        # title 
        PlotTools.titleFrictParam(Sim, simu_number)
        
        # save data
        if Save: 
            avaDir = pathlib.Path(avalancheDir)
            outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Velocity altitude and Thalweg time diagram'
            fU.makeADir(outDir)
            name = 'VelAltTT_'+'_'+avaDict[simu_number][0]['simName']
            plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
            fig.savefig(plotName)
            plotPath = pathlib.Path.cwd() / plotName 
            print("Plot for "+str(name)+" successfully saved at"+str(plotPath))

        if Show:
            plt.show()
    
    
# %% Plot the peak flow velocity, thalweg time diagram and energy line 

def plotPeakQuantThalTimeEnergyLine(avalancheDir,number_ava,Sim,plotDict,rasterTransfo,dictRaster,dictNodes,Save,Show,AvaNodes,modName):
   

    # preparing simulation data
    simDF = cfgUtils.createConfigurationInfo(
        avalancheDir, standardCfg='', writeCSV=False, specDir='')

    for simu_number in range(0, number_ava):

        # Load all infos from the peak files
        avaDir = pathlib.Path(avalancheDir)
        inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
        peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

        sim = Sim.index[simu_number]
        simu_numb = np.where(peakFilesDF.simID == sim)
        
        # preparing range time data for the simulation
        mtiInfo, cfgRangeTime = PlotTools.ThalwegTimeDiagram(avalancheDir, sim, simDF)
        # preparing the beta point data 
        x = mtiInfo['betaPoint'][0]
        y = mtiInfo['betaPoint'][1]
        
        for k in range(0, 3):

                #fig, ax = plt.subplots(1, 2, figsize=(pU.figW+10, pU.figH+3))
                fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
                gs = fig.add_gridspec(2,2)
                ax1 = fig.add_subplot(gs[0, 1])
                ax2 = fig.add_subplot(gs[:, 0])
                ax3 = fig.add_subplot(gs[1, 1])
            
                # ax1 Energy line

                # TO MAKE A PLOT SIMILAR TO MATTHIAS WORK, WITH VELOCITY PER CENTRES OF MASS
                # plot mass averaged center of mass 
                # energyLineDict = PostProcessingTools.energyLinePostProcessing() 
                # avaProfileMass = energyLineDict['avaProfileMass']
                # particlesIni = energyLineDict['particlesIni']
                # coefExt = energyLineDict['coefExt']
                # slopeExt = energyLineDict['slopeExt']
                # zEne = energyLineDict['zEne']
                # runOutAngleDeg = energyLineDict['runOutAngleDeg']
                # u2Path = energyLineDict['u2Path']
                # sGeomL = energyLineDict['sGeomL']
                # zGeomL = energyLineDict['zGeomL']
                # alphaDeg = energyLineDict['alphaDeg']
                # energyLineTestCfg = energyLineDict['energyLineTestCfg']
                # mu = energyLineDict['mu']
                # sIntersection = energyLineDict['sIntersection']
                # zIntersection = energyLineDict['zIntersection']
                # ax1.plot(avaProfileMass['s'], avaProfileMass['z'], '-y.', label='Center of mass altitude',
                #         lw=1, path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()])
                # extend this curve towards the bottom using a linear regression on the last x points
                # ax1.plot(avaProfileMass['s'][-1]*np.array([1, 1+coefExt]), avaProfileMass['z'][-1] +
                #         slopeExt*avaProfileMass['s'][-1]*np.array([0, coefExt]), ':k',
                #         label='Center of mass altitude extrapolation')
                # add center of mass velocity points and runout line
                # ax1.plot(avaProfileMass['s'][[0, -1]], zEne[[0, -1]], '-r', label='com1dfa energy line (%.2f°)' % runOutAngleDeg) 
                # scat = ax1.scatter(avaProfileMass['s'], zEne, marker='s', cmap=cmap, s=8*pU.ms, c=u2Path/(2*g),
                #                   label='Center of mass velocity altitude')
                # ax.plot(rasterTransfo['s'],pd.Series(maxPeakFlowVelocity).values[0], color=cmap(0.3))
                
                # Generate data for the simulation peakFile
                plotDict = PlotTools.PeakFields(
                    avalancheDir, peakFilesDF, simu_numb[0][3+k], demData='')
                # choose colormap
                cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[plotDict['resType']], np.amin(
                    plotDict['data']), np.amax(plotDict['data']), continuous=pU.contCmap)
                cmap.set_bad(alpha=0)
                
                g = 9.8 # gravitation constant 
                EnergyLineRaster = rasterTransfo['z'] + pd.Series(dictRaster['maxPeakFlowVelocity']).values[0]**2 / (2*g) # energy line
                velocity = pd.Series(dictRaster['maxPeakFlowVelocity']).values[0]
                scat = ax1.scatter(rasterTransfo['s'], EnergyLineRaster, marker='s', cmap=pU.cmapRangeTime, s=8*pU.ms, c=velocity)
                ax1.plot(rasterTransfo['s'], rasterTransfo['z'], '-y.', zorder=20, linewidth=0.1, markersize=0.8,label='Thalweg altitude')
                
                cbar2 = ax1.figure.colorbar(scat, ax=ax1, use_gridspec=True)
                cbar2.ax.set_title('[m/s]', pad=10)
                cbar2.ax.set_ylabel('Max peak flow velocity')
            
                # draw the horizontal and vertical bars 
                zLim = ax1.get_ylim()
                sLim = ax1.get_xlim()
                ax1.vlines(x=rasterTransfo['s'][0], ymin=rasterTransfo['z'][-1], ymax=rasterTransfo['z'][0],
                           color='r', linestyle='--')
                #ax1.vlines(x=sIntersection, color='b', ymin=zMin, ymax=zIntersection, linestyle='--')        
                ax1.hlines(y=rasterTransfo['z'][-1], xmin=0, xmax=rasterTransfo['s'][-1],
                           color='r', linestyle='--')
                deltaz = rasterTransfo['z'][0] - rasterTransfo['z'][-1]
                deltas = rasterTransfo['s'][-1] - rasterTransfo['s'][0]
                alpha = np.arctan(deltaz/deltas)*(180/math.pi) 
                ax1.text(750,800,'$\Delta z$='+str(round(deltaz,1))+'m', fontsize=10)
                ax1.text(750,750,'$\Delta s_{xy}$='+str(round(deltas,1))+'m', fontsize=10)
                ax1.text(750,700,r'$\alpha$='+str(round(alpha,2))+'°', fontsize=10)
                X = [0,rasterTransfo['s'][-1]]
                Y = [rasterTransfo['z'][0],rasterTransfo['z'][-1]]
                ax1.plot(X,Y,color='black', linestyle='dashdot', linewidth=0.8)
                #ax1.hlines(y=zIntersection, color='b', xmin=0, xmax=sIntersection, linestyle='--')
            
                # Labels 
                ax1.set_xlabel('$s_{xy}$ [m]', fontsize = 20)
                ax1.set_ylabel('z [m]', fontsize = 20)
                ax2.tick_params(axis='both', labelsize=15)
                ax1.set_xlim(sLim)
                ax1.set_ylim(zLim)
                ax1.legend(loc='upper right', fontsize = 10)
                ax1.set_title('Thalweg-Altitude')
                
                
                # ax2 Peak flow velocity 
                # choose colormap
                cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[plotDict['resType']], np.amin(
                    plotDict['data']), np.amax(plotDict['data']), continuous=pU.contCmap)
                cmap.set_bad(alpha=0)
                # uncomment this to set the under value for discrete cmap transparent
                # cmap.set_under(alpha=0)
                xllcenter = plotDict['raster']['header']['xllcenter']
                yllcenter = plotDict['raster']['header']['yllcenter']
                rowsMinPlot = plotDict['rowsMin']*plotDict['cellSize'] + yllcenter
                rowsMaxPlot = (plotDict['rowsMax']+1)*plotDict['cellSize'] + yllcenter
                colsMinPlot = plotDict['colsMin']*plotDict['cellSize'] + xllcenter
                colsMaxPlot = (plotDict['colsMax']+1)*plotDict['cellSize'] + xllcenter  
            
                #extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]
                extent = [0, colsMaxPlot-colsMinPlot, 0, rowsMaxPlot-rowsMinPlot]
            
                # add DEM hillshade with contour lines
                ls, CS = pU.addHillShadeContours(
                    ax2, plotDict['demConstrained'], plotDict['cellSize'], extent)
            
                # add peak field data
                im1 = ax2.imshow(plotDict['data'], cmap=cmap, norm=norm,
                                   extent=extent, origin='lower', aspect='equal', zorder=2)
                pU.addColorBar(im1, ax2, ticks, plotDict['unit'])
                # add the thalweg
                ax2.plot(rasterTransfo['x']-colsMinPlot, rasterTransfo['y']-rowsMinPlot, '-y.', zorder=20, linewidth=0.1, markersize=0.8, label='Avalanche thalweg')
                #ax2.set_xlim(0, colsMaxPlot-colsMinPlot)
                #ax2.set_ylim(0, rowsMaxPlot-rowsMinPlot)
                
                # add AvaNode data
                if AvaNodes: 
                    ax2.plot(dictNodes['Coordinates']['e10']-colsMinPlot, dictNodes['Coordinates']['n10']-rowsMinPlot, color='orange', label='AvaNode C10')
                    ax2.plot(dictNodes['Coordinates']['e09']-colsMinPlot, dictNodes['Coordinates']['n09']-rowsMinPlot, color='green', label='AvaNode C09')
                    ax2.plot(dictNodes['Coordinates']['e07']-colsMinPlot, dictNodes['Coordinates']['n07']-rowsMinPlot, color='brown', label='AvaNode C07')
            
                # # add labels
                ax2.set_xlabel('x [m] \n\n', fontsize=20)
                ax2.set_ylabel('y [m] \n\n', fontsize=20)
                ax2.tick_params(axis='both', labelsize=15)
                
                
                # ax3 thalweg time diagram    
                # Generate data for the simulation peakFile
                plotDict = PlotTools.PeakFields(avalancheDir, peakFilesDF, simu_numb[0][3+k], demData='')
                xllcenter = plotDict['raster']['header']['xllcenter']
                yllcenter = plotDict['raster']['header']['yllcenter']
                distance = np.sqrt((x+xllcenter-rasterTransfo['gridx'])**2 + (y+yllcenter-rasterTransfo['gridy'])**2)  
                sbetaPoint,lbetaPoint = np.unravel_index(np.argmin(distance, axis=None), distance.shape)
                # fetch required input info
                mti = mtiInfo['mti']
                rangeGates = mtiInfo['rangeGates'] + sbetaPoint
                timeList = mtiInfo['timeList']
                rangeList = mtiInfo['rangeList'] 
                rangeTimeResType = cfgRangeTime['GENERAL']['rangeTimeResType']
                maxVel, rangeVel, timeVel = dtAna.approachVelocity(mtiInfo)
                # in case time steps are not ordered - the colormesh x and y need to be ordered
                timeIndex = np.argsort(np.array(timeList))
                timeListNew = np.array(timeList)[timeIndex]
                mti = mti[:, timeIndex]
                # fetch velocity legend style info
                width = cfgRangeTime['PLOTS'].getfloat('width')
                height = cfgRangeTime['PLOTS'].getfloat('height')
                lw = cfgRangeTime['PLOTS'].getfloat('lw')
                textsize = cfgRangeTime['PLOTS'].getfloat('textsize')
            
                # plotting
                pc = ax3.pcolormesh(timeListNew, rangeGates,
                                      mti, cmap=pU.cmapRangeTime)
                ax3.plot(timeList, rangeList+sbetaPoint, '.', color='black',
                            markersize=4, label='avalanche front')
                ax3.set_xlabel('Time [s]\n\n', fontsize=20)
                ax3.set_title('Thalweg-Time')
                # add y label axis
                mtiInfo['type'] == 'thalwegTime'
                ax3.set_ylabel('$s_{xy}$ [m]', fontsize=20)
            
                # add colorbar and infobox
                unit = pU.cfgPlotUtils['unit' + rangeTimeResType]
                if mtiInfo['type'] == 'thalwegTime' and cfgRangeTime['GENERAL']['maxOrMean'].lower() == 'max':
                    avgType = 'max'
                else:
                    avgType = 'avg.'
                cName = '%s ' % avgType + pU.cfgPlotUtils['name' + rangeTimeResType]
                pU.addColorBar(pc, ax3, None, unit, title=cName)
                pU.putAvaNameOnPlot(ax3, cfgRangeTime['GENERAL']['avalancheDir'])
                # add range time velocity legend
                dtAnaPlots.rangeTimeVelocityLegend(
                    ax3, maxVel, width, height, lw, textsize)
            
                # add max velocity location
                ax3.plot(timeVel, rangeVel+sbetaPoint, 'r*', label='max velocity location')
                
                # add info on avalanche front in legend
                ax3.legend(facecolor='grey', framealpha=0.2,
                              loc='lower right', fontsize=8)
            
                # add experimental AvaNode data    
                if AvaNodes: 
                    ax3.plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['thalwegDomain']['sC07AfterBlasting'], color='brown', label='AvaNode C07')  # AvaNode velocity
                    ax3.plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['thalwegDomain']['sC09AfterBlasting'], color='green', label='AvaNode C09')  # AvaNode velocity
                    ax3.plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['thalwegDomain']['sC10AfterBlasting'], color='orange', label='AvaNode C10')  # AvaNode velocity
            
                
                # if tt-diagram add beta point info
                #if mtiInfo['type'] == 'thalwegTime':
                # invert y axis as ava flow starts from minus distance to beta point
                ax3.invert_yaxis()
                ax3.axhline(y=0.0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
                ax3.set_ylim(500)
                ax3.set_xlim(0,80)
                ax3.tick_params(axis='both', labelsize=15)
               
                # title 
                PlotTools.titleFrictParam(Sim, simu_number)
                
                # save data
                if Save: 
                    avaDir = pathlib.Path(avalancheDir)
                    outDir = avaDir / 'Outputs' / modName / 'Plots' / 'Energy line, Peak flow velocity and Thalweg time diagram'
                    fU.makeADir(outDir)
                    name = 'EnVelAltTT_'+peakFilesDF['resType'][7*simu_number+k+3]
                    plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
                    fig.savefig(plotName)
                    plotPath = pathlib.Path.cwd() / plotName 
                    print("Plot for "+str(name)+" successfully saved at"+str(plotPath))
                
                if Show:
                    plt.show()
            
            
# %% Plot the peak flow quantities with tracked particles and velocities 

def plotPeakQuantTrackedPartVel(avalancheDir,number_ava,Sim,avaDict,dictVelAltThalweg,dictVelAltThalwegPart,trackedPartProp,trackedPartPropAdapted,dictVelEnvelope,dictNodes,Save,AvaNodes,Show,modName):

    for simu_number in range(0, number_ava):

        # Load all infos from the peak files
        avaDir = pathlib.Path(avalancheDir)
        inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
        peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

        sim = Sim.index[simu_number]
        simu_numb = np.where(peakFilesDF.simID == sim)


        #fig, ax = plt.subplots(1, 2, figsize=(pU.figW+10, pU.figH+3))
        fig = plt.figure(figsize=(pU.figW+10, pU.figH+3))
        gs = fig.add_gridspec(3,3)
        ax1 = fig.add_subplot(gs[:, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 1])
        ax4 = fig.add_subplot(gs[2, 1])
        ax5 = fig.add_subplot(gs[0, 2])
        ax6 = fig.add_subplot(gs[1, 2])
        ax7 = fig.add_subplot(gs[2, 2])

        # ax1
        # Generate data for the simulation peakFile
        plotDict = PlotTools.PeakFields(
            avalancheDir, peakFilesDF, simu_numb[0][3], demData='')
        # choose colormap
        cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[plotDict['resType']], np.amin(
            plotDict['data']), np.amax(plotDict['data']), continuous=pU.contCmap)
        cmap.set_bad(alpha=0)
        # uncomment this to set the under value for discrete cmap transparent
        cmap.set_under(alpha=0)
        xllcenter = plotDict['raster']['header']['xllcenter']
        yllcenter = plotDict['raster']['header']['yllcenter']
        rowsMinPlot = plotDict['rowsMin']*plotDict['cellSize'] + yllcenter
        rowsMaxPlot = (plotDict['rowsMax']+1)*plotDict['cellSize'] + yllcenter
        colsMinPlot = plotDict['colsMin']*plotDict['cellSize'] + xllcenter
        colsMaxPlot = (plotDict['colsMax']+1)*plotDict['cellSize'] + xllcenter

        extent = [0, colsMaxPlot-colsMinPlot, 0, rowsMaxPlot-rowsMinPlot]

        # add DEM hillshade with contour lines
        ls, CS = pU.addHillShadeContours(
            ax1, plotDict['demConstrained'], plotDict['cellSize'], extent)

        # add peak field data
        import matplotlib.colors as colors 
        im1 = ax1.imshow(plotDict['data'], cmap='Blues',
                          extent=extent, origin='lower', aspect='equal', norm= colors.Normalize(vmin=-1, vmax=0.1), 
                          zorder=10, alpha=0.15)
        #pU.addColorBar(im1, ax1, ticks, plotDict['unit'])

        # add AvaNode data
        if AvaNodes: 
            ax1.plot(dictNodes['Coordinates']['e07']-colsMinPlot, dictNodes['Coordinates']['n07']-rowsMinPlot, color='brown', zorder=25, label='AvaNode C07')
            ax1.plot(dictNodes['Coordinates']['e09']-colsMinPlot, dictNodes['Coordinates']['n09']-rowsMinPlot, color='green', zorder=25, label='AvaNode C09')
            ax1.plot(dictNodes['Coordinates']['e10']-colsMinPlot, dictNodes['Coordinates']['n10']-rowsMinPlot, color='orange', zorder=25, label='AvaNode C10')
        
        # add center of mass path 
        cmap = cm.vik  # colormap
        #ax[0].plot(rasterTransfo['x']-colsMinPlot, rasterTransfo['y']-rowsMinPlot, '-y.', zorder=20, linewidth=0.3, markersize=2.5, label='Avalanche thalweg')
        ax1.plot(trackedPartProp[0]['x'][:,0]+xllcenter-colsMinPlot,trackedPartProp[0]['y'][:,0]+yllcenter-rowsMinPlot, color=cmap(0.25), zorder=0, linewidth=1.0, label='tracked particles')
        ax1.plot(trackedPartProp[0]['x']+xllcenter-colsMinPlot,trackedPartProp[0]['y']+yllcenter-rowsMinPlot, color=cmap(0.25), zorder=20, linewidth=1.0)
        # labels and ticks 
        ax1.set_xlabel('x [m] \n\n', fontsize=22)
        ax1.set_ylabel('y [m]', fontsize=22) 
        ax1.tick_params(axis='both', labelsize=13)
        fig.legend(loc='lower center', ncol=4, fancybox=True, shadow=True, fontsize=13)
        
        if AvaNodes: 
            circle = plt.Circle((79410.1-colsMinPlot,241997-rowsMinPlot), 7, color='r', fill=False, linewidth=5.0)
            #79410.1|241997 
            ax1.add_patch(circle)
        circle = plt.Circle((79410.1-colsMinPlot,241997-rowsMinPlot), 7, color='r', fill=False, linewidth=5.0)
        #79410.1|241997 
        ax1.add_patch(circle)
        #ax[0].set_xlim(200,300)
        #ax[0].set_ylim(500,600)
        
        
        # ax2
        # plot the travel length 
        cols = len(trackedPartPropAdapted[0][1]['s']) 
        rows = len(dictVelEnvelope['Time'][0]) 
        partTravelLength = np.array([[0 for i in range(cols)] for j in range(rows)])
        for i in range(1,len(dictVelEnvelope['Time'][0])):
            partTravelLength[i,:] = trackedPartPropAdapted[0][i]['s']
        
        ax2.plot(trackedPartProp[0]['time'],partTravelLength, zorder=1, linewidth=1.0, color=cmap(0.25))
        ax2.fill_between(dictVelEnvelope['Time'][0], dictVelEnvelope['SxyzMin'][0],
                        dictVelEnvelope['SxyzMax'][0], color=cmap(0.2), alpha=0.2, zorder=0)
        # add AvaNode data
        if AvaNodes: 
            ax2.plot(np.array(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'])-dictNodes['temporalDomain']['C07']['TimeAfterBlasting'][0], dictNodes['thalwegDomain']['sC07xyzAfterBlasting']-dictNodes['thalwegDomain']['sC07xyzAfterBlasting'][0], color='brown', zorder=25, label='AvaNode C07')
            ax2.plot(np.array(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'])-dictNodes['temporalDomain']['C09']['TimeAfterBlasting'][0], dictNodes['thalwegDomain']['sC09xyzAfterBlasting']-dictNodes['thalwegDomain']['sC09xyzAfterBlasting'][0], color='green', zorder=25, label='AvaNode C09')
            ax2.plot(np.array(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'])-dictNodes['temporalDomain']['C10']['TimeAfterBlasting'][0], dictNodes['thalwegDomain']['sC10xyzAfterBlasting']-dictNodes['thalwegDomain']['sC10xyzAfterBlasting'][0], color='orange', zorder=25, label='AvaNode C10')
        # labels and ticks 
        ax2.set_xlabel('Time [s] \n\n', fontsize=15)
        ax2.set_ylabel('$S_{xyz}$ particles [m]', fontsize=15) 
        ax2.tick_params(axis='both', labelsize=13)
        
        
        # ax3
        # plot velocities
        cmap = cm.vik  # colormap
        VelMagnitudePart = np.sqrt(trackedPartProp[0]['ux']**2 + trackedPartProp[0]['uy']**2 + trackedPartProp[0]['uz']**2)
        ax3.plot(trackedPartProp[0]['time'],VelMagnitudePart, zorder=1, linewidth=1.0, color=cmap(0.25))
        ax3.fill_between(dictVelEnvelope['Time'][0], dictVelEnvelope['Min'][0],
                        dictVelEnvelope['Max'][0], color=cmap(0.2), alpha=0.2, zorder=0)
        # labels and ticks 
        ax3.set_xlabel('Time [s] \n\n', fontsize=15)
        ax3.set_ylabel('Velocity [m/s]', fontsize=15) 
        ax3.tick_params(axis='both', labelsize=13)
        
        # Experiment (AvaRange)
        if AvaNodes: 
            ax3.plot(dictNodes['temporalDomain']['C07']['TimeAfterBlasting'], dictNodes['temporalDomain']['C07']
                          ['VelocityAfterBlasting'], color='brown', label='AvaNode C07', zorder=10)  
            ax3.plot(dictNodes['temporalDomain']['C09']['TimeAfterBlasting'], dictNodes['temporalDomain']['C09']
                          ['VelocityAfterBlasting'], color='green', label='AvaNode C09', zorder=10)  
            ax3.plot(dictNodes['temporalDomain']['C10']['TimeAfterBlasting'], dictNodes['temporalDomain']['C10']
                          ['VelocityAfterBlasting'], color='orange', label='AvaNode C10', zorder=10)
        # add legend
        #fig.legend(loc='lower center', ncol=4, fancybox=True, shadow=True, fontsize=15)
        
        # ax4
        # plot acceleration 
        rowsPart = np.shape(VelMagnitudePart)[0]
        colsPart = np.shape(VelMagnitudePart)[1]
        rows = np.shape(dictVelEnvelope['Velocity'])[1]
        cols = np.shape(dictVelEnvelope['Velocity'])[2]
        AccMagnitudePart = np.array([[0 for i in range(colsPart)] for j in range(rowsPart)])
        simAcc = np.array([[0 for i in range(cols)] for j in range(rows)])
        #simMinAcc = np.array([0 for i in range(len(simTime[0]))]) 
        #AccMagnitudePart = np.array(AccMagnitudePart)
        for i in range(len(VelMagnitudePart)-1):
            AccMagnitudePart[i+1] = (VelMagnitudePart[i+1] - VelMagnitudePart[i]) / (trackedPartProp[0]['time'][i+1] - trackedPartProp[0]['time'][i])
        for i in range(len(dictVelEnvelope['Time'][0])-1):
            simAcc[i+1] = (dictVelEnvelope['Velocity'][0][i+1] - dictVelEnvelope['Velocity'][0][i]) / (dictVelEnvelope['Time'][0][i+1] - dictVelEnvelope['Time'][0][i])
        #accelaration = [0]*np.shape(VelMagnitudePart)
        #acceleration[j] = AccMagnitudePart[i+1]
        ax4.plot(trackedPartProp[0]['time'], AccMagnitudePart, zorder=1, linewidth=1.0, color=cmap(0.25))
        ax4.fill_between(dictVelEnvelope['Time'][0], np.min(simAcc, axis=1),
                np.max(simAcc, axis=1), color=cmap(0.2), alpha=0.2, zorder=0)
        # labels and ticks 
        ax4.set_xlabel('Time [s] \n\n', fontsize=15)
        ax4.set_ylabel('Acceleration [m/s²]', fontsize=15) 
        ax4.tick_params(axis='both', labelsize=13)
        
        
        # ax5
        # plot travel length along thalweg 
        sPlot = [0]*len(avaDict[0])
        sAimecPlot = [0]*len(avaDict[0])
        for l in range(len(avaDict[0])):
            sPlot[l] = avaDict[0][l]['s'][:]
            sAimecPlot[l] = avaDict[0][l]['sAimec'][:]
        #ax5.plot(sPlot,sAimecPlot, zorder=1, linewidth=1.0, color=cmap(0.25))
        ax5.fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minSxyz'][simu_number],
                         dictVelAltThalweg['maxSxyz'][simu_number], color=cmap(0.2), alpha=0.2, zorder=0) 
        ax5.plot(dictVelAltThalwegPart['sXYThalweg'][simu_number], dictVelAltThalwegPart['maxSxyz'][simu_number], color=cmap(0.01))  # max
        ax5.plot(dictVelAltThalwegPart['sXYThalweg'][simu_number], dictVelAltThalwegPart['minSxyz'][simu_number], color=cmap(0.01))  # min
        # add AvaNode data
        if AvaNodes: 
            ax5.plot(dictNodes['thalwegDomain']['sC07'], dictNodes['thalwegDomain']['sC07xyz'], color='brown', zorder=25, label='AvaNode C07')
            ax5.plot(dictNodes['thalwegDomain']['sC09'], dictNodes['thalwegDomain']['sC09xyz'], color='green', zorder=25, label='AvaNode C09')
            ax5.plot(dictNodes['thalwegDomain']['sC10'], dictNodes['thalwegDomain']['sC10xyz'], color='orange', zorder=25, label='AvaNode C10')
        # labels and ticks 
        ax5.set_xlabel('$S_{xy}$ (thalweg) [m]', fontsize=15)
        ax5.set_ylabel('$S_{xyz}$ particles [m]', fontsize=15) 
        ax5.tick_params(axis='both', labelsize=13)
        
        
        # ax6
        # plot velocity along thalweg 
        ax6.fill_between(dictVelAltThalweg['sXYThalweg'][simu_number], dictVelAltThalweg['minVelocity'][simu_number],
                         dictVelAltThalweg['maxVelocity'][simu_number], color=cmap(0.2), alpha=0.2, zorder=0) 
        ax6.plot(dictVelAltThalwegPart['sXYThalweg'][simu_number], dictVelAltThalwegPart['maxVelocity'][simu_number], color=cmap(0.01))  # max
        ax6.plot(dictVelAltThalwegPart['sXYThalweg'][simu_number], dictVelAltThalwegPart['minVelocity'][simu_number], color=cmap(0.01))  # min
        # labels and ticks 
        ax6.set_xlabel('$S_{xy}$ [m]', fontsize=15)
        ax6.set_ylabel('Velocity [m/s]', fontsize=15) 
        ax6.tick_params(axis='both', labelsize=13)
        # AvaNodes data
        if AvaNodes: 
            ax6.plot(dictNodes['thalwegDomain']['sC07'], dictNodes['temporalDomain']['C07']['Velocity'], color='brown', label='AvaNode C07')
            ax6.plot(dictNodes['thalwegDomain']['sC09'], dictNodes['temporalDomain']['C09']['Velocity'], color='green', label='AvaNode C09')
            ax6.plot(dictNodes['thalwegDomain']['sC10'], dictNodes['temporalDomain']['C10']['Velocity'], color='orange', label='AvaNode C10')
             
        
        # ax7 
        # plot acceleration along the thalweg 
        # labels and ticks 
        ax7.set_xlabel('$S_{xy}$ [m]', fontsize=15)
        ax7.set_ylabel('Acceleration [m/s²]', fontsize=15)  
        ax7.tick_params(axis='both', labelsize=13)
        
        
        # title 
        PlotTools.titleFrictParam(Sim, simu_number)
     
        # saving the plot
        if Save:
            outDir = avaDir / 'Outputs' / modName / 'Plots'/'Peak flow quantities and Velocity of tracked particles'
            fU.makeADir(outDir)
            name = 'PfVelTrackedParticles_'+'_'+avaDict[simu_number][0]['simName']
            plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
            fig.savefig(plotName)
            plotPath = pathlib.Path.cwd() / plotName
            print("Plot for "+str(name)+" successfully saved at"+str(plotPath))
        
        if Show:
            plt.show()
     

            
            
    
    
    
    
    
    
    
    
    
    
    
    




