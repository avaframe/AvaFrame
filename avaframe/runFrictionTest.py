# -*- coding: utf-8 -*-
"""
Created on Oct Mon 12 2022

@author: dicko

Plot velocity envelopes and several comparison tools for different friction model parameters
Warning : can only be used to compare simulations run with the same friction 
model  
If the friction model used for the simulations is Coulomb, compares plots for different mu values.
If the friction model is Voellmy, compares plots for different mu ans xsi values.
If the friction model is SamosAT, compares plots for different mu and tau0. 

"""
# Python imports 
import numpy as np
from matplotlib import pyplot as plt
import statistics as stat
import os 
from cmcrameri import cm
from scipy.interpolate import interp1d 
import pathlib

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.com1DFA import particleTools 
from avaframe import runFindAvalancheInfo
from avaframe import NodeTools 


# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Load configuration info of all com1DFA simulations
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# # Extracting the different avalanche simulations of the output file 
Sim = runFindAvalancheInfo.postProcess(avalancheDir, cfgMain, simDF)

# Creating a dictionary containing information on each avalanche 
for i in range(0, len(Sim.index)):
    print('Avalanche found in the output file: '+Sim.index[i])
    dict[i], _ = particleTools.readPartFromPickle(
        avalancheDir, simName=Sim.index[i], flagAvaDir=True, comModule='com1DFA')
    
# Some general and useful parameters 
number_ava = len(Sim.index)
friction_model = Sim.frictModel[0]
stopping_model = Sim.stopCrit[0] 
avaDir = pathlib.Path(avalancheDir)
modName = 'com1DFA'

#%%Calculating the velocity envelope (time, max, min, mean, median) 

# Preparing time, maximum values, minimum values, mean, median
Time = [None]*number_ava
Max = [None]*number_ava
Min = [None]*number_ava
Mean = [None]*number_ava
Median = [None]*number_ava

for i in range(0, number_ava):
    number_time_steps = len(dict[i])
    # Reading and calculating the velocity magnitude for avalanche i  
    df = [None]*number_time_steps
    for k in range(0, number_time_steps):
        df[k] = np.sqrt(np.array(dict[i][k]['ux']**2 + dict[i][k]['uy']**2 + dict[i][k]['uz']**2)) 

    # Preparing time, maximum values, minimum values, mean, median for avalanche i 
    Time[i] = [None]*number_time_steps
    Max[i] = [None]*number_time_steps
    Min[i] = [None]*number_time_steps
    Mean[i] = [None]*number_time_steps
    Median[i] = [None]*number_time_steps

    for j in range(0, number_time_steps):
        Time[i][j] = dict[i][j]['t']
        Max[i][j] = max(df[j])
        Min[i][j] = min(df[j])
        Mean[i][j] = stat.mean(df[j])
        Median[i][j] = stat.median(df[j]) 
        
#%% Dissipation energy 

# Kinetic = [None]*number_files
# Potential = [None]*number_files
# Mechanic = [None]*number_files
# Dissipation = [None]*number_files
# Mass = [None]*number_files


# for i in range(0, number_files):
#     number_time_steps = len(L[i])
#     # Reading and calculating the velocity magnitude for avalanche i  
#     Kinetic[i] = [None]*number_time_steps
#     Potential[i] = [None]*number_time_steps
#     Mechanic[i] = [None]*number_time_steps
#     Dissipation[i] = [None]*number_time_steps
#     Mass[i] = [None]*number_time_steps
#     for k in range(0, number_time_steps):
#         Kinetic[i][k] = L[i][k]['kineticEne']
#         Potential[i][k] = L[i][k]['potentialEne']
#         Mass[i][k] = L[i][k]['m']
#         Mechanic[i][k] = Kinetic[i][k] + Potential[i][k] 

# Potential[0] = np.array(Potential[0]) - Potential[0][-1]
# Dissipation[0] = Potential[0][0] - np.array(Potential[0]+Kinetic[0]) 

# Mean_mass = np.mean(Mass[0][:],1)
# # Calculating the energy per unit mass 
# Kinetic[0] = np.array(Kinetic[0])/Mean_mass
# Potential[0] = np.array(Potential[0])/Mean_mass
# Mechanic[0] = np.array(Mechanic[0])/Mean_mass
# Dissipation[0] = np.array(Dissipation[0])/Mean_mass


# fig = plt.figure() 
# ax = fig.add_subplot(111) 
# plt.axhline(y=Potential[0][0], linestyle='dashed', color='black') 
# plt.plot(Time[0],Kinetic[0], label='Kin')
# plt.plot(Time[0],Potential[0], label='Pot')
# #plt.plot(Time[0],Mechanic[0], label='Mech')
# plt.plot(Time[0],Dissipation[0], label='Diss')
# ax.set_title("Avaframe, energy balance", fontsize=20)
# ax.set_ylabel("Energy per unit mass [m²/s²]", fontsize=15)
# ax.set_xlabel("Time [s]\n\n", fontsize=15)
# plt.legend()
# plt.show()


#%% Experiment 
    
# Calculating the Nodes velocity
NumberNodes = [7,9,10]  # Numbers of the nodes you want to plot 
Nodes = NodeTools.FindNodeVelocity(NumberNodes,cfgMain,avalancheDir)
if len(Nodes)==len(NumberNodes):
    key=1
else:
    key=0
    print("Not all experimental data have been found!")
    

#%% Boxplot 

# Node to make the comparisons 
Node = 'C10'

# Calculating the difference between the mean avaframe velocity and the sensor velocity
boxplot_data = [None]*number_ava
for i in range(0, number_ava):
    f = interp1d(Nodes[Node]['Time'] , Nodes[Node]['Velocity']) # interpolating the node velocity 
    if max(Time[i]) > max(Nodes[Node]['Time'] ):
        index = Time[i].index(next((j for j in Time[i] if int(j)>max(Nodes[Node]['Time'] )),None))
        vel_experimental_interpolated = f(Time[i][0:index-1]) 
        boxplot_data[i] = Mean[i][0:index-1] - vel_experimental_interpolated
    else:
        vel_experimental_interpolated = f(Time[i]) 
        boxplot_data[i] = Mean[i] - vel_experimental_interpolated

# Making the boxplot plot 
fig = plt.figure() 
ax = fig.add_subplot(111)
plt.boxplot(boxplot_data, whis=200)
# x-axis labels
if friction_model=='Coulomb':
    labels = ["mu ="+str(Sim.mu[i]) for i in range(0,number_ava)]
elif friction_model=='Voellmy':
    labels = ["mu ="+str(Sim.mu[i])+"\nxsi="+str(Sim.xsi[i]) for i in range(0,number_ava)]
elif friction_model=='samosAT':
    labels = ["mu ="+str(Sim.mu[i])+"\ntau0="+str(Sim.tau0[i]) for i in range(0,number_ava)]
elif friction_model=='VoellmyUpgraded':
    labels = ["mu ="+str(Sim.mu[i])+"\ntau0="+str(Sim.tau0[i])+"\nxsi="+str(Sim.xsi[i]) for i in range(0,number_ava)]
else:
    labels = ["No friction model found!" for i in range(0,number_ava)]

fig.suptitle('Velocity difference between the mean avaframe velocity and the Node velocity\nSeilbahnrinne, explicit '+friction_model+' model with stopCrit='+str(stopping_model), fontsize= 20)
ax.set_xticklabels(labels, color='black', fontsize=15)
ax.tick_params(axis='y', colors='black', labelsize=15) 
ax.set_ylabel("Velocity difference[m/s]", color='black', fontsize=15)
plt.show() 
    
#%% Plotting velocity envelope 
import avaframe.out3Plot.plotUtils as pU
# Preparing the subplots regarding the number of avalanches found 
if key == 1:  
    if number_ava == 1:
        fig, ax = plt.subplots(1,1,figsize=(20, 20),sharex=True) 
        nrow = 1
        ncol = 1 
    elif number_ava == 2:
        fig, ax = plt.subplots(2,1,figsize=(pU.figW+10, pU.figH+3),sharex=True,sharey=True)
        nrow = 2
        ncol = 1 
    elif number_ava == 3: 
        fig, ax = plt.subplots(1,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 1
        ncol = 3
    elif number_ava == 4: 
        fig, ax = plt.subplots(2,2,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 2
        ncol = 2 
    elif number_ava == 5: 
        fig, ax = plt.subplots(2,3,figsize=(20, 20),sharex=True,sharey=True)
        nrow = 2
        ncol = 3 
    elif number_ava == 6: 
        fig, ax = plt.subplots(2,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 2
        ncol = 3
    elif number_ava == 7: 
        fig, ax = plt.subplots(3,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 3
        ncol = 3 
    elif number_ava == 8: 
        fig, ax = plt.subplots(2,4,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 2
        ncol = 4
    elif number_ava == 9: 
        fig, ax = plt.subplots(3,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 3
        ncol = 3 
    else:
        print('Too many simulations for one plot!')
        
    max_time = max(max(Time[:]))
    simu_number = 0 
    
    ## To make comparisons between implicit and explicit methods 
    # friction_formulation = [None]*len(F.explicitFriction)
    # for i in range (0,len(F.explicitFriction)):
    #     if F.explicitFriction[i]==0:
    #         friction_formulation[i] = 'implict'
    #     else:
    #         friction_formulation[i] = 'explicit'
    
    # Colormap 
    cmap = cm.vik
    rgba = cmap(0.2) 
        
    # For 1D plots 
    if nrow == 1 and ncol == 1:
        ax.plot(Time[simu_number], Max[simu_number], color=cmap(0.01))  # max
        ax.plot(Time[simu_number], Min[simu_number], color=cmap(0.01))  # min
        ax.plot(Time[simu_number], Mean[simu_number], linestyle='dashed', color=cmap(0.7), markersize=1.8)  # mean
        ax.plot(Time[simu_number], Median[simu_number], linestyle='dotted', color=cmap(0.7), markersize=1.8)  # median
        # filling the space between the max and min
        ax.fill_between(Time[simu_number], Min[simu_number], Max[simu_number], color=cmap(0.2), alpha=0.2)
        # Experiment (AvaRange)
        ax.plot(Nodes['C07']['Time'], Nodes['C07']['Velocity'], color=cmap(0.80))  # AvaNode velocity
        ax.plot(Nodes['C09']['Time'], Nodes['C09']['Velocity'], color=cmap(0.90))  # AvaNode velocity
        ax.plot(Nodes['C10']['Time'], Nodes['C10']['Velocity'], color=cmap(0.95))  # AvaNode velocity
        # Plotting parameters
        if friction_model=='Coulomb':
            ax.set_title("mu ="+str(Sim.mu[simu_number]), fontsize=18)
        elif friction_model=='Voellmy':
            ax.set_title("mu ="+str(Sim.mu[simu_number])+", xsi="+str(Sim.xsi[simu_number]), fontsize=18)
        elif friction_model=='samosAT':
            ax.set_title("mu ="+str(Sim.mu[simu_number])+", tau0="+str(Sim.tau0[simu_number]), fontsize=18)
        elif friction_model=='VoellmyUpgraded':
            ax.set_title("mu ="+str(Sim.mu[simu_number])+", tau0="+str(Sim.tau0[simu_number])+", xsi="+str(Sim.xsi[i]), fontsize=15)
        else:
            ax.set_title("No friction model found!", fontsize=18)
        #ax[i].set_xlim(right=max_time) 
        
    elif nrow == 1 or ncol == 1:
        for i in range(0,max(nrow,ncol)):
            ax[i].plot(Time[simu_number], Max[simu_number], color=cmap(0.01))  # max
            ax[i].plot(Time[simu_number], Min[simu_number], color=cmap(0.01))  # min
            ax[i].plot(Time[simu_number], Mean[simu_number], linestyle='dashed', color=cmap(0.7), markersize=1.8)  # mean
            ax[i].plot(Time[simu_number], Median[simu_number], linestyle='dotted', color=cmap(0.7), markersize=1.8)  # median
            # filling the space between the max and min
            ax[i].fill_between(Time[simu_number], Min[simu_number], Max[simu_number], color=cmap(0.2), alpha=0.2)
            # Experiment (AvaRange)
            ax[i].plot(Nodes['C07']['Time'], Nodes['C07']['Velocity'], color=cmap(0.80))  # AvaNode velocity
            ax[i].plot(Nodes['C09']['Time'], Nodes['C09']['Velocity'], color=cmap(0.90))  # AvaNode velocity
            ax[i].plot(Nodes['C10']['Time'], Nodes['C10']['Velocity'], color=cmap(0.95))  # AvaNode velocity
            # Plotting parameters
            if friction_model=='Coulomb':
                ax[i].set_title("mu ="+str(Sim.mu[simu_number]), fontsize=18)
            elif friction_model=='Voellmy':
                ax[i].set_title("mu ="+str(Sim.mu[simu_number])+", xsi="+str(Sim.xsi[simu_number]), fontsize=18)
            elif friction_model=='samosAT':
                ax[i].set_title("mu ="+str(Sim.mu[simu_number])+", tau0="+str(Sim.tau0[simu_number]), fontsize=18)
            elif friction_model=='VoellmyUpgraded':
                ax[i].set_title("mu ="+str(Sim.mu[simu_number])+", tau0="+str(Sim.tau0[simu_number])+", xsi="+str(Sim.xsi[i]), fontsize=15)
            else:
                ax[i].set_title("No friction model found!", fontsize=18)
            #ax[i].set_xlim(right=max_time)
            simu_number += 1 
            
    # For 2D plots       
    else: 
        for i in range (0,nrow):
            for j in range (0,ncol):
                ax[i,j].plot(Time[simu_number], Max[simu_number], color=cmap(0.01))  # max
                ax[i,j].plot(Time[simu_number], Min[simu_number], color=cmap(0.01))  # min
                ax[i,j].plot(Time[simu_number], Mean[simu_number], color=cmap(0.65),linestyle='dashed', markersize=1.8)  # mean
                ax[i,j].plot(Time[simu_number], Median[simu_number], linestyle='dotted', color=cmap(0.65), markersize=1.8)  # median
                # filling the space between the max and min
                ax[i,j].fill_between(Time[simu_number], Min[simu_number], Max[simu_number], color=cmap(0.2), alpha=0.2)
                # Experiment (AvaRange)
                ax[i,j].plot(Nodes['C07']['Time'], Nodes['C07']['Velocity'], color=cmap(0.80))  # AvaNode velocity
                ax[i,j].plot(Nodes['C09']['Time'], Nodes['C09']['Velocity'], color=cmap(0.90))  # AvaNode velocity
                ax[i,j].plot(Nodes['C10']['Time'], Nodes['C10']['Velocity'], color=cmap(0.99))  # AvaNode velocity   
                # Plotting parameters
                if friction_model=='Coulomb':
                    ax[i,j].set_title("mu ="+str(Sim.mu[simu_number]), fontsize=18)
                elif friction_model=='Voellmy':
                    ax[i,j].set_title("mu ="+str(Sim.mu[simu_number])+", xsi="+str(Sim.xsi[simu_number]), fontsize=18)
                elif friction_model=='samosAT':
                    ax[i,j].set_title("mu ="+str(Sim.mu[simu_number])+", tau0="+str(Sim.tau0[simu_number]), fontsize=18)
                elif friction_model=='VoellmyUpgraded':
                    ax[i,j].set_title("mu ="+str(Sim.mu[simu_number])+", tau0="+str(Sim.tau0[simu_number])+", xsi="+str(Sim.xsi[i]), fontsize=15)
                else:
                    ax[i,j].set_title("No friction model found!", fontsize=18)
                #ax[i,j].set_xlim(right=max_time)
                simu_number += 1 
            
    # set legend position 
    fig.legend(['Maximum values','Minimum values','Mean','Median','Velocity envelope','AvaNode C07 Velocity','AvaNode C09 Velocity','AvaNode C10 Velocity'],loc='lower center', ncol=5, fancybox = True, shadow=False)
    fig.suptitle('Velocity envelope - Seilbahnrinne,\n explicit '+friction_model+' model with stopCrit='+str(stopping_model))

    for axs in ax.flat:
        axs.set_ylabel("Velocity[m/s]", fontsize=15)
        axs.set_xlabel("Time[s]\n\n", fontsize=15)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for axs in ax.flat:
        axs.label_outer()
    
    # set spacing to subplots
    fig.tight_layout()
    
    # saving the plot
    outDir = avaDir / 'Outputs' / modName /'ComparisonPlots'/'VelocityEnvelopes'
    fU.makeADir(outDir)
    plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
    fig.savefig(plotName)
    plotPath = pathlib.Path.cwd() / plotName
    
    plt.show() 


#%% function to plot a range-time diagram

#%% function to plot peak flow velocity 

import os
import numpy as np
from matplotlib import pyplot as plt
import pathlib 

import avaframe.out3Plot.plotUtils as pU
import avaframe.in1Data.getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf

from avaframe.data.avaSeilbahn.AvaNode_data import gps_imu_tools as git
from avaframe.data.avaSeilbahn.AvaNode_data import GPS_Class 

def PeakFields(avaDir, cfgFLAGS, modName, n, demData=''):
    """ Plot all peak fields and return dictionary with paths to plots
        with DEM in background

        Parameters
        ----------
        avaDir : str
            path to avalanche directoy
        cfgFLAGS : str
            general configuration, required to define if plots saved to reports directoy
        modName : str
            name of module that has been used to produce data to be plotted
        demData: dictionary
            optional - if not the dem in the avaDir/Inputs folder has been used but a different one

        Returns
        -------
        plotDict : dict
            dictionary with info on plots, like path to plot
        """

    # Load all infos on simulations
    avaDir = pathlib.Path(avaDir)
    inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    if demData == '':
        demFile = gI.getDEMPath(avaDir)
        demData = IOf.readRaster(demFile, noDataToNan=True)
        demDataField = demData['rasterData']
    else:
        # check if noDataValue is found and if replace with nans for plotting
        demDataField = np.where(demData['rasterData'] == demData['header']['noDataValue'], np.nan, demData['rasterData'])
    demField = demDataField


    # Initialise plot dictionary with simulation names
    # plotDict = {}
    # for sName in peakFilesDF['simName']:
    #     plotDict[sName] = {}

    #print(peakFilesDF['simName']) 
    # Loop through peakFiles 
    #for m in range(len(peakFilesDF['names'])):

    # Load names and paths of peakFiles
    name = peakFilesDF['names'][n]
    fileName = peakFilesDF['files'][n] 
    avaName = peakFilesDF['avaName'][n]
    resType = peakFilesDF['resType'][n]

    # Load data
    raster = IOf.readRaster(fileName, noDataToNan=True)
    data = raster['rasterData']

    # constrain data to where there is data
    cellSize = peakFilesDF['cellSize'][simu_number]
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(data, cellSize)
    dataConstrained = data[rowsMin:rowsMax+1, colsMin:colsMax+1]
    demConstrained = demField[rowsMin:rowsMax+1, colsMin:colsMax+1]

    data = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)
    unit = pU.cfgPlotUtils['unit%s' % resType]

    # Set extent of peak file
    ny = data.shape[0]
    nx = data.shape[1]
    Ly = ny*cellSize
    Lx = nx*cellSize
    

    return resType, data, raster, rowsMin, rowsMax, colsMin, colsMax, cellSize, demConstrained, unit, name, peakFilesDF

def change_coordinate():
    # Load information on AvAnodes data and changing coordinate system
    gps_c10 = GPS_Class.GPSData() 
    gps_c09 = GPS_Class.GPSData() 
    gps_c07 = GPS_Class.GPSData() 
    gps_c10.path = "/home/dick/Documents/AvaFrame/avaframe/data/avaSeilbahn/AvaNode_data/220222_C10_avalanche_GPS.txt"
    gps_c09.path = "/home/dick/Documents/AvaFrame/avaframe/data/avaSeilbahn/AvaNode_data/220222_C09_avalanche_GPS.txt"
    gps_c07.path = "/home/dick/Documents/AvaFrame/avaframe/data/avaSeilbahn/AvaNode_data/220222_C07_avalanche_GPS.txt"
    gps_c10.read_data() 
    gps_c09.read_data() 
    gps_c07.read_data() 
    n10,e10,z10 = git.gps_to_mercator(gps_c10) 
    n09,e09,z09 = git.gps_to_mercator(gps_c09)
    n07,e07,z07 = git.gps_to_mercator(gps_c07)
    return e07, n07, e09, n09, e10, n10 

#%% Plotting different tools together

friction_model = F.frictModel[0]   
stopping_model = F.stopCrit[0] 

for simu_number in range(0,number_files): 
    
    # Load all infos on simulations
    avaDir = pathlib.Path(avalancheDir)
    modName = 'com1DFA'
    inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    sim = F.index[simu_number] 
    simu_numb = np.where(peakFilesDF.simID== sim)
        
    for k in range (0,3):
        fig, ax = plt.subplots(1,3,figsize=(pU.figW+10, pU.figH+3))  
        # ax[0]
        
        # ax[1]
        # Generate plots for all peakFiles
        n = simu_numb[0][2+k]

        resType, data, raster, rowsMin, rowsMax, colsMin, colsMax, cellSize, demConstrained, unit, name, peakFilesDF = PeakFields(avalancheDir, cfgMain['FLAGS'], modName, n, demData='')
        e07, n07, e09, n09, e10, n10 = change_coordinate()
        
        # Figure  shows the result parameter data
        #fig, ax = plt.subplots(figsize=(pU.figW, pU.figH))
        # choose colormap
        cmap, col, ticks, norm = pU.makeColorMap(pU.colorMaps[resType], np.amin(data), np.amax(data), continuous=pU.contCmap)
        cmap.set_bad(alpha=0)
        # uncomment this to set the under value for discrete cmap transparent
        # cmap.set_under(alpha=0)
        xllcenter = raster['header']['xllcenter']
        yllcenter = raster['header']['yllcenter']
        rowsMinPlot = rowsMin*cellSize + yllcenter
        rowsMaxPlot = (rowsMax+1)*cellSize + yllcenter
        colsMinPlot = colsMin*cellSize + xllcenter
        colsMaxPlot = (colsMax+1)*cellSize + xllcenter
        
        # rowsMinPlot = rowsMin*cellSize
        # rowsMaxPlot = (rowsMax+1)*cellSize
        # colsMinPlot = colsMin*cellSize
        # colsMaxPlot = (colsMax+1)*cellSize
        extent = [colsMinPlot, colsMaxPlot, rowsMinPlot, rowsMaxPlot]
        
        # add DEM hillshade with contour lines
        ls, CS = pU.addHillShadeContours(ax[1], demConstrained, cellSize, extent)
        
        # add peak field data
        im1 = ax[1].imshow(data, cmap=cmap, norm=norm, extent=extent, origin='lower', aspect='equal', zorder=2)
        #im1 = ax[1].imshow(data, cmap=cmap, norm=norm, aspect='equal', zorder=2)
        pU.addColorBar(im1, ax[1], ticks, unit)
        
        # # add AvaNode data 
        ax[1].plot(e10,n10, color='orange', label='AvaNode C10')
        ax[1].plot(e09,n09, color='green', label='AvaNode C09')
        ax[1].plot(e07,n07, color='brown', label='AvaNode C07')
        
        # # add title, labels and ava Info
        # #title = str('%s' % name)
        # #ax.set_title(title +'\n')
        ax[1].set_xlabel('x [m] \n\n')
        # ax[1].set_ylabel('y [m]')
        
        
        
        # # ax[2] 
        # Colormap 
        cmap = cm.vik
        rgba = cmap(0.2) 
        ax[2].plot(Time[simu_number], Max[simu_number], color=cmap(0.01), label='Maximum and Minimum values')  # max
        ax[2].plot(Time[simu_number], Min[simu_number], color=cmap(0.01))  # min
        ax[2].plot(Time[simu_number], Mean[simu_number], linestyle='dashed', color=cmap(0.7), markersize=1.8, label='Mean')  # mean
        ax[2].plot(Time[simu_number], Median[simu_number], linestyle='dotted', color=cmap(0.7), markersize=1.8, label='Median')  # median
        # filling the space between the max and min
        ax[2].fill_between(Time[simu_number], Min[simu_number], Max[simu_number], color=cmap(0.2), alpha=0.2)
        # Experiment (AvaRange)
        ax[2].plot(time_experimental_C07, vel_experimental_C07, color='brown')  # AvaNode velocity
        ax[2].plot(time_experimental_C09, vel_experimental_C09, color='green')  # AvaNode velocity
        ax[2].plot(time_experimental_C10, vel_experimental_C10, color='orange')  # AvaNode velocity
        ax[2].set_ylabel("Velocity[m/s]", fontsize=15)
        ax[2].set_xlabel("Time[s]\n\n", fontsize=15)
        
        
        # Plot title 
        if friction_model=='Coulomb':
            fig.suptitle('Velocity envelope - Seilbahnrinne,\n explicit '+friction_model+' model with stopCrit='+str(stopping_model)+'\n mu ='+str(F.mu[simu_number]), fontsize=25)
        elif friction_model=='Voellmy':
            fig.suptitle('Velocity envelope - Seilbahnrinne,\n explicit '+friction_model+' model with stopCrit='+str(stopping_model)+'\n mu ='+str(F.mu[simu_number])+', xsi='+str(F.xsi[i]), fontsize=25 )
        elif friction_model=='samosAT':
            fig.suptitle('Velocity envelope - Seilbahnrinne,\n explicit '+friction_model+' model with stopCrit='+str(stopping_model)+'\n mu ='+str(F.mu[simu_number])+', tau0='+str(F.tau0[simu_number]), fontsize=25 )
        elif friction_model=='VoellmyUpgraded':
            fig.suptitle('Seilbahnrinne, explicit '+friction_model+' model with stopCrit='+str(stopping_model)+', mu ='+str(F.mu[simu_number])+', tau0='+str(F.tau0[simu_number])+', xsi='+str(F.xsi[i])+'\n', fontsize=23 )
        else:
            fig.suptitle("No friction model found!", fontsize=25)
            
            
        # Add legend 
        fig.legend(loc='lower center', ncol=3, fancybox=True, shadow=True, fontsize=10)
        
        
        # if cfgFLAGS.getboolean('showPlot'):
        #     plt.show()
        avaDir = pathlib.Path(avalancheDir)
        outDir = avaDir / 'Outputs' / modName /'ComparisonPlots'/'GraphicInterface'
        fU.makeADir(outDir)
        plotName = outDir / ('%s.%s' % (name, pU.outputFormat))
        fig.savefig(plotName)
        plotPath = pathlib.Path.cwd() / plotName
        plt.close(fig)
