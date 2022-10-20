# -*- coding: utf-8 -*-
"""
Created on Oct Mon 12 2022

@author: dicko

Plot velocity envelopes for different friction model parameters
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

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.com1DFA import particleTools 
import runFindAvalancheInfo


# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Load configuration info of all com1DFA simulations
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# Extracting the different avalanche simulations of the output file 
F = runFindAvalancheInfo.postProcess(avalancheDir, cfgMain, simDF)

# Creating a dictionary containing information on each avalanche 
L = [None]*len(F.index)
for i in range(0, len(F.index)):
    print('Avalanche found in the output file: '+F.index[i])
    L[i], _ = particleTools.readPartFromPickle(
        avalancheDir, simName=F.index[i], flagAvaDir=True, comModule='com1DFA')

#%%Calculating the velocity envelope

number_files = len(F.index)

# Preparing time, maximum values, minimum values, mean, median
Time = [None]*number_files
Max = [None]*number_files
Min = [None]*number_files
Mean = [None]*number_files
Median = [None]*number_files


for i in range(0, number_files):
    number_time_steps = len(L[i])
    # Reading and calculating the velocity magnitude for avalanche i  
    df = [None]*number_time_steps
    for k in range(0, number_time_steps):
        df[k] = L[i][k]['ux']**2 + L[i][k]['uy']**2 + L[i][k]['uz']**2
        for l in range (0,len(df[k])):
            df[k][l] = np.sqrt(df[k][l])

    # Preparing time, maximum values, minimum values, mean, median for avalanche i 
    Time[i] = [None]*number_time_steps
    Max[i] = [None]*number_time_steps
    Min[i] = [None]*number_time_steps
    Mean[i] = [None]*number_time_steps
    Median[i] = [None]*number_time_steps

    for j in range(0, number_time_steps):
        Time[i][j] = L[i][j]['t']
        Max[i][j] = max(df[j])
        Min[i][j] = min(df[j])
        Mean[i][j] = stat.mean(df[j])
        Median[i][j] = stat.median(df[j])
        
#%% Experiment 

# Reading the experiment files 
file = avalancheDir + "/220222_C10_avalanche_GPS.txt"
isExist = os.path.exists(file)
if isExist== False:
    print("No experiment data found: the AvaNode data should be in "+avalancheDir+"/220222_C10_avalanche_GPS.txt")
    key = 0 
else:      
    Time_experiment = [x.split(',')[0] for x in open(file).readlines()]
    velN = [x.split(',')[7] for x in open(file).readlines()]
    velE = [x.split(',')[8]for x in open(file).readlines()]
    velD = [x.split(',')[9] for x in open(file).readlines()]
    
    # Calculating the velocity norm
    vel = [None]*(len(velN)-1)
    time = [None]*(len(velN)-1)
    for i in range(1, len(velN)):
        vel[i-1] = np.sqrt(int(velN[i])**2 + int(velE[i])**2 + int(velD[i])**2) * 10**-3
        time[i-1] = int(Time_experiment[i])*10**-6
    
    # Deleting the null velocity before the blasting  
    threshold = 1             # Warning: arbitrary 
    index = vel.index(next((i for i in vel if int(i)>threshold),None))
    vel_experimental= vel[index:]
    time_experimental = [i-time[index] for i in time[index:]] 
    
    # Key to plot only if the experimental data have been found 
    key = 1 
    
    
    
    
#%% Plotting 

# Preparing the subplots regarding the number of avalanches found 
if key == 1:  
    if number_files == 1:
        fig, ax = plt.subplots(1,1,figsize=(20, 20),sharex=True) 
        nrow = 1
        ncol = 1 
    elif number_files == 2:
        fig, ax = plt.subplots(2,1,figsize=(20, 20),sharex=True,sharey=True)
        nrow = 2
        ncol = 1 
    elif number_files == 3: 
        fig, ax = plt.subplots(2,2,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 2
        ncol = 2 
    elif number_files == 4: 
        fig, ax = plt.subplots(2,2,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 2
        ncol = 2 
    elif number_files == 5: 
        fig, ax = plt.subplots(2,3,figsize=(20, 20),sharex=True,sharey=True)
        nrow = 2
        ncol = 3 
    elif number_files == 6: 
        fig, ax = plt.subplots(2,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 2
        ncol = 3
    elif number_files == 7: 
        fig, ax = plt.subplots(3,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 3
        ncol = 3 
    elif number_files == 8: 
        fig, ax = plt.subplots(2,4,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 2
        ncol = 4
    elif number_files == 9: 
        fig, ax = plt.subplots(3,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 3
        ncol = 3 
    else:
        print('Too many simulations for one plot!')
        
    max_time = max(max(Time[:]))
    simu_number = 0 
    friction_model = F.frictModel[0]   
    stopping_model = F.stopCrit[0] 
    
    ## To make comparisons between implicit and explicit methods 
    # friction_formulation = [None]*len(F.explicitFriction)
    # for i in range (0,len(F.explicitFriction)):
    #     if F.explicitFriction[i]==0:
    #         friction_formulation[i] = 'implict'
    #     else:
    #         friction_formulation[i] = 'explicit'
    
    # For 1D plots 
    if nrow == 1 or ncol == 1:
        for i in range(0,max(nrow,ncol)):
            ax[i].plot(Time[simu_number], Max[simu_number], color='blue')  # max
            ax[i].plot(Time[simu_number], Min[simu_number], color='blue')  # min
            ax[i].plot(Time[simu_number], Mean[simu_number], color='black',
                     linestyle='dashed', markersize=1.8)  # mean
            ax[i].plot(Time[simu_number], Median[simu_number], linestyle='dotted',
                     color='black', markersize=1.8)  # median
            # filling the space between the max and min
            ax[i].fill_between(Time[simu_number], Min[simu_number], Max[simu_number], color='blue', alpha=0.2)
            # Experiment (AvaRange)
            ax[i].plot(time_experimental, vel_experimental, color = 'red')  # AvaNode velocity
            # Plotting parameters
            if friction_model=='Coulomb':
                ax[i].set_title("mu ="+str(F.mu[simu_number]), fontsize=18)
            elif friction_model=='Voellmy':
                ax[i].set_title("mu ="+str(F.mu[simu_number])+", xsi="+str(F.xsi[simu_number]), fontsize=18)
            elif friction_model=='samosAT':
                ax[i].set_title("mu ="+str(F.mu[simu_number])+", tau0="+str(F.tau0[simu_number]), fontsize=18)
            else:
                ax[i].set_title("No friction model found!", fontsize=18)
            ax[i].set_xlim(right=max_time)
            simu_number += 1 
            
    # For 2D plots         
    else: 
        for i in range (0,nrow):
            for j in range (0,ncol):
                ax[i,j].plot(Time[simu_number], Max[simu_number], color='blue')  # max
                ax[i,j].plot(Time[simu_number], Min[simu_number], color='blue')  # min
                ax[i,j].plot(Time[simu_number], Mean[simu_number], color='black',
                         linestyle='dashed', markersize=1.8)  # mean
                ax[i,j].plot(Time[simu_number], Median[simu_number], linestyle='dotted',
                         color='black', markersize=1.8)  # median
                # filling the space between the max and min
                ax[i,j].fill_between(Time[simu_number], Min[simu_number], Max[simu_number], color='blue', alpha=0.2)
                # Experiment (AvaRange)
                ax[i,j].plot(time_experimental, vel_experimental, color = 'red')  # AvaNode velocity
                # Plotting parameters
                if friction_model=='Coulomb':
                    ax[i,j].set_title("mu ="+str(F.mu[simu_number]), fontsize=18)
                elif friction_model=='Voellmy':
                    ax[i,j].set_title("mu ="+str(F.mu[simu_number])+", xsi="+str(F.xsi[simu_number]), fontsize=18)
                elif friction_model=='samosAT':
                    ax[i,j].set_title("mu ="+str(F.mu[simu_number])+", tau0="+str(F.tau0[simu_number]), fontsize=18)
                else:
                    ax[i,j].set_title("No friction model found!", fontsize=18)
                ax[i,j].set_xlim(right=max_time)
                simu_number += 1 
        # If the number of subplots does not match with the number of simulations, delete the last subplot
        if simu_number!=number_files:
            fig.delaxes(ax[nrow,ncol])
            
    # set legend position 
    fig.legend(['Maximum values','Minimum values','Mean','Median','Velocity envelope','AvaNode Velocity'],loc='lower center', ncol=6, fancybox = True, shadow=False)
    fig.suptitle('Velocity envelope - Seilbahnrinne,\n '+friction_model+' model with stopCrit='+str(stopping_model))

    for axs in ax.flat:
        axs.set_ylabel("Velocity[m/s]", fontsize=15)
        axs.set_xlabel("Time[s]\n\n", fontsize=15)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for axs in ax.flat:
        axs.label_outer()
    
    # set spacing to subplots
    fig.tight_layout()
    
    plt.show() 


    
   