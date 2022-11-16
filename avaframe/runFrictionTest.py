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
from scipy.interpolate import interp1d 

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.com1DFA import particleTools 
import runFindAvalancheInfo


# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Load configuration info of all com1DFA simulations
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# # Extracting the different avalanche simulations of the output file 
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
        
#%% Dissipation energy 

Kinetic = [None]*number_files
Potential = [None]*number_files
Mechanic = [None]*number_files
Dissipation = [None]*number_files
Mass = [None]*number_files


for i in range(0, number_files):
    number_time_steps = len(L[i])
    # Reading and calculating the velocity magnitude for avalanche i  
    Kinetic[i] = [None]*number_time_steps
    Potential[i] = [None]*number_time_steps
    Mechanic[i] = [None]*number_time_steps
    Dissipation[i] = [None]*number_time_steps
    Mass[i] = [None]*number_time_steps
    for k in range(0, number_time_steps):
        Kinetic[i][k] = L[i][k]['kineticEne']
        Potential[i][k] = L[i][k]['potentialEne']
        Mass[i][k] = L[i][k]['m']
        Mechanic[i][k] = Kinetic[i][k] + Potential[i][k] 

Potential[0] = np.array(Potential[0]) - Potential[0][-1]
Dissipation[0] = Potential[0][0] - np.array(Potential[0]+Kinetic[0]) 

Mean_mass = np.mean(Mass[0][:],1)
# Calculating the energy per unit mass 
Kinetic[0] = np.array(Kinetic[0])/Mean_mass
Potential[0] = np.array(Potential[0])/Mean_mass
Mechanic[0] = np.array(Mechanic[0])/Mean_mass
Dissipation[0] = np.array(Dissipation[0])/Mean_mass


fig = plt.figure() 
ax = fig.add_subplot(111) 
plt.axhline(y=Potential[0][0], linestyle='dashed', color='black') 
plt.plot(Time[0],Kinetic[0], label='Kin')
plt.plot(Time[0],Potential[0], label='Pot')
#plt.plot(Time[0],Mechanic[0], label='Mech')
plt.plot(Time[0],Dissipation[0], label='Diss')
ax.set_title("Avaframe, energy balance", fontsize=20)
ax.set_ylabel("Energy per unit mass [m²/s²]", fontsize=15)
ax.set_xlabel("Time [s]\n\n", fontsize=15)
plt.legend()
plt.show()


#%% Experiment 

# Reading the experiment files for AvaNode C07 
file = avalancheDir + "/AvaNode_data/220222_C07_avalanche_GPS.txt"
isExist = os.path.exists(file)
if isExist== False:
    print("No experiment data found: the AvaNode data should be in "+avalancheDir+"/AvaNode_data/220222_C10_avalanche_GPS.txt")
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
    vel_experimental_C07= vel[index:]
    time_experimental_C07 = [i-time[index] for i in time[index:]] 
    
    # Key to plot only if the experimental data have been found 
    key = 1 

# Reading the experiment files for AvaNode C09
file = avalancheDir + "/AvaNode_data/220222_C09_avalanche_GPS.txt"
isExist = os.path.exists(file)
if isExist== False:
    print("No experiment data found: the AvaNode data should be in "+avalancheDir+"/AvaNode_data/220222_C09_avalanche_GPS.txt")
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
    vel_experimental_C09 = vel[index:]
    time_experimental_C09 = [i-time[index] for i in time[index:]] 
    
    # Key to plot only if the experimental data have been found 
    key = 1 
    
    
# Reading the experiment files for AvaNode C10  
file = avalancheDir + "/AvaNode_data/220222_C10_avalanche_GPS.txt"
isExist = os.path.exists(file)
if isExist== False:
    print("No experiment data found: the AvaNode data should be in "+avalancheDir+"/AvaNode_data/220222_C10_avalanche_GPS.txt")
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
    vel_experimental_C10= vel[index:]
    time_experimental_C10 = [i-time[index] for i in time[index:]] 
    
    # Key to plot only if the experimental data have been found 
    key = 1 

#%% Boxplot 

# Calculating the difference between the avaframe velocity and the sensor velocity
boxplot_data = [None]*number_files
for i in range(0, number_files):
    f = interp1d(time_experimental_C10, vel_experimental_C10)
    if max(Time[i]) > max(time_experimental_C10):
        index = Time[i].index(next((j for j in Time[i] if int(j)>max(time_experimental_C10)),None))
        vel_experimental_C10_interpolated = f(Time[i][0:index-1]) 
        boxplot_data[i] = Mean[i][0:index-1] - vel_experimental_C10_interpolated
    else:
        vel_experimental_C10_interpolated = f(Time[i]) 
        boxplot_data[i] = Mean[i] - vel_experimental_C10_interpolated

# Making the boxplot plot 
fig = plt.figure() 
ax = fig.add_subplot(111)
plt.boxplot(boxplot_data, whis=200)
# x-axis labels
friction_model = F.frictModel[0]
stopping_model = F.stopCrit[0]  
if friction_model=='Coulomb':
    labels = ["mu ="+str(F.mu[i]) for i in range(0,number_files)]
elif friction_model=='Voellmy':
    labels = ["mu ="+str(F.mu[i])+"\nxsi="+str(F.xsi[i]) for i in range(0,number_files)]
elif friction_model=='samosAT':
    labels = ["mu ="+str(F.mu[i])+"\ntau0="+str(F.tau0[i]) for i in range(0,number_files)]
elif friction_model=='VoellmyUpgraded':
        labels = ["mu ="+str(F.mu[i])+"\ntau0="+str(F.tau0[i])+"\nxsi="+str(F.xsi[i]) for i in range(0,number_files)]
else:
    labels = ["No friction model found!" for i in range(0,number_files)]

fig.suptitle('Velocity difference between the mean avaframe velocity and the Node velocity\nSeilbahnrinne, explicit '+friction_model+' model with stopCrit='+str(stopping_model), fontsize= 20)
ax.set_xticklabels(labels, color='black', fontsize=15)
ax.tick_params(axis='y', colors='black', labelsize=15) 
ax.set_ylabel("Velocity difference[m/s]", color='black', fontsize=15)
plt.show() 
    
#%% Plotting velocity envelope 

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
        fig, ax = plt.subplots(1,3,figsize=(20, 20),sharex=True,sharey=True) 
        nrow = 1
        ncol = 3
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
    
    # Colormap 
    cmap = cm.vik
    rgba = cmap(0.2) 
        
    # For 1D plots 
    if nrow == 1 or ncol == 1:
        for i in range(0,max(nrow,ncol)):
            ax[i].plot(Time[simu_number], Max[simu_number], color=cmap(0.01))  # max
            ax[i].plot(Time[simu_number], Min[simu_number], color=cmap(0.01))  # min
            ax[i].plot(Time[simu_number], Mean[simu_number], linestyle='dashed', color=cmap(0.7), markersize=1.8)  # mean
            ax[i].plot(Time[simu_number], Median[simu_number], linestyle='dotted', color=cmap(0.7), markersize=1.8)  # median
            # filling the space between the max and min
            ax[i].fill_between(Time[simu_number], Min[simu_number], Max[simu_number], color=cmap(0.2), alpha=0.2)
            # Experiment (AvaRange)
            ax[i].plot(time_experimental_C07, vel_experimental_C07, color=cmap(0.80))  # AvaNode velocity
            ax[i].plot(time_experimental_C09, vel_experimental_C09, color=cmap(0.90))  # AvaNode velocity
            ax[i].plot(time_experimental_C10, vel_experimental_C10, color=cmap(0.95))  # AvaNode velocity
            # Plotting parameters
            if friction_model=='Coulomb':
                ax[i].set_title("mu ="+str(F.mu[simu_number]), fontsize=18)
            elif friction_model=='Voellmy':
                ax[i].set_title("mu ="+str(F.mu[simu_number])+", xsi="+str(F.xsi[simu_number]), fontsize=18)
            elif friction_model=='samosAT':
                ax[i].set_title("mu ="+str(F.mu[simu_number])+", tau0="+str(F.tau0[simu_number]), fontsize=18)
            elif friction_model=='VoellmyUpgraded':
                ax[i].set_title("mu ="+str(F.mu[simu_number])+", tau0="+str(F.tau0[simu_number])+", xsi="+str(F.xsi[i]), fontsize=15)
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
                ax[i,j].plot(time_experimental_C07, vel_experimental_C07, color=cmap(0.80))  # AvaNode velocity
                ax[i,j].plot(time_experimental_C09, vel_experimental_C09, color=cmap(0.90))  # AvaNode velocity
                ax[i,j].plot(time_experimental_C10, vel_experimental_C10, color=cmap(0.99))  # AvaNode velocity   
                # Plotting parameters
                if friction_model=='Coulomb':
                    ax[i,j].set_title("mu ="+str(F.mu[simu_number]), fontsize=18)
                elif friction_model=='Voellmy':
                    ax[i,j].set_title("mu ="+str(F.mu[simu_number])+", xsi="+str(F.xsi[simu_number]), fontsize=18)
                elif friction_model=='samosAT':
                    ax[i,j].set_title("mu ="+str(F.mu[simu_number])+", tau0="+str(F.tau0[simu_number]), fontsize=18)
                elif friction_model=='VoellmyUpgraded':
                    ax[i,j].set_title("mu ="+str(F.mu[simu_number])+", tau0="+str(F.tau0[simu_number])+", xsi="+str(F.xsi[i]), fontsize=15)
                else:
                    ax[i,j].set_title("No friction model found!", fontsize=18)
                ax[i,j].set_xlim(right=max_time)
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
    
    plt.show() 


    
   