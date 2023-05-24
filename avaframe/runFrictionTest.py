# -*- coding: utf-8 -*-
"""
Created on Oct Mon 12 2022

@author: Oscar Dick

Plot velocity envelopes and several comparison tools for different simulations.
One specific simulation can also be investigated.
The present script has been used to produce the plots in my master thesis as well.

Last change: 10/05/2023

"""
# Python imports
import pathlib

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.com1DFA import particleTools
from avaframe.Tools import PostProcessingTools
# from avaframe.Tools import NodeTools
from avaframe.Tools import PlotFunctions

###############################################################################################################################################################################################################################################################################################################################################################################
##### README ##################################################################################################################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################################################################################################################
# The present script calls each of the plot functions made so far.

# The plotting functions are in PlotFunctions. Useful tools for the plotting are in PlotTools. Useful tools to post-process the
# data and prepare them for the plots are in PostProcessingTools. NodeTools contain useful tools to process the AvaNodes data

# For all: the three most interesting plots are produced by:
#   * Peak flow quantities and thalweg-altitude diagram: PlotFunctions.plotPeakVelVelThalwegEnvelope
#   * Thalweg-Altitude energy line, peak flow quantities and thalweg time diagram: PlotFunctions.plotPeakQuantThalTimeEnergyLine
#   * Comparative plots (vel,acc,s_{xyz}..): PlotFunctions.plotPeakQuantTrackedPartVel

# For Michi: NodeTools.produceAvaNodesDict(avalancheDir) (line 95) produces a dictionary containing all the relevant information on the
# AvaNodes (velocity, trajectory in the same coordinate system than the simulations...). And PlotFunctions.plotPeakQuantTrackedPartVel
# (line 205) produces the 7 comparative plots (vel,acc,s_{xyz}..) you are interested in, even though it is  not finished unfortunately...

#%%############################################################################################################################################################################################################################################################################################################################################################################
##### SECTION TO PREPARE THE DATA FOR THE PLOTS ###############################################################################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################################################################################################################
# WARNING: THE BELOW SECTON IS SPECIFICALLY DEDICATED TO EXTRACT POST-PROCESSING INFORMATION. IT IS LIKELY TO BE IMPROVED


# TO run this script:
# run com1DFA with particles, trackParticles and thalwegTimeDiagram
# run runAna3AIMECTransform.py to add thalweg coordinate system to particlesDict from com1DFA sim

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Load configuration info of all com1DFA simulations and extract the different avalanche simulations of the output file
print('Extracting the different avalanche simulations from the ini file')
Sim, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# Some general and useful parameters
number_ava = len(Sim.index)
avaDir = pathlib.Path(avalancheDir)
modName = 'com1DFA'

# Creating a dictionary containing information on each avalanche
avaDict = {}      # creating a dictionary containing information on the different avalanche simulations
trackedPartProp = {}      # creating a dictionary containing information on the tracked particles for the different simulations
for i in range(0, len(Sim.index)):
    # log simulation number
    print('Avalanche found in the output file: %s' % Sim.index[i])
    # filling the avalanches dictionary with information on avalanche Sim.index[i]
    avaDict[i], _ = particleTools.readPartFromPickle(
        avalancheDir, simName=Sim.index[i], flagAvaDir=True, comModule='Aimec')
    # filling the tracked particles dictionary with information on avalanche Sim.index[i]
    print('Extracting information on the tracked particles for simulation %s' % Sim.index[i])
    pathFile = avalancheDir + '/Outputs/' + modName + '/configurationFiles/' + Sim.simName[i] +'.ini'
    cfgParticles = cfgUtils.readCfgFile(avaDir, module='', fileName=pathFile)
    trackedPartProp[i] = PostProcessingTools.trackedParticlesDictionary(cfgParticles,avaDict[i])


# Running Aimec post processing
print('Running Aimec post processing...')
rasterTransfo, resAnalysisDF, plotDict, newRasters = PostProcessingTools.aimecPostProcess(avalancheDir)

# Calculating the velocity envelope of the simulations (time, max, min, mean, median)
print('Calculating velocity time envelopes...')
dictVelEnvelope = PostProcessingTools.velocityEnvelope(avaDict)

# Extracting experiments data
dictNodes ={}

# Calculating the velocity envelope in the thalweg coordinate system (velocity-altitude-thalweg envelope)
print('Calculating velocity altitude thalweg envelopes...')
dictVelAltThalweg = PostProcessingTools.velocityEnvelopeThalweg(avaDict)

# Calculating the velocity-altitude-thalweg information for the particles
print('Calculating velocity altitude thalweg envelopes for the tracked numerical particles...')
# Calculating the velocity envelope in the thalweg coordinate system (velocity-altitude-thalweg envelope) for the tracked particles
# changing the dictionary structure to make it similar to the particles dictionary
trackedPartPropAdapted = PostProcessingTools.dictChangeTrackedPart(trackedPartProp)
# calculating and adding for each particle sAimec and lAimec (s and l along the thalweg in the projected xy plan)
trackedPartPropAdapted = PostProcessingTools.sAimeclAimecCalculation(trackedPartPropAdapted,avalancheDir)
# calculating and adding for each particle s (trajectory length in the xyz plan)
trackedPartPropAdapted = PostProcessingTools.sCalculation(trackedPartPropAdapted)
# calculating the envelope along the thalweg
#sAimecPart,maxSxyzPart,minSxyzPart, sSortedDeduplicatedPart, sBetaPointPart, simVelThalwegPart, simMeanVelThalwegPart, simMedianVelThalwegPart, simMinVelThalwegPart, simMaxVelThalwegPart, simAltitudeThalwegPart,  simMaxAltThalwegPart, simMinAltThalwegPart, maxAccPart, minAccPart =   PostProcessingTools.velocityEnvelopeThalweg(trackedPartPropAdapted)
dictVelAltThalwegPart = PostProcessingTools.velocityEnvelopeThalweg(trackedPartPropAdapted)

# Calculating the max and min pfv pft and ppr envelopes
print('Calculating the max and min peak flow quantities...')
# Choose your simulation
simNum = 0
# The thresholds pfvMinThreshold pftMinThreshold pprMinThreshold are used to define the masked array
dictRaster = PostProcessingTools.rasterVelField(simNum,avalancheDir,pfvMinThreshold=1,pftMinThreshold=0.1,pprMinThreshold=10)

#%%############################################################################################################################################################################################################################################################################################################################################################################
##### PLOTS ###################################################################################################################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################################################################################################################

AvaNodes = False
# %% Plotting velocity altitude thalweg diagram
# Choose your simulation
simu_number = 0
# # Save the plot the Velocity Altitude Thalweg diagram
# PlotFunctions.plotVelocityAltitudeThalweg(simu_number,dictVelAltThalweg,dictNodes,avaDict,avaDir,rasterTransfo,dictRaster,Title=False,Save=True,AvaNodes=AvaNodes,TrackedPart=True,Raster=False,EnergyLine=False,modName='com1DFA')


# %% Plotting thalweg time diagram for each particle
# Choose your simulation
# simu_number = 0
# # Save the plot of the Thalweg time diagram for each particle
# PlotFunctions.plotThalwegTimeParticles(simu_number,dictVelEnvelope,dictVelAltThalweg,avaDict,avaDir,Save=True,Show=False,modName='com1DFA')

# # %% Plotting velocity time envelopes for all the simulations in the outputs file
# PlotFunctions.plotVelocityTimeEnvelope(number_ava,dictVelEnvelope,dictNodes,Sim,avaDir,AvaNodes=AvaNodes,Show=False,modName='com1DFA')


# %% Plotting range time diagram with peak flow velocity and the velocity time envelope
# if avalancheDir == 'data/avaSeilbahn':
#     PlotFunctions.plotRangeTimePeakVelVelTimeEnvelope(Sim,avalancheDir,number_ava,dictVelEnvelope,dictNodes,Save=True,Show=False,AvaNodes=AvaNodes,modName='com1DFA')


# OK
# %% Plotting peak flow quantities and the velocity thalweg diagram
PlotFunctions.plotPeakVelVelThalwegEnvelope(Sim, avalancheDir, number_ava, rasterTransfo,
    dictVelAltThalweg, Save=True, modName='com1DFA')


# # %% Plotting peak flow quantities and velocity time envelope
# PlotFunctions.plotPeakVelVelEnvelope(Sim,avalancheDir,number_ava,dictNodes,dictVelEnvelope,Save=True,AvaNodes=False,Show=False,modName='com1DFA')
# # PlotFunctions.plotBoxplot(avalancheDir,Sim,number_ava,dictNodes,dictVelEnvelope,TrackedPart=False,Save=True,Show=False,modName='com1DFA')


# %% Plotting velocity time envelopes for all the simulations in the outputs file
# PlotFunctions.plotVelocityTimeEnvelope(number_ava,dictVelEnvelope,dictNodes,Sim,avaDir,AvaNodes=AvaNodes,Show=False,modName='com1DFA')


# %% Plotting range time diagram with peak flow velocity and the velocity time envelope
# PlotFunctions.plotRangeTimePeakVelVelTimeEnvelope(Sim,avalancheDir,number_ava,dictVelEnvelope,dictNodes,Save=True,Show=False,AvaNodes=AvaNodes,modName='com1DFA')


# OK
# %% Plotting peak flow quantities and the velocity thalweg diagram
# PlotFunctions.plotPeakVelVelThalwegEnvelope(Sim,avalancheDir,number_ava,dictNodes,rasterTransfo,dictVelAltThalweg,Save=True,AvaNodes=AvaNodes,modName='com1DFA')


# %% Plotting peak flow quantities and velocity time envelope
# PlotFunctions.plotPeakVelVelEnvelope(Sim,avalancheDir,number_ava,dictNodes,dictVelEnvelope,Save=True,AvaNodes=AvaNodes,Show=False,modName='com1DFA')


# # %%Plotting velocity altitude thalweg diagram for the tracked particles
# # Choose your simulation
# simu_number = 0
# PlotFunctions.plotVelocityAltitudeThalweg(simu_number,dictVelAltThalwegPart,dictNodes,avaDict,avaDir,rasterTransfo,dictRaster,Title=False,Save=True,AvaNodes=AvaNodes,TrackedPart=True,Raster=False,EnergyLine=False,modName='com1DFA')
#
#
# # %% Plotting velocity-altitude-thalweg envelope of the whole avalanche flow for all the different friction parameters
# PlotFunctions.plotVelocityAltitudeThalwegAllSim(dictVelAltThalweg,Sim,dictNodes,avaDir,number_ava,TrackedPart=False,Title=False,Save=True,Show=False,AvaNodes=AvaNodes,modName='com1DFA')

# %% Plotting velocity-altitude-thalweg envelope of the tracked numerical particles for all the different friction parameters
# PlotFunctions.plotVelocityAltitudeThalwegAllSim(dictVelAltThalwegPart,Sim,dictNodes,avaDir,number_ava,TrackedPart=True,Title=False,Save=True,Show=False,AvaNodes=AvaNodes,modName='com1DFA')


#PlotFunctions.plotVelocityAltitudeThalweg(simu_number,dictVelAltThalwegPart,dictNodes,avaDict,avaDir,rasterTransfo,dictRaster,Title=False,Save=True,AvaNodes=AvaNodes,TrackedPart=True,Raster=False,EnergyLine=False,modName='com1DFA')


# %% Plotting velocity-altitude-thalweg envelope of the whole avalanche flow for all the different friction parameters
# PlotFunctions.plotVelocityAltitudeThalwegAllSim(dictVelAltThalweg,Sim,dictNodes,avaDir,number_ava,TrackedPart=False,Title=False,Save=True,Show=False,AvaNodes=AvaNodes,modName='com1DFA')

# %% Plotting velocity-altitude-thalweg envelope of the tracked numerical particles for all the different friction parameters
# PlotFunctions.plotVelocityAltitudeThalwegAllSim(dictVelAltThalwegPart,Sim,dictNodes,avaDir,number_ava,TrackedPart=True,Title=False,Save=True,Show=False,AvaNodes=AvaNodes,modName='com1DFA')


# %%Plotting velocity altitude thalweg diagram for the whole avalanche flow with the raster data and the energy line
# Choose your simulation
# simu_number = 0
# PlotFunctions.plotVelocityAltitudeThalweg(simu_number,dictVelAltThalweg,dictNodes,avaDict,avaDir,rasterTransfo,dictRaster,Title=False,Save=True,AvaNodes=AvaNodes,TrackedPart=False,Raster=True,EnergyLine=True,modName='com1DFA')


# %% Plot the energy line
# Choose your simulation
# simu_number = 0
# PlotFunctions.plotEnergyLine(avalancheDir,avaDict,simu_number,Sim,dictVelAltThalweg,rasterTransfo,dictRaster,Show=False,Save=True,modName='com1DFA')


# %% Plotting velocity altitude thalweg and thalweg time diagram
# Choose your simulation
# simu_number = 0
# PlotFunctions.plotVelAltThalTimeDiag(avalancheDir,number_ava,simu_number,avaDict,Sim,dictVelAltThalweg,rasterTransfo,dictNodes,Show=False,Save=True,AvaNodes=AvaNodes,modName='com1DFA')


# OK
# %% Plotting peak flow velocity, thalweg time diagram and energy line
from avaframe.Tools import PlotFunctions
PlotFunctions.plotPeakQuantThalTimeEnergyLine(avalancheDir, number_ava, Sim, plotDict,
    rasterTransfo, dictRaster, Save=True, Show=False, modName='com1DFA')

# OK
# %% Plotting peak flow quantities with tracked particles and velocities
PlotFunctions.plotPeakQuantTrackedPartVel(avalancheDir, number_ava, Sim, avaDict, dictVelAltThalweg,
    dictVelAltThalwegPart, trackedPartProp, trackedPartPropAdapted, dictVelEnvelope,
    Save=True, Show=False, modName='com1DFA')
