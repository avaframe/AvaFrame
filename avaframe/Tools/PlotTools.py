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


#%% function to plot peak flow velocity

def PeakFields(avaDir, peakFilesDF, simNum, demData=''):
    """ Prepare the data for the plot of peak file number simNum and return dictionary with data for the plot
        with DEM in background

        Parameters
        ----------
        avaDir : str
            path to avalanche directoy
        peakFilesDF : dict
            dictionary containing info on the peak files
        simNum : int
            simulation number
        demData: dictionary
            optional - if not the dem in the avaDir/Inputs folder has been used but a different one

        Returns
        -------
        plotDict : dict
            dictionary with info for the plot
        """

    # getting the dem data
    if demData == '':
        demFile = gI.getDEMPath(avaDir)
        demData = IOf.readRaster(demFile, noDataToNan=True)
        demDataField = demData['rasterData']
    else:
        # check if noDataValue is found and if replace with nans for plotting
        demDataField = np.where(demData['rasterData'] == demData['header']['noDataValue'], np.nan, demData['rasterData'])
    demField = demDataField

    # Load names and paths of peakFiles
    print(peakFilesDF['avaName'][simNum])
    name = peakFilesDF['names'][simNum]
    fileName = peakFilesDF['files'][simNum]
    avaName = peakFilesDF['avaName'][simNum]
    resType = peakFilesDF['resType'][simNum]

    # Load data
    raster = IOf.readRaster(fileName, noDataToNan=True)
    data = raster['rasterData']

    # constrain data to where there is data
    cellSize = peakFilesDF['cellSize'][simNum]
    rowsMin, rowsMax, colsMin, colsMax = pU.constrainPlotsToData(data, cellSize)
    dataConstrained = data[rowsMin:rowsMax+1, colsMin:colsMax+1]
    demConstrained = demField[rowsMin:rowsMax+1, colsMin:colsMax+1]

    data = np.ma.masked_where(dataConstrained == 0.0, dataConstrained)
    dataNotConstrained = np.ma.masked_where(data != 0.0, data)
    unit = pU.cfgPlotUtils['unit%s' % resType]

    # Preparing output dictionary
    plotDict = {}
    plotDict['resType'] = resType
    plotDict['data'] = data
    plotDict['dataNotConstrained'] = dataNotConstrained
    plotDict['raster'] = raster
    plotDict['rowsMin'] = rowsMin
    plotDict['rowsMax'] = rowsMax
    plotDict['colsMin'] = colsMin
    plotDict['colsMax'] = colsMax
    plotDict['cellSize'] = cellSize
    plotDict['demConstrained'] = demConstrained
    plotDict['dem'] = demField
    plotDict['unit'] = unit
    plotDict['name'] = name

    return  plotDict


#%% function to plot a range-time diagram

def RangeTimeDiagram(avalancheDir, index,simDF):

    # Load all input Parameters from config file
    # get the configuration of an already imported module
    cfgRangeTime = cfgUtils.getModuleConfig(dtAna)
    cfgRangeTime['GENERAL']['avalancheDir'] = avalancheDir

    # fetch dem data
    demInputs = gI.readDEM(avalancheDir)

    # fetch all flow parameter result fields
    configDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'configurationFiles')
    if (configDir.is_dir() is False) or ( len(list(configDir.glob('*.ini'))) == 0):
        fU.fileNotFoundMessage(('No configuration files found in %s - consider first running avalanche simulations' %
            configDir))

    # fetch info on available simulations
    #simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False, specDir='')

    #for index, simDFrow in simDF.iterrows():

    # add simHash info
    cfgRangeTime['GENERAL']['simHash'] = index

    # check if DEM in Inputs has been used, if not lead simulation DEM
    #if 'DEMremeshed' in simDFrow['DEM']:
    if 'DEMremeshed' in simDF.DEM[index]:
        # get dem dictionary - already read DEM with correct mesh cell size
        dem = gI.initializeDEM(simDF.avalancheDir[index], simDF.DEM[index])
        #log.info('Remeshed DEM read from: %s/%s' % (simDF.avalancheDir[index], simDF.DEM[index]))
    else:
        dem = demInputs.copy()

    # setup required data
    mtiInfo = dtAna.setupRangeTimeDiagram(dem, cfgRangeTime)
    # create plot of range distance already masked with radar field of view
    #dtAnaPlots.radarFieldOfViewPlot(mtiInfo['radarFov'], mtiInfo['rangeMasked'], cfgRangeTime,
   #     mtiInfo['rangeGates'], dem)
    #np.savetxt('rangeMasked.txt', mtiInfo['rangeMasked'])
    mtiInfo['plotTitle'] = 'range-time diagram %s' % index

    # fetch all flow parameter result fields
    flowFieldsDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'peakFiles', 'timeSteps')
    simNameSuffix = simDF.simName[index] + '_' + cfgRangeTime['GENERAL']['rangeTimeResType']
    flowFields = fU.fetchFlowFields(flowFieldsDir, suffix=simNameSuffix)

    # check if simulation results are available
    if len(flowFields) == 0:
        fU.fileNotFoundMessage(('No flow variable results found in %s - consider first running avalanche simulations' %
            flowFieldsDir))

    for flowField in flowFields:

        # read flow field data
        flowFieldDict = IOf.readRaster(flowField)
        flowF = flowFieldDict['rasterData']

        # extract avalanche front distance to radar and average values of range gates for mti plot
        mtiInfo = dtAna.extractFrontAndMeanValuesRadar(cfgRangeTime, flowF, mtiInfo)
        timeStep, _ = dtAna.fetchTimeStepFromName(flowField)
        mtiInfo['timeList'].append(timeStep[0])

        # create plot
        #dtAnaPlots.plotRangeTime(mtiInfo, cfgRangeTime, F, index )
    return(mtiInfo,cfgRangeTime)



#%% function to plot a thalweg-time diagram
def ThalwegTimeDiagram(avalancheDir, index, simDF):
    # Load all input Parameters from config file
    # get the configuration of an already imported module
    cfgRangeTime = cfgUtils.getModuleConfig(dtAna)
    cfgRangeTime['GENERAL']['avalancheDir'] = avalancheDir

    # fetch dem data
    demInputs = gI.readDEM(avalancheDir)

    # fetch all flow parameter result fields
    configDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'configurationFiles')
    if (configDir.is_dir() is False) or ( len(list(configDir.glob('*.ini'))) == 0):
        fU.fileNotFoundMessage(('No configuration files found in %s - consider first running avalanche simulations' %
            configDir))

    # fetch info on available simulations
    #simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False, specDir='')

    #for index, simDFrow in simDF.iterrows():

    # add simHash info
    cfgRangeTime['GENERAL']['simHash'] = index

    # check if DEM in Inputs has been used, if not lead simulation DEM
    #if 'DEMremeshed' in simDFrow['DEM']:
    if 'DEMremeshed' in simDF.DEM[index]:
        # get dem dictionary - already read DEM with correct mesh cell size
        dem = gI.initializeDEM(simDF.avalancheDir[index], simDF.DEM[index])
        #log.info('Remeshed DEM read from: %s/%s' % (simDF.avalancheDir[index], simDF.DEM[index]))
    else:
        dem = demInputs.copy()

    # setup required data
    mtiInfo = dtAna.setupThalwegTimeDiagram(dem, cfgRangeTime)
    # create plot of range distance already masked with radar field of view
    #dtAnaPlots.radarFieldOfViewPlot(mtiInfo['radarFov'], mtiInfo['rangeMasked'], cfgRangeTime,
   #     mtiInfo['rangeGates'], dem)
    #np.savetxt('rangeMasked.txt', mtiInfo['rangeMasked'])
    mtiInfo['plotTitle'] = mtiInfo['plotTitle'] = 'thalweg-time diagram %s' % index

    # fetch all flow parameter result fields
    flowFieldsDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'peakFiles', 'timeSteps')
    simNameSuffix = simDF.simName[index] + '_' + cfgRangeTime['GENERAL']['rangeTimeResType']
    flowFields = fU.fetchFlowFields(flowFieldsDir, suffix=simNameSuffix)

    # check if simulation results are available
    if len(flowFields) == 0:
        fU.fileNotFoundMessage(('No flow variable results found in %s - consider first running avalanche simulations' %
            flowFieldsDir))

    for flowField in flowFields:
        # read flow field data
        flowFieldDict = IOf.readRaster(flowField)
        flowF = flowFieldDict['rasterData']

        # extract avalanche front distance to radar and average values of range gates for mti plot
        mtiInfo = dtAna.extractFrontAndMeanValuesTT(cfgRangeTime, flowF, dem['header'], mtiInfo)
        timeStep, _ = dtAna.fetchTimeStepFromName(flowField)
        mtiInfo['timeList'].append(timeStep[0])

        # create plot
        #dtAnaPlots.plotRangeTime(mtiInfo, cfgRangeTime, F, index )
    return(mtiInfo,cfgRangeTime)


#%% function to generate the data for particles plot on the dem

def demParticles(simu_number):

    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']

    avaDir = pathlib.Path(avalancheDir)
    modName = 'com1DFA'

    # Load all infos from the peak files
    inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    plotDict = PeakFields(avalancheDir, peakFilesDF, simu_number, demData='')

    return plotDict



#%% function to generate the correct title for the plot

def titleFrictParam(Sim,simIndex):
    friction_model = Sim.frictModel[0]
    if friction_model == 'Coulomb':
        title="$\mu$ ="+str(Sim['mucoulomb'].loc[simIndex])
    elif friction_model == 'Voellmy':
        title = "$\mu$ ="+str(Sim['muvoellmy'].loc[simIndex]) +r", $\xi$="+str(Sim['xsivoellmy'].loc[simIndex])
    elif friction_model == 'samosAT':
        title = "$\mu$ ="+str(Sim['musamosat'].loc[simIndex]) +r", $\tau_0$="+str(Sim['tau0samosat'].loc[simIndex])
    elif friction_model == 'VoellmyUpgraded':
        title = "$\mu$ ="+str(Sim['muvoellmyupgraded'].loc[simIndex])+r", $\tau_0$="+str(
            Sim['tau0voellmyupgraded'].loc[simIndex])+r", $\xi$="+str(Sim['xsivoellmyupgraded'].loc[simIndex])
    else:
        title = "No friction model found!"
    return title


#%% function to generate the correct labels for the boxplots

def labelFrictParamBoxplot(Sim,simu_number):
    friction_model = Sim.frictModel[0]
    if friction_model == 'Coulomb':
        title="$\mu$ ="+str(Sim.mu[simu_number])
    elif friction_model == 'Voellmy':
        title = "$\mu$ ="+str(Sim.mu[simu_number])+"\n"+r"$\xi$="+str(Sim.xsi[simu_number])
    elif friction_model == 'samosAT':
        title = "$\mu$ ="+str(Sim.mu[simu_number])+"\n"+r"$\tau_0$="+str(Sim.tau0[simu_number])
    elif friction_model == 'VoellmyUpgraded':
        title = "$\mu$ ="+str(Sim.mu[simu_number])+"\n"+r"$\tau_0$="+str(
            Sim.tau0[simu_number])+"\n"+r"$\xi$="+str(Sim.xsi[simu_number])
    else:
        title = "No friction model found!"
    return title


#%% function to prepare the subplots depending on the number of simulations

def subplot(number_ava):
    if number_ava == 1:
        fig, ax = plt.subplots(1, 1, figsize=(20, 20), sharex=True)
        nrow = 1
        ncol = 1
    elif number_ava == 2:
        fig, ax = plt.subplots(2, 1, figsize=(pU.figW+10, pU.figH+3), sharex=True, sharey=True)
        nrow = 2
        ncol = 1
    elif number_ava == 3:
        fig, ax = plt.subplots(3, 1, figsize=(20, 20), sharex=True, sharey=True)
        nrow = 3
        ncol = 1
    elif number_ava == 4:
        fig, ax = plt.subplots(2, 2, figsize=(20, 20), sharex=True, sharey=True)
        nrow = 2
        ncol = 2
    elif number_ava == 5:
        fig, ax = plt.subplots(2, 3, figsize=(20, 20), sharex=True, sharey=True)
        nrow = 2
        ncol = 3
    elif number_ava == 6:
        fig, ax = plt.subplots(2, 3, figsize=(20, 20), sharex=True, sharey=True)
        nrow = 2
        ncol = 3
    elif number_ava == 7:
        fig, ax = plt.subplots(3, 3, figsize=(20, 20), sharex=True, sharey=True)
        nrow = 3
        ncol = 3
    elif number_ava == 8:
        fig, ax = plt.subplots(2, 4, figsize=(20, 20), sharex=True, sharey=True)
        nrow = 2
        ncol = 4
    elif number_ava == 9:
        fig, ax = plt.subplots(3, 3, figsize=(pU.figW+10, pU.figH+3), sharex=True, sharey=True)
        nrow = 3
        ncol = 3
    else:
        print('Too many simulations for one plot!')
    return (nrow,ncol,fig,ax)
