"""
    Run script for creating a range-time diagram from avalanche simulations
"""

import matplotlib.pyplot as plt
from cmcrameri import cm
import pathlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import configparser

# Local imports
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana5Utils.distanceTimeAnalysis as dtAna
import avaframe.in1Data.getInput as gI
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.outDistanceTimeAnalysis as dtAnaPlots


#+++++++++++++USER INPUT++++++++
# if True use preprocessedData
preProcessedData = False
#++++++++++++++++++++++++++

# log file name; leave empty to use default runLog.log
logName = 'runRangeTimeDiagram'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
# write config to log file
cfgRangeTime = cfgUtils.getModuleConfig(dtAna)
cfgRangeTime['GENERAL']['avalancheDir'] = avalancheDir

if preProcessedData:

    log.info('Distance -time analysis performed with preprocessed data')
    # fetch preprocesses mtiInfo dicts
    mtiInfoDicts = dtAna.importMTIData(avalancheDir, 'com1DFA', inputDir='', simHash='')

    # loop over all mtiInfo dicts found
    for mtiInfo in mtiInfoDicts:
        # log file path of info dict
        log.info('Process mtiInfo: %s' % mtiInfo['name'].stem)

        # read configuration from mtiInfo and convert to configparser
        cfgRangeTimeRead = cfgUtils.convertDictToConfigParser(mtiInfo['configurationSettings'])
        # set section GENEREAL but keep section PLOTS from module ini file to be adjustable
        cfgRangeTime['GENERAL'] = cfgRangeTimeRead['GENERAL']

        # TODO: not so nice, logged settings twice now but required here as overwrites other settings
        log.info('Settings used for analysis')
        for section in cfgRangeTime.sections():
            for key in cfgRangeTime.items(section):
                log.info('section: %s - values %s' % (section, key))

        # produce radar field of view plot if type is range time
        if mtiInfo['type'] != 'rangeTime':
            log.warning('Found %s - but type is thalweg-time cannot produce simulated range-time diagram' % mtiInfo['name'].stem)
        else:
            # create plot of range distance already masked with radar field of view
            dtAnaPlots.radarFieldOfViewPlot(mtiInfo['radarFov'], mtiInfo['rangeMasked'], cfgRangeTime,
                mtiInfo['rangeGates'], mtiInfo['dem'])
            # create plot
            dtAnaPlots.plotRangeTime(mtiInfo, cfgRangeTime)

else:
    # fetch dem data
    demInputs = gI.readDEM(avalancheDir)

    # fetch all flow parameter result fields
    configDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'configurationFiles')
    if (configDir.is_dir() is False) or ( len(list(configDir.glob('*.ini'))) == 0):
        fU.fileNotFoundMessage(('No configuration files found in %s - consider first running avalanche simulations' %
            configDir))

    # fetch info on available simulations
    simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False, specDir='')
    for index, simDFrow in simDF.iterrows():

        log.info('Perform range-time diagram for simulation: %s' % index)

        # add simHash info
        cfgRangeTime['GENERAL']['simHash'] = index

        # check if DEM in Inputs has been used, if not lead simulation DEM
        if 'DEMremeshed' in simDFrow['DEM']:
            # get dem dictionary - already read DEM with correct mesh cell size
            dem = gI.initializeDEM(simDFrow['avalancheDir'], simDFrow['DEM'])
            log.info('Remeshed DEM read from: %s/%s' % (simDFrow['avalancheDir'], simDFrow['DEM']))
        else:
            dem = demInputs.copy()

        # setup required data
        mtiInfo = dtAna.setupRangeTimeDiagram(dem, cfgRangeTime)
        # create plot of range distance already masked with radar field of view
        dtAnaPlots.radarFieldOfViewPlot(mtiInfo['radarFov'], mtiInfo['rangeMasked'], cfgRangeTime,
            mtiInfo['rangeGates'], dem)
        np.savetxt('rangeMasked.txt', mtiInfo['rangeMasked'])
        mtiInfo['plotTitle'] = 'range-time diagram %s' % index

        # fetch all flow parameter result fields
        flowFieldsDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'peakFiles', 'timeSteps')
        simNameSuffix = index + '_' + cfgRangeTime['GENERAL']['rangeTimeResType']
        flowFields =  fU.fetchFlowFields(flowFieldsDir, suffix=simNameSuffix)

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
        dtAnaPlots.plotRangeTime(mtiInfo, cfgRangeTime)
