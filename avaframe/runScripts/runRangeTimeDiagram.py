"""
    Run script for module com2AB
"""

import matplotlib.pyplot as plt
from cmcrameri import cm
import pathlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3Plot import outAB
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana5Utils.distanceTimeAnalysis as dtAna
import avaframe.in1Data.getInput as gI
import avaframe.in2Trans.ascUtils as IOf
import avaframe.out3Plot.outDistanceTimeAnalysis as dtAnaPlots


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

# fetch dem data
demInputs = gI.readDEM(avalancheDir)

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
    np.savetxt('rangeMasked.txt', mtiInfo['rangeMasked'])
    mtiInfo['plotTitle'] = 'range-time diagram %s' % index

    # fetch all flow parameter result fields
    flowFieldsDir = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'peakFiles', 'timeSteps')
    simNameSuffix = index + '_' + cfgRangeTime['GENERAL']['rangeTimeResType']
    flowFields =  dtAna.fetchFlowFields(flowFieldsDir, suffix=simNameSuffix)

    for flowField in flowFields:

        # read flow field data
        flowFieldDict = IOf.readRaster(flowField)
        flowF = flowFieldDict['rasterData']

        # extract avalanche front distance to radar
        mtiInfo = dtAna.extractFrontInSim(flowF, mtiInfo,
            cfgRangeTime['GENERAL'].getfloat('thresholdResult'))

        # extract average values of range gates for mti plot
        mtiInfo = dtAna.extractMeanValuesAtRangeGates(cfgRangeTime, flowF, mtiInfo)
        timeStep, _ = dtAna.fetchTimeStepFromName(flowField)
        mtiInfo['timeList'].append(timeStep[0])

    # create plot
    dtAnaPlots.plotRangeTime(mtiInfo, cfgRangeTime['GENERAL'])
