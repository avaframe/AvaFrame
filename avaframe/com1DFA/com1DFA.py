"""
    Main functions for python DFA kernel
"""

import logging
import os
import sys
import time
import pathlib
import pandas as pd
from functools import partial
import platform
if os.name == 'nt':
    # from multiprocessing.pool import ThreadPool as Pool
    from multiprocessing import Pool
elif platform.system() == 'Darwin':
    from multiprocessing.pool import ThreadPool as Pool
else:
    from multiprocessing import Pool

# Local imports
import avaframe.com1DFA.com1DFATools as com1DFATools
from avaframe.in3Utils import cfgUtils
import avaframe.com1DFA.deriveParameterSet as dP
import avaframe.com1DFA.com1DFACore as com1DFACore
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.log2Report import generateReport as gR
import avaframe.in3Utils.fileHandlerUtils as fU

#######################################
# Set flags here
#######################################
# create local logger
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')


def com1DFAPreprocess(cfgMain, typeCfgInfo, cfgInfo):
    """ preprocess information from configuration, read input data and gather into inputSimFiles,
        create one config object for each of all desired simulations,
        create dataFrame with one line per simulations of already existing sims in avalancheDir

        Parameters
        ------------
        cfgMain: configparser object
            main configuration of AvaFrame
        typeCfgInfo: str
            name of type of cfgInfo (cfgFromFile or cfgFromObject)
        cfgInfo: str or pathlib Path or configparser object
            path to configuration file if overwrite is desired - optional
            if not local (if available) or default configuration will be loaded
            if cfgInfo is a configparser object take this as initial config

        Returns
        --------
        simDict: dict
            dictionary with one key per simulation to perform including its config object
        inputSimFiles: dict
            dictionary with input files info
        outDir: str
            path to store outputs
    """

    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    # Get info about the current module
    cuModule = sys.modules[__name__]

    # read initial configuration
    if typeCfgInfo in ['cfgFromFile', 'cfgFromDefault']:
        cfgStart = cfgUtils.getModuleConfig(cuModule, fileOverride=cfgInfo, toPrint=False)

        # cfgStart = cfgUtils.getModuleConfig(com1DFA, fileOverride=cfgInfo, toPrint=False)
    elif typeCfgInfo == 'cfgFromObject':
        cfgStart = cfgInfo

    # fetch input data and create work and output directories
    inputSimFilesAll, outDir, simDFExisting, simNameExisting = com1DFATools.initializeInputs(avalancheDir,
                                                                    cfgStart['GENERAL'].getboolean('cleanDEMremeshed'))

    # create dictionary with one key for each simulation that shall be performed
    simDict = dP.createSimDict(avalancheDir, cuModule, cfgStart, inputSimFilesAll, simNameExisting)

    return simDict, outDir, inputSimFilesAll, simDFExisting


def com1DFAMain(cfgMain, cfgInfo=''):
    """ preprocess information from ini and run all desired simulations, create outputs and reports

        Parameters
        ------------
        cfgMain: configparser object
            main configuration of AvaFrame
        cfgInfo: str or pathlib Path or configparser object
            path to configuration file if overwrite is desired - optional
            if not local (if available) or default configuration will be loaded
            if cfgInfo is a configparser object take this as initial config
            if path to a directory is provided - load one or multiple override configuration files

        Returns
        --------
        dem: dict
            dictionary with dem header and raster data (that has been used for the final run)
        plotDict: dict
            information on result plot paths
        reportDictList: list
            list of report dictionaries for all performed simulations
        simDF: pandas dataFrame
            configuration dataFrame of the simulations computed (if no simulation computed, configuration dataFrame
            of the already existing ones)
    """

    avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    # fetch type of cfgInfo
    typeCfgInfo = com1DFATools.checkCfgInfoType(cfgInfo)
    if typeCfgInfo == 'cfgFromDir':
        # preprocessing to create configuration objects for all simulations to run by reading multiple cfg files
        simDict, inputSimFiles, simDFExisting, outDir = com1DFATools.createSimDictFromCfgs(cfgMain, cfgInfo)
    else:
        # preprocessing to create configuration objects for all simulations to run
        simDict, outDir, inputSimFiles, simDFExisting = com1DFAPreprocess(cfgMain, typeCfgInfo, cfgInfo)

    log.info('The following simulations will be performed')
    for key in simDict:
        log.info('Simulation: %s' % key)

    # initialize reportDict list
    reportDictList = list()

    # is there any simulation to run?
    if bool(simDict):
        # reset simDF and timing
        simDF = pd.DataFrame()
        tCPUDF = pd.DataFrame()
        dem = dict()

        startTime = time.time()

        # create a few necessary output dirs:
        outDirData = outDir / 'particles'
        fU.makeADir(outDirData)
        outDirPeakAll = outDir / 'peakFiles'
        fU.makeADir(outDirPeakAll)
        outDirPeak = outDir / 'peakFiles' / 'timeSteps'
        fU.makeADir(outDirPeak)

        log.info("--- STARTING (potential) PARALLEL PART ----")
        # Get number of CPU Cores wanted
        nCPU = cfgUtils.getNumberOfProcesses(cfgMain, len(simDict))

        results = {}
        pool = Pool(processes=nCPU)
        start_time = time.time()

        for keyValue in simDict.items():
            key, myVal = keyValue
            log.info(f'Submitting {keyValue} to worker process')
            results[key] = pool.apply_async(com1DFACore.com1DFACoreTask, args=(keyValue,
                                                                         inputSimFiles,
                                                                         avalancheDir,
                                                                         outDir))
            # results.append(result)

        resultChecked = []
        # This code block can be (re)run whenever you want to check on the progress of the pool.
        while 1:
            time.sleep(10)
            running, successful, error = [], [], []
            for key, result in results.items():
                try:
                    if result.successful():
                        successful.append(key)
                        resultChecked.append(key)
                        res = result.get()
                        simDF = pd.concat([simDF, res[0]], axis=0)
                        tCPUDF = pd.concat([tCPUDF, res[1]], axis=0)
                        dem = res[2]  # only last dem is used
                        reportDictList.append(res[3])
                    else:
                        error.append(key)
                except ValueError:
                    running.append(key)
            for key in successful:
                del results[key]

            log.debug(f'Running: {sorted(running)}')
            log.info(f'Currently {len(running)} sims running, '
                     f'{len(resultChecked)} are done, '
                     f'{len(error)} with errors')
            # print('Successful:', sorted(successful))
            # print('Error:', sorted(error))
            # print('ResultChecked:', sorted(resultChecked))
            if len(running) < 1:
                break

        pool.close()
        pool.join()

        # for result in results:
        #     res = result.get()
        #     simDF = pd.concat([simDF, res[0]], axis=0)
        #     tCPUDF = pd.concat([tCPUDF, res[1]], axis=0)
        #     dem = res[2]  #only last dem is used
        #     reportDictList.append(res[3])



        # # Supply compute task with inputs
        # com1DFACoreTaskWithInput = partial(com1DFACore.com1DFACoreTask, simDict, inputSimFiles, avalancheDir, outDir)

        # # Create parallel pool and run
        # # with multiprocessing.Pool(processes=nCPU) as pool:
        # with Pool(processes=nCPU) as pool:
        #     results = pool.map(com1DFACoreTaskWithInput, simDict)
        #     pool.close()
        #     pool.join()

        # # Split results to according structures
        # for result in results:
        #     simDF = pd.concat([simDF, result[0]], axis=0)
        #     tCPUDF = pd.concat([tCPUDF, result[1]], axis=0)
        #     dem = result[2]  #only last dem is used
        #     reportDictList.append(result[3])

        timeNeeded = '%.2f' % (time.time() - startTime)
        log.info('Overall (parallel) com1DFA computation took: %s s ' % timeNeeded)
        log.info("--- ENDING (potential) PARALLEL PART ----")

        # postprocessing: writing report, creating plots
        dem, plotDict, reportDictList, simDFNew = com1DFAPostprocess(
            simDF, tCPUDF, simDFExisting, cfgMain, dem, reportDictList
        )

        return dem, plotDict, reportDictList, simDFNew

    else:

        log.warning('There is no simulation to be performed for releaseScenario')

        return 0, {}, [], ''

def com1DFAPostprocess(simDF, tCPUDF, simDFExisting, cfgMain, dem, reportDictList):
    """ postprocessing of simulation results: save configuration to csv, create plots and report

        Parameters
        -----------
        simDF: pandas DataFrame
            dataframe with one line per simulation and info on parameters used
        tCPUDF:
            computation time
        simDFExisting: pandas DataFrame
            dataframe with one line per simulation and info on parameters used before
            simulations have been performed
        cfgMain: configparser object
            global avaframe config
        dem: dict
            dem dictionary
        reportDictList: list
            list of dictionaries for each simulation with info for report creation

        Returns
        --------
        dem: dict
            dictionary with dem header and raster data (that has been used for final sim)
        plotDict: dict
            information on result plot paths
        reportDictList: list
            list of report dictionaries for all performed simulations
        simDFNew: pandas dataFrame
            configuration dataFrame of the simulations computed and the ones that have been already
            in the Outputs folder (if no simulation computed, configuration dataFrame
            of the already existing ones)
    """

    modName = 'com1DFA'
    avalancheDir = cfgMain['MAIN']['avalancheDir']

    # prepare for writing configuration info
    simDF = cfgUtils.convertDF2numerics(simDF)
    # add cpu time info to the dataframe
    simDF = simDF.join(tCPUDF)

    # write the actually simulated sims to a separate csv file,
    # this is used for the qgis connector
    cfgUtils.writeAllConfigurationInfo(avalancheDir, simDF, specDir='',
                                       csvName='latestSims.csv')

    # append new simulations configuration to old ones (if they exist),
    # return total dataFrame and write it to csv
    simDFNew = pd.concat([simDF, simDFExisting], axis=0)
    cfgUtils.writeAllConfigurationInfo(avalancheDir, simDFNew, specDir='')

    # Set directory for report
    reportDir = pathlib.Path(avalancheDir, 'Outputs', modName, 'reports')
    # Generate plots for all peakFiles
    plotDict = oP.plotAllPeakFields(avalancheDir, cfgMain['FLAGS'], modName, demData=dem)

    if cfgMain['FLAGS'].getboolean('createReport'):
        # write report
        gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'].getboolean('reportOneFile'), plotDict)

    return dem, plotDict, reportDictList, simDFNew
