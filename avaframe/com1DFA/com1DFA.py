"""
    Main functions for python DFA kernel
"""

import logging
import os
import time
import pathlib
import numpy as np
import pandas as pd
import math
import copy
import pickle
from datetime import datetime
import matplotlib.path as mpltPath
from itertools import product
import matplotlib.tri as tri
from shapely.geometry import Polygon as sPolygon
from functools import partial
import platform
if os.name == 'nt':
    from multiprocessing.pool import ThreadPool as Pool
elif platform.system() == 'Darwin':
    from multiprocessing.pool import ThreadPool as Pool
else:
    from multiprocessing import Pool

# Local imports
from avaframe.version import getVersion
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.geoTrans as geoTrans
from avaframe.in3Utils import initializeProject as iP
import avaframe.com1DFA.timeDiscretizations as tD
import avaframe.out3Plot.outCom1DFA as outCom1DFA
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.com1DFA.com1DFATools as com1DFATools
import avaframe.com1DFA.particleTools as particleTools
import avaframe.com1DFA.DFAfunctionsCython as DFAfunC
import avaframe.com1DFA.DFAToolsCython as DFAtllsC
import avaframe.com1DFA.damCom1DFA as damCom1DFA
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
import avaframe.in3Utils.geoTrans as gT
import avaframe.out3Plot.outDebugPlots as debPlot
import avaframe.com1DFA.deriveParameterSet as dP
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in1Data import getInput as gI
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.log2Report import generateReport as gR
from avaframe.com1DFA import particleInitialisation as pI
from avaframe.com1DFA import checkCfg
from avaframe.ana5Utils import distanceTimeAnalysis as dtAna
import avaframe.out3Plot.outDistanceTimeAnalysis as dtAnaPlots
import avaframe.out3Plot.outAIMEC as oA
import avaframe.out3Plot.outQuickPlot as oQ
import threading

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

    # read initial configuration
    if typeCfgInfo in ['cfgFromFile', 'cfgFromDefault']:
        cfgStart = cfgUtils.getModuleConfig(com1DFA, fileOverride=cfgInfo, toPrint=False)
    elif typeCfgInfo == 'cfgFromObject':
        cfgStart = cfgInfo

    # fetch input data and create work and output directories
    inputSimFilesAll, outDir, simDFExisting, simNameExisting = com1DFATools.initializeInputs(avalancheDir, cfgStart['GENERAL'].getboolean('cleanDEMremeshed'))

    # create dictionary with one key for each simulation that shall be performed
    simDict = dP.createSimDict(avalancheDir, com1DFA, cfgStart, inputSimFilesAll, simNameExisting)

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
        exportFlag = simDict[key]['cfgSim']['EXPORTS'].getboolean('exportData')

    # initialize reportDict list
    reportDictList = list()

    # is there any simulation to run?
    if bool(simDict):
        # reset simDF and timing
        simDF = pd.DataFrame()
        tCPUDF = pd.DataFrame()
        dem = dict()

        startTime = time.time()

        log.info("--- STARTING (potential) PARALLEL PART ----")
        # Get number of CPU Cores wanted
        nCPU = cfgUtils.getNumberOfProcesses(cfgMain, len(simDict))

        # Supply compute task with inputs
        com1DFACoreTaskWithInput = partial(com1DFACoreTask, simDict, inputSimFiles, avalancheDir, outDir)

        # Create parallel pool and run
        # with multiprocessing.Pool(processes=nCPU) as pool:
        with Pool(processes=nCPU) as pool:
            results = pool.map(com1DFACoreTaskWithInput, simDict)
            pool.close()
            pool.join()

        # Split results to according structures
        for result in results:
            simDF = pd.concat([simDF, result[0]], axis=0)
            tCPUDF = pd.concat([tCPUDF, result[1]], axis=0)
            dem = result[2] #only last dem is used
            reportDictList.append(result[3])

        timeNeeded = '%.2f' % (time.time() - startTime)
        log.info('Overall (parallel) com1DFA computation took: %s s ' %timeNeeded)
        log.info("--- ENDING (potential) PARALLEL PART ----")

        # postprocessing: writing report, creating plots
        dem, plotDict, reportDictList, simDFNew = com1DFAPostprocess(
            simDF, tCPUDF, simDFExisting, cfgMain, dem, reportDictList, exportData=exportFlag
        )

        return dem, plotDict, reportDictList, simDFNew

    else:

        log.warning('There is no simulation to be performed for releaseScenario')

        return 0, {}, [], ''


def com1DFACoreTask(simDict, inputSimFiles, avalancheDir, outDir, cuSim):
    """ This is a subdivision of com1DFAMain to allow for parallel execution.
        Please read this in the context of the com1DFAMain function.
    """


    simDF = pd.DataFrame()
    tCPUDF = pd.DataFrame()

    # load configuration object for current sim
    cfg = simDict[cuSim]["cfgSim"]

    # check configuraton for consistency
    checkCfg.checkCfgConsistency(cfg)

    # fetch simHash for current sim
    simHash = simDict[cuSim]["simHash"]

    log.info('%s runs as process: %s, %s' % (cuSim,  os.getpid(), threading.current_thread().ident))

    # append configuration to dataframe
    simDF = cfgUtils.appendCgf2DF(simHash, cuSim, cfg, simDF)

    # log simulation name
    log.info("Run simulation: %s" % cuSim)

    # ++++++++++PERFORM com1DFA SIMULAITON++++++++++++++++
    (
        dem,
        reportDict,
        cfgFinal,
        tCPU,
        particlesList,
    ) = com1DFA.com1DFACore(cfg, avalancheDir, cuSim, inputSimFiles, outDir, simHash=simHash)

    simDF.at[simHash, "nPart"] = str(int(particlesList[0]["nPart"]))

    # append time to data frame
    tCPUDF = cfgUtils.appendTcpu2DF(simHash, tCPU, tCPUDF)

    # create hash to check if configuration didn't change
    simHashFinal = cfgUtils.cfgHash(cfgFinal)
    if simHashFinal != simHash:
        cfgUtils.writeCfgFile(avalancheDir, com1DFA, cfg, fileName="%s_butModified" % simHash)
        message = "Simulation configuration has been changed since start"
        log.error(message)
        raise AssertionError(message)

    # return simDF, tCPUDF, simDFExisting, cfg, cfgMain, dem, reportDictList
    return simDF, tCPUDF, dem, reportDict


def com1DFAPostprocess(simDF, tCPUDF, simDFExisting, cfgMain, dem, reportDictList, exportData):
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
        exportData: bool
            if True result fields are exported and plots generated

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

    # create plots and report
    reportDir = pathlib.Path(avalancheDir, 'Outputs', modName, 'reports')
    fU.makeADir(reportDir)
    # Generate plots for all peakFiles
    if exportData:
        plotDict = oP.plotAllPeakFields(avalancheDir, cfgMain['FLAGS'], modName, demData=dem)
    else:
        plotDict = ''
        # create contour line plot
        reportDictList, _ = outCom1DFA.createContourPlot(reportDictList, avalancheDir, simDF)


    if cfgMain['FLAGS'].getboolean('createReport'):
        # write report
        gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'].getboolean('reportOneFile'), plotDict)

    return dem, plotDict, reportDictList, simDFNew


def com1DFACore(cfg, avaDir, cuSimName, inputSimFiles, outDir, simHash=''):
    """ Run main com1DFA model

    This will compute a dense flow avalanche with the settings specified in cfg and the name cuSimName

    Parameters
    ----------
    cfg : configparser object
        configuration object for simulation to be performed
    cuSimName: str
        name of simulation
    inputSimFiles: dict
        dictionary with input files, release scenario chosen according to inputSimFiles['releaseScenario']
    avaDir : str or pathlib object
        path to avalanche directory
    outDir: str or pathlib object
        path to Outputs
    simHash: str
        unique sim ID

    Returns
    -------
    reportDictList : list
        list of dictionaries that contain information on simulations that can be used for report generation
    dem: dict
        dictionary with info on header and dem data
    reportDict: dict
        dictionary that contains information on simulation that can be used for report generation
    cfg: configparser object
        configuration object for simulation to be performed
    infoDict['tCPU']: dict
        info on cpu timing
    particlesList: list
        list of particle dictionaries for all saving time steps
    """

    # select release area input data according to chosen release scenario
    inputSimFiles = gI.selectReleaseFile(inputSimFiles, cfg['INPUT']['releaseScenario'])

    # create required input from input files
    demOri, inputSimLines = prepareInputData(inputSimFiles, cfg)

    if cfg['GENERAL'].getboolean('iniStep'):
        # append buffered release Area
        inputSimLines = pI.createReleaseBuffer(cfg, inputSimLines)

    # set thickness values for the release area, entrainment and secondary release areas
    relName, inputSimLines, badName = prepareReleaseEntrainment(cfg, inputSimFiles['releaseScenario'], inputSimLines)

    log.debug('Perform %s simulation' % cuSimName)

    # +++++++++PERFORM SIMULAITON++++++++++++++++++++++
    # for timing the sims
    startTime = time.time()

    # initialize particles, fields, dem
    particles, fields, dem, reportAreaInfo = initializeSimulation(cfg, outDir, demOri, inputSimLines, cuSimName)

    # ------------------------
    #  Start time step computation
    Tsave, particlesList, fieldsList, infoDict = DFAIterate(cfg, particles, fields, dem, inputSimLines, simHash=simHash)

    # write mass balance to File
    writeMBFile(infoDict, avaDir, cuSimName)

    tCPUDFA = '%.2f' % (time.time() - startTime)
    log.info(('cpu time DFA = %s s' % (tCPUDFA)))

    # export particles dictionaries of saving time steps
    # (if particles is not in resType, only first and last time step are saved)
    outDirData = outDir / 'particles'
    fU.makeADir(outDirData)
    savePartToPickle(particlesList, outDirData, cuSimName)

    # export particles properties for visulation
    if cfg['VISUALISATION'].getboolean('writePartToCSV'):
        particleTools.savePartToCsv(cfg['VISUALISATION']['visuParticleProperties'], particlesList, outDir)

    # write report dictionary
    reportDict = createReportDict(avaDir, cuSimName, relName, inputSimLines, cfg, reportAreaInfo)
    # add time and mass info to report
    reportDict = reportAddTimeMassInfo(reportDict, tCPUDFA, infoDict)

    # Result parameters to be exported
    if cfg['EXPORTS'].getboolean('exportData'):
        exportFields(cfg, Tsave, fieldsList, dem, outDir, cuSimName)
    else:
        # fetch contourline info
        contDictXY = outCom1DFA.fetchContCoors(dem['header'],
            fieldsList[-1][cfg['VISUALISATION']['contourResType']], cfg['VISUALISATION'], cuSimName)
        reportDict['contours'] = contDictXY

    if particlesList[-1]['nExitedParticles'] != 0.:
        log.warning('%d particles have been removed during simulation because they exited the domain' %  particlesList[-1]['nExitedParticles'])

    cfgTrackPart = cfg['TRACKPARTICLES']
    # track particles
    if cfgTrackPart.getboolean('trackParticles'):
        particlesList, trackedPartProp, track = trackParticles(cfgTrackPart, dem, particlesList)
        if track:
            outDirData = outDir / 'particles'
            fU.makeADir(outDirData)
            outCom1DFA.plotTrackParticle(outDirData, particlesList, trackedPartProp, cfg, dem, cuSimName)

    return dem, reportDict, cfg, infoDict['tCPU'], particlesList


def prepareReleaseEntrainment(cfg, rel, inputSimLines):
    """ set thickness values for release, secondary release and entrainment
        set flag to append _AF to release scenario name if it includes an underscore

    Parameters
    ----------
    cfg : dict
        configuration parameters - keys: relTh, secRelArea, secondaryRelTh
    rel : str
        path to release file
    inputSimLines: dict
        dictionary with dictionaries with input data infos (releaseLine, entLine, ...)

    Returns
    -------
    relName : str
        release name
    inputSimLines : dict
        dictionary with dictionaries with input data infos now updated with thickness values (releaseLine, entLine, ...)
    badName : boolean
        changed release name
    """

    # Set release areas and release thickness
    relName = rel.stem
    badName = False
    if '_' in relName:
        badName = True
        log.warning('Release area scenario file name includes an underscore \
        the suffix _AF will be added for the simulation name')

    # set release thickness
    if cfg['GENERAL'].getboolean('relThFromFile') is False:
        releaseLine = setThickness(cfg, inputSimLines['releaseLine'], 'relTh')
        inputSimLines['releaseLine'] = releaseLine
    log.debug('Release area scenario: %s - perform simulations' % (relName))

    if cfg['GENERAL'].getboolean('iniStep'):
        # set release thickness for buffer
        releaseLineBuffer = setThickness(cfg, inputSimLines['releaseLineBuffer'], 'relTh')
        inputSimLines['releaseLineBuffer'] = releaseLineBuffer

    if cfg.getboolean('GENERAL', 'secRelArea'):
        secondaryReleaseLine = setThickness(cfg, inputSimLines['secondaryReleaseLine'], 'secondaryRelTh')
    else:
        inputSimLines['entResInfo']['flagSecondaryRelease'] = 'No'
        secondaryReleaseLine = None
    inputSimLines['secondaryReleaseLine'] = secondaryReleaseLine

    if cfg['GENERAL']['simTypeActual'] in ['ent', 'entres']:
        # set entrainment thickness
        entLine = setThickness(cfg, inputSimLines['entLine'], 'entTh')
        inputSimLines['entLine'] = entLine

    return relName, inputSimLines, badName


def setThickness(cfg, lineTh, typeTh):
    """ set thickness in line dictionary of release, entrainment, second. release

    Parameters
    -----------
    cfg: config parser
        configuration settings
    lineTh: dict
        dictionary with info on line (e.g. release area line)
    typeTh: str
        type of thickness to be set (e.g. relTh for release thickness)

    Returns
    --------
    lineTh: dict
        updated dictionary with new key: thickness and thicknessSource
    """

    # create thickness flag name
    thFlag = typeTh + 'FromShp'
    # set thickness source info
    if cfg['GENERAL'].getboolean(thFlag):
        if cfg['INPUT']['thFromIni'] != '' and typeTh in cfg['INPUT']['thFromIni']:
            lineTh['thicknessSource'] = ['ini file'] * len(lineTh['thickness'])
        else:
            lineTh['thicknessSource'] = ['shp file'] * len(lineTh['thickness'])
    else:
        lineTh['thicknessSource'] = ['ini file'] * len(lineTh['thickness'])

    # set thickness value info read from ini file that has been updated from shp or ini previously
    for count, id in enumerate(lineTh['id']):
        if cfg['GENERAL'].getboolean(thFlag):
            thName = typeTh + id
            lineTh['thickness'][count] = cfg['GENERAL'].getfloat(thName)

        else:
            thName = typeTh
            lineTh['thickness'][count] = cfg['GENERAL'].getfloat(thName)

    return lineTh


def prepareInputData(inputSimFiles, cfg):
    """ Fetch input data

    Parameters
    ----------
    inputSimFiles : dict
        dictionary containing

        - relFile : str, path to release area file
        - demFile : str, path to dem file in Inputs/
        - secondaryReleaseFile : str, path to secondaryRelease file
        - entFiles : str, path to entrainment file
        - resFile : str, path to resistance file
        - entResInfo : flag dict
        flag if Yes entrainment and/or resistance areas found and used for simulation
        flag True if a Secondary Release file found and activated

    cfg: configparser object
        configuration for simType and secondary rel

    Returns
    -------
    demOri : dict
        dictionary with dem info (header original origin), raster data correct mesh cell size
        this dem has been remeshed/read from remeshed if chosen cell size is not equal to cell size
        of DEM in Inputs/

    inputSimLines : dict
        dictionary containing

        - releaseLine : dict, release line dictionary
        - secondaryReleaseLine : dict, secondaryRelease line dictionary
        - entLine : dict, entrainment line dictionary
        - resLine : dict, resistance line dictionary
        - entrainmentArea : str, entrainment file name
        - resistanceArea : str, resistance file name
        - entResInfo : flag dict
        flag if Yes entrainment and/or resistance areas found and used for simulation
        flag True if a Secondary Release file found and activated

    """

    # load data
    entResInfo = inputSimFiles['entResInfo'].copy()
    relFile = inputSimFiles['releaseScenario']
    relThFile = inputSimFiles['relThFile']

    # get dem dictionary - already read DEM with correct mesh cell size
    demOri = gI.initializeDEM(cfg['GENERAL']['avalancheDir'], demPath=cfg['INPUT']['DEM'])

    # read data from relThFile
    if relThFile != None and cfg['GENERAL'].getboolean('relThFromFile'):
        relThField = IOf.readRaster(relThFile)
        relThFieldData = relThField['rasterData']
        relThFieldDataOrig = relThFieldData.copy()
        if demOri['header']['ncols'] != relThField['header']['ncols'] or \
            demOri['header']['nrows'] != relThField['header']['nrows']:

            message = ('Release thickness field read from %s does not match the number of rows and columns of the dem'
                % inputSimFiles['relThFile'])
            log.error(message)
            raise AssertionError(message)
        elif np.isnan(relThFieldData).any() == True:
            message = ('Release thickness field contains nans - not allowed no release thickness must be set to 0')
            log.error(message)
            raise AssertionError(message)
    else:
        relThFieldData = ''

    # get line from release area polygon
    releaseLine = shpConv.readLine(relFile, 'release1', demOri)
    releaseLine['file'] = relFile
    releaseLine['type'] = 'Release'

    # get line from secondary release area polygon
    if cfg['GENERAL'].getboolean('secRelArea'):
        if entResInfo['flagSecondaryRelease'] == 'Yes':
            secondaryReleaseFile = inputSimFiles['secondaryReleaseFile']
            secondaryReleaseLine = shpConv.readLine(secondaryReleaseFile, '', demOri)
            secondaryReleaseLine['fileName'] = secondaryReleaseFile
            secondaryReleaseLine['type'] = 'Secondary release'
        else:
            message = 'No secondary release file found'
            log.error(message)
            raise FileNotFoundError(message)
    else:
        secondaryReleaseLine = None
        # set False
        entResInfo['flagSecondaryRelease'] = 'No'

    # get line from entrainement area polygon
    if cfg['GENERAL']['simTypeActual'] in ['ent', 'entres']:
        entFile = inputSimFiles['entFile']
        entLine = shpConv.readLine(entFile, '', demOri)
        entrainmentArea = entFile.name
        entLine['fileName'] = entFile
        entLine['type'] = 'Entrainment'
    else:
        entLine = None
        entrainmentArea = ''

    # get line from resistance area polygon
    if cfg['GENERAL']['simTypeActual'] in ['entres', 'res']:
        resFile = inputSimFiles['resFile']
        resLine = shpConv.readLine(resFile, '', demOri)
        resistanceArea = resFile.name
        resLine['fileName'] = resFile
        resLine['type'] = 'Resistance'
    else:
        resLine = None
        resistanceArea = ''

    # get line from dam
    if cfg['GENERAL'].getboolean('dam'):
        if entResInfo['dam'] == 'Yes':
            damFile = inputSimFiles['damFile']
            damLine = shpConv.readLine(damFile, '', demOri)
            damLine['fileName'] = [damFile]
            damLine['type'] = 'Dam'
        else:
            message = 'Take dam is set, but no dam file found, skipping it'
            log.debug(message)
            damLine = None
    else:
        damLine = None

    inputSimLines = {'releaseLine': releaseLine, 'secondaryReleaseLine': secondaryReleaseLine,
                     'entLine': entLine, 'resLine': resLine, 'damLine': damLine, 'entrainmentArea': entrainmentArea,
                     'resistanceArea': resistanceArea, 'entResInfo': entResInfo,
                     'relThField': relThFieldData}

    return demOri, inputSimLines


def createReportDict(avaDir, logName, relName, inputSimLines, cfg, reportAreaInfo):
    """ create simulaton report dictionary

    Parameters
    ----------
    logName : str
        simulation scenario name
    relName : str
        release name
    relDict : dict
        release dictionary
    cfg : configparser
        simulation configuration
    entrainmentArea : str
        entrainment file name
    resistanceArea : str
        resistance file name

    Returns
    -------
    reportST : dict
        simulation scenario dictionary
    """

    # load parameters
    dateTimeInfo = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    entInfo = reportAreaInfo['entrainment']
    resInfo = reportAreaInfo['resistance']
    entrainmentArea = inputSimLines['entrainmentArea']
    resistanceArea = inputSimLines['resistanceArea']
    relDict = inputSimLines['releaseLine']
    secRelAreaInfo = reportAreaInfo['secRelArea']
    if secRelAreaInfo == 'No':
        secRelAreaFlag = 'No'
    else:
        secRelAreaFlag = 'Yes'

    # Get default cfg and convert to dict for comparison
    cfgGen = cfg['GENERAL']

    # Create dictionary
    reportST = {}
    reportST = {'headerLine': {'type': 'title', 'title': 'com1DFA Simulation'},
                'avaName': {'type': 'avaName', 'name': str(avaDir)},
                'simName': {'type': 'simName', 'name': logName},
                'time': {'type': 'time', 'time': dateTimeInfo},
                'Simulation Parameters': {
                'type': 'list',
                'Program version': getVersion(),
                'Parameter set': 'Default',
                'Release Area Scenario': relName,
                'Entrainment': entInfo,
                'Resistance': resInfo,
                'Secondary release area': secRelAreaFlag,
                'Density [kgm-3]': cfgGen['rho'],
                'Friction model': cfgGen['frictModel']},
                'Release Area': {'type': 'columns', 'Release area scenario': relName, 'Release Area': relDict['Name'],
                                 'Release thickness [m]': relDict['thickness']}}

    # add frict parameters
    if cfgGen['frictModel'].lower() == 'samosat':
        reportST['Friction model'] = {'type': 'columns', 'model': 'samosAT', 'mu': cfgGen['musamosat'], 'tau0': cfgGen['tau0samosat'],
            'Rs0': cfgGen['Rs0samosat'], 'kappa': cfgGen['kappasamosat'], 'R':  cfgGen['Rsamosat'], 'B' :  cfgGen['Bsamosat']}
    elif cfgGen['frictModel'].lower() == 'samosatsmall':
        reportST['Friction model'] = {'type': 'columns', 'model': 'samosATSmall', 'mu': cfgGen['musamosatsmall'], 'tau0': cfgGen['tau0samosatsmall'],
            'Rs0': cfgGen['Rs0samosatsmall'], 'kappa': cfgGen['kappasamosatsmall'], 'R':  cfgGen['Rsamosatsmall'], 'B' :  cfgGen['Bsamosatsmall']}
    elif cfgGen['frictModel'].lower() == 'samosatmedium':
        reportST['Friction model'] = {'type': 'columns', 'model': 'samosATMedium', 'mu': cfgGen['musamosatmedium'], 'tau0': cfgGen['tau0samosatmedium'],
            'Rs0': cfgGen['Rs0samosatmedium'], 'kappa': cfgGen['kappasamosatmedium'], 'R':  cfgGen['Rsamosatmedium'], 'B' :  cfgGen['Bsamosatmedium']}
    elif cfgGen['frictModel'].lower() == 'voellmy':
        reportST['Friction model'] = {'type': 'columns', 'model': 'Voellmy', 'mu': cfgGen['muvoellmy'], 'xsi': cfgGen['xsivoellmy']}
    elif cfgGen['frictModel'].lower() == 'coulomb':
        reportST['Friction model'] = {'type': 'columns', 'model': 'Coulomb', 'mu': cfgGen['mucoulomb']}
    elif cfgGen['frictModel'].lower() == 'wetsnow':
        reportST['Friction model'] = {'type': 'columns', 'model': 'wetsnow', 'mu': cfgGen['mu0wetsnow'], 'xsi':  cfgGen['mu0wetsnow']}

    # check if secondary release area
    if secRelAreaFlag == 'Yes':
        reportST['Secondary release Area'] = secRelAreaInfo

    # Check if parameter set is modified from default, and add section to report
    if '_C_' in logName:
        reportST['Simulation Parameters']['Parameter set'] = 'Changed'
        reportST.update({'Parameters changed from default': {
                        'type': 'list'
                         }})
        cfgDict = cfgUtils.convertConfigParserToDict(cfg)
        _, changedVals = com1DFATools.compareSimCfgToDefaultCfgCom1DFA(cfgDict)

        for key, val in changedVals.items():
            # Format key string for better readability (without loosing the parameter path)
            keyStr = key.replace('root','')
            keyStr = keyStr.replace("\'][\'","->")
            keyStr = keyStr.replace("[\'","")
            keyStr = keyStr.replace("\']","")

            valStr = val['new_value'] + ' (default is ' + val['old_value'] +')'

            reportST['Parameters changed from default'][keyStr] = valStr

    if entInfo == 'Yes':
        entDict = inputSimLines['entLine']
        reportST.update({'Entrainment area':
                         {'type': 'columns',
                             'Entrainment area scenario': entrainmentArea,
                             'Entrainment thickness [m]': entDict['thickness'],
                             'Entrainment density [kgm-3]': cfgGen['rhoEnt']}})
    if resInfo == 'Yes':
        reportST.update({'Resistance area': {'type': 'columns', 'Resistance area scenario': resistanceArea}})

    reportST['Release Area'].update(reportAreaInfo['Release area info'])

    return reportST


def reportAddTimeMassInfo(reportDict, tCPUDFA, infoDict):
    """ Add time and mass info to report """

    # add mass info
    reportDict['Simulation Parameters'].update({'Initial mass [kg]': ('%.2f' % infoDict['initial mass'])})
    reportDict['Simulation Parameters'].update({'Final mass [kg]': ('%.2f' % infoDict['final mass'])})
    reportDict['Simulation Parameters'].update({'Entrained mass [kg]': ('%.2f' % infoDict['entrained mass'])})
    reportDict['Simulation Parameters'].update({'Entrained volume [m3]': ('%.2f' % infoDict['entrained volume'])})

    # add stop info
    reportDict['Simulation Parameters'].update(infoDict['stopInfo'])

    # add computation time to report dict
    reportDict['Simulation Parameters'].update({'Computation time [s]': tCPUDFA})

    return reportDict


def initializeMesh(cfg, demOri, num):
    """ Create rectangular mesh

    Reads the DEM information, computes the normal vector field and
    boundries to the DEM. Also generates the grid for the neighbour search

    Parameters
    ----------
    demOri : dict
        dictionary with initial dem information
    num : int
        chose between 4, 6 or 8 (using then 4, 6 or 8 triangles) or
        1 to use the simple cross product method

    Returns
    -------
    dem : dict
        dictionary relocated in (0,0) and completed with normal field and
        boundaries as well as neighbour search grid information
    """

    # set origin to 0, 0 for computations, store original origin
    dem = setDEMoriginToZero(demOri)
    dem['originalHeader'] = demOri['header'].copy()

    # read dem header
    headerDEM = dem['header']
    nColsDEM = headerDEM['ncols']
    nRowsDEM = headerDEM['nrows']
    cszDEM = headerDEM['cellsize']

    # get normal vector of the grid mesh
    dem = DFAtls.getNormalMesh(dem, num)

    # Prepare SPH grid
    headerNeighbourGrid = {}
    cszNeighbourGrid = cfg.getfloat('sphKernelRadius')
    headerNeighbourGrid['cellsize'] = cszNeighbourGrid
    headerNeighbourGrid['ncols'] = np.ceil(nColsDEM * cszDEM / cszNeighbourGrid)
    headerNeighbourGrid['nrows'] = np.ceil(nRowsDEM * cszDEM / cszNeighbourGrid)
    headerNeighbourGrid['xllcenter'] = 0
    headerNeighbourGrid['yllcenter'] = 0
    dem['headerNeighbourGrid'] = headerNeighbourGrid

    # get real Area
    dem = DFAtls.getAreaMesh(dem, num)
    projArea = nColsDEM * nRowsDEM * cszDEM * cszDEM
    areaRaster = dem['areaRaster']
    log.debug('Largest cell area: %.2f m²' % (np.nanmax(areaRaster)))
    log.debug('Projected Area : %.2f' % projArea)
    log.debug('Total Area : %.2f' % np.nansum(areaRaster))

    return dem


def setDEMoriginToZero(demOri):
    """ set origin of DEM to 0,0 """

    dem = copy.deepcopy(demOri)
    dem['header']['xllcenter'] = 0
    dem['header']['yllcenter'] = 0

    return dem


def initializeSimulation(cfg, outDir, demOri, inputSimLines, logName):
    """ create simulaton report dictionary

    Parameters
    ----------
    cfg : str
        simulation scenario name
    outDir : pathlib path
        path to output directory (to save the dam foot line if a dam is used)
    demOri : dict
        dictionary with original dem
    inputSimLines : dict
        releaseLine : dict
            release line dictionary
        releaseLineBuffer: dict
            release line buffer dictionary - optional if iniStep True
        secondaryReleaseLine : dict
            secondary release line dictionary
        entLine : dict
            entrainment line dictionary
        resLine : dict
            resistance line dictionary
    logName : str
        simulation scenario name

    Returns
    -------
    particles : dict
        particles dictionary at initial time step
        list of secondary release particles to be used
    fields : dict
        fields dictionary at initial time step
    dem : dict
        dictionary with new dem (lower left center at origin)
    """
    cfgGen = cfg['GENERAL']
    methodMeshNormal = cfg.getfloat('GENERAL', 'methodMeshNormal')
    thresholdPointInPoly = cfgGen.getfloat('thresholdPointInPoly')
    relThField = inputSimLines['relThField']

    # -----------------------
    # Initialize mesh
    log.debug('Initializing Mesh')
    dem = initializeMesh(cfgGen, demOri, methodMeshNormal)

    # ------------------------
    log.debug('Initializing main release area')
    # process release info to get it as a raster
    if cfg['GENERAL'].getboolean('iniStep'):
        releaseLine = inputSimLines['releaseLineBuffer']
        releaseLineReal = inputSimLines['releaseLine']
        # check if release features overlap between features
        prepareArea(releaseLineReal, dem, thresholdPointInPoly, combine=True, checkOverlap=True)
        buffer1 = (cfg['GENERAL'].getfloat('sphKernelRadius') * cfg['GENERAL'].getfloat('additionallyFixedFactor')
                   * cfg['GENERAL'].getfloat('bufferZoneFactor'))
        if len(relThField) == 0:
            # if no release thickness field or function - set release according to shapefile or ini file
            # this is a list of release rasters that we want to combine
            releaseLineReal = prepareArea(releaseLineReal, dem, buffer1, thList=releaseLineReal['thickness'],
                                          combine=True, checkOverlap=False)
        else:
            # if relTh provided - set release thickness with field or function
            releaseLineReal = prepareArea(releaseLineReal, dem, buffer1, combine=True, checkOverlap=False)

    else:
        releaseLine = inputSimLines['releaseLine']
        # check if release features overlap between features
        prepareArea(releaseLine, dem, thresholdPointInPoly, combine=True, checkOverlap=True)

    if len(relThField) == 0:
        # if no release thickness field or function - set release according to shapefile or ini file
        # this is a list of release rasters that we want to combine
        releaseLine = prepareArea(releaseLine, dem, np.sqrt(2), thList=releaseLine['thickness'],
                                  combine=True, checkOverlap=False)
    else:
        # if relTh provided - set release thickness with field or function
        releaseLine = prepareArea(releaseLine, dem, np.sqrt(2), combine=True, checkOverlap=False)

    # plot release area scenario
    outCom1DFA.plotReleaseScenarioView(cfgGen['avalancheDir'], releaseLine, dem,
                                       ('Release Scenario %s' % inputSimLines['releaseLine']['file'].stem), logName)
    # compute release area
    header = dem['header']
    csz = header['cellsize']
    relRaster = releaseLine['rasterData']
    # for area computation use smaller threshold to identify raster cells that lie within release line
    # as for creating particles a bigger radius is chosen as particles that lie outside are removed afterwards
    releaseLineArea = releaseLine.copy()
    relAreaActualList, relAreaProjectedList, _ = gI.computeAreasFromRasterAndLine(releaseLineArea, dem)
    relAreaProjected = np.sum(relAreaProjectedList)
    relAreaActual = np.sum(relAreaActualList)
    reportAreaInfo = {'Release area info': {'Projected Area [m2]': '%.2f' % (relAreaProjected),
                                            'Actual Area [m2]': '%.2f' % (relAreaActual)}}

    # ------------------------
    # initialize simulation
    # create primary release area particles and fields
    releaseLine['header'] = dem['originalHeader']
    inputSimLines['releaseLine']['header'] = dem['originalHeader']
    particles = initializeParticles(cfgGen, releaseLine, dem, inputSimLines=inputSimLines,
                                    logName=logName, relThField=relThField)
    particles, fields = initializeFields(cfg, dem, particles, releaseLine)

    # initialize Dam
    damLine = inputSimLines['damLine']
    # FSO: disabled writing of damFootLine inside the initializeWallLines.
    # This is not threadsafe and leads to errors on multiprocessing
    # damFootLinePath = outDir / 'dam' / 'damFootLine.shp'
    # damLine = damCom1DFA.initializeWallLines(cfgGen, dem, damLine, damFootLinePath)
    damLine = damCom1DFA.initializeWallLines(cfgGen, dem, damLine, '')
    dem['damLine'] = damLine

    # perform initialisation step for redistributing particles
    if cfg['GENERAL'].getboolean('iniStep'):
        startTimeIni = time.time()
        particles, fields = pI.getIniPosition(cfg, particles, dem, fields, inputSimLines, relThField)
        tIni = time.time() - startTimeIni
        log.info('Ini step for initialising particles finalized, total mass: %.2f, number of particles: %d' %
                 (np.sum(particles['m']), particles['nPart']))
        log.debug('Time needed for ini step: %.2f s' % (tIni))
    # ------------------------
    # process secondary release info to get it as a list of rasters
    secondaryReleaseInfo, reportAreaInfo = initializeSecRelease(inputSimLines, dem, relRaster, reportAreaInfo)

    particles['secondaryReleaseInfo'] = secondaryReleaseInfo

    # initialize entrainment and resistance
    # get info of simType and whether or not to initialize resistance and entrainment
    simTypeActual = cfgGen['simTypeActual']
    entrMassRaster, entrEnthRaster, reportAreaInfo = initializeMassEnt(dem, simTypeActual, inputSimLines['entLine'],
                                                                       reportAreaInfo, thresholdPointInPoly,
                                                                       cfgGen)

    # check if entrainment and release overlap
    entrMassRaster = geoTrans.checkOverlap(entrMassRaster, relRaster, 'Entrainment', 'Release', crop=True)
    entrEnthRaster = geoTrans.checkOverlap(entrEnthRaster, relRaster, 'Entrainment', 'Release', crop=True)
    # check for overlap with the secondary release area
    if secondaryReleaseInfo['flagSecondaryRelease'] == 'Yes':
        for secRelRaster in secondaryReleaseInfo['rasterData']:
            entrMassRaster = geoTrans.checkOverlap(entrMassRaster, secRelRaster, 'Entrainment', 'Secondary release ',
                                                   crop=True)
            entrEnthRaster = geoTrans.checkOverlap(entrEnthRaster, secRelRaster, 'Entrainment', 'Secondary release ',
                                                   crop=True)
    # surfacic entrainment mass available (unit kg/m²)
    fields['entrMassRaster'] = entrMassRaster
    fields['entrEnthRaster'] = entrEnthRaster
    entreainableMass = np.nansum(fields['entrMassRaster']*dem['areaRaster'])
    log.info('Mass available for entrainment: %.2f kg' % (entreainableMass))

    log.debug('Initializing resistance area')
    cResRaster, reportAreaInfo = initializeResistance(cfgGen, dem, simTypeActual, inputSimLines['resLine'],
                                                      reportAreaInfo, thresholdPointInPoly)
    fields['cResRaster'] = cResRaster
    return particles, fields, dem, reportAreaInfo


def initializeParticles(cfg, releaseLine, dem, inputSimLines='', logName='', relThField=''):
    """ Initialize DFA simulation

    Create particles and fields dictionary according to config parameters
    release raster and dem

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    releaseLine: dict
        dictionary with info on release
    dem : dict
        dictionary with dem information
    inputSimLines: dictionary
        info on input files; real releaseline info required for iniStep
    relThField: 2D numpy array
        if the release thickness is not uniform, give here the releaseRaster

    Returns
    -------
    particles : dict
        particles dictionary at initial time step
    fields : dict
        fields dictionary at initial time step
    """

    # get simulation parameters
    rho = cfg.getfloat('rho')
    gravAcc = cfg.getfloat('gravAcc')
    cpIce = cfg.getfloat('cpIce')
    TIni = cfg.getfloat('TIni')
    avaDir = cfg['avalancheDir']
    massPerParticleDeterminationMethod = cfg['massPerParticleDeterminationMethod']
    interpOption = cfg.getfloat('interpOption')

    # read dem header
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    csz = header['cellsize']
    # if the release is not constant but given by a varying function, we need both the mask giving the cells
    # to be initialized and the raster giving the flow thickness value
    relRasterMask = releaseLine['rasterData']
    if len(relThField) == 0:
        relRaster = releaseLine['rasterData']
    else:
        log.info('Release thickness read from relThFile')
        relRaster = relThField
    areaRaster = dem['areaRaster']

    # get the initialization method used
    relThForPart = getRelThFromPart(cfg, releaseLine, relThField)
    massPerPart, nPPK = com1DFATools.getPartInitMethod(cfg, csz, relThForPart)

    # initialize arrays
    partPerCell = np.zeros(np.shape(relRaster), dtype=np.int64)
    # find all non empty cells (meaning release area)
    indRelY, indRelX = np.nonzero(relRasterMask)
    if inputSimLines != '':
        indRelYReal, indRelXReal = np.nonzero(inputSimLines['releaseLine']['rasterData'])
    else:
        indRelYReal, indRelXReal = np.nonzero(relRaster)
    iReal = list(zip(indRelYReal, indRelXReal))

    # get approximate ratio between projected and real release area
    # because relRasterMask has a none 0 value where the release is but we want a 1 there
    realArea = np.sum(areaRaster * np.where(relRasterMask > 0, 1, 0))
    projectedArea = csz*csz*np.size(indRelY)
    ratioArea = projectedArea / realArea

    # make option available to read initial particle distribution from file
    if cfg.getboolean('initialiseParticlesFromFile'):
        if cfg.getboolean('iniStep'):
            message = 'If initialiseParticlesFromFile is used, iniStep cannot be performed - chose only one option'
            log.error(message)
            raise AssertionError(message)
        particles, hPartArray = particleTools.initialiseParticlesFromFile(cfg, avaDir, releaseLine['file'].stem)
    else:
        # initialize random generator
        rng = np.random.default_rng(cfg.getint('seed'))

        nPart = 0
        xPartArray = np.empty(0)
        yPartArray = np.empty(0)
        mPartArray = np.empty(0)
        aPartArray = np.empty(0)
        hPartArray = np.empty(0)
        idFixed = np.empty(0)
        if len(relThField) != 0 and cfg.getboolean('iniStep'):
            # set release thickness to a constant value for initialisation
            relRaster = np.where(relRaster > 0., cfg.getfloat('relTh'), 0.)
            log.warning('relThField!= 0, but relRaster set to relTh value (from ini)')
        # loop on non empty cells
        for indRelx, indRely in zip(indRelX, indRelY):
            # compute number of particles for this cell
            hCell = relRaster[indRely, indRelx]
            aCell = areaRaster[indRely, indRelx]
            xPart, yPart, mPart, n, aPart = particleTools.placeParticles(hCell, aCell, indRelx, indRely, csz,
                                                                         massPerPart, nPPK, rng, cfg, ratioArea)
            nPart = nPart + n
            partPerCell[indRely, indRelx] = n
            # initialize particles position, mass, height...
            xPartArray = np.append(xPartArray, xPart)
            yPartArray = np.append(yPartArray, yPart)
            mPartArray = np.append(mPartArray, mPart * np.ones(n))
            aPartArray = np.append(aPartArray, aPart * np.ones(n))
            if (indRely, indRelx) in iReal:
                idFixed = np.append(idFixed, np.zeros(n))
            else:
                idFixed = np.append(idFixed, np.ones(n))

        hPartArray = DFAtllsC.projOnRaster(xPartArray, yPartArray, relRaster, csz, ncols, nrows, interpOption)
        hPartArray = np.asarray(hPartArray)
        # for the MPPKR option use hPart and aPart to define the mass of the particle (this means, within a cell
        # partticles have the same area but may have different flow thickness which means a different mass)
        if massPerParticleDeterminationMethod == 'MPPKR':
            mPartArray = rho * aPartArray * hPartArray
        # create dictionnary to store particles properties
        particles = {}
        particles['nPart'] = nPart
        particles['x'] = xPartArray
        particles['y'] = yPartArray
        # adding z component
        particles, _ = geoTrans.projectOnRaster(dem, particles, interp='bilinear')
        particles['m'] = mPartArray
        particles['idFixed'] = idFixed
    # initialize enthalpy
    particles['totalEnthalpy'] = TIni * cpIce + gravAcc * particles['z']

    particles['massPerPart'] = massPerPart
    particles['mTot'] = np.sum(particles['m'])
    particles['h'] = hPartArray
    particles['ux'] = np.zeros(np.shape(hPartArray))
    particles['uy'] = np.zeros(np.shape(hPartArray))
    particles['uz'] = np.zeros(np.shape(hPartArray))
    particles['uAcc'] = np.zeros(np.shape(hPartArray))
    particles['velocityMag'] = np.zeros(np.shape(hPartArray))
    particles['trajectoryLengthXY'] = np.zeros(np.shape(hPartArray))
    particles['trajectoryLengthXYCor'] = np.zeros(np.shape(hPartArray))
    particles['trajectoryLengthXYZ'] = np.zeros(np.shape(hPartArray))
    particles['trajectoryAngle'] = np.zeros(np.shape(hPartArray))
    particles['stoppCriteria'] = False
    mPartArray = particles['m']
    kineticEne = np.sum(0.5 * mPartArray * DFAtls.norm2(particles['ux'], particles['uy'], particles['uz']))
    particles['kineticEne'] = kineticEne
    particles['potentialEne'] = np.sum(gravAcc * mPartArray * particles['z'])
    particles['peakKinEne'] = kineticEne
    particles['peakForceSPH'] = 0.0
    particles['forceSPHIni'] = 0.0
    particles['peakMassFlowing'] = 0
    particles['simName'] = logName
    particles['xllcenter'] = dem['originalHeader']['xllcenter']
    particles['yllcenter'] = dem['originalHeader']['yllcenter']
    particles['nExitedParticles'] = 0.

    # remove particles that might lay outside of the release polygon
    if not cfg.getboolean('iniStep') and not cfg.getboolean('initialiseParticlesFromFile'):
        particles = checkParticlesInRelease(particles, releaseLine, cfg.getfloat('thresholdPointInPoly'))

    # add a particles ID:
    # integer ranging from 0 to nPart in the initialisation.
    # Everytime that a new particle is created, it gets a new ID that is > nID
    # where nID is the number of already used IDs
    # (enable tracking of particles even if particles are added or removed)
    # unique identifier for each particle
    particles['ID'] = np.arange(particles['nPart'])
    # keep track of the identifier (usefull to add identifier to newparticles)
    particles['nID'] = particles['nPart']
    # keep track of parents (usefull for new particles created after splitting)
    particles['parentID'] = np.arange(particles['nPart'])

    # initialize time
    t = 0
    particles['t'] = t

    relCells = np.size(indRelY)
    partPerCell = particles['nPart']/relCells

    if massPerParticleDeterminationMethod != 'MPPKR':
        # we need to set the nPPK
        aTot = np.sum(particles['m'] / (rho * particles['h']))
        # average number of particles per kernel radius
        nPPK = particles['nPart'] * math.pi * csz**2 / aTot
    particles['nPPK'] = nPPK

    log.info('Initialized particles. MTot = %.2f kg, %s particles in %.2f cells.' %
             (particles['mTot'], particles['nPart'], relCells))
    log.info('Mass per particle = %.2f kg and particles per cell = %.2f.' %
             (particles['mTot']/particles['nPart'], partPerCell))

    if debugPlot:
        debPlot.plotPartIni(particles, dem)

    return particles


def getRelThFromPart(cfg, releaseLine, relThField):
    """ get release thickness for initialising particles - use max value

        Parameters
        -----------
        cfg: configparser object
            configuration settings
        releaseLine: dict
            info on releaseline (thickness)
        relThField: numpy array or str
            release thickness field if used, else empty string

        Returns
        --------
        relThForPart: float
            max value of release thickness
    """

    if len(relThField) != 0:
        relThForPart = np.amax(relThField)
    elif cfg.getboolean('relThFromShp'):
        relThForPart = np.amax(np.asarray(releaseLine['thickness'], dtype=float))
    else:
        relThForPart = cfg.getfloat('relTh')

    return relThForPart


def initializeFields(cfg, dem, particles, releaseLine):
    """Initialize fields and update particles flow thickness. Eventually build bond array if snowSlide is activated

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    dem : dict
        dictionary with dem information
    particles : dict
        particles dictionary at initial time step

    Returns
    -------
    particles : dict
        particles dictionary at initial time step updated with the flow thickness
    fields : dict
        fields dictionary at initial time step
    releaseLine : dict
        release line dictionary
    """
    # read config
    cfgGen = cfg['GENERAL']
    # what result types are desired as output (we need this to decide which fields we actually need to compute)
    resTypes = fU.splitIniValueToArraySteps(cfgGen['resType'])
    resTypesReport = fU.splitIniValueToArraySteps(cfg['REPORT']['plotFields'])
    resTypesLast = list(set(resTypes + resTypesReport))
    # read dem header
    header = dem['header']
    ncols = header['ncols']
    nrows = header['nrows']
    # initialize fields
    fields = {}
    fields['pfv'] = np.zeros((nrows, ncols))
    fields['pft'] = np.zeros((nrows, ncols))
    fields['FV'] = np.zeros((nrows, ncols))
    fields['FT'] = np.zeros((nrows, ncols))
    fields['FM'] = np.zeros((nrows, ncols))
    fields['Vx'] = np.zeros((nrows, ncols))
    fields['Vy'] = np.zeros((nrows, ncols))
    fields['Vz'] = np.zeros((nrows, ncols))
    # for optional fields, initialize with dummys (minimum size array). The cython functions then need something
    # even if it is empty to run properly
    if ('TA' in resTypesLast) or ('pta' in resTypesLast):
        fields['pta'] = np.zeros((nrows, ncols))
        fields['TA'] = np.zeros((nrows, ncols))
        fields['computeTA'] = True
        log.debug('Computing Travel Angle')
    else:
        fields['pta'] = np.zeros((1, 1))
        fields['TA'] = np.zeros((1, 1))
        fields['computeTA'] = False
    if ('pke' in resTypesLast):
        fields['pke'] = np.zeros((nrows, ncols))
        fields['computeKE'] = True
        log.debug('Computing Kinetic energy')
    else:
        fields['pke'] = np.zeros((1, 1))
        fields['computeKE'] = False
    if ('P' in resTypesLast) or ('ppr' in resTypesLast):
        fields['P'] = np.zeros((nrows, ncols))
        fields['ppr'] = np.zeros((nrows, ncols))
        fields['computeP'] = True
        log.debug('Computing Pressure')
    else:
        fields['P'] = np.zeros((1, 1))
        fields['ppr'] = np.zeros((1, 1))
        fields['computeP'] = False

    particles = DFAfunC.getNeighborsC(particles, dem)
    particles, fields = DFAfunC.updateFieldsC(cfgGen, particles, dem, fields)

    if cfgGen.getint('snowSlide') == 1:
        # Initialize the bonds between particles if the snowSlide is activated
        # get particles
        x = particles['x']
        y = particles['y']
        # get release outline (in order to remove spurious bonds)
        xOutline = releaseLine['x']-dem['originalHeader']['xllcenter']
        yOutline = releaseLine['y']-dem['originalHeader']['yllcenter']
        # original triangulation (make delaunay triangulation on points)
        triangles = tri.Triangulation(x, y)

        # Cleaning up the triangles (remove unwanted bonds)
        # masking triangles exiting the release (plan small buffer to be sure to keep all inner edges)
        # this happends on non-convex release areas
        outline = sPolygon(zip(xOutline, yOutline)).buffer(cfgGen.getfloat('thresholdPointInPoly'))
        mask = [not outline.contains(sPolygon(zip(x[tri], y[tri])))
                for tri in triangles.get_masked_triangles()]
        triangles.set_mask(mask)
        # masking triangles with sidelength bigger than some threshold
        # (this removes spurious bonds that appear on the sides of the release polygon)
        triang = triangles.triangles
        # Mask off unwanted triangles.
        xtri = x[triang] - np.roll(x[triang], 1, axis=1)
        ytri = y[triang] - np.roll(y[triang], 1, axis=1)
        maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
        if cfgGen['initPartDistType'] == 'triangular':
            # if it is a regular triangular mesh, remove all bonds longer than the triangle side size (np.nanmin(maxi))
            maxRadius = np.nanmin(maxi)*1.1
        else:
            # otherwise remove bonds that are longer than sqrt(5)csz (worst case senario with one particle per cell)
            maxRadius = dem['header']['cellsize']*math.sqrt(5)
        triangles.set_mask(maxi > maxRadius)

        # build the bond array from the triagular mesh (put it to a format that suits us for cython)
        particles = DFAfunC.initializeBondsC(particles, triangles)
        # debugg plot
        if debugPlot:
            debPlot.plotBondsSnowSlideFinal(cfg, particles, dem)
    return particles, fields


def initializeSecRelease(inputSimLines, dem, relRaster, reportAreaInfo):
    """ Initialize secondary release area

    Parameters
    ----------
    inputSimLines : dict
        dict with

        - entResInfo : dict, with the flagSecondaryRelease
        - secondaryReleaseLine : dict, secondary release line dictionary

    dem: dict
        dem dictionary
    relRaster: 2D numpy array
        release Raster (to check overlap)
    reportAreaInfo: dict
        simulation area information dictionary

    Returns
    -------
    secondaryReleaseInfo: dict
        inputSimLines['secondaryReleaseLine'] dictionary completed with

        - header: the dem original header
        - rasterData: list of secondary release rasters (without the overlapping part with the release)
        - flagSecondaryRelease: 'Yes' if a secondary release is there

    reportAreaInfo: dict
        updated simulation area information dictionary

    """
    if inputSimLines['entResInfo']['flagSecondaryRelease'] == 'Yes':
        secondaryReleaseInfo = inputSimLines['secondaryReleaseLine']
        log.info('Initializing secondary release area: %s' % secondaryReleaseInfo['fileName'])
        log.info('Secondary release area features: %s' % (secondaryReleaseInfo['Name']))
        secondaryReleaseInfo['header'] = dem['originalHeader']

        # fetch secondary release areas
        secondaryReleaseInfo = prepareArea(secondaryReleaseInfo, dem, np.sqrt(2),
                                           thList=secondaryReleaseInfo['thickness'], combine=False)
        # remove overlaping parts of the secondary release area with the main release areas
        noOverlaprasterList = []
        for secRelRatser, secRelName in zip(secondaryReleaseInfo['rasterData'], secondaryReleaseInfo['Name']):
            noOverlaprasterList.append(geoTrans.checkOverlap(secRelRatser, relRaster, 'Secondary release ' + secRelName,
                                                             'Release', crop=True))

        secondaryReleaseInfo['flagSecondaryRelease'] = 'Yes'
        # replace the rasterData with noOverlaprasterList (which is the list of rasterData without the overlapping
        # part with the release)
        secondaryReleaseInfo['rasterData'] = noOverlaprasterList
        reportAreaInfo['secRelArea'] = {'type': 'columns',
                                                    'Secondary release area scenario': secondaryReleaseInfo['fileName'].stem,
                                                    'features': secondaryReleaseInfo['Name'].copy(),
                                                    'thickness [m]': secondaryReleaseInfo['thickness'].copy()}
    else:
        secondaryReleaseInfo = {}
        secondaryReleaseInfo['flagSecondaryRelease'] = 'No'
        reportAreaInfo['secRelArea'] = 'No'

    return secondaryReleaseInfo, reportAreaInfo


def initializeMassEnt(dem, simTypeActual, entLine, reportAreaInfo, thresholdPointInPoly, cfgGen):
    """ Initialize mass for entrainment

    Parameters
    ----------
    dem: dict
        dem dictionary
    simTypeActual: str
        simulation type
    entLine: dict
        entrainment line dictionary
    reportAreaInfo: dict
        simulation area information dictionary
    thresholdPointInPoly: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside
    cfgGen: config parser
        General configuration

    Returns
    -------
    entrMassRaster : 2D numpy array
        raster of available mass for entrainment
    reportAreaInfo: dict
        simulation area information dictionary completed with entrainment area info
    """
    # read dem header
    header = dem['originalHeader']
    ncols = header['ncols']
    nrows = header['nrows']
    if 'ent' in simTypeActual:
        entrainmentArea = entLine['fileName']
        log.info('Initializing entrainment area: %s' % (entrainmentArea))
        log.info('Entrainment area features: %s' % (entLine['Name']))
        entLine = prepareArea(entLine, dem, thresholdPointInPoly, thList=entLine['thickness'])
        entrMassRaster = entLine['rasterData']
        # ToDo: not used in samos but implemented
        # tempRaster = cfgGen.getfloat('entTempRef') + (dem['rasterData'] - cfgGen.getfloat('entMinZ')) * cfgGen.getfloat('entTempGrad')
        # entrEnthRaster = np.where(tempRaster < 0, tempRaster*cfgGen.getfloat('cpIce'),
        #                           tempRaster*cfgGen.getfloat('cpWtr')/cfgGen.getfloat('hFusion'))
        entrEnthRaster = np.where(entrMassRaster > 0, cfgGen.getfloat('entTempRef')*cfgGen.getfloat('cpIce'), 0)
        reportAreaInfo['entrainment'] = 'Yes'
    else:
        entrMassRaster = np.zeros((nrows, ncols))
        entrEnthRaster = np.zeros((nrows, ncols))
        reportAreaInfo['entrainment'] = 'No'

    entrMassRaster = entrMassRaster * cfgGen.getfloat('rhoEnt')

    return entrMassRaster, entrEnthRaster, reportAreaInfo


def initializeResistance(cfg, dem, simTypeActual, resLine, reportAreaInfo, thresholdPointInPoly):
    """ Initialize resistance matrix

    Parameters
    ----------
    dem: dict
        dem dictionary
    simTypeActual: str
        simulation type
    resLine: dict
        resistance line dictionary
    reportAreaInfo: dict
        simulation area information dictionary
    thresholdPointInPoly: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside

    Returns
    -------
    cResRaster : 2D numpy array
        raster of resistance coefficients
    reportAreaInfo: dict
        simulation area information dictionary completed with entrainment area info
    """
    d = cfg.getfloat('dRes')
    cw = cfg.getfloat('cw')
    sres = cfg.getfloat('sres')
    # read dem header
    header = dem['originalHeader']
    ncols = header['ncols']
    nrows = header['nrows']
    if simTypeActual in ['entres', 'res']:
        resistanceArea = resLine['fileName']
        log.info('Initializing resistance area: %s' % (resistanceArea))
        log.info('Resistance area features: %s' % (resLine['Name']))
        resLine = prepareArea(resLine, dem, thresholdPointInPoly)
        mask = resLine['rasterData']
        cResRaster = 0.5 * d * cw / (sres*sres) * mask
        reportAreaInfo['resistance'] = 'Yes'
    else:
        cResRaster = np.zeros((nrows, ncols))
        reportAreaInfo['resistance'] = 'No'

    return cResRaster, reportAreaInfo


def DFAIterate(cfg, particles, fields, dem, inputSimLines, simHash=''):
    """ Perform time loop for DFA simulation
     Save results at desired intervals

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at initial time step

        - secondaryReleaseParticles : list, of secondary release area particles dictionaries at initial time step

    fields : dict
        fields dictionary at initial time step
    dem : dict
        dictionary with dem information
    inputSimLines : dict
        dictionary with input data dictionaries (releaseLine, entLine, ...)

    Returns
    -------
    particlesList : list
        list of particles dictionary
    fieldsList : list
        list of fields dictionary (for each time step saved)
    tCPU : dict
        computation time dictionary
    infoDict : dict
        Dictionary of all simulations carried out

    """

    cfgGen = cfg['GENERAL']
    # Initialise cpu timing
    tCPU = {'timeLoop': 0, 'timeForce': 0., 'timeForceSPH': 0., 'timePos': 0., 'timeNeigh': 0., 'timeField': 0.}

    # Load configuration settings
    tEnd = cfgGen.getfloat('tEnd')
    dtSave = fU.splitTimeValueToArrayInterval(cfgGen['tSteps'], tEnd)
    sphOption = cfgGen.getint('sphOption')
    log.debug('using sphOption %s:' % sphOption)
    # desired output fields
    resTypes = fU.splitIniValueToArraySteps(cfgGen['resType'])
    # add particles to the results type if trackParticles option is activated
    if cfg.getboolean('TRACKPARTICLES', 'trackParticles'):
        resTypes = list(set(resTypes + ['particles']))
    # make sure to save all desiered resuts for first and last time step for
    # the report
    resTypesReport = fU.splitIniValueToArraySteps(cfg['REPORT']['plotFields'])
    # always add particles to first and last time step
    resTypesLast = list(set(resTypes + resTypesReport + ['particles']))
    # derive friction type
    # turn friction model into integer
    frictModelsList = ['samosat', 'coulomb', 'voellmy', 'wetsnow', 'samosatsmall', 'samosatmedium', 'voellmyminshear']
    frictModel = cfgGen['frictModel'].lower()
    frictType = frictModelsList.index(frictModel) + 1
    log.debug('Friction Model used: %s, %s' % (frictModelsList[frictType-1], frictType))

    # Initialise Lists to save fields and add initial time step
    particlesList = []
    fieldsList = []
    timeM = []
    massEntrained = []
    massTotal = []

    # setup a result fields info data frame to save max values of fields and avalanche front
    resultsDF = setupresultsDF(resTypesLast, cfg['VISUALISATION'].getboolean('createRangeTimeDiagram'))

    # TODO: add here different time stepping options
    log.debug('Use standard time stepping')
    # Initialize time and counters
    nSave = 1
    tCPU['nSave'] = nSave
    nIter = 1
    nIter0 = 1
    particles['iterate'] = True
    t = particles['t']
    log.debug('Saving results for time step t = %f s', t)
    fieldsList, particlesList = appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypesLast)
    zPartArray0 = copy.deepcopy(particles['z'])

    # create range time diagram
    # check if range-time diagram should be performed, if yes - initialize
    if cfg['VISUALISATION'].getboolean('createRangeTimeDiagram'):
        demRT = dtAna.setDemOrigin(dem)
        mtiInfo, dtRangeTime, cfgRangeTime = dtAna.initializeRangeTime(dtAna, cfg, demRT, simHash)
        # fetch initial time step too
        mtiInfo, dtRangeTime = dtAna.fetchRangeTimeInfo(cfgRangeTime, cfg, dtRangeTime, t,
                                                        demRT['header'], fields, mtiInfo)
        cfgRangeTime['GENERAL']['simHash'] = simHash

    # add initial time step to Tsave array
    Tsave = [0]
    # derive time step for first iteration
    if cfgGen.getboolean('sphKernelRadiusTimeStepping'):
        dtSPHKR = tD.getSphKernelRadiusTimeStep(dem, cfgGen)
        dt = dtSPHKR
    else:
        # get time step
        dt = cfgGen.getfloat('dt')
    particles['dt'] = dt
    t = t + dt

    # Start time step computation
    while t <= tEnd*(1.+1.e-13) and particles['iterate']:
        startTime = time.time()
        log.debug('Computing time step t = %f s, dt = %f s' % (t, dt))
        # Perform computations
        particles, fields, zPartArray0, tCPU = computeEulerTimeStep(cfgGen, particles, fields, zPartArray0, dem, tCPU,
                                                                    frictType)
        # set max values of fields to dataframe
        if cfg['VISUALISATION'].getboolean('createRangeTimeDiagram'):
            rangeValue = mtiInfo['rangeList'][-1]
        else:
            rangeValue = ''
        resultsDF = addMaxValuesToDF(resultsDF, fields, t, resTypesLast, rangeValue=rangeValue)

        tCPU['nSave'] = nSave
        particles['t'] = t

        # write mass balance info
        massEntrained.append(particles['massEntrained'])
        massTotal.append(particles['mTot'])
        timeM.append(t)
        # print progress to terminal
        print("time step t = %f s\r" % t, end="")

        # create range time diagram
        # determine avalanche front and flow characteristics in respective coodrinate system
        if cfg['VISUALISATION'].getboolean('createRangeTimeDiagram') and t >= dtRangeTime[0]:
            mtiInfo, dtRangeTime = dtAna.fetchRangeTimeInfo(cfgRangeTime, cfg, dtRangeTime, t, demRT['header'],
                                                            fields, mtiInfo)

            # create plots for tt diagram animation
            if cfgRangeTime['PLOTS'].getboolean('animate') and cfg['VISUALISATION'].getboolean('TTdiagram'):
                TTResType = cfgRangeTime['GENERAL']['rangeTimeResType']
                dtAnaPlots.animationPlot(demRT, fields[TTResType], demRT['header']['cellsize'], TTResType, cfgRangeTime, mtiInfo, t)

        # make sure the array is not empty
        if t >= (dtSave[0] - 1.e-8):
            Tsave.append(t)
            log.debug('Saving results for time step t = %f s', t)
            log.debug('MTot = %f kg, %s particles' % (particles['mTot'], particles['nPart']))
            log.debug(('cpu time Force = %s s' % (tCPU['timeForce'] / nIter)))
            log.debug(('cpu time ForceSPH = %s s' % (tCPU['timeForceSPH'] / nIter)))
            log.debug(('cpu time Position = %s s' % (tCPU['timePos'] / nIter)))
            log.debug(('cpu time Neighbour = %s s' % (tCPU['timeNeigh'] / nIter)))
            log.debug(('cpu time Fields = %s s' % (tCPU['timeField'] / nIter)))
            fieldsList, particlesList = appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypes)

            # remove saving time steps that have already been saved
            dtSave = updateSavingTimeStep(dtSave, cfg['GENERAL'], t)

            # debugg plot
            if debugPlot:
                debPlot.plotBondsSnowSlideFinal(cfg, particles, dem, inputSimLines)

        # derive time step
        if cfgGen.getboolean('sphKernelRadiusTimeStepping'):
            dt = dtSPHKR
        else:
            # get time step
            dt = cfgGen.getfloat('dt')
        particles['dt'] = dt

        t = t + dt
        nIter = nIter + 1
        nIter0 = nIter0 + 1
        tCPUtimeLoop = time.time() - startTime
        tCPU['timeLoop'] = tCPU['timeLoop'] + tCPUtimeLoop

    tCPU['nIter'] = nIter
    log.info('Ending computation at time t = %f s', t-dt)
    log.debug('Saving results for time step t = %f s', t-dt)
    log.info('MTot = %f kg, %s particles' % (particles['mTot'], particles['nPart']))
    log.debug('Computational performances:')
    log.debug(('cpu time Force = %s s' % (tCPU['timeForce'] / nIter)))
    log.debug(('cpu time ForceSPH = %s s' % (tCPU['timeForceSPH'] / nIter)))
    log.debug(('cpu time Position = %s s' % (tCPU['timePos'] / nIter)))
    log.debug(('cpu time Neighbour = %s s' % (tCPU['timeNeigh'] / nIter)))
    log.debug(('cpu time Fields = %s s' % (tCPU['timeField'] / nIter)))
    log.debug(('cpu time timeLoop = %s s' % (tCPU['timeLoop'] / nIter)))
    log.debug(('cpu time total other = %s s' % ((tCPU['timeForce'] + tCPU['timeForceSPH']
                                               + tCPU['timePos'] + tCPU['timeNeigh']
                                               + tCPU['timeField']) / nIter)))
    Tsave.append(t-dt)

    fieldsList, particlesList = appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypesLast)
    # debugg plot
    if debugPlot:
        debPlot.plotBondsSnowSlideFinal(cfg, particles, dem, inputSimLines)

    # create infoDict for report and mass log file
    infoDict = {'massEntrained': massEntrained, 'timeStep': timeM, 'massTotal': massTotal, 'tCPU': tCPU,
                'final mass': massTotal[-1], 'initial mass': massTotal[0], 'entrained mass': np.sum(massEntrained),
                'entrained volume': (np.sum(massEntrained)/cfgGen.getfloat('rhoEnt'))}

    # determine if stop criterion is reached or end time
    stopCritNotReached = particles['iterate']
    avaTime = particles['t']
    stopCritPer = cfgGen.getfloat('stopCrit') * 100.
    # update info dict with stopping info for report
    if stopCritNotReached:
        infoDict.update({'stopInfo': {'Stop criterion': 'end Time reached: %.2f' % avaTime,
                                      'Avalanche run time [s]': '%.2f' % avaTime}})
    else:
        infoDict.update({'stopInfo': {'Stop criterion': '< %.2f percent of PKE' % stopCritPer,
                                      'Avalanche run time [s]': '%.2f' % avaTime}})

    # create range time diagram
    # export data for range-time diagram
    if cfg['VISUALISATION'].getboolean('createRangeTimeDiagram'):
        lastTimeStep = t - dt
        # first append final time step
        mtiInfo, dtRangeTime = dtAna.fetchRangeTimeInfo(cfgRangeTime, cfg, dtRangeTime, lastTimeStep, demRT['header'],
                                                        fields, mtiInfo)
        dtAna.exportData(mtiInfo, cfgRangeTime, 'com1DFA')

    # save resultsDF to file
    resultsDFPath = pathlib.Path(cfgGen['avalancheDir'], 'Outputs', 'com1DFA', 'resultsDF_%s.csv' % simHash)
    resultsDF.to_csv(resultsDFPath)

    return Tsave, particlesList, fieldsList, infoDict


def setupresultsDF(resTypes, cfgRangeTime):
    """ setup result fields max values dataframe for initial time step
        for all resTypes used and optional for avalanche front

        Parameters
        -----------
        resTypes: list
            list of all resultTypes
        cfgRangeTime: bool
            config info if range time diagram should be performed and rangeList is available

        Returns
        --------
        resultsDF: dataframe
            data frame with on line for iniital time step and max and mean values of fields
    """

    resDict = {'timeStep': [0.0]}
    for resT in resTypes:
        if resT != 'particles':
            resDict['max' + resT] = [0.0]
    if cfgRangeTime:
        resDict['rangeList'] = [0.0]
    resultsDF = pd.DataFrame.from_dict(resDict)
    resultsDF = resultsDF.set_index('timeStep')

    return resultsDF


def addMaxValuesToDF(resultsDF, fields, timeStep, resTypes, rangeValue=''):
    """ add max values of peakFields to dataframe and optionally rangeValue

        Parameters
        -----------
        fields: dict
            dict with all result type fields
        resultsDF: dataframe
            data frame with on line for each time step and max and mean values of fields
        timeStep: float
            computation time step
        resTypes: list
            list of all resultTypes
        rangeValue: float
            avalanche front location -optional

        Returns
        --------
        resultsDF: data frame
            updated data frame
    """

    newLine = []
    for resT in resTypes:
        if resT != 'particles':
            newLine.append(np.nanmax(fields[resT]))

    if rangeValue != '':
        newLine.append(rangeValue)
    resultsDF.loc[timeStep] = newLine

    return resultsDF


def updateSavingTimeStep(dtSave, cfg, t):
    """ update saving time step list

        Parameters
        -----------
        dtSave: list
            list of time steps that shall be saved
        cfg: configparser object
            configuration settings, end time step
        t: float
            actual time step

        Returns
        --------
        dtSave: list
            updated list of saving time steps

    """

    if dtSave.size == 1:
        dtSave = np.asarray([2*cfg.getfloat('tEnd')])
    else:
        indSave = np.where(dtSave > (t+1.e-8))
        dtSave = dtSave[indSave]

    return dtSave


def appendFieldsParticles(fieldsList, particlesList, particles, fields, resTypes):
    """ append fields and optionally particle dictionaries to list for export

        Parameters
        ------------
        particles: dict
            dictionary with particle properties
        fields: dict
            dictionary with all result type fields
        resTypes: list
            list with all result types that shall be exported

        Returns
        -------
        Fields: list
            updated list with desired result type fields dictionary
        Particles: list
            updated list with particles dicionaries
    """

    fieldAppend = {}
    for resType in resTypes:
        if resType == 'particles':
            particlesList.append(copy.deepcopy(particles))
        elif resType != '':
            fieldAppend[resType] = copy.deepcopy(fields[resType])
    fieldsList.append(fieldAppend)

    return fieldsList, particlesList


def writeMBFile(infoDict, avaDir, logName):
    """ write mass balance info to file

        Parameters
        -----------
        infoDict: dict
            info on mass
        avaDir: str or pathlib path
            path to avalanche directory
        logName: str
            simulation name
    """

    t = infoDict['timeStep']
    massEntrained = infoDict['massEntrained']
    massTotal = infoDict['massTotal']

    # write mass balance info to log file
    massDir = pathlib.Path(avaDir, 'Outputs', 'com1DFA')
    fU.makeADir(massDir)
    with open(massDir / ('mass_%s.txt' % logName), 'w') as mFile:
        mFile.write('time, current, entrained\n')
        for m in range(len(t)):
            mFile.write('%.02f,    %.06f,    %.06f\n' % (t[m], massTotal[m], massEntrained[m]))


def computeEulerTimeStep(cfg, particles, fields, zPartArray0, dem, tCPU, frictType):
    """ compute next time step using an euler forward scheme

    Parameters
    ----------
    cfg: configparser
        configuration for DFA simulation
    particles : dict
        particles dictionary at t
    fields : dict
        fields dictionary at t
    zPartArray0 : dict
        z coordinate of particles at t=0s
    dem : dict
        dictionary with dem information
    tCPU : dict
        computation time dictionary
    frictType: int
        indicator for chosen type of friction model

    Returns
    -------
    particles : dict
        particles dictionary at t + dt
    fields : dict
        fields dictionary at t + dt
    tCPU : dict
        computation time dictionary
    """
    # get forces
    startTime = time.time()

    # loop version of the compute force
    log.debug('Compute Force C')
    particles, force, fields = DFAfunC.computeForceC(cfg, particles, fields, dem, frictType)
    tCPUForce = time.time() - startTime
    tCPU['timeForce'] = tCPU['timeForce'] + tCPUForce

    # compute lateral force (SPH component of the calculation)
    startTime = time.time()
    if cfg.getint('sphOption') == 0:
        force['forceSPHX'] = np.zeros(np.shape(force['forceX']))
        force['forceSPHY'] = np.zeros(np.shape(force['forceY']))
        force['forceSPHZ'] = np.zeros(np.shape(force['forceZ']))
    else:
        log.debug('Compute Force SPH C')
        particles, force = DFAfunC.computeForceSPHC(cfg, particles, force, dem, cfg.getint('sphOption'), gradient=0)
    tCPUForceSPH = time.time() - startTime
    tCPU['timeForceSPH'] = tCPU['timeForceSPH'] + tCPUForceSPH

    # add bonding force if required (if snowSlide is activated)
    if cfg.getint('snowSlide') == 1:
        force, particles = DFAfunC.computeCohesionForceC(cfg, particles, force)

    # update velocity and particle position
    startTime = time.time()
    # particles = updatePosition(cfg, particles, dem, force)
    log.debug('Update position C')
    particles = DFAfunC.updatePositionC(cfg, particles, dem, force, fields, typeStop=0)
    tCPUPos = time.time() - startTime
    tCPU['timePos'] = tCPU['timePos'] + tCPUPos

    # Split particles
    if cfg.getint('splitOption') == 0:
        # split particles with too much mass
        # this only splits particles that grew because of entrainment
        log.debug('Split particles')
        particles = particleTools.splitPartMass(particles, cfg)
    elif cfg.getint('splitOption') == 1:
        # split merge operation
        # first update fields (compute grid values) because we need the h of the particles to get the aPart
        # ToDo: we could skip the update field and directly do the split merge. This means we would use the old h
        startTime = time.time()
        log.debug('update Fields C')
        particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)
        tcpuField = time.time() - startTime
        tCPU['timeField'] = tCPU['timeField'] + tcpuField
        # Then split merge particles
        particles = particleTools.splitPartArea(particles, cfg, dem)
        particles = particleTools.mergePartArea(particles, cfg, dem)

    # release secondary release area?
    if particles['secondaryReleaseInfo']['flagSecondaryRelease'] == 'Yes':
        particles, zPartArray0 = releaseSecRelArea(cfg, particles, fields, dem, zPartArray0)

    # get particles location (neighbours for sph)
    startTime = time.time()
    log.debug('get Neighbours C')
    particles = DFAfunC.getNeighborsC(particles, dem)

    tCPUNeigh = time.time() - startTime
    tCPU['timeNeigh'] = tCPU['timeNeigh'] + tCPUNeigh

    # update fields (compute grid values)
    startTime = time.time()
    log.debug('update Fields C')
    if fields['computeTA']:
        particles = DFAfunC.computeTrajectoryAngleC(particles, zPartArray0)
    particles, fields = DFAfunC.updateFieldsC(cfg, particles, dem, fields)
    tCPUField = time.time() - startTime
    tCPU['timeField'] = tCPU['timeField'] + tCPUField

    return particles, fields, zPartArray0, tCPU


def prepareArea(line, dem, radius, thList='', combine=True, checkOverlap=True):
    """ convert shape file polygon to raster

    Parameters
    ----------
    line: dict
        line dictionary
    dem : dict
        dictionary with dem information
    radius : float
        include all cells which center is in the polygon or close enough
    thList: list
        thickness values for all features in the line dictionary
    combine : Boolean
        if True sum up the rasters in the area list to return only 1 raster
        if False return the list of distinct area rasters
        this option works only if thList is not empty
    checkOverlap : Boolean
        if True check if features are overlaping and return an error if it is the case
        if False check if features are overlaping and average the value for overlaping areas
        (Attention: if combine is set to False, you do not see the result of the averaging
        since the list of raters was not affected by the averaging step)

    Returns
    -------
    updates the line dictionary with the rasterData: Either
        contains either

        -  Raster: 2D numpy array, raster of the area (returned if relRHlist is empty OR if combine is set
        to True)
        - RasterList: list, list of 2D numpy array rasters (returned if relRHlist is not empty AND
        if combine is set to False)

    """
    NameRel = line['Name']
    StartRel = line['Start']
    LengthRel = line['Length']
    RasterList = []

    for i in range(len(NameRel)):
        name = NameRel[i]
        start = StartRel[i]
        end = start + LengthRel[i]
        avapath = {}
        avapath['x'] = line['x'][int(start):int(end)]
        avapath['y'] = line['y'][int(start):int(end)]
        avapath['Name'] = name
        # if relTh is given - set relTh
        if thList != '':
            log.info('%s feature %s, thickness: %.2f - read from %s' % (line['type'], name, thList[i],
                     line['thicknessSource'][i]))
            Raster = polygon2Raster(dem['originalHeader'], avapath, radius, th=thList[i])
        else:
            Raster = polygon2Raster(dem['originalHeader'], avapath, radius)
        RasterList.append(Raster)

    # if RasterList not empty check for overlap between features
    Raster = np.zeros(np.shape(dem['rasterData']))
    for rast in RasterList:
        ind1 = Raster > 0
        ind2 = rast > 0
        indMatch = np.logical_and(ind1, ind2)
        if indMatch.any():
            # if there is an overlap, raise error
            if checkOverlap:
                message = 'Features are overlaping - this is not allowed'
                log.error(message)
                raise AssertionError(message)
            else:
                # if there is an overlap, take average of values for the overlapping cells
                Raster = np.where(((Raster > 0) & (rast > 0)), (Raster + rast)/2, Raster + rast)
        else:
            Raster = Raster + rast
    if debugPlot:
        debPlot.plotAreaDebug(dem, avapath, Raster)
    if combine:
        line['rasterData'] = Raster
        return line
    else:
        line['rasterData'] = RasterList
        return line


def polygon2Raster(demHeader, Line, radius, th=''):
    """ convert line to raster

    Parameters
    ----------
    demHeader: dict
        dem header dictionary
    Line : dict
        line dictionary
    radius : float
        include all cells which center is in the polygon or close enough
    th: float
        thickness value ot the line feature

    Returns
    -------
    Mask : 2D numpy array
        updated raster
    """
    # adim and center dem and polygon
    ncols = demHeader['ncols']
    nrows = demHeader['nrows']
    xllc = demHeader['xllcenter']
    yllc = demHeader['yllcenter']
    csz = demHeader['cellsize']
    xCoord0 = (Line['x'] - xllc) / csz
    yCoord0 = (Line['y'] - yllc) / csz
    if (xCoord0[0] == xCoord0[-1]) and (yCoord0[0] == yCoord0[-1]):
        xCoord = np.delete(xCoord0, -1)
        yCoord = np.delete(yCoord0, -1)
    else:
        xCoord = copy.deepcopy(xCoord0)
        yCoord = copy.deepcopy(yCoord0)
        xCoord0 = np.append(xCoord0, xCoord0[0])
        yCoord0 = np.append(yCoord0, yCoord0[0])

    # get the raster corresponding to the polygon
    polygon = np.stack((xCoord, yCoord), axis=-1)
    path = mpltPath.Path(polygon)
    # add a tolerance to include cells for which the center is on the lines
    # for this we need to know if the path is clockwise or counter clockwise
    # to decide if the radius should be positif or negatif in contains_points
    is_ccw = geoTrans.isCounterClockWise(path)
    r = (radius*is_ccw - radius*(1-is_ccw))
    x = np.linspace(0, ncols-1, ncols)
    y = np.linspace(0, nrows-1, nrows)
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    points = np.stack((X, Y), axis=-1)
    mask = path.contains_points(points, radius=r)
    Mask = mask.reshape((nrows, ncols)).astype(int)
    # thickness field is provided, then return array with ones
    if th != '':
        log.debug('REL set from dict, %.2f' % th)
        Mask = np.where(Mask > 0, th, 0.)
    else:
        Mask = np.where(Mask > 0, 1., 0.)

    if debugPlot:
        debPlot.plotRemovePart(xCoord0, yCoord0, demHeader, X, Y, Mask, mask)

    return Mask


def checkParticlesInRelease(particles, line, radius):
    """ remove particles laying outside the polygon

    Parameters
    ----------
    particles : dict
        particles dictionary
    line: dict
        line dictionary
    radius: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside

    Returns
    -------
    particles : dict
        particles dictionary where particles outside of the polygon have been removed
    """
    NameRel = line['Name']
    StartRel = line['Start']
    LengthRel = line['Length']
    Mask = np.full(np.size(particles['x']), False)
    for i in range(len(NameRel)):
        name = NameRel[i]
        start = StartRel[i]
        end = start + LengthRel[i]
        avapath = {}
        avapath['x'] = line['x'][int(start):int(end)]
        avapath['y'] = line['y'][int(start):int(end)]
        avapath['Name'] = name
        mask = pointInPolygon(line['header'], particles, avapath, radius)
        Mask = np.logical_or(Mask, mask)

    # also remove particles with negative mass
    mask = np.where(particles['m'] <= 0, False, True)
    Mask = np.logical_and(Mask, mask)
    nRemove = len(Mask)-np.sum(Mask)
    if nRemove > 0:
        particles = particleTools.removePart(particles, Mask, nRemove, '')
        log.debug('removed %s particles because they are not within the release polygon' % (nRemove))

    return particles


def pointInPolygon(demHeader, points, Line, radius):
    """ find particles within a polygon

    Parameters
    ----------
    demHeader: dict
        dem header dictionary
    points: dict
        points to check
    Line : dict
        line dictionary
    radius: float
        threshold val that decides if a point is in the polygon, on the line or
        very close but outside

    Returns
    -------
    Mask : 1D numpy array
        Mask of particles to keep
    """
    xllc = demHeader['xllcenter']
    yllc = demHeader['yllcenter']
    xCoord0 = (Line['x'] - xllc)
    yCoord0 = (Line['y'] - yllc)
    if (xCoord0[0] == xCoord0[-1]) and (yCoord0[0] == yCoord0[-1]):
        xCoord = np.delete(xCoord0, -1)
        yCoord = np.delete(yCoord0, -1)
    else:
        xCoord = copy.deepcopy(xCoord0)
        yCoord = copy.deepcopy(yCoord0)
        xCoord0 = np.append(xCoord0, xCoord0[0])
        yCoord0 = np.append(yCoord0, yCoord0[0])

    # get the raster corresponding to the polygon
    polygon = np.stack((xCoord, yCoord), axis=-1)
    path = mpltPath.Path(polygon)
    # add a tolerance to include cells for which the center is on the lines
    # for this we need to know if the path is clockwise or counter clockwise
    # to decide if the radius should be positif or negatif in contains_points
    is_ccw = geoTrans.isCounterClockWise(path)
    r = (radius*is_ccw - radius*(1-is_ccw))
    points2Check = np.stack((points['x'], points['y']), axis=-1)
    mask = path.contains_points(points2Check, radius=r)
    mask = np.where(mask > 0, True, False)

    if debugPlot:
        debPlot.plotPartAfterRemove(points, xCoord0, yCoord0, mask)

    return mask


def releaseSecRelArea(cfg, particles, fields, dem, zPartArray0):
    """ Release secondary release area if trigered
    Initialize particles of the trigured secondary release area and add them
    to the simulation (particles dictionary)
    """

    secondaryReleaseInfo = particles['secondaryReleaseInfo']
    flowThicknessField = fields['FT']
    secRelRasterList = secondaryReleaseInfo['rasterData']
    secRelRasterNameList = secondaryReleaseInfo['Name']
    count = 0
    indexRel = []
    for secRelRaster, secRelRasterName in zip(secRelRasterList, secRelRasterNameList):
        # do the two arrays intersect (meaning a flowing particle entered the
        # secondary release area)
        mask = (secRelRaster > 0) & (flowThicknessField > 0)
        if mask.any():
            # create secondary release area particles
            log.info('Initializing secondary release area feature %s' % secRelRasterName)
            secRelInfo = shpConv.extractFeature(secondaryReleaseInfo, count)
            secRelInfo['rasterData'] = secRelRaster
            secRelParticles = initializeParticles(cfg, secRelInfo, dem)
            # release secondary release area by just appending the particles
            log.info('Releasing secondary release area %s at t = %.2f s' % (secRelRasterName, particles['t']))
            particles = particleTools.mergeParticleDict(particles, secRelParticles)
            # save index of secRel feature
            indexRel.append(secRelRasterName)
            # save initial z position for travel angle computation
            zPartArray0 = np.append(zPartArray0, copy.deepcopy(secRelParticles['z']))
        count = count + 1

    secondaryReleaseInfo['rasterData'] = secRelRasterList
    particles['secondaryReleaseInfo'] = secondaryReleaseInfo
    for item in indexRel:
        iR = secRelRasterNameList.index(item)
        # remove it from the secondary release area list
        secRelRasterList.pop(iR)
        secondaryReleaseInfo = shpConv.removeFeature(secondaryReleaseInfo, iR)
        secRelRasterNameList.pop(iR)

    # update secondaryReleaseInfo
    secondaryReleaseInfo['rasterData'] = secRelRasterList
    particles['secondaryReleaseInfo'] = secondaryReleaseInfo

    return particles, zPartArray0


def savePartToPickle(dictList, outDir, logName):
    """ Save each dictionary from a list to a pickle in outDir; works also for one dictionary instead of list
        Note: particle coordinates are still in com1DFA reference system with origin 0,0

        Parameters
        ---------
        dictList: list or dict
            list of dictionaries or single dictionary
        outDir: str
            path to output directory
        logName : str
            simulation Id
    """

    if isinstance(dictList, list):
        for dict in dictList:
            fi = open(outDir / ("particles_%s_%09.4f.p" % (logName, dict['t'])), "wb")
            pickle.dump(dict, fi )
            fi.close()
    else:
        fi = open(outDir / ("particles_%s_%09.4f.p" % (logName, dictList['t'])), "wb")
        pickle.dump(dictList, fi)
        fi.close()


def trackParticles(cfgTrackPart, dem, particlesList):
    """ track particles from initial area

        Find all particles in an initial area. Find the same particles in
        the other time steps (+ the children if they were splitted).
        Extract time series of given properties of the tracked particles

        Parameters
        -----------
        cfgTrackPart: configParser
            centerTrackPartPoint : str
                centerTrackPartPoint of the location of the particles to
                track (x|y coordinates)
            radius : str
                radius of the circle around point
            particleProperties: str
                list of particles properties to extract ('x', 'y', 'ux', 'm'...)
        dem: dict
            dem dictionary
        particlesList: list
            list of particles dictionary

        Returns
        -------
        particlesList : list
            Particles list of dict updated with the 'trackedParticles' array
            (in the array, ones for particles that are tracked, zeros otherwise)
        trackedPartProp: dict
            dictionary with time series of the wanted properties for tracked
            particles
        track: boolean
            False if no particles are tracked
    """

    # read particle properties to be extracted
    particleProperties = cfgTrackPart['particleProperties']
    if particleProperties == '':
        particleProperties = ['x', 'y', 'z', 'ux', 'uy', 'uz', 'm', 'h']
    else:
        particleProperties = set(['x', 'y', 'z', 'ux', 'uy', 'uz', 'm', 'h'] + particleProperties.split('|'))
    # read location of particle to be tracked
    radius = cfgTrackPart.getfloat('radius')
    centerList = cfgTrackPart['centerTrackPartPoint']
    centerList = centerList.split('|')
    centerTrackPartPoint = {'x': np.array([float(centerList[0])]),
                            'y': np.array([float(centerList[1])])}
    centerTrackPartPoint['x'] = (centerTrackPartPoint['x']
                                 - dem['originalHeader']['xllcenter'])
    centerTrackPartPoint['y'] = (centerTrackPartPoint['y']
                                 - dem['originalHeader']['yllcenter'])

    # start by finding the particles to be tracked
    particles2Track, track = particleTools.findParticles2Track(particlesList[0], centerTrackPartPoint, radius)
    if track:
        # find those same particles and their children in the particlesList
        particlesList, nPartTracked = particleTools.getTrackedParticles(particlesList, particles2Track)

        # extract the wanted properties for the tracked particles
        trackedPartProp = particleTools.getTrackedParticlesProperties(particlesList, nPartTracked, particleProperties)
    else:
        trackedPartProp = None

    return particlesList, trackedPartProp, track


def readFields(inDir, resType, simName='', flagAvaDir=True, comModule='com1DFA', timeStep='', atol=1.e-6):
    """ Read ascii files within a directory and return List of dictionaries

    Parameters
    -----------
    inDir: str
        path to input directory
    resType: list
        list of desired result types, if string converted to list
    simName : str
        simulation name
    flagAvaDir: bool
        if True inDir corresponds to an avalanche directory and pickles are
        read from avaDir/Outputs/com1DFA/particles
    comModule: str
        module that computed the particles
    timeStep: float or list of floats
        desired time step if difference to time step of field file is smaller than atol
        field is found - optional
    atol: float
        look for matching time steps with atol tolerance - default is atol=1.e-6

    Returns
    -------
    fieldsList : list
        list of fields dictionaries
    fieldHeader: dict
        raster header corresponding to first element in fieldsList
    timeList: list
        tme corresponding to elements in fieldsList

    """

    if isinstance(resType, list) is False:
        resType = [resType]

    if flagAvaDir:
        inDir = pathlib.Path(inDir, 'Outputs', comModule, 'peakFiles', 'timeSteps')

    # initialise list of fields dictionaries
    fieldsList = []
    timeList = []
    first = True
    for r in resType:
        # search for all files within directory
        if simName:
            name = '*' + simName + '*' + r + '*.asc'
        else:
            name = '*' + r + '*.asc'
        FieldsNameList = list(inDir.glob(name))
        timeListTemp = [float(element.stem.split('_t')[-1]) for element in FieldsNameList]
        FieldsNameList = [x for _, x in sorted(zip(timeListTemp, FieldsNameList))]
        count = 0
        for fieldsName in FieldsNameList:
            t = float(fieldsName.stem.split('_t')[-1])
            if timeStep == '' or np.isclose(timeStep, t, atol=atol).any():
                # initialize field Dict
                if first:
                    fieldsList.append({})
                    timeList.append(t)
                field = IOf.readRaster(fieldsName)
                fieldsList[count][r] = field['rasterData']
                count = count + 1
        first = False

    if count == 0:
        log.warning('No matching fields found in %s' % inDir)
        fieldHeader = None
        fieldsList = []
    else:
        fieldHeader = field['header']

    return fieldsList, fieldHeader, timeList


def exportFields(cfg, Tsave, fieldsList, dem, outDir, logName):
    """ export result fields to Outputs directory according to result parameters and time step
        that can be specified in the configuration file

        Parameters
        -----------
        cfg: dict
            configurations
        Tsave: list
            list of time step that corresponds to each dict in Fields
        Fields: list
            list of Fields for each dtSave
        outDir: str
            outputs Directory

        Returns
        --------
        exported peak fields are saved in Outputs/com1DFA/peakFiles
    """

    resTypesGen = fU.splitIniValueToArraySteps(cfg['GENERAL']['resType'])
    resTypesReport = fU.splitIniValueToArraySteps(cfg['REPORT']['plotFields'])
    if 'particles' in resTypesGen:
        resTypesGen.remove('particles')
    if 'particles' in resTypesReport:
        resTypesReport.remove('particles')
    numberTimes = len(Tsave)-1
    countTime = 0
    for timeStep in Tsave:
        if (countTime == numberTimes) or (countTime == 0):
            # for last time step we need to add the report fields
            resTypes = list(set(resTypesGen + resTypesReport))
        else:
            resTypes = resTypesGen
        for resType in resTypes:
            resField = fieldsList[countTime][resType]
            if resType == 'ppr':
                # convert from Pa to kPa
                resField = resField * 0.001
            if resType == 'pke':
                # convert from J/cell to kJ/m²
                # (by dividing the peak kinetic energy per cell by the real area of the cell)
                resField = resField * 0.001 / dem['areaRaster']
            dataName = (logName + '_' + resType + '_' + 't%.2f' % (Tsave[countTime]) + '.asc')
            # create directory
            outDirPeak = outDir / 'peakFiles' / 'timeSteps'
            fU.makeADir(outDirPeak)
            outFile = outDirPeak / dataName
            IOf.writeResultToAsc(
                dem['originalHeader'], resField, outFile, flip=True)
            if countTime == numberTimes:
                log.debug('Results parameter: %s exported to Outputs/peakFiles for time step: %.2f - FINAL time step ' %
                          (resType, Tsave[countTime]))
                dataName = logName + '_' + resType + '.asc'
                # create directory
                outDirPeakAll = outDir / 'peakFiles'
                fU.makeADir(outDirPeakAll)
                outFile = outDirPeakAll / dataName
                IOf.writeResultToAsc(dem['originalHeader'], resField, outFile, flip=True)
            else:
                log.debug('Results parameter: %s has been exported to Outputs/peakFiles for time step: %.2f ' %
                          (resType, Tsave[countTime]))
        countTime = countTime + 1


def prepareVarSimDict(standardCfg, inputSimFiles, variationDict, simNameExisting=''):
    """ Prepare a dictionary with simulations that shall be run with varying parameters following the variation dict

        Parameters
        -----------
        standardCfg : configParser object
            default configuration or local configuration
        inputSimFiles: dict
            info dict on available input data
        variationDict: dict
            dictionary with parameter to be varied as key and list of it's values
        simNameExisting: list
            list of simulation names that already exist (optional). If provided,
            only carry on simulations that do not exist

        Returns
        -------
        simDict: dict
            dicionary with info on simHash, releaseScenario, release area file path,
            simType and contains full configuration configparser object for simulation run
    """

    # get list of simulation types that are desired
    if 'simTypeList' in variationDict:
        simTypeList = variationDict['simTypeList']
        del variationDict['simTypeList']
    else:
        simTypeList = standardCfg['GENERAL']['simTypeList'].split('|')
    # get a list of simulation types that are desired AND available
    standardCfg, simTypeList = getSimTypeList(standardCfg, simTypeList, inputSimFiles)

    # set simTypeList (that has been checked if available) as parameter in variationDict
    variationDict['simTypeList'] = simTypeList
    # create a dataFrame with all possible combinations of the variationDict values
    variationDF = pd.DataFrame(product(*variationDict.values()), columns=variationDict.keys())

    # generate a dictionary of full simulation info for all simulations to be performed
    # simulation info must contain: simName, releaseScenario, relFile, configuration as dictionary
    simDict = {}

    # loop over all simulations that shall be performed according to variationDF
    # one row per simulation
    for row in variationDF.itertuples():
        # convert full configuration to dict
        cfgSim = cfgUtils.convertConfigParserToDict(standardCfg)
        # create release scenario name for simulation name
        rel, cfgSim = gI.fetchReleaseFile(inputSimFiles, row._asdict()['releaseScenario'], cfgSim, variationDict['releaseScenario'])
        relName = rel.stem
        if '_' in relName:
            relNameSim = relName + '_AF'
        else:
            relNameSim = relName

        # update info for parameters that are given in variationDF
        for parameter in variationDict:
            # add simType
            cfgSim['GENERAL']['simTypeActual'] = row._asdict()['simTypeList']
            # update parameter value - now only single value for each parameter
            keyList = ['relThPercentVariation', 'entThPercentVariation',
                       'secondaryRelThPercentVariation', 'relThRangeVariation',
                       'entThRangeVariation', 'secondaryRelThRangeVariation',
                       'relThRangeFromCiVariation', 'entThRangeFromCiVariation',
                       'secondaryRelThRangeFromCiVariation', 'relThDistVariation',
                       'entThDistVariation', 'secondaryRelThDistVariation']
            if parameter in keyList:
                # set thickness value according to percent variation info
                cfgSim = dP.setThicknessValueFromVariation(parameter, cfgSim, cfgSim['GENERAL']['simTypeActual'], row)
            elif parameter == 'releaseScenario':
                cfgSim['INPUT'][parameter] = row._asdict()[parameter]
            else:
                cfgSim['GENERAL'][parameter] = row._asdict()[parameter]

        # update INPUT section - delete non relevant parameters
        if cfgSim['GENERAL']['simTypeActual'] not in ['ent', 'entres']:
            cfgSim['INPUT'].pop('entrainmentScenario', None)
            cfgSim['INPUT'].pop('entThId', None)
            cfgSim['INPUT'].pop('entThThickness', None)
            cfgSim['INPUT'].pop('entThCi95', None)
        if cfgSim['GENERAL']['secRelArea'] == 'False':
            cfgSim['INPUT'].pop('secondaryReleaseScenario', None)
            cfgSim['INPUT'].pop('secondaryRelThId', None)
            cfgSim['INPUT'].pop('secondaryRelThThickness', None)
            cfgSim['INPUT'].pop('secondaryRelThCi95', None)

        # check if DEM in Inputs has desired mesh size
        pathToDem = dP.checkDEM(cfgSim, inputSimFiles['demFile'])
        cfgSim['INPUT']['DEM'] = pathToDem

        # add thickness values if read from shp and not varied
        cfgSim = dP.appendShpThickness(cfgSim)

        # check differences to default and add indicator to name
        defID, _ = com1DFATools.compareSimCfgToDefaultCfgCom1DFA(cfgSim)

        # if frictModel is samosATAuto compute release vol
        if cfgSim['GENERAL']['frictModel'].lower() == 'samosatauto':
            pathToDemFull = pathlib.Path(cfgSim['GENERAL']['avalancheDir'], 'Inputs', pathToDem)
            relVolume = fetchRelVolume(rel, cfgSim, pathToDemFull)
        else:
            relVolume = ''

        # check sphKernelRadius setting
        cfgSim = checkCfg.checkCellSizeKernelRadius(cfgSim)

        # only keep friction model parameters that are used
        cfgSim = checkCfg.checkCfgFrictionModel(cfgSim, relVolume=relVolume)

        # set frictModelIndicator, this needs to happen AFTER checkCfgFrictModel
        frictIndi = com1DFATools.setFrictTypeIndicator(cfgSim)

        # convert back to configParser object
        cfgSimObject = cfgUtils.convertDictToConfigParser(cfgSim)
        # create unique hash for simulation configuration
        simHash = cfgUtils.cfgHash(cfgSimObject)

        simName = '_'.join(filter(None, [relNameSim, simHash,
                            defID, frictIndi, row._asdict()['simTypeList'],
                            cfgSim['GENERAL']['modelType']]))

        # check if simulation exists. If yes do not append it
        if simName not in simNameExisting:
            simDict[simName] = {'simHash': simHash, 'releaseScenario': relName,
                                'simType': row._asdict()['simTypeList'], 'relFile': rel,
                                'cfgSim': cfgSimObject}
            # write configuration file
            cfgUtils.writeCfgFile(cfgSimObject['GENERAL']['avalancheDir'], com1DFA, cfgSimObject, fileName=simName)
        else:
            log.warning('Simulation %s already exists, not repeating it' % simName)

    inputSimFiles.pop('demFile')

    return simDict


def getSimTypeList(standardCfg, simTypeList, inputSimFiles):
    """ Define available simulation types of requested types

        Parameters
        -----------
        standardCfg : configParser object
            default configuration or local configuration
        simTypeList: List
            list of simTypes to conpute (ent, null...)
        inputSimFiles: dict
            info dict on available input data

        Returns
        --------
        standardCfg : configParser object
            configuration with updated 'secRelArea' depending on if a secondary release file is available or not
        simTypeList: list
            list of requested simTypes where also the required input data is available
    """

    # read entrainment resistance info
    entResInfo = inputSimFiles['entResInfo']

    # check if set simType is a valid option
    validSimTypesStr = 'available|null|ent|entres|res'
    validSimTypes = validSimTypesStr.split('|')
    validArray = [True if item in validSimTypes else False for item in simTypeList]
    if False in validArray:
        message = 'A non-valid entry found in simType, valid Types are %s' % validSimTypesStr
        log.error(message)
        raise AssertionError(message)

    # define simulation type
    if 'available' in simTypeList:
        if entResInfo['flagEnt'] == 'Yes' and entResInfo['flagRes'] == 'Yes':
            simTypeList.append('entres')
        elif entResInfo['flagEnt'] == 'Yes' and entResInfo['flagRes'] == 'No':
            simTypeList.append('ent')
        elif entResInfo['flagEnt'] == 'No' and entResInfo['flagRes'] == 'Yes':
            simTypeList.append('res')
        # always add null simulation
        simTypeList.append('null')
        simTypeList.remove('available')

    # remove duplicate entries
    simTypeList = set(simTypeList)
    simTypeList = sorted(list(simTypeList), reverse=False)

    if 'ent' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagEnt'] == 'No':
            message = 'No entrainment file found'
            log.error(message)
            raise FileNotFoundError(message)
    if 'res' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagRes'] == 'No':
            message = 'No resistance file found'
            log.error(message)
            raise FileNotFoundError(message)
    if standardCfg['GENERAL'].getboolean('secRelArea'):
        if entResInfo['flagSecondaryRelease'] == 'No':
            standardCfg['GENERAL']['secRelArea'] = 'False'
        else:
            log.info('Using the secondary release area file: %s' % inputSimFiles['secondaryReleaseFile'])

    return standardCfg, simTypeList


def runOrLoadCom1DFA(avalancheDir, cfgMain, runDFAModule=True, cfgFile='', deleteOutput=True):
    """ Run or load DFA results depending on runDFAModule=True or False

        Parameters
        -----------
        avalancheDir: pathlib path
            avalanche directory path
        cfgMain : configParser object
            main avaframe configuration
        runDFAModule: bool
            True to run the DFA simulation Falso to load the results dataFrame and dem
        cfgFile: str or pathlib path
            path to cfgFile to read overall configuration - optional if not provided the local or default config is used
        deleteOutput: Boolean
            True to delete the com1DFA output dir before running com1DFA (used only if runDFAModule=True)

        Returns
        --------
        dem: dict
            dem dictionary
        simDF: dataFrame
            simulation results dataframe
        resTypeList: list
            list of output files resTypes available (ppr, pft...)
    """
    if runDFAModule:
        # clean avalanche directory
        iP.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=deleteOutput)
        # Run the DFA simulation
        dem, _, _, simDF = com1DFA.com1DFAMain(cfgMain, cfgInfo=cfgFile)
    else:
        # read simulation dem
        demOri = gI.readDEM(avalancheDir)
        dem = com1DFA.setDEMoriginToZero(demOri)
        dem['originalHeader'] = demOri['header'].copy()
        # load DFA results
        simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)
        if simDF is None:
            message = 'Did not find any com1DFA simulations in %s/Outputs/com1DFA/' % avalancheDir
            log.error(message)
            raise FileExistsError(message)

    dataDF, resTypeList = fU.makeSimFromResDF(avalancheDir, 'com1DFA', inputDir='', simName='')
    simDF = simDF.reset_index().merge(dataDF, on='simName').set_index('index')
    return dem, simDF, resTypeList


def fetchRelVolume(releaseFile, cfg, pathToDem, radius=0.01):
    """ compute release area volume using release line and thickness info and dem

        Parameters
        -----------
        releaseFile: pathlib path
            path to release area shp file
        cfg: dict
            config settings of current sim
        pathToDem: pathlib path
            path to dem file used for current sim
        radius : float
            include all cells which center is in the release line or close enough

        Returns
        ---------
        relVolume: float
            volume of release area in m3


    """

    # convert back to configParser object
    cfg = cfgUtils.convertDictToConfigParser(cfg)

    # read simulation dem
    demVol = IOf.readRaster(pathToDem, noDataToNan=True)
    demVol['originalHeader'] = demVol['header'].copy()
    methodMeshNormal = cfg['GENERAL'].getfloat('methodMeshNormal')
    # get normal vector of the grid mesh
    demVol = DFAtls.getNormalMesh(demVol, methodMeshNormal)
    demVol = DFAtls.getAreaMesh(demVol, methodMeshNormal)

    # create release line
    releaseLine = {}
    releaseLine = shpConv.readLine(releaseFile, 'release1', demVol)
    # check if release features overlap between features
    thresholdPointInPoly = cfg['GENERAL'].getfloat('thresholdPointInPoly')
    prepareArea(releaseLine, demVol, thresholdPointInPoly, combine=True, checkOverlap=True)
    releaseLine['type'] = 'Release'

    # check if release thickness provided as field or constant value
    if cfg['GENERAL']['relThFromFile'] == 'True':

        # read relThField from file
        relThFilePath = pathlib.Path(cfg['GENERAL']['avalancheDir'], 'Inputs', cfg['INPUT']['relThFile'])
        relThFieldFull = IOf.readRaster(relThFilePath)
        relThField = relThFieldFull['rasterData']

        # create raster from polygon
        releaseLine = prepareArea(releaseLine, demVol, radius, combine=True, checkOverlap=False)

        # mask the relThField with raster from polygon
        releaseLineMask = np.ma.masked_where(releaseLine['rasterData']==0., releaseLine['rasterData'])
        releaseLineField = np.ma.masked_where(np.ma.getmask(releaseLineMask),relThField)
        relVolumeField = np.ma.masked_where(np.ma.getmask(releaseLineMask),relThField) * demVol['areaRaster']
        relVolume = np.nansum(relVolumeField)

        if debugPlot:
            debPlot.plotVolumeRelease(releaseLine, relThField, releaseLineField)
    else:
        relThField = ''

        # set thickness values on releaseLine
        releaseLine = setThickness(cfg, releaseLine, 'relTh')
        # when creating raster from polygon apply release thickness
        releaseLine = prepareArea(releaseLine, demVol, radius, thList=releaseLine['thickness'],
            combine=True, checkOverlap=False)

        # compute release volume using raster and dem area
        relVolume = np.nansum(releaseLine['rasterData'] * demVol['areaRaster'])

    log.info('release volume is: %.2f m3' % relVolume)

    return relVolume
