"""
    Python wrapper to execute the compiled com1Exe file and set desired simulation options
"""

# Load modules
import os
import subprocess
import shutil
import logging
import numpy as np
import pickle
import time
from datetime import datetime

# Local imports
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initialiseDirs as iD
from avaframe.in1Data import getInput as gI
from avaframe.in2Trans import shpConversion as sP
import avaframe.in2Trans.ascUtils as IOf

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def execCom1Exe(com1Exe, cintFile, avaDir, fullOut=False, simName=''):
    """ Execute compiled com1Exe file using cintFile to set configuration
        and run options

        Parameters
        ----------
        com1Exe : str
            path to com1Exe
        cintFile : str
            path to cint file
        avaDir : str
            path to avalanche directoy
        fullOut : bool
            flag if True print full output from com1Exe to terminal and log, default False
        simName : str
            optional - name of simulation, will be used to create a log
    """

    # define command line
    runCommand = com1Exe + ' -offscreen ' + cintFile

    # initialise log file to save stoudt
    if simName != '':
        f_log = open(os.path.join(avaDir, 'Outputs', 'com1DFAOrig', 'start%s.log' % (simName)), 'w')

    # Call command
    proc = subprocess.Popen(runCommand, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)

    # loop through output and save to logFile if file provided
    for line in proc.stdout:
        if simName != '':
            f_log.write(line)
        if fullOut:
            log.info((line.rstrip()))
        elif 'BatchCom1DFA' in line:
            if 'Computing' in line:
                log.info(line.rstrip())
            else:
                log.debug(line.rstrip())
        elif 'error' in line:
            log.info(line.rstrip())

    # make process wait for previous process to finish
    reVal = proc.wait()


def copyReplace(origFile, workFile, searchString, replString):
    """ Modifiy cintFiles to be used to set simulation configuration"""

    # Check if input files match, if not save origFile to workFile
    try:
        shutil.copy2(origFile, workFile)
    except shutil.SameFileError:
        pass

    # Read file
    fileData = None
    with open(workFile, 'r') as file:
        fileData = file.read()

    # Replace target string
    fileData = fileData.replace(searchString, str(replString))

    # Write new info to file
    with open(workFile, 'w') as file:
        file.write(fileData)


def getSimulation(cfg, rel, entResInfo):
    """ Get a list of all simulations that shall be performed according to simTypeList in configuration file;
        and a dictionary with information on release area


    Parameters
    ----------
    cfg : dict
        general configuration
    entResInfo: dict
        info about available entrainment and resistance info

    Returns
    -------
    simTypeList: list
        list of simulation types to be performed
    """

    # Set release areas and simulation name
    relName = os.path.splitext(os.path.basename(rel))[0]
    simName = relName
    badName = False
    if '_' in relName:
        badName = True
        log.warning('Release area scenario file name includes an underscore \
        the suffix _AF will be added')
        simName = relName + '_AF'
    relDict = sP.SHP2Array(rel)
    for k in range(len(relDict['thickness'])):
        if relDict['thickness'][k] == 'None':
            relDict['thickness'][k] = cfg['DEFVALUES'].getfloat('RelTh')
        else:
            relDict['thickness'][k] = float(relDict['thickness'][k])

    log.info('Release area scenario: %s - perform simulations' % (relName))

    # read list of desired simulation types
    simTypeList = cfg['GENERAL']['simTypeList'].split('|')

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

    if 'ent' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagEnt'] == 'No':
            log.error('No entrainment file found')
            raise FileNotFoundError
    if 'res' in simTypeList or 'entres' in simTypeList:
        if entResInfo['flagRes'] == 'No':
            log.error('No resistance file found')
            raise FileNotFoundError

    return simTypeList, relName, relDict,  badName


def com1DFAOrigMain(cfg, avaDir):
    """ Run main model

    This will compute a dense flow avalanche

    Parameters
    ----------
    cfg : dict
        configuration read from ini file
    avaDir : str
        path to avalanche directory

    Returns
    -------
    reportDictList : list
        list of dictionaries that contain information on simulations that can be used for report generation
    """

    # Setup configuration
    cfgGen = cfg['GENERAL']
    com1Exe = cfgGen['com1Exe']
    modName = 'com1DFAOrig'
    fullOut = cfgGen.getboolean('fullOut')
    cfgPar = cfg['PARAMETERVAR']
    resDir = os.path.join(avaDir, 'Work', 'com1DFAOrig')
    # Get path of module
    modPath = os.path.dirname(__file__)
    # Standard values for parameters that can be varied
    defValues = cfg['DEFVALUES']
    entrainmentTH = defValues['defaultEntH']

    # Log current avalanche directory
    log.debug('Your current avalanche name: %s' % avaDir)

    # Create output and work directories
    workDir, outDir = iD.initialiseRunDirs(avaDir, modName, True)

    # Load input data
    dem, rels, ent, res, entResInfo = gI.getInputData(avaDir, cfgGen)

    # Parameter variation
    if cfgPar.getboolean('parameterVar'):
        varPar = cfgPar['varPar']
    else:
        varPar = 'Mu'

    # Initialise full experiment log file
    dateTimeInfo = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(workDir, 'ExpLog.txt'), 'w') as logFile:
        logFile.write("NoOfSimulation,SimulationRunName,%s\n" % varPar)

    # Counter for release area loop
    countRel = 0

    # Setup simulation dictionaries for report genereation
    reportDictList = []

    # Loop through release areas
    for rel in rels:
        startTime = time.time()

        # load release area and simulation type
        simTypeList, relName, relDict, badName = getSimulation(cfg, rel, entResInfo)
        entrainmentArea = ''
        resistanceArea = ''
        if 'ent' in simTypeList:
            entrainmentArea = os.path.splitext(os.path.basename(ent))[0]
            if cfg['ENTRAINMENT'].getboolean('setEntThickness'):
                entrainmentTH = cfg['ENTRAINMENT']['entH']
                log.info('Entrainment thickness is changed! set to %s' % entrainmentTH)
        if 'res' in simTypeList:
            resistanceArea = os.path.splitext(os.path.basename(res))[0]

        # Initialise CreateProject cint file
        templateFile = os.path.join(modPath, 'CreateProject.cint')
        workFile = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'CreateProject.cint')
        projDir = os.path.join(avaDir, 'Work', 'com1DFAOrig', relName)
        demName = os.path.splitext(os.path.basename(dem))[0]

        # Set Parameters in cint file
        copyReplace(templateFile, workFile, '##PROJECTDIR##', projDir)
        copyReplace(workFile, workFile, '##DHMFILE##', dem)
        copyReplace(workFile, workFile, '##DHMNAME##', demName)
        copyReplace(workFile, workFile, '##RELFILE##', rel)
        copyReplace(workFile, workFile, '##ENTFILE##', ent)
        copyReplace(workFile, workFile, '##RESFILE##', res)
        # Setup Project
        execCom1Exe(com1Exe, workFile, avaDir, fullOut)

        # loop over simulations
        for simTypeActual in simTypeList:

            simName = relName + '_' + simTypeActual + '_dfa'
            log.info('Perform %s simulation ' % simName)

            # set entrainment and resistance info
            resInfo = 'No'
            entInfo = 'No'
            if 'ent' in simTypeActual:
                entInfo = 'Yes'
            if 'res' in simTypeActual:
                resInfo = 'Yes'

            # Initialise CreateSimulations cint file and set parameters
            templateFile = os.path.join(modPath, 'Create%sSimulation.cint' % simTypeActual)
            workFile = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'Create%sSimulation.cint' % simTypeActual)

            # Write required info to cint file
            copyReplace(templateFile, workFile, '##PROJECTDIR##', projDir)
            copyReplace(workFile, workFile, '##BASESIMNAME##', relName)
            execCom1Exe(com1Exe, workFile, avaDir, fullOut)

            # If parameter shall be varied
            if cfgPar.getboolean('parameterVar'):

                log.info('Parameter variation used varying: %s' % cfgPar['varPar'])

                # read values of parameter variation in config file
                itemsRaw = fU.splitIniValueToArraySteps(cfg['PARAMETERVAR']['varParValues'])
                items = []
                for itemR in itemsRaw:
                    items.append('%.5f' % float(itemR))
                for item in items:
                    startTime = time.time()
                    logName = simName + '_' + item
                    log.info('Perform simulation with %s = %s: logName = %s' % (cfgPar['varPar'], item, logName))
                    templateFile = os.path.join(modPath, '%s%s.cint' % (cfgPar['varRunCint'], cfgPar['varPar']))
                    workFile = os.path.join(avaDir, 'Work', 'com1DFAOrig',
                                            '%s%sBasic.cint' % (cfgPar['varRunCint'], cfgPar['varPar']))
                    copyReplace(templateFile, workFile, '##PROJECTDIR##', projDir)
                    copyReplace(workFile, workFile, '##RESDIR##', resDir)
                    copyReplace(workFile, workFile, '##NAME##', simName)
                    copyReplace(workFile, workFile, '##COUNTREL##', countRel)
                    copyReplace(workFile, workFile, '##VALUE##', item)
                    copyReplace(workFile, workFile, '##ENTH##', entrainmentTH)
                    execCom1Exe(com1Exe, workFile, avaDir, fullOut, logName)

                    # Create dictionary
                    reportVar = {}
                    reportVar = {'headerLine': {'type': 'title', 'title': 'com1DFAOrig Simulation'},
                    'avaName': {'type': 'avaName', 'name': avaDir},
                    'simName': {'type': 'simName', 'name': logName},
                    'time': {'type': 'time', 'time': dateTimeInfo},
                        'Simulation Parameters': {
                            'type': 'list',
                            'Release Area Scenario': relName,
                            'Release Area': relDict['Name'],
                            'Entrainment': entInfo,
                            'Resistance': resInfo,
                            'Parameter variation on': cfgPar['varPar'],
                            'Parameter value': item},
                        'Release area': {'type': 'columns', 'Release area scenario': relName,  'Release features': relDict['Name']}}

                    if cfgPar['varPar'] == 'RelTh':
                        reportVar['Simulation Parameters'].update({'Mu': defValues['Mu']})
                        reportVar['Simulation Parameters'].update({'Release thickness [m]': item})
                        reportVar['Simulation Parameters'].update({'Entrainment thickness [m]': float(entrainmentTH)})
                        reportVar['Release area'].update({'Release thickness [m]': item})
                    elif cfgPar['varPar'] == 'Mu':
                        reportVar['Simulation Parameters'].update({'Release thickness [m]': relDict['thickness']})
                        reportVar['Simulation Parameters'].update({'Mu': item})
                        reportVar['Simulation Parameters'].update({'Entrainment thickness [m]': float(entrainmentTH)})
                        reportVar['Release area'].update({'Release thickness [m]': relDict['thickness']})

                    if entInfo == 'Yes':
                        reportVar['Entrainment area'] = {'type': 'columns', 'Entrainment area scenario': entrainmentArea, 'Entrainment thickness [m]': float(entrainmentTH)}
                    if resInfo == 'Yes':
                        reportVar['Resistance area'] = {'type': 'columns', 'Resistance area scenario': resistanceArea}

                    endTime = time.time()
                    timeNeeded =  '%.2f' % (endTime - startTime)
                    log.info(('Took %s seconds to calculate.' % (timeNeeded)))
                    reportVar['Simulation Parameters'].update({'Computation time [s]': timeNeeded})

                    # Add to report dictionary list
                    reportDictList.append(reportVar)
                    startTime = time.time()

                    # Count total number of simulations
                    countRel = countRel + 1

            else:

                # set entrainment and resistance
                startTime = time.time()
                templateFile = os.path.join(modPath, 'runStandardSimulation.cint')
                workFile = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'runStandardSimulation.cint')
                logName = simName + '_' + defValues['Mu']
                # Write required info to cint file
                copyReplace(templateFile, workFile, '##PROJECTDIR##', projDir)
                copyReplace(workFile, workFile, '##RESDIR##', resDir)
                copyReplace(workFile, workFile, '##NAME##', simName)
                copyReplace(workFile, workFile, '##COUNTREL##', countRel)
                copyReplace(workFile, workFile, '##VARPAR##', 'Mu')
                copyReplace(workFile, workFile, '##VALUE##', defValues['Mu'])
                copyReplace(workFile, workFile, '##ENTH##', entrainmentTH)
                execCom1Exe(com1Exe, workFile, avaDir, fullOut, logName)
                # save initial particle distribution
                saveInitialParticleDistribution(avaDir, logName, dem)

                # Create dictionary
                reportST = {}
                reportST = {}
                reportST = {'headerLine': {'type': 'title', 'title': 'com1DFAOrig Simulation'},
                'avaName': {'type': 'avaName', 'name': avaDir},
                'simName': {'type': 'simName', 'name': logName},
                'time': {'type': 'time', 'time': dateTimeInfo},
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': relName,
                        'Release Area': relDict['Name'],
                        'Entrainment': entInfo,
                        'Resistance': resInfo,
                        'Parameter variation on': '',
                        'Parameter value': '',
                        'Mu': defValues['Mu'],
                        'Release thickness [m]': relDict['thickness'],
                        'Entrainment thickness [m]': float(entrainmentTH)},
                    'Release Area': {'type': 'columns', 'Release area scenario': relName, 'Release features': relDict['Name'], 'Release thickness [m]': relDict['thickness']}}

                if entInfo == 'Yes':
                    reportST.update({'Entrainment area': {'type': 'columns', 'Entrainment area scenario': entrainmentArea, 'Entrainment thickness [m]': float(entrainmentTH)}})
                if resInfo == 'Yes':
                    reportST.update({'Resistance area': {'type': 'columns', 'Resistance area scenario': resistanceArea}})

                endTime = time.time()
                timeNeeded =  '%.2f' % (endTime - startTime)
                log.info(('Took %s seconds to calculate.' % (timeNeeded)))

                reportST['Simulation Parameters'].update({'Computation time [s]': timeNeeded})

                # Add to report dictionary list
                reportDictList.append(reportST)

                # Count total number of simulations
                countRel = countRel + 1

    # If parameter shall be varied
    if cfgPar.getboolean('parameterVar'):

        # Initialise CreateSimulations cint file and set parameters
        templateFile = os.path.join(modPath, 'CreatenullSimulation.cint')
        workFile = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'CreatenullSimulation.cint')

        # Write required info to cint file
        copyReplace(templateFile, workFile, '##PROJECTDIR##', projDir)
        copyReplace(workFile, workFile, '##BASESIMNAME##', relName)
        execCom1Exe(com1Exe, workFile, avaDir, fullOut)

        # Also perform one standard simulation
        simST = relName + '_null_dfa'
        logName = simST + '_' + defValues[cfgPar['varPar']]
        log.info('Also perform one standard simulation: %s' % simST)
        templateFile = os.path.join(modPath, 'runStandardSimulation.cint')
        workFile = os.path.join(avaDir, 'Work', 'com1DFAOrig', 'runStandardSimulation.cint')
        # Write required info to cint file
        copyReplace(templateFile, workFile, '##PROJECTDIR##', projDir)
        copyReplace(workFile, workFile, '##RESDIR##', resDir)
        copyReplace(workFile, workFile, '##NAME##', simST)
        copyReplace(workFile, workFile, '##COUNTREL##', countRel)
        copyReplace(workFile, workFile, '##VARPAR##', cfgPar['varPar'])
        copyReplace(workFile, workFile, '##VALUE##', defValues[cfgPar['varPar']])
        copyReplace(workFile, workFile, '##ENTH##', entrainmentTH)
        execCom1Exe(com1Exe, workFile, avaDir, fullOut, logName)
        # Create dictionary
        reportNull = {}
        reportNull = {'headerLine': {'type': 'title', 'title': 'com1DFAOrig Simulation'},
        'avaName': {'type': 'avaName', 'name': avaDir},
        'simName': {'type': 'simName', 'name': logName},
        'time': {'type': 'time', 'time': dateTimeInfo},
            'Simulation Parameters': {
                'type': 'list',
                'Release Area Scenario': relName,
                'Release Area': relDict['Name'],
                'Entrainment': 'No',
                'Resistance': 'No',
                'Parameter variation on': '',
                'Parameter value': '',
                'Mu': defValues['Mu'],
                'Release thickness [m]': relDict['thickness'],
                'Entrainment thickness [m]': float(entrainmentTH)},
            'Release area': {'type': 'columns', 'Release area scenario': relName, 'Release features': relDict['Name'], 'Release thickness [m]': relDict['thickness']}}

        endTime = time.time()
        timeNeeded = '%.2f' % (endTime - startTime)
        log.info(('Took %s seconds to calculate.' % (timeNeeded)))
        reportNull['Simulation Parameters'].update({'Computation time [s]': timeNeeded})

        # Add to report dictionary list
        reportDictList.append(reportNull)

        # Count total number of simulations
        countRel = countRel + 1


    log.debug('Avalanche Simulations performed')

    # Setup input from com1DFAOrig and exort to Outputs/com1DFAOrig
    if cfgPar.getboolean('parameterVar'):
        fU.exportcom1DFAOrigOutput(avaDir, cfgPar, cfgGen.getboolean('addTSteps'))
    else:
        cfgParHere = ''
        fU.exportcom1DFAOrigOutput(avaDir, cfgParHere, cfgGen.getboolean('addTSteps'))

    log.info('Exported results to Outputs/com1DFAOrig')

    return reportDictList


def saveInitialParticleDistribution(avaDir, simName, dem):
    x = np.empty(0)
    y = np.empty(0)
    z = np.empty(0)
    m = np.empty(0)
    DEM = IOf.readRaster(dem)
    header = DEM['header']
    # Read log file
    fileName = os.path.join(os.getcwd(), avaDir, 'Outputs', 'com1DFAOrig', 'start%s.log' % (simName))
    with open(fileName, 'r') as file:
        for line in file:
            if "IPD" in line:
                ltime = line.split(', ')
                x = np.append(x, float(ltime[1]))
                y = np.append(y, float(ltime[2]))
                z = np.append(z, float(ltime[3]))
                m = np.append(m, float(ltime[4]))

    x = x - header['xllcenter']
    y = y - header['yllcenter']
    particles = {'t': 0.0, 'x': x, 'y': y, 'z': z, 'm': m}

    partDit = os.path.join(os.getcwd(), avaDir, 'Outputs', 'com1DFAOrig', 'particles', simName)
    fU.makeADir(partDit)
    savePartToPickle(particles, partDit)


def savePartToPickle(dictList, outDir):
    """ Save each dictionary from a list to a pickle in outDir; works also for one dictionary instead of list
        Parameters
        ---------
        dictList: list or dict
            list of dictionaries or single dictionary
        outDir: str
            path to output directory
    """

    if isinstance(dictList, list):
        for dict in dictList:
            pickle.dump(dict, open(os.path.join(outDir, "particles%09.4f.p" % (dict['t'])), "wb"))
    else:
        pickle.dump(dictList, open(os.path.join(outDir, "particles%09.4f.p" % (dictList['t'])), "wb"))
