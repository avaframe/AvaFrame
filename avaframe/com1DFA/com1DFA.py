"""
    This is a python wrapper to execute a compiled SamosAT file

    This file is part of Avaframe.
"""

# Load modules
import os
import sys
import glob
import subprocess
import shutil
import numpy as np
import logging
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import ascUtils as aU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def execSamos(samosAT, cintFile, avaDir, fullOut=False, simName=''):
    """ Execute compiled SamosAT file using cintFile to set configuration
        and run options """

    # define command line
    runCommand = samosAT + ' -offscreen ' + cintFile

    # initialise log file to save stoudt
    if simName != '':
        f_log = open(os.path.join(avaDir, 'Outputs', 'com1DFA', 'start%s.log' % (simName)), 'w')

    # Call command
    proc = subprocess.Popen(runCommand, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)

    # loop through output and save to logFile if file provided
    for line in proc.stdout:
        if simName != '':
            f_log.write(line)
        if fullOut:
            log.info((line.rstrip()))
        elif 'BatchSamos' in line:
            log.info(line.rstrip())
        elif 'error' in line:
            log.info(line.rstrip())

    # make process wait for previous process to finish
    reVal = proc.wait()


def initialiseRun(avaDir, flagEnt, flagRes, cfgPar, inputf='shp'):
    """ Initialise Simulation Run with input data """

    # Set directories for inputs, outputs and current work
    inputDir = os.path.join(avaDir, 'Inputs')
    outputDir = os.path.join(avaDir, 'Outputs', 'com1DFA')
    fU.makeADir(outputDir)
    workDir = os.path.join(avaDir, 'Work', 'com1DFA')
    # If Work directory already exists - error
    if os.path.isdir(workDir):
        log.error('Work directory %s already exists - delete first!' % (workDir))
    else:
        os.makedirs(workDir)
    log.info('Directory: %s created' % workDir)

    # Set flag if there is an entrainment or resistance area
    flagEntRes = False

    # Initialise release areas, default is to look for shapefiles
    if inputf == 'nxyz':
        relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.nxyz')
    else:
        relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.shp')
        log.info('Release area files are: %s' % relFiles)

    # Initialise resistance areas
    if flagRes:
        resFiles = glob.glob(inputDir+os.sep + 'RES' + os.sep+'*.shp')
        if len(resFiles) < 1:
            log.warning('No resistance file')
            resFiles.append('')  # Kept this for future enhancements
        else:
            flagEntRes = True
    else:
        resFiles = []
        resFiles.append('')

    # Initialise entrainment areas
    if flagEnt:
        entFiles = glob.glob(inputDir+os.sep + 'ENT' + os.sep+'*.shp')
        if len(entFiles) < 1:
            log.warning('No entrainment file')
            entFiles.append('')  # Kept this for future enhancements
        else:
            flagEntRes = True
    else:
        entFiles = []
        entFiles.append('')

    # Initialise DEM
    demFile = glob.glob(inputDir+os.sep+'*.asc')

    # Parameter variation
    if cfgPar.getboolean('flagVarPar'):
        varPar = cfgPar['varPar']
    else:
        varPar = 'Mu'

    # Initialise full experiment log file
    with open(os.path.join(workDir, 'ExpLog.txt'), 'w') as logFile:
        logFile.write("NoOfSimulation,SimulationRunName,%s\n" % varPar)

    if flagEntRes:
        log.info('Simulations are performed using entrainment and resistance')
    else:
        log.info('Standard simulation is performed without entrainment and resistance')

    # return DEM, first item of release, entrainment and resistance areas
    return demFile[0], relFiles, entFiles[0], resFiles[0], flagEntRes


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


def runSamos(cfg, avaDir):
    """ Run main model"""

    # Setup configuration
    cfgGen = cfg['GENERAL']
    samosAT = cfgGen['samosAT']
    flagEnt = cfgGen.getboolean('flagEnt')
    flagRes = cfgGen.getboolean('flagRes')
    inputf = cfgGen['inputf']
    fullOut = cfgGen.getboolean('flagOut')
    cfgPar = cfg['PARAMETERVAR']
    resDir = os.path.join(avaDir, 'Work', 'com1DFA')
    # Get path of module
    modPath = os.path.dirname(__file__)
    # Standard values for parameters that can be varied
    defValues = cfg['DEFVALUES']

    # Log chosen settings
    log.info('The chosen settings: entrainment - %s , resistance - %s ' % (flagEnt, flagRes))

    # Log current avalanche directory
    log.info('Your current avalanche test name: %s' % avaDir)

    # Load input data
    dem, rels, ent, res, flagEntRes = initialiseRun(avaDir, flagEnt, flagRes, cfgPar, inputf)
    entrainmentArea = ''
    resistanceArea = ''
    if flagEntRes:
        entrainmentArea = os.path.splitext(os.path.basename(ent))[0]
        resistanceArea = os.path.splitext(os.path.basename(res))[0]

    # Get cell size from DEM header
    demData = aU.readASCheader(dem)
    cellSize = demData.cellsize

    # Counter for release area loop
    countRel = 0

    # Setup simulation dictionaries for report genereation
    reportDictList = []

    # Loop through release areas
    for rel in rels:

        # Set release areas and simulation name
        relName = os.path.splitext(os.path.basename(rel))[0]
        simName = relName
        log.info('Release area: %s - perform simulations' % (relName))
        if flagEntRes:
            log.info('Entrainment area: %s and resistance area: %s' % (entrainmentArea, resistanceArea))

        # Initialise CreateProject cint file
        templateFile = os.path.join(modPath, 'CreateProject.cint')
        workFile = os.path.join(avaDir, 'Work', 'com1DFA', 'CreateProject.cint')
        projDir = os.path.join(avaDir, 'Work', 'com1DFA', simName)
        demName = os.path.splitext(os.path.basename(dem))[0]

        # Set Parameters in cint file
        copyReplace(templateFile, workFile, '##BASEPATH##', avaDir)
        copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
        copyReplace(workFile, workFile, '##DHMFILE##', dem)
        copyReplace(workFile, workFile, '##DHMNAME##', demName)
        copyReplace(workFile, workFile, '##CELLSIZE##', cellSize)
        copyReplace(workFile, workFile, '##RELFILE##', rel)
        copyReplace(workFile, workFile, '##ENTFILE##', ent)
        copyReplace(workFile, workFile, '##RESFILE##', res)
        # Setup Project
        execSamos(samosAT, workFile, avaDir, fullOut)

        if flagEntRes:
            # Initialise CreateSimulations cint file and set parameters
            templateFile = os.path.join(modPath, 'CreateSimulations.cint')
            workFile = os.path.join(avaDir, 'Work', 'com1DFA', 'CreateSimulations.cint')
            cuSim = [simName + '_entres_dfa', simName + '_null_dfa']
        else:
            # Initialise CreateSimulations cint file and set parameters
            templateFile = os.path.join(modPath, 'CreateBasicSimulation.cint')
            workFile = os.path.join(avaDir, 'Work', 'com1DFA', 'CreateBasicSimulation.cint')
            cuSim = [simName + '_null_dfa']

        # Write required info to cint file
        copyReplace(templateFile, workFile, '##BASEPATH##', os.getcwd())
        copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
        copyReplace(workFile, workFile, '##BASESIMNAME##', simName)
        execSamos(samosAT, workFile, avaDir, fullOut)

        # If parameter shall be varied
        if cfgPar.getboolean('flagVarPar'):

            # Also perform one standard simulation
            simST = simName + '_null_dfa'
            logName = simST + '_' + defValues[cfgPar['varPar']]
            log.info('Also perform one standard simulation: %s' % simST)
            templateFile = os.path.join(modPath, 'runBasicST.cint')
            workFile = os.path.join(avaDir, 'Work', 'com1DFA', 'runBasicST.cint')
            # Write required info to cint file
            copyReplace(templateFile, workFile, '##BASEPATH##', os.getcwd())
            copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
            copyReplace(workFile, workFile, '##RESDIR##', resDir)
            copyReplace(workFile, workFile, '##NAME##', simST)
            copyReplace(workFile, workFile, '##COUNTREL##', countRel)
            copyReplace(workFile, workFile, '##VARPAR##', cfgPar['varPar'])
            copyReplace(workFile, workFile, '##VALUE##', defValues[cfgPar['varPar']])
            execSamos(samosAT, workFile, avaDir, fullOut, logName)

            # Create dictionary
            reportNull = {}
            reportNull = {'headerLine': {'type': 'title', 'title': 'com1DFA Simulation'},
            'simName': {'type': 'simName', 'name': logName},
                'Simulation Parameters': {
                    'type': 'list',
                    'Release Area': relName,
                    'Entrainment Area': '',
                    'Resistance Area': '',
                    'Parameter variation on': '',
                    'Parameter value': '',
                    'Mu': defValues['Mu'],
                    'Release thickness [m]': defValues['RelTh']},
                'Release area': {'type': 'columns', 'Release area scenario': relName}}

            # Add to report dictionary list
            reportDictList.append(reportNull)

            # Count total number of simulations
            countRel = countRel + 1

            if cfgPar.getboolean('flagVarEnt') and (simName + '_entres_dfa') in cuSim:
                sim = simName + '_entres_dfa'
            else:
                sim = simName + '_null_dfa'

            log.info('Parameter variation used, varying: %s' % cfgPar['varPar'])

            # Values of parameter variations in config file as string
            varParValues = cfgPar['varParValues']
            itemsRaw = varParValues.split('_')
            items = []
            for itemR in itemsRaw:
                items.append('%.3f' % float(itemR))
            for item in items:
                logName = sim + '_' + item
                log.info('Perform simulation with %s = %s: logName = %s' % (cfgPar['varPar'], item, logName))
                templateFile = os.path.join(modPath, '%s%sBasic.cint' % (cfgPar['varRunCint'], cfgPar['varPar']))
                workFile = os.path.join(avaDir, 'Work', 'com1DFA',
                                        '%s%sBasic.cint' % (cfgPar['varRunCint'], cfgPar['varPar']))
                copyReplace(templateFile, workFile, '##BASEPATH##', os.getcwd())
                copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
                copyReplace(workFile, workFile, '##RESDIR##', resDir)
                copyReplace(workFile, workFile, '##NAME##', sim)
                copyReplace(workFile, workFile, '##COUNTREL##', countRel)
                copyReplace(workFile, workFile, '##VALUE##', item)
                execSamos(samosAT, workFile, avaDir, fullOut, logName)

                # Create dictionary
                reportVar = {}
                reportVar = {'headerLine': {'type': 'title', 'title': 'com1DFA Simulation'},
                'simName' : {'type': 'simName', 'name': logName},
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': relName,
                        'Entrainment Area': entrainmentArea,
                        'Resistance Area': resistanceArea,
                        'Parameter variation on': cfgPar['varPar'],
                        'Parameter value': item},
                    'Release area': {'type': 'columns', 'Release area scenario': relName},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': entrainmentArea},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': resistanceArea}}

                if cfgPar['varPar'] == 'RelTh':
                    reportVar['Simulation Parameters'].update({'Mu': defValues['Mu']})
                    reportVar['Simulation Parameters'].update({'Release thickness [m]': item})
                elif cfgPar['varPar'] == 'Mu':
                    reportVar['Simulation Parameters'].update({'Release thickness [m]': defValues['RelTh']})
                    reportVar['Simulation Parameters'].update({'Mu': item})

                # Add to report dictionary list
                reportDictList.append(reportVar)

                # Count total number of simulations
                countRel = countRel + 1

        else:
            for sim in cuSim:
                templateFile = os.path.join(modPath, 'runBasicST.cint')
                workFile = os.path.join(avaDir, 'Work', 'com1DFA', 'runBaiscST.cint')
                logName = sim + '_' + defValues['Mu']
                # Write required info to cint file
                copyReplace(templateFile, workFile, '##BASEPATH##', os.getcwd())
                copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
                copyReplace(workFile, workFile, '##RESDIR##', resDir)
                copyReplace(workFile, workFile, '##NAME##', sim)
                copyReplace(workFile, workFile, '##COUNTREL##', countRel)
                copyReplace(workFile, workFile, '##VARPAR##', 'Mu')
                copyReplace(workFile, workFile, '##VALUE##', defValues['Mu'])
                execSamos(samosAT, workFile, avaDir, fullOut, logName)

                # Create dictionary
                reportST = {}
                reportST = {}
                reportST = {'headerLine': {'type': 'title', 'title': 'com1DFA Simulation'},
                'simName': {'type': 'simName', 'name': logName},
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': relName,
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Parameter variation on': '',
                        'Parameter value': '',
                        'Mu': defValues['Mu'],
                        'Release thickness [m]': defValues['RelTh']},
                    'Release Area': {'type': 'columns', 'Release area scenario': relName}}

                if 'entres' in sim:
                    reportST['Simulation Parameters'].update({'Entrainment Area': entrainmentArea})
                    reportST['Simulation Parameters'].update({'Resistance Area': resistanceArea})
                    reportST.update({'Entrainment area': {'type': 'columns', 'Entrainment area scenario': entrainmentArea}})
                    reportST.update({'Resistance area': {'type': 'columns', 'Resistance area scenario': resistanceArea}})

                # Add to report dictionary list
                reportDictList.append(reportST)

                # Count total number of simulations
                countRel = countRel + 1

    log.info('Avalanche Simulations performed')

    # Setup input from com1DFA and exort to Outputs/com1DFA
    if cfgPar.getboolean('flagVarPar'):
        fU.exportcom1DFAOutput(avaDir, cfgPar)
    else:
        fU.exportcom1DFAOutput(avaDir)

    log.info('Exported results to Outputs/com1DFA')

    return reportDictList
