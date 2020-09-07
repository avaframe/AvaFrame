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

    runCommand = samosAT + ' -offscreen ' + cintFile

    proc = subprocess.Popen(runCommand, shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True)

    # initialise log file to save stoudt
    if simName != '':
        f_log = open(os.path.join(avaDir, 'Outputs', 'com1DFA', 'start%s.log' % (simName)), 'w')
    # FSO--- loop through output
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
        resFiles = []
        resFiles.append('')

    # Initialise entrainment areas
    if flagEnt:
        entFiles = glob.glob(inputDir+os.sep + 'ENT' + os.sep+'*.shp')
        if len(entFiles) < 1:
            log.warning('No entrainment file')
            entFiles.append('')  # Kept this for future enhancements
    else:
        entFiles = []
        entFiles.append('')

    # Initialise DEM
    demFile = glob.glob(inputDir+os.sep+'*.asc')

    # Parameter variation
    if cfgPar.getboolean('flagVarPar'):
        varPar = cfgPar['varPar']
        print('flag set varpar', varPar)
    else:
        varPar = 'Mu'

    # Initialise full experiment log file
    with open(os.path.join(workDir, 'ExpLog.txt'), 'w') as logFile:
            logFile.write("NoOfSimulation,SimulationRunName,%s\n" % varPar)

    # return DEM, first item of release, entrainment and resistance areas
    return demFile[0], relFiles, entFiles[0], resFiles[0]


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

def writeParFile(avaDir, cfgPar, simName):
    """ Write text file to read parameters if parameter variation is used """

    # Values of parameter variations in config file as string
    varPar = cfgPar['varParValues']
    varParList = []
    items = varPar.split('_')
    for item in items:
        varParList.append(float(item))

    # Write parameter values to text file
    varFile = os.path.join(avaDir, 'Work' ,'com1DFA', simName+'_VarPar.txt')
    varF = open(varFile, 'w')
    for vals in varParList:
        # Important write Par in correct sequence from small to big
        varF.write('%.3f\n' % vals)
    varF.close()

    log.info('Parameter variation input file written')


def runSamos(cfg, avaDir):
    """ Run main model"""

    # Setup configuration
    cfgGen = cfg['GENERAL']
    samosAT = cfgGen['samosAT']
    flagEnt = cfgGen.getboolean('flagEnt')
    flagRes = cfgGen.getboolean('flagRes')
    inputf = cfgGen['inputf']
    fullOut = cfgGen.getboolean('flagOut')
    cfgPar = cfg['PARAMATERVAR']
    resDir = os.path.join(avaDir, 'Work', 'com1DFA')
    # Get path of module
    modPath = os.path.dirname(__file__)

    # Log chosen settings
    log.info('The chosen settings: entrainment - %s , resistance - %s ' % (flagEnt, flagRes))

    # Log current avalanche directory
    log.info('Your current avalanche test name: %s' % avaDir)

    # Load input data
    dem, rels, res, ent = initialiseRun(avaDir, flagEnt, flagRes, cfgPar, inputf)

    # Get cell size from DEM header
    demData = aU.readASCheader(dem)
    cellSize = demData.cellsize

    # Counter for release area loop
    countRel = 0

    # Loop through release areas
    for rel in rels:

        # Set release areas and simulation name
        relName = os.path.splitext(os.path.basename(rel))[0]
        simName = relName
        log.info('Release area: %s - perform simulations' % (relName))

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
        copyReplace(workFile, workFile, '##RESFILE##', ent)
        copyReplace(workFile, workFile, '##ENTFILE##', res)
        # Setup Project
        execSamos(samosAT, workFile, avaDir, fullOut)

        # Initialise CreateSimulations cint file and set parameters
        templateFile = os.path.join(modPath, 'CreateSimulations.cint')
        workFile = os.path.join(avaDir, 'Work', 'com1DFA', 'CreateSimulations.cint')
        copyReplace(templateFile, workFile, '##BASEPATH##', os.getcwd())
        copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
        copyReplace(workFile, workFile, '##BASESIMNAME##', simName)
        execSamos(samosAT, workFile, avaDir, fullOut)

        # If mu shall be varied
        if cfgPar.getboolean('flagVarPar'):
            log.info('Parameter variation used, varying: %s' % cfgPar['varPar'])
            writeParFile(avaDir, cfgPar, simName)
            templateFile = os.path.join(modPath, '%s%s.cint' % (cfgPar['varRunCint'], cfgPar['varPar']))
            workFile = os.path.join(avaDir, 'Work', 'com1DFA', '%s%s.cint' % (cfgPar['varRunCint'], cfgPar['varPar']))
            copyReplace(templateFile, workFile, '##BASEPATH##', os.getcwd())
            copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
            copyReplace(workFile, workFile, '##RESDIR##', resDir)
            copyReplace(workFile, workFile, '##COUNTREL##', countRel)
            # Count total number of simulations
            countRel = countRel + 3

        else:
            templateFile = os.path.join(modPath, '%s.cint' % (cfgGen['RunCint']))
            workFile = os.path.join(avaDir, 'Work', 'com1DFA', '%s.cint' % (cfgGen['RunCint']))
            copyReplace(templateFile, workFile, '##BASEPATH##', os.getcwd())
            copyReplace(workFile, workFile, '##PROJECTDIR##', projDir)
            copyReplace(workFile, workFile, '##RESDIR##', resDir)
            copyReplace(workFile, workFile, '##COUNTREL##', countRel)
            # Count total number of simulations
            countRel = countRel + 2

        execSamos(samosAT, workFile, avaDir, fullOut, relName)

    log.info('Avalanche Simulations performed')

    # Setup input from com1DFA and exort to Outputs/com1DFA
    if cfgPar.getboolean('flagVarPar'):
        fU.exportcom1DFAOutput(avaDir, cfgPar)
    else:
        fU.exportcom1DFAOutput(avaDir)

    log.info('Exported results to Outputs/com1DFA')
