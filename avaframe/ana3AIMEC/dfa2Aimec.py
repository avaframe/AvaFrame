"""
    This is a helper function to export required data from com1DFA to be used by Aimec.

    This file is part of Avaframe.
"""

# Load modules
import os
import glob
import logging
import numpy as np
import shutil

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def makeAimecDirs(avaDir):
    """ Make directories where Aimec reads input data from """

    # Set directories for Aimec inputs
    workDir = os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA')
    flowDepthDir = os.path.join(workDir, 'dfa_depth')
    pressureDir = os.path.join(workDir, 'dfa_pressure')
    massDir = os.path.join(workDir, 'dfa_mass_balance')
    massDirTemp = os.path.join(workDir, 'dfa_mass_balance_temp')


    if os.path.isdir(workDir):
        log.warning('Be careful directories in %s already existed - but got now deleted' % (workDir))
        shutil.rmtree(workDir, ignore_errors=True)

    os.makedirs(workDir)
    os.makedirs(flowDepthDir)
    os.makedirs(pressureDir)
    os.makedirs(massDir)
    os.makedirs(massDirTemp)

    log.info('Aimec Work folders created to start postprocessing com1DFA data')


def writeAimecPathsFile(cfgSetup, avaDir):
    """ Write a pathFile to inform Aimec where its input data is located """

    # Initialise DEM
    inputDir = os.path.join(avaDir, 'Inputs')
    dem = glob.glob(inputDir+os.sep+'*.asc')

    # Load parameters for Aimec postprocessing
    pressureLimit = float(cfgSetup['pressureLimit'])
    domainWidth = float(cfgSetup['domainWidth'])

    # Path to com1DFA output in Aimec format
    workDir = os.path.join(avaDir, 'Work', 'ana3AIMEC')

    # Create empty variable
    emptyVar = ""

    with open(os.path.join(workDir, 'aimecPathFile.txt'), 'w') as pfile:
        pfile.write('pathPressure=%s,\n' % (os.path.join(workDir, 'dfa_pressure')))
        pfile.write('pathMass=%s,\n' % (os.path.join(workDir, 'dfa_mass_balance')))
        pfile.write('pathDocDamage=%s,\n' % (emptyVar))
        pfile.write('pathDocRadar=%s,\n' % (emptyVar))
        pfile.write('pathNumInfo=%s,\n' % (emptyVar))
        pfile.write('pathFlowHeight=%s,\n' % (os.path.join(workDir, 'dfa_depth')))
        pfile.write('pathAvalanchePath=%s,\n' % (os.path.join(inputDir, 'avalanche_path.xyz')))
        pfile.write('calcPressureLimit=%s,\n' % (pressureLimit))
        pfile.write('domainWidth=%s,\n' % (domainWidth))
        pfile.write('pathDepoArea=%s,\n' % (emptyVar))
        pfile.write('pathAOI=%s,\n' % (emptyVar))
        pfile.write('pathDHM=%s,\n' % (dem[0]))
        pfile.write('pathEnergy=%s,\n' % (emptyVar))
        pfile.write('pathVelocity=%s,\n' % (emptyVar))
        pfile.write('pathResult=%s,\n' % (os.path.join(workDir, 'AimecResults')))


def extractMBInfo(avaDir):
    """ Extract the mass balance info from the log file """

    # Get release area names
    inputDir = os.path.join(avaDir, 'Inputs')
    relFiles = glob.glob(inputDir+os.sep + 'REL'+os.sep + '*.shp')
    relNames = []
    for rels in relFiles:
        relNames.append(os.path.splitext(os.path.basename(rels))[0])

    # Get logFile
    [logDictExp, indSims] = readLogFile(avaDir)
    simName = []
    for name in logDictExp['simName']:
        simName.append(name.split('_')[0])
    relNames = set(simName)


    # Read mass data from log and save to file for each simulation run
    countFile = 0
    for relName in relNames:
        log.info('These are the release areas: %s ' % (relName))

        # Initialise fields
        time = []
        mass = []
        entrMass = []
        indRun = [0]
        countMass = 0
        flagStop = 0

        # Read log file
        with open(os.path.join(os.getcwd(), avaDir, 'Outputs', 'com1DFA', 'start%s.log' % (relName)), 'r') as file:
            for line in file:
                if "computing time step" in line:
                    ltime = line.split()[3]
                    timeNum = ltime.split('...')[0]
                    time.append(float(timeNum))
                elif "entrained DFA mass" in line:
                    entrMass.append(float(line.split()[3]))
                elif "total DFA mass" in line:
                    mass.append(float(line.split()[3]))
                    countMass = countMass + 1
                elif "terminated" in line:
                    indRun.append(countMass)
        # Save to dictionary
        logDict = {'time' : np.asarray(time), 'mass' : np.asarray(mass), 'entrMass' : np.asarray(entrMass)}

        # Write mass balance info files
        for k in range(len(indRun)-1):
            with open(os.path.join(os.getcwd(), avaDir, 'Work','ana3AIMEC', 'com1DFA', 'dfa_mass_balance_temp', '%06d.txt' % (countFile + 1)), 'w') as MBFile:
                MBFile.write('time, current, entrained\n')
                for m in range(indRun[k], indRun[k] + indRun[k+1] - indRun[k]-1):
                    MBFile.write('%.02f,    %.06f,    %.06f\n' % (logDict['time'][m], logDict['mass'][m], logDict['entrMass'][m]))
            countFile = countFile + 1

    # Delete the files that are not in the local Exp Log
    countSims = 0
    for l in range(len(logDictExp['simName'])):
        if l in indSims:
            fname = ('%06d.txt' % (l+1))
            fnameNew = ('%06d.txt' % (countSims+1))
            shutil.copyfile(os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA', 'dfa_mass_balance_temp', fname),
             os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA', 'dfa_mass_balance', fnameNew))
            countSims = countSims + 1
            log.info('NEW files: %s new is %s' % (fname, fnameNew))

    shutil.rmtree(os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA', 'dfa_mass_balance_temp'))

def readLogFile(avaDir):
    """ Read experiment log file """

    # Initialise directories
    inputDir = os.path.join(avaDir, 'Outputs', 'com1DFA')
    workDirMain = os.path.join(avaDir, 'Work', 'ana3AIMEC')

    logFile = open(os.path.join(inputDir, 'ExpLog.txt'), 'r')
    log.info('Take com1DFA full experiment log')

    noSim = []      # number of Simulation
    simName = []    # name of Simulation
    Mu = []         # Mu parameter value

    lines = logFile.readlines()[1:]
    for line in lines:
        vals = line.strip().split()
        noSim.append(float(vals[0]))
        simName.append(vals[1])
        Mu.append(float(vals[2]))

    # Save info to dictionary
    suffix = ['pfd', 'ppr', 'pv', 'fd']
    logDict = {'noAva' : noSim, 'simName' : simName, 'Mu' : Mu, 'suffix' : suffix}

    # Read the experiment log - if copied to local_ExpLog take this!
    if os.path.isfile(os.path.join(workDirMain, 'local_ExpLog.txt')):
        logFileLocal = open(os.path.join(workDirMain, 'local_ExpLog.txt'), 'r')
        log.info('Take local (potentially modified) experiment log')
    else:
        logFileLocal = open(os.path.join(inputDir, 'ExpLog.txt'), 'r')
        log.warning('There is no file local_ExpLog - using all simulations')

    # Read simulation names from local exp Log
    simName = []
    lines = logFileLocal.readlines()[1:]
    for line in lines:
        vals = line.strip().split()
        simName.append(vals[1])

    # Identify common simulations
    setSim = set(simName)
    indSims = [i for i, item in enumerate(logDict['simName']) if item in setSim]

    return logDict, indSims


def getDFAData(avaDir, cfgDFA):
    """ Export the required input data from com1DFA output """

    # Initialise directories
    inputDir = os.path.join(avaDir, 'Outputs', 'com1DFA')
    workDir = os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA')
    workDirMain = os.path.join(avaDir, 'Work', 'ana3AIMEC')

    # Read log file information
    [logDict, indSims] = readLogFile(avaDir)
    # Get number of values
    sNo = len(logDict['noAva'])
    sufNo = len(logDict['suffix'])

    # Path to com1DFA results
    resPath = os.path.join(inputDir, cfgDFA['filesDir'])

    countpfd = 0
    countppr = 0

    for m in range(sufNo):
        if logDict['suffix'][m] == 'pfd':
            for k in range(sNo):
                if k in indSims:
                    shutil.copy('%s%.03f/%s/raster/%s_%s.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                                          logDict['simName'][k], logDict['suffix'][m]),
                                          '%s/dfa_depth/%06d.txt' % (workDir, countpfd+1))
                    # log.info('%s%f/%s/raster/%s_%s.asc to the new file %s/dfa_depth/%06d.txt' % (resPath,
                    #                     logDict['Mu'][k], logDict['simName'][k],
                    #                     logDict['simName'][k], logDict['suffix'][m], outputDir, countpfd))
                    countpfd = countpfd + 1

        elif logDict['suffix'][m] == 'ppr':
            for k in range(sNo):
                if k in indSims:
                    shutil.copy('%s%.03f/%s/raster/%s_%s.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                                          logDict['simName'][k], logDict['suffix'][m]),
                                          '%s/dfa_pressure/%06d.txt' % (workDir, countppr+1))
                    log.info('Simulation %s is copied to ana3AIMEC' % logDict['simName'][k])
                    countppr = countppr + 1


def mainDfa2Aimec(avalancheDir, cfgDFA, cfgSetup):
    """ Exports the required data from com1DFA to be used by Aimec """
    # Setup input from com1DFA
    makeAimecDirs(avalancheDir)
    getDFAData(avalancheDir, cfgDFA)
    writeAimecPathsFile(cfgSetup, avalancheDir)
    extractMBInfo(avalancheDir)
