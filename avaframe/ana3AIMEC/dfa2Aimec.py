"""
    Helper function to export required data from com1DFA to be used by Aimec.
"""

# Load modules
import os
import glob
import logging
import numpy as np
import shutil

# local modules
from avaframe.in3Utils import fileHandlerUtils as fU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def writeAimecPathsFile(cfgSetup, avaDir, dirName='com1DFA'):
    """ Write a pathFile to inform Aimec where its input data is located """

    # Initialise DEM
    inputDir = os.path.join(avaDir, 'Inputs')
    dem = glob.glob(inputDir+os.sep+'*.asc')

    # Load parameters for Aimec postprocessing
    pressureLimit = float(cfgSetup['pressureLimit'])
    domainWidth = float(cfgSetup['domainWidth'])

    # Path to com1DFA output in Aimec format
    workDir = os.path.join(avaDir, 'Work', 'ana3AIMEC', dirName)

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


def extractMBInfo(avaDir, simName=''):
    """ Extract the mass balance info from the log file """

    if simName != '':
        simNames = [simName]
        nameDir = simName
    else:
        # Get info from ExpLog
        nameDir = 'com1DFA'
        logName = os.path.join(avaDir, 'Outputs', 'com1DFA', 'ExpLog.txt')
        logDictExp = fU.readLogFile(logName)
        names = logDictExp['fullName']
        simNames = sorted(set(names), key=lambda s: (s.split("_")[0], s.split("_")[1], s.split("_")[3]))

    # Read mass data from log and save to file for each simulation run
    countFile = 0
    for simName in simNames:
        log.info('This is the simulation name %s for mod com1DFA ' % (simName))

        # Initialise fields
        time = []
        mass = []
        entrMass = []
        indRun = [0]
        countMass = 0
        flagStop = 0

        # Read log file
        fileName = os.path.join(os.getcwd(), avaDir, 'Outputs', 'com1DFA', 'start%s.log' % (simName))
        with open(fileName, 'r') as file:
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
        logDict = {'time': np.asarray(time), 'mass': np.asarray(mass),
                   'entrMass': np.asarray(entrMass)}

        # Write mass balance info files
        for k in range(len(indRun)-1):
            saveName = os.path.join(os.getcwd(), avaDir, 'Work', 'ana3AIMEC', nameDir, 'dfa_mass_balance', '%06d.txt' % (countFile + 1))
            with open(saveName, 'w') as MBFile:
                MBFile.write('time, current, entrained\n')
                for m in range(indRun[k], indRun[k] + indRun[k+1] - indRun[k]-1):
                    MBFile.write('%.02f,    %.06f,    %.06f\n' %
                                 (logDict['time'][m], logDict['mass'][m], logDict['entrMass'][m]))
            log.info('Saved to dfa_mass_balance/%s ' % (os.path.basename(saveName)))
            countFile = countFile + 1


def getMBInfo(avaDir, simName=''):
    """ Get MB info """

    # Get info from ExpLog
    if simName != '':
        mbFile = os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'mass_%s.txt' % simName)
        saveName = os.path.join(avaDir, 'Work', 'ana3AIMEC', simName, 'dfa_mass_balance', '%06d.txt' % (2))
        shutil.copy(mbFile, saveName)
        log.info('%s copied to %s' % (mbFile, saveName))
        countFile = 0

    else:
        mbFiles = glob.glob(os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'mass*.txt'))
        mbNames = sorted(set(mbFiles), key=lambda s: (s.split("_")[1], s.split("_")[2], s.split("_")[4]))

        countFile = 0
        for mFile in mbNames:
            saveName = os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA', 'dfa_mass_balance', '%06d.txt' % (countFile + 1))
            shutil.copy(mFile, saveName)
            log.info('%s copied to %s' % (mFile, saveName))
            countFile = countFile + 1

    return countFile


def creatAimecDirs(avaDir, cfgSetup, dirName='com1DFA'):

    # Create required directories
    workDir = os.path.join(avaDir, 'Work', 'ana3AIMEC', dirName)
    fU.makeADir(workDir)
    flowDepthDir = os.path.join(workDir, 'dfa_depth')
    fU.makeADir(flowDepthDir)
    pressureDir = os.path.join(workDir, 'dfa_pressure')
    fU.makeADir(pressureDir)
    speedDir = os.path.join(workDir, 'dfa_speed')
    fU.makeADir(speedDir)
    massDir = os.path.join(workDir, 'dfa_mass_balance')
    fU.makeADir(massDir)
    log.debug('Aimec Work folders created to start postprocessing data')

    # Write the paths to this files to a file
    writeAimecPathsFile(cfgSetup, avaDir, dirName=dirName)

    return workDir, flowDepthDir, pressureDir, speedDir, massDir



def dfaComp2Aimec(avaDir, cfgSetup):
    """ Exports the required data from com1DFA to be used by Aimec """

    # look for matching simulations
    comModules = cfgSetup['comModules'].split('|')
    refModule = comModules[0]
    compModule = comModules[1]
    log.info('Reference data is from module: %s' % refModule)
    log.info('Comparison data is from module: %s' % compModule)

    # Lead all infos on refernce simulations
    inputDirRef = os.path.join(avaDir, 'Outputs', refModule, 'peakFiles')
    refData = fU.makeSimDict(inputDirRef)

    # Load all infos on comparison module simulations
    inputDirComp = os.path.join(avaDir, 'Outputs', compModule, 'peakFiles')
    compData = fU.makeSimDict(inputDirComp)
    allSims = []

    for countRef, simNameRef in enumerate(refData['simName']):
        for countComp, simNameComp in enumerate(compData['simName']):
            if simNameRef == simNameComp:
                allSims.append(simNameRef)
                # Create required directories
                workDir, flowDepthDir, pressureDir, speedDir, massDir = creatAimecDirs(avaDir, cfgSetup, dirName=simNameRef)
                suffix = {'type' : ['pfd', 'ppr', 'pfv'], 'directory' : ['dfa_depth', 'dfa_pressure', 'dfa_speed']}
                countM = 0
                for suf, dir in zip(suffix['type'], suffix['directory']):
                    if refData['resType'][countRef] == suf and compData['resType'][countComp] == suf:
                        shutil.copy(refData['files'][countRef], '%s/%s/%06d.txt' % (workDir, dir, 1))
                        log.info('%s   : %s/%s/%06d.txt' % (refData['files'][countRef], workDir, dir, 1))
                        shutil.copy(compData['files'][countComp], '%s/%s/%06d.txt' % (workDir, dir, 2))
                        log.info('%s   : %s/%s/%06d.txt' % (compData['files'][countComp], workDir, dir, 2))
                        if countM == 0:
                            extractMBInfo(avaDir, simName=simNameRef)
                            getMBInfo(avaDir, simName=simNameRef)
                            countM = countM + 1

    allSims = sorted(set(allSims))
    return allSims


def mainDfa2Aimec(avaDir, cfgSetup):
    """ Exports the required data from com1DFA to be used by Aimec """

    # Create required directories
    workDir, flowDepthDir, pressureDir, speedDir, massDir= creatAimecDirs(avaDir, cfgSetup)

    # Setup input from com1DFA and export to Work ana3AIMEC
    suffix = {'type' : ['pfd', 'ppr', 'pfv'], 'directory' : ['dfa_depth', 'dfa_pressure', 'dfa_speed']}
    for suf, dir in zip(suffix['type'], suffix['directory']):
        fU.getDFAData(avaDir, workDir, suf, nameDir=dir)

    # Extract mb info
    extractMBInfo(avaDir)
