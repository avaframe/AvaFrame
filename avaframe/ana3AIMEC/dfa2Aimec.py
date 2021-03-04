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


def extractMBInfo(avaDir, countStart=0):
    """ Extract the mass balance info from the log file """

    # Get info from ExpLog
    logName = os.path.join(avaDir, 'Outputs', 'com1DFA', 'ExpLog.txt')
    logDictExp = fU.readLogFile(logName)
    names = logDictExp['fullName']
    simNames = sorted(set(names), key=lambda s: (s.split("_")[0], s.split("_")[1], s.split("_")[3]))

    # Read mass data from log and save to file for each simulation run
    countFile = countStart
    for simName in simNames:
        log.info('This is the simulation name: %s ' % (simName))

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
            saveName = os.path.join(os.getcwd(), avaDir, 'Work', 'ana3AIMEC', 'com1DFA', 'dfa_mass_balance', '%06d.txt' % (countFile + 1))
            with open(saveName, 'w') as MBFile:
                MBFile.write('time, current, entrained\n')
                for m in range(indRun[k], indRun[k] + indRun[k+1] - indRun[k]-1):
                    MBFile.write('%.02f,    %.06f,    %.06f\n' %
                                 (logDict['time'][m], logDict['mass'][m], logDict['entrMass'][m]))
            log.info('Saved to dfa_mass_balance/%s ' % (os.path.basename(saveName)))
            countFile = countFile + 1

    return countFile

def getMBInfo(avaDir, countStart=0):
    """ Get MB info """

    # Get info from ExpLog
    mbFiles = glob.glob(os.path.join(avaDir, 'Outputs', 'com1DFAPy', 'mass*.txt'))
    mbNames = sorted(set(mbFiles), key=lambda s: (s.split("_")[1], s.split("_")[2], s.split("_")[4]))

    countFile = countStart
    for mFile in mbNames:
        saveName = os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA', 'dfa_mass_balance', '%06d.txt' % (countFile + 1))
        shutil.copy(mFile, saveName)
        log.info('%s copied to %s' % (mFile, saveName))
        countFile = countFile + 1

    return countFile


def mainDfa2Aimec(avaDir, cfgSetup):
    """ Exports the required data from com1DFA to be used by Aimec """

    # Create required directories
    workDir = os.path.join(avaDir, 'Work', 'ana3AIMEC', 'com1DFA')
    fU.makeADir(workDir)
    flowDepthDir = os.path.join(workDir, 'dfa_depth')
    fU.makeADir(flowDepthDir)
    pressureDir = os.path.join(workDir, 'dfa_pressure')
    fU.makeADir(pressureDir)
    speedDir = os.path.join(workDir, 'dfa_speed')
    fU.makeADir(speedDir)
    massDir = os.path.join(workDir, 'dfa_mass_balance')
    fU.makeADir(massDir)
    log.info('Aimec Work folders created to start postprocessing com1DFA data')

    # Write the paths to this files to a file
    writeAimecPathsFile(cfgSetup, avaDir)

    comModules = cfgSetup['comModules'].split('|')
    countStart = 0
    countFile = 0
    for comMod in comModules:

        # Setup input from com1DFA and export to Work ana3AIMEC
        suffix = {'type' : ['pfd', 'ppr', 'pfv'], 'directory' : ['dfa_depth', 'dfa_pressure', 'dfa_speed']}
        for suf, dir in zip(suffix['type'], suffix['directory']):
            nFiles = fU.getDFAData(avaDir, workDir, suf, comModule=comMod, nameDir=dir, countStart=countStart)
        countStart = nFiles / len(suffix['type'])

        # Extract the MB info
        if comMod == 'com1DFA':
            countFile = extractMBInfo(avaDir, countStart=countFile)
        elif comMod == 'com1DFAPy':
            countFile = getMBInfo(avaDir, countStart=countFile)
