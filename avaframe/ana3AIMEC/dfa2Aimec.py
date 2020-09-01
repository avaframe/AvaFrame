"""
    This is a helper function to export required data from com1DFA to be used by Aimec.

    This file is part of Avaframe.
"""

# Load modules
import os
import glob
import logging
import numpy as np
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


def extractMBInfo(avaDir):
    """ Extract the mass balance info from the log file """

    # Get info from ExpLog
    logName = os.path.join(avaDir, 'Outputs', 'com1DFA', 'ExpLog.txt')
    logDictExp = fU.readLogFile(logName)
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
        logDict = {'time': np.asarray(time), 'mass': np.asarray(mass),
                   'entrMass': np.asarray(entrMass)}

        # Write mass balance info files
        for k in range(len(indRun)-1):
            with open(os.path.join(os.getcwd(), avaDir, 'Work', 'ana3AIMEC', 'com1DFA', 'dfa_mass_balance', '%06d.txt' % (countFile + 1)), 'w') as MBFile:
                MBFile.write('time, current, entrained\n')
                for m in range(indRun[k], indRun[k] + indRun[k+1] - indRun[k]-1):
                    MBFile.write('%.02f,    %.06f,    %.06f\n' %
                                 (logDict['time'][m], logDict['mass'][m], logDict['entrMass'][m]))
            countFile = countFile + 1


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

    # Setup input from com1DFA and export to Work ana3AIMEC
    suffix = {'type' : ['pfd', 'ppr', 'pv'], 'directory' : ['dfa_depth', 'dfa_pressure', 'dfa_speed']}
    countsuf = 0
    for suf in suffix['type']:
        fU.getDFAData(avaDir, workDir, suf, suffix['directory'][countsuf])
        countsuf = countsuf + 1

    # Write the paths to this files to a file
    writeAimecPathsFile(cfgSetup, avaDir)

    # Extract the MB info
    extractMBInfo(avaDir)
