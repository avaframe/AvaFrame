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
    workDir = os.path.join(avaDir, 'Work/ana3AIMEC/com1DFA')
    flowDepthDir = os.path.join(workDir, 'dfa_depth')
    pressureDir = os.path.join(workDir, 'dfa_pressure')
    massDir = os.path.join(workDir, 'dfa_mass_balance')

    if os.path.isdir(workDir):
        log.warning('Be careful Work/ana3AIMEC/com1DFA directories already existed - but got now deleted')
        shutil.rmtree(workDir, ignore_errors=True)

    os.makedirs(workDir)
    os.makedirs(flowDepthDir)
    os.makedirs(pressureDir)
    os.makedirs(massDir)
    log.info('Aimec Work folders created to start postprocessing com1DFA data')

<<<<<<< HEAD

=======
>>>>>>> Add function to export output data for postprocessing with aimec
def writeAimecPathsFile(cfgSetup, avaDir):
    """ Write a pathFile to inform Aimec where its input data is located """

    # Initialise DEM
    inputDir = os.path.join(avaDir, 'Inputs')
    dem = glob.glob(inputDir+os.sep+'*.asc')

    # Load parameters for Aimec postprocessing
    pressureLimit = float(cfgSetup['pressureLimit'])
    domainWidth = float(cfgSetup['domainWidth'])

    # Path to com1DFA output in Aimec format
    workDir = os.path.join(avaDir, 'Work/ana3AIMEC')

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


def getDFAData(avaDir, cfgDFA):
    """ Export the required input data from com1DFA output """

    # Initialise directories
    inputDir = os.path.join(avaDir, 'Outputs/com1DFA/')
    workDir = os.path.join(avaDir, 'Work/ana3AIMEC/com1DFA')
    workDirMain = os.path.join(avaDir, 'Work/ana3AIMEC')

    noSim = []      # number of Simulation
    simName = []    # name of Simulation
    Mu = []         # Mu parameter value
    # Read the experiment log - if copied to local_ExpLog take this!
    if os.path.isfile(os.path.join(workDirMain, 'local_ExpLog.txt')):
        logFile = open(os.path.join(workDirMain, 'local_ExpLog.txt'), 'r')
        log.info('Take local (potentially modified) experiment log')
    else:
        logFile = open(os.path.join(inputDir, 'ExpLog.txt'), 'r')
        log.info('Take com1DFA full experiment log')

    lines = logFile.readlines()[1:]
    for line in lines:
        vals = line.strip().split()
        noSim.append(float(vals[0]))
        simName.append(vals[1])
        Mu.append(float(vals[2]))

    # Save info to dictionary
    suffix = ['pfd', 'ppr', 'pv', 'fd']
    logDict = {'noAva' : noSim, 'simName' : simName, 'Mu' : Mu, 'suffix' : suffix}
    sNo = len(noSim)
    sufNo = len(suffix)

    # Path to com1DFA results
    resPath = os.path.join(inputDir, cfgDFA['filesDir'])

    countpfd = 0
    countppr = 0

    for m in range(sufNo):
        if logDict['suffix'][m] == 'pfd':
            for k in range(sNo):
                shutil.copy('%s%.03f/%s/raster/%s_%s.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                                      logDict['simName'][k], logDict['suffix'][m]),
                                      '%s/dfa_depth/%06d.txt' % (workDir, countpfd))
                # log.info('%s%f/%s/raster/%s_%s.asc to the new file %s/dfa_depth/%06d.txt' % (resPath,
                #                     logDict['Mu'][k], logDict['simName'][k],
                #                     logDict['simName'][k], logDict['suffix'][m], outputDir, countpfd))
                countpfd = countpfd + 1

        elif logDict['suffix'][m] == 'ppr':
            for k in range(sNo):
                shutil.copy('%s%.03f/%s/raster/%s_%s.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                                      logDict['simName'][k], logDict['suffix'][m]),
                                      '%s/dfa_pressure/%06d.txt' % (workDir, countppr))
                # log.info('%s%f/%s/raster/%s_%s.asc to the new file %s/dfa_pressure/%06d.txt' % (resPath,
                #                     logDict['Mu'][k], logDict['simName'][k],
                #                     logDict['simName'][k], logDict['suffix'][m], outputDir, countppr))
                countppr = countppr + 1
