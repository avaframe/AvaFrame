"""
    Directory and file handling helper functions

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


def makeADir(dirName, flagRemDir=False):
    """ Make directories """

    # If directory already exists - Delete directory first is default
    if os.path.isdir(dirName):
        if flagRemDir:
            log.warning('Be careful directories in %s already existed - but got now deleted' % (dirName))
            shutil.rmtree(dirName, ignore_errors=True)
            os.makedirs(dirName)
        else:
            log.warning('Be careful directory %s already existed - data saved on top of it' % (dirName))
    else:
        os.makedirs(dirName)
    log.info('Directory: %s created' % dirName)


def readLogFile(avaDir):
    """ Read experiment log file and return an index for all simulations
        and a dictionary that contains all simulations
    """

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

    # Save info to dictionary, add all result parameters that are saved in com1DFA Outputs
    suffix = ['pfd', 'ppr', 'pv', 'fd']
    logDict = {'noAva': noSim, 'simName': simName, 'Mu': Mu, 'suffix': suffix}

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


def getDFAData(avaDir, com1DFAOutput, workDir, suffix, nameDir=''):
    """ Export the required input data from com1DFA output to desired location and option for renaming """

    # Initialise directories
    inputDir = os.path.join(avaDir, 'Outputs', 'com1DFA')

    # Read log file information
    [logDict, indSims] = readLogFile(avaDir)

    # Get number of values
    sNo = len(logDict['noAva'])
    sufNo = len(logDict['suffix'])

    # Path to com1DFA results
    resPath = os.path.join(inputDir, com1DFAOutput)

    countsuf = 0
    for m in range(sufNo):
        if logDict['suffix'][m] == suffix:
            for k in range(sNo):
                if k in indSims:
                    if nameDir != '':
                        shutil.copy('%s%.03f/%s/raster/%s_%s.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                                    logDict['simName'][k], logDict['suffix'][m]),
                                    '%s/%s/%06d.txt' % (workDir, nameDir, countsuf+1))
                    else:
                        shutil.copy2('%s%.03f/%s/raster/%s_%s.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                                    logDict['simName'][k], logDict['suffix'][m]), workDir)
                    countsuf = countsuf + 1


def getRefData(avaDir, outputDir, suffix, nameDir=''):
    """ Grab reference data and save to outputDir

        Inputs:
        avaDir          avalanche directory
        suffix          result parameter abbreviation (e.g. 'ppr')
        outputDir       folder where files should be copied to
    """

    # Input directory and load input datasets
    ava = avaDir.split(os.sep)[1]
    refDir = os.path.join('..', 'benchmarks', ava)

    dataRefFiles = glob.glob(refDir+os.sep + '*%s.asc' % suffix)

    # copy these files to desired working directory for outQuickPlot
    for files in dataRefFiles:
        if nameDir != '':
            shutil.copy(files, '%s/%s/000000.txt' % (outputDir, nameDir))
        else:
            shutil.copy2(files, outputDir)

    # Give status information
    if os.path.isdir(refDir) == False:
        log.error('%s does not exist - no files for reference found' % refDir)
    elif dataRefFiles == []:
        log.error('No files found in %s' % refDir)
    else:
        log.info('Reference files copied from directory: %s' % refDir)
