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

# Local imports
import avaframe.in3Utils.ascUtils as IOf


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def makeADir(dirName):
    """ Make directories """

    # If directory already exists - Delete directory first is default
    if os.path.isdir(dirName):
        log.warning('Be careful directory %s already existed - data saved on top of it' % (dirName))
    else:
        os.makedirs(dirName)
    log.info('Directory: %s created' % dirName)


def readLogFile(logName):
    """ Read experiment log file and make dictionary that contains all simulations """

    # Read log file
    logFile = open(logName, 'r')
    log.info('Take com1DFA full experiment log')

    # Save info to dictionary, add all result parameters that are saved in com1DFA Outputs
    suffix = ['pfd', 'ppr', 'pv', 'fd']
    logDict = {'noSim': [], 'simName': [], 'Mu': [], 'suffix': suffix}

    lines = logFile.readlines()[1:]
    countSims = 1
    for line in lines:
        vals = line.strip().split()
        logDict['noSim'].append(countSims)
        logDict['simName'].append(vals[1])
        logDict['Mu'].append(float(vals[2]))
        countSims = countSims + 1

    return logDict

#
def checkCommonSims(logName, localLogName):
    """ Check which files are common between local and full ExpLog """

    if os.path.isfile(localLogName) == False:
        localLogName = logName

    # Read log files and extract info
    logDict = readLogFile(logName)
    logDictLocal = readLogFile

    # Identify common simulations
    setSim = set(logDictLocal['simName'])
    indSims = [i for i, item in enumerate(logDict['simName']) if item in setSim]

    log.info('Common simulations are: %d' % indSims)

    return indSims


def getDFAData(avaDir, workDir, suffix, nameDir=''):
    """ Export the required data from com1DFA output to Aimec Work directory and rename  """

    # Lead all infos on simulations
    inputDir = os.path.join(avaDir, 'Outputs', 'com1DFA', 'peakFiles')
    data = makeSimDict(inputDir)

    countsuf = 0
    for m in range(len(data['files'])):
        if data['resType'][m] == suffix:
            if nameDir == '':
                shutil.copy(data['files'][m], workDir)
            else:
                shutil.copy(data['files'][m], '%s/%s/%06d.txt' % (workDir, nameDir, countsuf+1))
                print(data['files'][m], '%s/%s/%06d.txt' % (workDir, nameDir, countsuf+1))
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


def exportcom1DFAOutput(avaDir):
    """ Export the simulation results from com1DFA output to desired location

        Inputs:     avaDir:     name of avalanche
                    workDir:    directory where data shall be exported to

        Outputs:    simulation result files saved to Outputs/com1DFA
    """

    # Initialise directories
    inputDir = os.path.join(avaDir, 'Work', 'com1DFA')
    outDir = os.path.join(avaDir, 'Outputs', 'com1DFA')
    outDirPF = os.path.join(outDir, 'peakFiles')
    outDirRep = os.path.join(outDir, 'reports')
    makeADir(outDir)
    makeADir(outDirPF)
    makeADir(outDirRep)

    # Read log file information
    logName = os.path.join(inputDir, 'ExpLog.txt')
    logDict = readLogFile(logName)

    # Get number of values
    sNo = len(logDict['noSim'])

    # Path to com1DFA results
    resPath = os.path.join(inputDir, 'FullOutput_mu_')

    # Export peak files and reports
    for k in range(sNo):
        shutil.copy('%s%.03f/%s/raster/%s_pfd.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                    logDict['simName'][k]),
                    '%s/%s_%s_pfd.asc' % (outDirPF, logDict['simName'][k], logDict['Mu'][k]))
        shutil.copy('%s%.03f/%s/raster/%s_ppr.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                    logDict['simName'][k]),
                    '%s/%s_%s_ppr.asc' % (outDirPF, logDict['simName'][k], logDict['Mu'][k]))
        shutil.copy('%s%.03f/%s/raster/%s_pv.asc' % (resPath, logDict['Mu'][k], logDict['simName'][k],
                    logDict['simName'][k]),
                    '%s/%s_%s_pv.asc' % (outDirPF, logDict['simName'][k], logDict['Mu'][k]))
        shutil.copy('%s%.03f/%s.html' % (resPath, logDict['Mu'][k], logDict['simName'][k]),
                    '%s/%s_%s.html' % (outDirRep, logDict['simName'][k], logDict['Mu'][k]))

    # Export ExpLog to Outputs/com1DFA
    shutil.copy2('%s/ExpLog.txt' % inputDir, outDir)


def makeSimDict(inputDir):
    """ Create a dictionary that contains all info on simulations:

            files:          full file path
            names:          file name
            simType:        entres or null (e.g. entres is simulation with entrainment and resistance)
            resType:        which result parameter (e.g. 'ppr' is peak pressure)
            releaseArea:    release area
            Mu:             value of Mu parameter
            cellSize:       cell size of raster file
    """

    # Load input datasets from input directory
    datafiles = glob.glob(inputDir+os.sep + '*.asc')

    # Sort datafiles by name
    datafiles = sorted(datafiles)

    # Make dictionary of input data info
    data = {'files': [], 'names': [], 'resType': [], 'simType': [],
            'releaseArea': [], 'cellSize' : [], 'Mu' : []}

    for m in range(len(datafiles)):
        data['files'].append(datafiles[m])
        name = os.path.splitext(os.path.basename(datafiles[m]))[0]
        data['names'].append(name)
        nameParts = name.split('_')
        data['releaseArea'].append(nameParts[0])
        data['simType'].append(nameParts[1])
        data['Mu'].append(nameParts[3])
        data['resType'].append(nameParts[4])
        header = IOf.readASCheader(datafiles[m])
        data['cellSize'].append(header.cellsize)

    return data
