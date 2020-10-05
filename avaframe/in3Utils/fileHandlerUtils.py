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


def readLogFile(logName, cfg=''):
    """ Read experiment log file and make dictionary that contains all simulations """

    # Read log file
    logFile = open(logName, 'r')
    log.info('Take com1DFA full experiment log')

    # Parameter variation
    if cfg != '':
        varPar = cfg['varPar']
    else:
        varPar = 'Mu'

    # Save info to dictionary, add all result parameters that are saved in com1DFA Outputs
    suffix = ['pfd', 'ppr', 'pv', 'fd']
    logDict = {'noSim' : [], 'simName' : [], varPar: [], 'fullName' : [], 'suffix ': suffix}

    lines = logFile.readlines()[1:]
    countSims = 1
    for line in lines:
        vals = line.strip().split()
        logDict['noSim'].append(countSims)
        logDict['simName'].append(vals[1])
        logDict[varPar].append(float(vals[2]))
        logDict['fullName'].append(vals[1] + '_' + '%.3f' % float(vals[2]))
        countSims = countSims + 1

    return logDict


def extractParameterInfo(avaDir, simName):
    """ Extract info about simulation parameters from the log file """

    # Get info from ExpLog
    logName = os.path.join(avaDir, 'Outputs', 'com1DFA', 'ExpLog.txt')
    logDictExp = readLogFile(logName)
    names = logDictExp['fullName']
    simNames = sorted(set(names), key=lambda s: s.split("_")[3])

    # Initialise parameter dictionary
    parameterDict = {}

    # Initialise fields
    time = []
    mass = []
    entrMass = []
    flagNoStop = True
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

    # Set values in dictionary
    parameterDict['release mass [kg]'] = np.asarray(mass)[0]
    parameterDict['final time step [s]'] = np.asarray(time)[-1]
    parameterDict['current mass [kg]'] = np.asarray(mass)[-1]

    return parameterDict


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
                log.info('%s   : %s/%s/%06d.txt' % (data['files'][m], workDir, nameDir, countsuf+1))
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

    # Create simulations dictionary
    data = makeSimDict(refDir)

    # copy these files to desired working directory for outQuickPlot
    for m in range(len(data['files'])):
        if data['resType'][m] == suffix:
            if nameDir != '':
                shutil.copy(data['files'][m], '%s/%s/000000.txt' % (outputDir, nameDir))
            else:
                shutil.copy2(data['files'][m], outputDir)

    # Give status information
    if os.path.isdir(refDir) == False:
        log.error('%s does not exist - no files for reference found' % refDir)
    elif data['files'] == []:
        log.error('No files found in %s' % refDir)
    else:
        log.info('Reference files copied from directory: %s' % refDir)


def exportcom1DFAOutput(avaDir, cfg=''):
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
    if cfg != '':
        logDict = readLogFile(logName, cfg)
        varPar = cfg['varPar']
    else:
        logDict = readLogFile(logName)
        varPar = 'Mu'

    # Get number of values
    sNo = len(logDict['noSim'])

    # Path to com1DFA results
    resPath = os.path.join(inputDir, 'FullOutput_%s_' % varPar)

    # Export peak files and reports
    for k in range(sNo):
        pathFrom = os.path.join('%s%.03f' % (resPath, logDict[varPar][k]),
                                logDict['simName'][k], 'raster',
                                '%s_pfd.asc' % logDict['simName'][k])
        pathTo = os.path.join(outDirPF, '%s_%.03f_pfd.asc' % (logDict['simName'][k], logDict[varPar][k]))
        shutil.copy(pathFrom, pathTo)
        pathFrom = os.path.join('%s%.03f' % (resPath, logDict[varPar][k]),
                                logDict['simName'][k], 'raster',
                                '%s_ppr.asc' % logDict['simName'][k])
        pathTo = os.path.join(outDirPF, '%s_%.03f_ppr.asc' % (logDict['simName'][k], logDict[varPar][k]))
        shutil.copy(pathFrom, pathTo)
        pathFrom = os.path.join('%s%.03f' % (resPath, logDict[varPar][k]),
                                logDict['simName'][k], 'raster',
                                '%s_pv.asc' % logDict['simName'][k])
        pathTo = os.path.join(outDirPF, '%s_%.03f_pv.asc' % (logDict['simName'][k], logDict[varPar][k]))
        shutil.copy(pathFrom, pathTo)
        pathFrom = os.path.join('%s%.03f' % (resPath, logDict[varPar][k]),
                                '%s.html' % logDict['simName'][k])
        pathTo = os.path.join(outDirRep, '%s_%.03f.html' % (logDict['simName'][k], logDict[varPar][k]))
        shutil.copy(pathFrom, pathTo)

    # Export ExpLog to Outputs/com1DFA
    shutil.copy2(os.path.join('%s' % inputDir, 'ExpLog.txt'), outDir)


def makeSimDict(inputDir, varPar=''):
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

    # Check if parameter variation other than Mu
    if varPar == '':
        varPar = 'Mu'

    # Make dictionary of input data info
    data = {'files': [], 'names': [], 'resType': [], 'simType': [], 'simName' : [],
            'modelType' : [], 'releaseArea': [], 'cellSize': [], varPar: []}

    for m in range(len(datafiles)):
        data['files'].append(datafiles[m])
        name = os.path.splitext(os.path.basename(datafiles[m]))[0]
        data['names'].append(name)
        nameParts = name.split('_')
        data['releaseArea'].append(nameParts[0])
        data['simType'].append(nameParts[1])
        data['modelType'].append(nameParts[2])
        data[varPar].append(nameParts[3])
        data['resType'].append(nameParts[4])
        data['simName'].append('_'.join(nameParts[0:4]))
        header = IOf.readASCheader(datafiles[m])
        data['cellSize'].append(header.cellsize)

    return data
