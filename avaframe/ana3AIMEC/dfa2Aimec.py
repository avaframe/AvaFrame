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
import avaframe.out3Plot.plotUtils as pU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def extractCom1DFAMBInfo(avaDir, pathDict, simNameInput=''):
    """ Extract the mass balance info from the log file """

    if simNameInput != '':
        simNames = [simNameInput]
        nameDir = simNameInput
    else:
        # Get info from ExpLog
        nameDir = 'com1DFA'
        logLoc = os.path.join(avaDir, 'Outputs', 'com1DFA')
        logName = os.path.join(logLoc, 'ExpLog.txt')
        logDictExp = fU.readLogFile(logName)
        names = logDictExp['fullName']
        simNames = sorted(set(names), key=lambda s: (s.split("_")[0], s.split("_")[1], s.split("_")[3]))

    # Read mass data from log and save to file for each simulation run
    countFile = 0
    for simName in simNames:
        log.debug('This is the simulation name %s for mod com1DFA ' % (simName))

        # Initialise fields
        time = []
        mass = []
        entrMass = []
        indRun = [0]
        countMass = 0
        flagStop = 0

        # Read log file
        locFiles = os.path.join(avaDir, 'Outputs', 'com1DFA')
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
            saveName = os.path.join(locFiles, 'mass_%s.txt' % (simName))
            with open(saveName, 'w') as MBFile:
                MBFile.write('time, current, entrained\n')
                for m in range(indRun[k], indRun[k] + indRun[k+1] - indRun[k]-1):
                    MBFile.write('%.02f,    %.06f,    %.06f\n' %
                                 (logDict['time'][m], logDict['mass'][m], logDict['entrMass'][m]))
            if simNameInput != '':
                pathDict[simName]['mb'].append(saveName)
            else:
                pathDict['mb'].append(saveName)
            log.debug('Mass file saved to %s ' % (saveName))
            log.debug('Added to pathDict[mb] %s ' % (saveName))
            countFile = countFile + 1

    return pathDict


def getMBInfo(avaDir, pathDict, comMod, simName=''):
    """ Get MB info """

    # Get info from ExpLog
    if simName != '':
        mbFile = os.path.join(avaDir, 'Outputs', comMod, 'mass_%s.txt' % simName)
        pathDict[simName]['mb'].append(mbFile)
        log.info('Added to pathDict[mb] %s' % (mbFile))

    else:
        mbFiles = glob.glob(os.path.join(avaDir, 'Outputs', comMod, 'mass*.txt'))
        mbNames = sorted(set(mbFiles), key=lambda s: (s.split("_")[1], s.split("_")[2], s.split("_")[4]))

        for mFile in mbNames:
            pathDict['mb'].append(mFile)
            log.debug('Added to pathDict[mb] %s' % (mFile))

    return pathDict

def getRefMB(testName, pathDict, simName):
    """ Get MB info """

    # Get info from ExpLog
    mbFile = os.path.join('..', 'benchmarks', testName, 'mass_%s.txt' % simName)
    pathDict[simName]['mb'].append(mbFile)
    log.info('Added to pathDict[mb] %s' % (mbFile))

    return pathDict


def dfaComp2Aimec(avaDir, cfgSetup):
    """ Exports the required data from com1DFA to be used by Aimec """

    # look for matching simulations
    comModules = cfgSetup['comModules'].split('|')
    refModule = comModules[0]
    compModule = comModules[1]
    log.info('Reference data is from module: %s' % refModule)
    log.info('Comparison data is from module: %s' % compModule)

    # Lead all infos on refernce simulations
    if refModule == 'benchmarkReference':
        testName = cfgSetup['testName']
        inputDirRef = os.path.join('..', 'benchmarks', testName)
        refData = fU.makeSimDict(inputDirRef)
        # this is required to find matching sim Names
        for m in range(len(refData['simName'])):
            refData['simName'][m] = refData['simName'][m].replace('ref', 'dfa')
    else:
        inputDirRef = os.path.join(avaDir, 'Outputs', refModule, 'peakFiles')
        refData = fU.makeSimDict(inputDirRef)

    # Load all infos on comparison module simulations
    inputDirComp = os.path.join(avaDir, 'Outputs', compModule, 'peakFiles')
    compData = fU.makeSimDict(inputDirComp)

    pathDict = {}
    simNamesMatch = {}
    count = 0
    # initialise path dicionary with subdictionary for each simulation
    for simNameRef in refData['simName']:
        pathDict.update({simNameRef: {'ppr': [], 'pfd': [], 'pfv': [], 'mb': []}})
        simNamesMatch[simNameRef] = False

    for countRef, simNameRef in enumerate(refData['simName']):
        for countComp, simNameComp in enumerate(compData['simName']):
            if simNameRef == simNameComp:
                suffix = ['pfd', 'ppr', 'pfv']
                for suf in suffix:
                    if refData['resType'][countRef] == suf and compData['resType'][countComp] == suf:
                        pathDict[simNameRef][suf].append(refData['files'][countRef])
                        log.info('Added to pathDict[%s] %s ' % (suf, refData['files'][countRef]))
                        pathDict[simNameRef][suf].append(compData['files'][countComp])
                        log.info('Added to pathDict[%s] %s' % (suf, compData['files'][countComp]))
                        if simNamesMatch[simNameRef] == False:
                            for comMod in comModules:
                                if comMod == 'ref':
                                    pathDict = getRefMB(testName, pathDict, simNameRef)
                                elif comMod == 'com1DFA':
                                    pathDict = extractCom1DFAMBInfo(avaDir, pathDict, simNameInput=simNameRef)
                                else:
                                    pathDict = getMBInfo(avaDir, pathDict, comMod, simName=simNameRef)
                        simNamesMatch[simNameRef] = True
                pathDict[simNameRef]['compType'] = ['comModules', refModule, compModule]
                pathDict[simNameRef]['referenceFile'] = 0

    for key in simNamesMatch:
        if simNamesMatch[key] == False:
            log.info('no matching files found for simulation: %s' % key)
            del pathDict[key]

    # info about colourmap
    pathDict['contCmap'] = True

    return pathDict


def mainDfa2Aimec(avaDir, comModule='com1DFA'):
    """ Exports the required data from com1DFA to be used by Aimec """

    # path dictionary for Aimec
    pathDict = {'ppr': [], 'pfd': [], 'pfv': [], 'mb': []}

    # Setup input from com1DFA and save file paths to dictionary
    suffix = ['pfd', 'ppr', 'pfv']
    for suf in suffix:
        pathDict = fU.getDFADataPaths(avaDir, pathDict, suf, comModule)

    # Extract mb info
    if comModule == 'com1DFA':
        pathDict = extractCom1DFAMBInfo(avaDir, pathDict)
    elif comModule == 'com1DFAPy':
        pathDict = getMBInfo(avaDir, pathDict, comModule)

    pathDict['compType'] = ['singleModule', comModule]

    # info about colormap
    pathDict['contCmap'] = pU.contCmap

    return pathDict


def indiDfa2Aimec(avaDir, suffix, comModule='com1DFA', inputDir=''):
    """ Exports the required data from com1DFA to be used by Aimec """

    # path dictionary for Aimec
    pathDict = {suffix: []}

    # Setup input from com1DFA and save file paths to dictionary
    pathDict = fU.getDFADataPaths(avaDir, pathDict, suffix, comModule, inputDir=inputDir)

    pathDict['compType'] = ['singleModule', comModule]

    # info about colormap
    pathDict['contCmap'] = pU.contCmap

    return pathDict
