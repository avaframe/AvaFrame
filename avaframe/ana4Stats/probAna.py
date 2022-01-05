"""

This is a simple function for computing a probability map of all peak files of one parameter that
exceed a particular threshold

"""

import numpy as np
import logging
import pathlib

import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def createComModConfig(cfgProb, avaDir, modName):
    """ create configuration file for performing sims with modName com module

        Parameters
        -----------
        cfgProb: configParser object
            configuration settings
        avaDir: pathlib path
            path to avalanche directory
        modName: module
            computational module

        Returns
        -------
        cfgFiles: dict
            dictionary of paths to newly generated configuration files for com module for all parameters

    """

    # setup where configuration file is saved
    outDir = avaDir / 'Outputs'
    fU.makeADir(outDir)

    # loop over all parameters for performing parameter variation
    varParList = cfgProb['PROBRUN']['varParList'].split('|')
    cfgFiles = {}
    for varName in varParList:
        # define configuration files
        # get filename of module
        modNameString = str(pathlib.Path(modName.__file__).stem)
        cfgFile = outDir / ('probRun%sCfg%s.ini' % (modNameString, varName))

        # use default com module settings or local settings
        if cfgProb['PROBRUN'].getboolean('defaultSetup'):
            modCfg = cfgUtils.getDefaultModuleConfig(modName)
            modCfg, refIn = updateCfgRange(modCfg, cfgProb, varName)
            with open(cfgFile, 'w') as configfile:
                modCfg.write(configfile)
        else:
            modCfg = cfgUtils.getModuleConfig(modName)
            modCfg, refIn = updateCfgRange(modCfg, cfgProb, varName)
            with open(cfgFile, 'w') as configfile:
                modCfg.write(configfile)
        # append cfgFiles to list
        cfgFiles[varName] = {'cfgFile': cfgFile, 'referenceIncluded': refIn}

    return cfgFiles


def updateCfgRange(cfg1, cfgProb, varName):
    """ update cfg with a range for parameters in cfgProb

        Parameters
        -----------
        cfg1: configparser object
            configuration object to update
        cfgProb: configParser object
            configparser object with info on update
        varName: str
            name of parameter used for variation

        Returns
        --------
        cfg1: configParser
            updated configuration object
        refIn: bool
            True if the reference value had to be added to parameter variation

    """

    # set reference values of parameters - override values in com module configurations
    varParList = cfgProb['PROBRUN']['varParList'].split('|')
    # also for the other parameters that are varied subsequently
    for varPar in varParList:
        if any(chars in cfg1['GENERAL'][varPar] for chars in ['|', '$', ':']):
            message = 'Only one reference values is allowed for %s: but %s is given' % (varPar, cfg1['GENERAL'][varPar])
            log.error(message)
            raise AssertionError(message)

    # get range, steps and reference value of parameter to perform variations
    valVariation = cfgProb['PROBRUN']['%sVariation' % varName]
    valSteps = cfgProb['PROBRUN']['%sSteps' % varName]
    valVal = cfg1['GENERAL'][varName]

    refIn = False

    # set variation in configuration
    if varName in ['relTh', 'entTh', 'secondaryRelTh']:
        if cfgProb['PROBRUN'].getboolean('percentVariation'):
            parName = varName + 'PercentVariation'
            cfg1['GENERAL'][parName] = valVariation + '$' + valSteps
        else:
            parName = varPar + 'rangeVariation'
            cfg1['GENERAL'][parName] = valVariation + '$' + valSteps
    else:
        # set variation
        if cfgProb['PROBRUN'].getboolean('percentVariation'):
            cfg1['GENERAL'][varName] = '%s$%s$%s' % (valVal, valVariation, valSteps)
        else:
            valStart = str(float(valVal) - float(valVariation))
            valStop = str(float(valVal) + float(valVariation))
            cfg1['GENERAL'][varName] = '%s:%s:%s' % (valStart, valStop, valSteps)

    return cfg1, refIn


def probAnalysis(avaDir, cfg, module, parametersDict='', inputDir=''):
    """ Compute propability map of a given set of simulation result exceeding a particular threshold and save to outDir

        Parameters
        ----------
        avaDir: str
            path to avalanche directory
        cfg : dict
            configuration read from ini file of probAna function
        module
            computational module that was used to run the simulations
        parametersDict: dict
            dictionary with simulation parameters to filter simulations
        inputDir : str
            optional - path to directory where data that should be analysed can be found in
            a subfolder called peakFiles and configurationFiles, required if not in module results

    """

    # get filename of module
    modName = pathlib.Path(module.__file__).stem
    avaDir = pathlib.Path(avaDir)

    # set output directory
    outDir = avaDir / 'Outputs' / 'ana4Stats'
    fU.makeADir(outDir)

    # fetch all result files and filter simulations according to parametersDict
    simNameList = cfgUtils.filterSims(avaDir, parametersDict, specDir=inputDir)

    # initialize flag if analysis has been performed or e.g. no matching files found
    analysisPerformed = False
    if simNameList != []:

        if inputDir == '':
            inputDir = avaDir / 'Outputs' / modName / 'peakFiles'
            flagStandard = True
            peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)
        else:
            inputDirPF = inputDir / 'peakFiles'
            peakFilesDF = fU.makeSimDF(inputDirPF, avaDir=avaDir)

        # get header info from peak files - this should be the same for all peakFiles
        header = IOf.readASCheader(peakFilesDF['files'][0])
        cellSize = header['cellsize']
        nRows = header['nrows']
        nCols = header['ncols']
        xllcenter = header['xllcenter']
        yllcenter = header['yllcenter']
        noDataValue = header['noDataValue']

        # Initialise array for computations
        probSum = np.zeros((nRows, nCols))
        count = 0

        # Loop through peakFiles and compute probability
        for m in range(len(peakFilesDF['names'])):

            # only take simulations that match filter criteria from parametersDict
            if peakFilesDF['simName'][m] in simNameList:
                # Load peak field for desired peak field parameter
                if peakFilesDF['resType'][m] == cfg['GENERAL']['peakVar']:

                    # Load data
                    fileName = peakFilesDF['files'][m]
                    data = np.loadtxt(fileName, skiprows=6)
                    dataLim = np.zeros((nRows, nCols))

                    log.info('File Name: %s , simulation parameter %s ' % (fileName, cfg['GENERAL']['peakVar']))

                    # Check if peak values exceed desired threshold
                    dataLim[data > float(cfg['GENERAL']['peakLim'])] = 1.0
                    probSum = probSum + dataLim
                    count = count + 1

        # Create probability map ranging from 0-1
        probMap = probSum / count
        unit = pU.cfgPlotUtils['unit%s' % cfg['GENERAL']['peakVar']]
        log.info('probability analysis performed for peak parameter: %s and a peak value threshold of: %s %s' % (cfg['GENERAL']['peakVar'], cfg['GENERAL']['peakLim'], unit))
        log.info('%s peak fields added to analysis' % count)

        # # Save to .asc file
        avaName = avaDir.name
        outFileName = '%s_probMap%s.asc' % (avaName, cfg['GENERAL']['peakLim'])
        outFile = outDir / outFileName
        IOf.writeResultToAsc(header, probMap, outFile)
        analysisPerformed = True

    return analysisPerformed
