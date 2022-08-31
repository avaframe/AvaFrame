"""

This is a simple function for computing a probability map of all peak files of one parameter that
exceed a particular threshold

"""

import numpy as np
import logging
import pathlib

import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in1Data.computeFromDistribution as cP
import avaframe.com1DFA.deriveParameterSet as dP


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def createComModConfig(cfgProb, avaDir, modName, cfgFileMod=''):
    """ create configuration file for performing sims with modName com module

        Parameters
        -----------
        cfgProb: configParser object
            configuration settings
        avaDir: pathlib path
            path to avalanche directory
        modName: module
            computational module
        cfgFileMod: str
            path to cfgFile for computational module name - optional

        Returns
        -------
        cfgFiles: dict
            dictionary of paths to newly generated configuration files for com module for all parameters

    """

    # setup where configuration file is saved
    outDir = avaDir / 'Outputs'
    fU.makeADir(outDir)

    # loop over all parameters for performing parameter variation
    variationsDict = makeDictFromVars(cfgProb['PROBRUN'])
    cfgFiles = {}
    for varName in variationsDict:
        # define configuration files
        # get filename of module
        modNameString = str(pathlib.Path(modName.__file__).stem)
        cfgFile = outDir / ('probRun%sCfg%s.ini' % (modNameString, varName))

        # use cfgFile, local com module settings or default settings if local not available
        modCfg = cfgUtils.getModuleConfig(modName, fileOverride=cfgFileMod)
        modCfg = updateCfgRange(modCfg, cfgProb, varName, variationsDict[varName])
        with open(cfgFile, 'w') as configfile:
            modCfg.write(configfile)
        # append cfgFiles to list
        cfgFiles[varName] = {'cfgFile': cfgFile}

    return cfgFiles


def updateCfgRange(cfg, cfgProb, varName, varDict):
    """ update cfg with a range for parameters in cfgProb

        Parameters
        -----------
        cfg: configparser object
            configuration object to update
        cfgProb: configParser object
            configparser object with info on update
        varName: str
            name of parameter used for variation
        varDict: dict
            dictionary with variationValue and numberOfSteps for varName

        Returns
        --------
        cfg: configParser
            updated configuration object

    """

    # set reference values of parameters - override values in com module configurations
    varParList = cfgProb['PROBRUN']['varParList'].split('|')
    # also for the other parameters that are varied subsequently
    # first check if no parameter variation in provided for these parameters in the com module ini
    # if so - errror
    for varPar in varParList:
        if any(chars in cfg['GENERAL'][varPar] for chars in ['|', '$', ':']):
            message = ('Only one reference value is allowed for %s: but %s is given' %
                (varPar, cfg['GENERAL'][varPar]))
            log.error(message)
            raise AssertionError(message)
        elif varPar in ['entTh', 'relTh', 'secondaryRelTh']:
            thPercentVariation = varPar + 'PercentVariation'
            thRangeVariation = varPar + 'RangeVariation'
            thDistVariation = varPar + 'DistVariation'
            if cfg['GENERAL'][thPercentVariation] != '' or cfg['GENERAL'][thRangeVariation] != '' or cfg['GENERAL'][thDistVariation] != '':
                message = ('Only one reference value is allowed for %s: but %s %s or %s %s is given' %
                    (varPar, thPercentVariation, cfg['GENERAL'][thPercentVariation],
                     thRangeVariation, cfg['GENERAL'][thRangeVariation]))
                log.error(message)
                raise AssertionError(message)

    # this is now done for parameter VARNAME from inputs
    # get range, steps and reference value of parameter to perform variations
    valVariation = varDict['variationValue']
    valSteps = varDict['numberOfSteps']
    valVal = cfg['GENERAL'][varName]

    if cfgProb['PROBRUN']['variationType'].lower() == 'normaldistribution':
        # get computeFromDistribution configuration and apply override
        cfgDist = cfgUtils.getModuleConfig(cP, fileOverride='', modInfo=False, toPrint=False,
                                              onlyDefault=cfgProb['computeFromDistribution_override']['defaultConfig'])
        cfgDist, cfgProb = cfgHandling.applyCfgOverride(cfgDist, cfgProb, cP, addModValues=False)

    # set variation in configuration
    if varName in ['relTh', 'entTh', 'secondaryRelTh']:
        # if variation using normal distribution
        if cfgProb['PROBRUN']['variationType'].lower() == 'normaldistribution':
            parName = varName + 'DistVariation'
            if valVariation == '':
                valVariation = '-'
            parValue = (cfgProb['PROBRUN']['variationType'] + '$'
                + valSteps + '$'  + valVariation + '$'
                + cfgDist['GENERAL']['minMaxInterval'] + '$'
                + cfgDist['GENERAL']['buildType'] + '$'
                + cfgDist['GENERAL']['support'])
        # if variation using percent
        elif cfgProb['PROBRUN']['variationType'].lower() == 'percent':
            parName = varName + 'PercentVariation'
            parValue = valVariation + '$' + valSteps
        # if variation using absolute range
        elif cfgProb['PROBRUN']['variationType'].lower() == 'range':
            parName = varName + 'RangeVariation'
            parValue = valVariation + '$' + valSteps
        else:
            message = ('Variation Type: %s - not a valid option, options are: percent, range, normaldistribution' % cfgProb['PROBRUN']['variationType'])
            log.error(message)
            raise AssertionError(message)
        # write parameter variation for varName in config file
        cfg['GENERAL'][parName] = parValue
        cfg['GENERAL']['addStandardConfig'] = cfgProb['PROBRUN']['addStandardConfig']
    else:
        # set variation
        if  cfgProb['PROBRUN']['variationType'].lower() == 'normaldistribution':
            cfgDist = {'sampleSize': valSteps, 'mean': valVal,
                'buildType': cfgProb['computeFromDistribution_override']['buildType'],
                'buildValue': valVariation,
                'minMaxInterval':  cfgDist['GENERAL']['minMaxInterval'],
                'support': cfgDist['GENERAL']['support']}
            _, valValues, _, _ = cP.extractNormalDist(cfgDist)
            cfg['GENERAL'][varName] = dP.writeToCfgLine(valValues)
        elif cfgProb['PROBRUN']['variationType'].lower() == 'percent':
            cfg['GENERAL'][varName] = '%s$%s$%s' % (valVal, valVariation, valSteps)
            valValues = fU.splitIniValueToArraySteps(cfg['GENERAL'][varName])

        elif cfgProb['PROBRUN']['variationType'].lower() == 'range':
            if '-' in valVariation or '+' in valVariation:
                valStart = str(float(valVal) + float(valVariation))
                valStop = float(valVal)
            else:
                valStart = str(float(valVal) - float(valVariation))
                valStop = str(float(valVal) + float(valVariation))
            cfg['GENERAL'][varName] = '%s:%s:%s' % (valStart, valStop, valSteps)
            valValues = np.linspace(float(valStart), float(valStop), int(valSteps))
        else:
            message = ('Variation Type: %s - not a valid option, options are: percent, range, normaldistribution' % cfgProb['PROBRUN']['variationType'])
            log.error(message)
            raise AssertionError(message)

        # if reference value is not in variation - add reference values
        if cfgProb['PROBRUN'].getboolean('addStandardConfig'):
            if np.isclose(valValues, float(valVal), atol=1.e-7, rtol=1.e-8).any():
                log.info('Reference value is in variation for %s' % (varName))
            else:
                log.info('Reference value of %s: %s is added' % (varName, str(valVal)))
                cfg['GENERAL'][varName] = cfg['GENERAL'][varName] + '&' + str(valVal)

    # add a scenario Name to VISUALISATION
    cfg['VISUALISATION']['scenario'] = varName

    return cfg


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
    simNameList = cfgHandling.filterSims(avaDir, parametersDict, specDir=inputDir)

    # initialize flag if analysis has been performed or e.g. no matching files found
    analysisPerformed = False
    if simNameList == []:
        # no matching sims found for filtering criteria
        log.warning('No matching simulations found for filtering criteria')
        return analysisPerformed

    # if matching sims found - perform analysis
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
    log.info('probability analysis performed for peak parameter: %s and a peak value '
             'threshold of: %s %s' % (cfg['GENERAL']['peakVar'], cfg['GENERAL']['peakLim'], unit))
    log.info('%s peak fields added to analysis' % count)

    # # Save to .asc file
    avaName = avaDir.name
    outFileName = '%s_probMap%s.asc' % (avaName, cfg['GENERAL']['peakLim'])
    outFile = outDir / outFileName
    IOf.writeResultToAsc(header, probMap, outFile)
    analysisPerformed = True

    return analysisPerformed


def makeDictFromVars(cfg):
    """ create a dictionary with info on parameter variation for all parameter in
        varParList

        Parameters
        -----------
        cfg: configparser object
            configuration settings, here varParList, variationValue, numberOfSteps

        Returns
        --------
        variationsDict: dict
            dictionary with for each varName, varVariation, varSteps

    """

    varParList = cfg['varParList'].split('|')
    varValues = cfg['variationValue'].split('|')
    varSteps = cfg['numberOfSteps'].split('|')

    # check if value is provided for each parameter
    if (len(varParList) == len(varValues) == len(varSteps)) is False:
        message = 'For every parameter in varParList a variationValue and numberOfSteps needs to be provided'
        log.error(message)
        raise AssertionError

    variationsDict = {}
    for idx, val in enumerate(varParList):
        variationsDict[val] = {'variationValue': varValues[idx], 'numberOfSteps': varSteps[idx]}

    return variationsDict
