"""

This is a simple function for computing a probability map of all peak files of one parameter that
exceed a particular threshold

"""

import numpy as np
import logging
import pathlib
from scipy.stats import qmc

import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
import avaframe.in1Data.computeFromDistribution as cP
import avaframe.com1DFA.deriveParameterSet as dP
from avaframe.in3Utils import geoTrans as gT
from avaframe.out3Plot import statsPlots as sP
from avaframe.in1Data import getInput as gI


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
        cfgFiles: list
            list of paths to newly generated configuration files for com module inlcuding
            parameter variations
    """

    # setup where configuration file is saved
    modNameString = str(pathlib.Path(modName.__file__).stem)
    outDir = avaDir / 'Work' / ('%sConfigFiles' % modNameString)
    fU.makeADir(outDir)

    if cfgProb['PROBRUN'].getint('samplingStrategy') == 2:
        log.info('Probability run performed by varying one parameter at a time - local approach.')
        cfgFiles = cfgFilesLocalApproach(cfgProb, modName, outDir, cfgFileMod)

    else:
        log.info('Probability run perfromed drawing parameter set from full sample.')
        cfgFiles = cfgFilesGlobalApproach(cfgProb, modName, outDir, cfgFileMod)

    return cfgFiles, outDir


def cfgFilesGlobalApproach(cfgProb, modName, outDir, cfgFileMod):
    """ create configuration files with all parameter variations - drawn from full sample
        for performing sims with modName comModule

        Parameters
        -----------
        cfgProb: configParser object
            configuration settings
        avaDir: pathlib path
            path to avalanche directory
        modName: module
            computational module
        cfgFileMod: str
            path to cfgFile for computational module - optional

        Returns
        -------
        cfgFiles: list
            list of paths to newly generated configuration files for com module inlcuding parameter
            variations
    """

    # create sample of all parameter variations
    paramValuesD = createSampleWithPercentVariation(cfgProb, modName, fileOverride=cfgFileMod)

    # create plot of parameter sample if variation of two parameters
    sP.plotSample(paramValuesD, outDir)

    # write cfg files one for each parameter set drawn from full sample
    cfgFiles = createCfgFiles(paramValuesD, modName, cfgProb, cfgPath=outDir,
                              fileOverride=cfgFileMod)

    return cfgFiles


def cfgFilesLocalApproach(cfgProb, modName, outDir, cfgFileMod):
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

    # loop over all parameters for performing parameter variation
    variationsDict = makeDictFromVars(cfgProb['PROBRUN'])
    cfgFiles = []
    for varName in variationsDict:
        # define configuration files
        # get filename of module
        modNameString = str(pathlib.Path(modName.__file__).stem)
        cfgFile = outDir / ('probRun%sCfg%s.ini' % (modNameString, varName))

        # use cfgFile, local com module settings or default settings if local not available
        modCfg = fetchStartCfg(modName, cfgProb['PROBRUN'].getboolean('defaultComModuleCfg'), cfgFileMod)
        modCfg = updateCfgRange(modCfg, cfgProb, varName, variationsDict[varName])

        with open(cfgFile, 'w') as configfile:
            modCfg.write(configfile)
        # append cfgFiles to list
        cfgFiles.append(cfgFile)

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
    _, _ = checkParameterSettings(cfg, varParList)

    # this is now done for parameter VARNAME from inputs
    # get range, steps and reference value of parameter to perform variations
    valVariation = varDict['variationValue']
    valSteps = varDict['numberOfSteps']
    valVal = cfg['GENERAL'][varName]

    if cfgProb['PROBRUN']['variationType'].lower() == 'normaldistribution':
        # get computeFromDistribution configuration and apply override
        cfgDist = cfgUtils.getModuleConfig(cP, fileOverride='', modInfo=False, toPrint=False,
                                              onlyDefault=cfgProb['computeFromDistribution_override'].getboolean('defaultConfig'))
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

def checkParameterSettings(cfg, varParList):
    """ check if parameter settings in comMod configuration do not inlcude variation for parameters to be varied

        Parameters
        -----------
        cfg: configparser object
            configuration settings
        varParList: list
            list of parameters (names) that shall be varied

    """

    # set a list of all thickness parameters that are set to be read from shp file
    thReadFromShp = []

    # loop over all parameters and check if no variation is set and if read from shp
    for varPar in varParList:
        if any(chars in cfg['GENERAL'][varPar] for chars in ['|', '$', ':']):
            message = ('Only one reference value is allowed for %s: but %s is given' %
                (varPar, cfg['GENERAL'][varPar]))
            log.error(message)
            raise AssertionError(message)
        elif varPar in ['entTh', 'relTh', 'secondaryRelTh']:
            thFromShp = varPar + 'FromShp'
            # check if reference settings have already variation of varPar
            _ = checkForNumberOfReferenceValues(cfg, varPar)
            # check if th read from shp file
            if cfg['GENERAL'].getboolean(thFromShp):
                thReadFromShp.append(varPar)

    return True, thReadFromShp


def checkForNumberOfReferenceValues(cfg, varPar):
    """ check if in reference configuration no variation option of varPar is set
        if set - throw error

        Parameters
        -----------
        cfg: configparser object
            reference configuration settings
        varPar: str
            name of parameter to be checked

    """

    thPercentVariation = varPar + 'PercentVariation'
    thRangeVariation = varPar + 'RangeVariation'
    thDistVariation = varPar + 'DistVariation'

    # check if variation is set
    if cfg['GENERAL'][thPercentVariation] != '' or cfg['GENERAL'][thRangeVariation] != '' or cfg['GENERAL'][thDistVariation] != '':
        message = ('Only one reference value is allowed for %s: but %s %s or %s %s is given' %
            (varPar, thPercentVariation, cfg['GENERAL'][thPercentVariation],
             thRangeVariation, cfg['GENERAL'][thRangeVariation]))
        log.error(message)
        raise AssertionError(message)

    return True


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
    contourDict = {}

    # Loop through peakFiles and compute probability
    for m in range(len(peakFilesDF['names'])):

        # only take simulations that match filter criteria from parametersDict
        if peakFilesDF['simName'][m] in simNameList:
            # Load peak field for desired peak field parameter
            if peakFilesDF['resType'][m] == cfg['GENERAL']['peakVar']:

                # Load data
                fileName = peakFilesDF['files'][m]
                dataLim = np.zeros((nRows, nCols))
                fileData = IOf.readRaster(fileName)
                data = np.flipud(fileData['rasterData'])

                # fetch contourline info
                xGrid, yGrid, _, _ = gT.makeCoordGridFromHeader(fileData['header'])
                x, y = pU.fetchContourCoords(xGrid, yGrid, data, float(cfg['GENERAL']['peakLim']))
                contourDict[fileName.stem] = {'x': x, 'y': y}

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

    return analysisPerformed, contourDict


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


def fetchThicknessInfo(avaDir):
    """ Fetch input data for avaDir and thickness info

        Parameters
        ------------
        cfg: configparser object
            configuration settings
        avaDir: pathlib path or str
            path to avalanche directory

        Returns
        -----------
    """

    # fetch input data - dem, release-, entrainment- and resistance areas (and secondary release areas)
    inputSimFilesAll = gI.getInputDataCom1DFA(avaDir)

    # get thickness of release and entrainment areas (and secondary release areas) -if thFromShp = True
    inputSimFilesAll = gI.getThicknessInputSimFiles(inputSimFilesAll, avaDir)

    return inputSimFilesAll


def createSampleWithPercentVariation(cfgProb, comMod, fileOverride=''):
    """ Create a sample of parameters for a desired parameter variation,
        and draw nSample sets of parameter values from full sample
        if thickness values read from shp for comMod, convert sample values for these
        to a percent variation

        Parameters
        ------------
        cfgProb: configparser object
            configuration settings for parameter variation
        comMod: computational module
            module to perform then sims for parameter variation
        fileOverride: pathlib path
            optional - path to configuration file for comMod settings
            (for all parameters except those to be varied)

        Returns
        --------
        paramValuesD: dict
         dictionary used to pass parameter variation values
            names: list
                list of parameter names (that are varied)
            values: numpy nd array
                as many rows as sets of parameter values and as many rows as parameters
            typeList: list
                list of types of parameters (float, ...)
            thReadFromShp: list
                list of parameter names where the base value is read from shape
        """

    # random generator initialized with seed
    randomGen = np.random.default_rng(cfgProb['PROBRUN'].getint('sampleSeed'))

    # get filename of module
    modName = str(pathlib.Path(comMod.__file__).stem)

    # read initial configuration
    cfgStart = fetchStartCfg(comMod, cfgProb['PROBRUN'].getboolean('defaultComModuleCfg'), fileOverride)

    # fetch parameter names for parameter variation and variation value
    varParList = cfgProb['PROBRUN']['varParList'].split('|')
    valVariationValue = cfgProb['PROBRUN']['variationValue'].split('|')
    # check if thickness parameters are actually read from shp file
    _, thReadFromShp = checkParameterSettings(cfgStart, varParList)

    if len(thReadFromShp) > 0:
        if cfgProb['PROBRUN']['variationType'].lower() != 'percent':
            message = 'If thickness values read from shp file - only percent variation is allowed!'
            log.error(message)
            raise AssertionError(message)

    # initialze lower and upper bounds required to get a sample for the parameter values
    lowerBounds = []
    upperBounds = []
    for idx, varPar in enumerate(varParList):
        # if thickness parameters are read from shapefile - convert to a percent variation value
        # can be used in ini for thPercentVariation
        if varPar in thReadFromShp:
            lowerBounds.append((100. - float(valVariationValue[idx])))
            upperBounds.append((100. + float(valVariationValue[idx])))
        else:
            # if parameter value directly set in configuration modify the value directly
            varVal = cfgStart['GENERAL'].getfloat(varPar)
            if cfgProb['PROBRUN']['variationType'].lower() == 'percent':
                lB = varVal - varVal * (float(valVariationValue[idx]) / 100.)
                uB = varVal + varVal * ( float(valVariationValue[idx]) / 100.)
            elif cfgProb['PROBRUN']['variationType'].lower() == 'range':
                lB = varVal - float(valVariationValue[idx])
                uB = varVal + float(valVariationValue[idx])
            # update bounds
            lowerBounds.append(lB)
            upperBounds.append(uB)

    # create a sample of parameter values using scipy latin hypercube sampling
    if cfgProb['PROBRUN']['sampleMethod'].lower() == 'latin':
        sampler = qmc.LatinHypercube(d=len(varParList), seed=randomGen)
        sample = sampler.random(n=int(cfgProb['PROBRUN']['nSample']))
        sample = qmc.scale(sample, lowerBounds, upperBounds)
        log.info('Parameter sample created using latin hypercube sampling')
    else:
        message = ('Sampling method: %s not a valid option' % cfgProb['PROBRUN']['sampleMethod'])
        log.error(message)
        raise AssertionError(message)

    # create dictionary with all the info
    paramValuesD = {'names': varParList,
                    'values': sample,
                    'typeList': cfgProb['PROBRUN']['varParType'].split('|'),
                    'thReadFromShp': thReadFromShp}

    return paramValuesD


def createCfgFiles(paramValuesD, comMod, cfg, cfgPath='', fileOverride=''):
    """ create all config files required to run com Module from parameter variations using paramValues

        Parameters
        -----------
        paramValuesD: dict
            dictionary with parameter names and values (array of all sets of parameter values, one row per value set)
        comMod: com module
            computational module
        cfg: configparser object
            configuration settings
        cfgPath: str
            path where cfg files should be saved to
        fileOverride: pathlib path
            if config of comMod should be overwritten provide file path

        Returns
        --------
        cfgFiles: list
            list of cfg file paths for comMod including the updated values of the parameters to vary

    """

    # get filename of module
    modName = str(pathlib.Path(comMod.__file__).stem)

    # read initial configuration
    cfgStart = fetchStartCfg(comMod, cfg['PROBRUN'].getboolean('defaultComModuleCfg'), fileOverride)

    # create one cfgFile with one line of the parameter values from the full parameter variation
    cfgFiles = []
    for count1, pVal in enumerate(paramValuesD['values']):
        for index, par in enumerate(paramValuesD['names']):
            # convert percent variation to a +- variation in percent and of step 1
            if par in paramValuesD['thReadFromShp']:
                thVal = pVal[index] - 100.
                signVal = ['+' if thVal > 0 else '']
                cfgStart['GENERAL'][par + 'PercentVariation'] = signVal[0] + str(thVal) + '$1'
            else:
                cfgStart['GENERAL'][par] = str(pVal[index])
        cfgStart['VISUALISATION']['scenario'] = str(count1)
        cfgF = pathlib.Path(cfgPath, ('%d_%sCfg.ini' % (count1, modName)))
        with open(cfgF, 'w') as configfile:
            cfgStart.write(configfile)
        # append file path to list of cfg files
        cfgFiles.append(cfgF)

    return cfgFiles


def fetchStartCfg(comMod, onlyDefault, fileOverride):
    """ fetch start configuration of comMod
        if onlyDefault True and fileOverride path provided throw error

        Parameters
        -----------
        comMod: computational module
            module where configuration is read from
        onlyDefault: bool
            if True - read default config and not local
        fileOverride: pathlib path
            path to optional override configuration file

        Returns
        --------
        cfgStart: configparser object
            configuration object of comMod
    """

    if fileOverride != '' and onlyDefault:
        message = 'fileOverride provided AND defaultComModuleCfg set to True, only one is allowed'
        log.error(message)
        raise AssertionError(message)
    else:
        cfgStart = cfgUtils.getModuleConfig(comMod, fileOverride=fileOverride, toPrint=False,
                                            onlyDefault=onlyDefault)

    return cfgStart
