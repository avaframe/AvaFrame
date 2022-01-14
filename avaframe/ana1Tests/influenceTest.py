"""
Influence Test module

This module objetive is to test the influence of some parameters on the runout
area computed by the com1DFA module.
The chosen refence solution will be the one computed by com2AB module.
"""

# imports
import math
import logging
import configparser
config = configparser.RawConfigParser()
config.optionxform = str
import pathlib
import pandas as pd
import os

# local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in1Data import getInput as gI
from avaframe.out3Plot import outAB
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.com2AB.com2AB as com2AB
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.com1DFA.deriveParameterSet as dP
from avaframe.out3Plot import outInfluenceTest

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)



def mainInfluenceTest(avalancheDir, influenceTestCfg, simDF, resAIMEC, pathDict):
    """ Get simulation data frame, and AIMEC results to compare the runout, mmpr and mmpfd of different simulations
    while varying inifiles parameters
    Parameters
    -----------
    avalancheDir: str or pathlib path
        avalanche directory
    influenceTestCfg: path to ana1Tests/influenceTestCfg.ini
        contains settings for samos and ata simulations
    simDF: Data-Frame
        contains results of DFA simulations
    resAIMEC: dict
        dict with AIMEC results
    pathDict: path dictionary
        needed to acces AIMEC results
    resAB : dict
        dict with com2AB results
    """

    # Collect data from config file
    log.info('Collect data from config file')
    resInfluenceTest = collectDataIniFile(influenceTestCfg)

    # Get each DFA simulation runout
    log.info('Get each DFA simulation runout')
    simDF = getAIMECResults(simDF, resAIMEC, pathDict)

    visc2StudyCfg = pathlib.Path('ana1Tests', resInfluenceTest['visc2Study'] + 'Cfg.ini')
    refCfg = pathlib.Path('ana1Tests', 'refCfg.ini')

    # Get ATA simulation label and, if needed, get those from SAMOS simulations
    IDList = getSimLabel(avalancheDir, visc2StudyCfg, simDF, param2Study = resInfluenceTest['param2Study'])
    refID = getSimLabel(avalancheDir, refCfg, simDF)[0]

    # analyze and plot the results
    resInfluenceTest = outInfluenceTest.analyzeResults(simDF, avalancheDir, IDList, refID, resInfluenceTest)
    outInfluenceTest.plotResults(avalancheDir, resInfluenceTest)

    # Check if last line is computed
    log.info('Main Function well computed')



def collectDataIniFile(cfg):
    """ Get some informations related to parameters to study in .ini files
    Parameters
    -----------
    cfg: path to .ini file
    Returns
    -------
    param2Study: string
        main parameter under study
    visc2Study: string (can be 'ata' or 'samos')
        chose if we analyze ata or samos viscosity results
    otherParamValueList: List of float List
        contains values taken by each parameter who is specified to vary in .ini files
    otherParamNameList: string List
        contains the name of each parameter who is specified to vary in .ini files
    xAxisArray: float List
        xAxis array for ploting results
    """
    config.read(cfg)
    param2Study = config['INFLUENCETEST']['param2Study']
    visc2Study = config['INFLUENCETEST']['visc2Study']
    otherParamValueList = []
    otherParamNameList = []
    for paramName in config[visc2Study]:
        if paramName != param2Study:
            otherParamNameList.append(paramName)
            otherParamValueList.append([])
    nParam = len(otherParamNameList)
    for i in range(nParam):
        paramName = otherParamNameList[i]
        otherParamValueList[i] = fromIni2List(config[visc2Study][paramName])
    xAxisArray = fromIni2List(config[visc2Study][param2Study])
    for index in range(len(xAxisArray)):
        xAxisArray[index] = float(xAxisArray[index])
    xAxisArray.sort()
    print('xAxisArray : ', xAxisArray)

    resInfluenceTest = {}
    resInfluenceTest['param2Study'] = param2Study
    resInfluenceTest['visc2Study'] = visc2Study
    resInfluenceTest['otherParamValueList'] = otherParamValueList
    resInfluenceTest['otherParamNameList'] = otherParamNameList
    resInfluenceTest['xAxisArray'] = xAxisArray

    return resInfluenceTest



def fromIni2List(iniExpression):
    """ Take a .ini suit of values and turns it into a list (for example: 1|2|3 becomes ['1', '2', '3'])
    Parameters
    ----------
    iniExpression: string
    Returns
    -------
    List: string List
    """
    List = ['']
    index=0
    for char in iniExpression:
        if char!='|':
            List[index] = List[index] + char
        elif char=='|':
            List.append('')
            index = index + 1

    return List



def getABRunout(resAB):
    """ Get alpha-beta runout
    Parameters
    ----------
    resAB : dict
        dict with com2AB results
    Returns
    -------
    abRunout: float
    """
    AvaPath = resAB['AvaPath']
    NameAva = AvaPath['Name']
    abRunout=[]
    for i in range(len(NameAva)):
        name = NameAva[i]
        resAB = outAB.processABresults(resAB, name)
        ids_alpha = resAB[name]['ids_alpha']
        sAB = resAB[name]['s'][ids_alpha]
        abRunout.append(sAB)

    return abRunout



def getAIMECResults(simDF, resAIMEC, pathDict):
    """ Get DFA runout
    Parameters
    ----------
    simDF: Data-Frame
        contains results of DFA simulations
    resAIMEC: dict
        dict with AIMEC results
    pathDict: path dictionary
        needed to acces AIMEC results
    Returns
    -------
    simDF: DataFrame
    """

    index = 0
    for pathDictID in pathDict['simID']:
        for simDFID in simDF.index:
            if simDFID==pathDictID:
                simDF.loc[simDFID, 'runout'] = resAIMEC['runout'][0][index]
                simDF.loc[simDFID, 'MMPPR'] = resAIMEC['MMPPR'][index]
                simDF.loc[simDFID, 'MMPFD'] = resAIMEC['MMPFD'][index]
                index = index + 1

    return simDF



def getSimLabel(avalancheDir, cfg, simDF, param2Study=0):
    """ Returns a list with every simulation hashes generated by an .ini file
    Parameters
    ----------
    avalancheDir: str or pathlib path
        avalanche directory
    cfg: path to an .ini file
        contains settings for samos and ata simulations
    param2Study: string
        main parameter under study
    simDF: Data-Frame
        contains results of DFA simulations
    """
    # get information on simulations that shall be performed according to parameter variation
    modCfg, variationDict = dP.getParameterVariationInfo(avalancheDir, com1DFA, cfg, variationDict='')
    # check if parameter variation on release or entrainment thickness is working - where thickness is read from
    dP.checkRelEntThVariation(modCfg, variationDict)
    # fetch input data - dem, release-, entrainment- and resistance areas
    inputSimFiles = gI.getInputDataCom1DFA(avalancheDir, modCfg['FLAGS'])
    # create a list of simulations and generate an individual configuration object for each simulation
    # if need to reproduce exactly the hash - need to be strings with exactely the same number of digits!!
    # first get already existing simulations
    simDFOld, simNameOld = cfgUtils.readAllConfigurationInfo(avalancheDir, specDir='')
    # Get each simulation label
    IDList, _ = com1DFA.prepareVarSimDict(modCfg, inputSimFiles, variationDict, simNameOld=simNameOld)

    # if param2Study != 0:
    #     tri(IDList, simDF, param2Study)

    return IDList



def generateAutoIniFile(cfg, viscName):
    """ generates .ini files starting from ana1Tests/influenceTestCfg.ini and combining datas with the section
    related to the viscosity to study.
    Parameters
    ----------
    cfg: path to an .ini file
    viscName: string
        has to be 'ata' or 'samos'
    """
    log.info("Write " + viscName + " simulation settings in '" + viscName + "Cfg.ini' file")
    config.read(cfg)
    paramNameList = []
    paramValueList = []
    for param in config[viscName]:
        print('param = ', param)
        paramNameList.append(param)
        paramValueList.append(config[viscName][param])

    if viscName == 'samos':
        paramNameList.append('viscOption')
        paramValueList.append('1')
    elif viscName == 'ata':
        paramNameList.append('viscOption')
        paramValueList.append('2')

    simTypeList = config['GENERAL']['simTypeList']
    if simTypeList != 'null':
        log.info('Please modify source code yourself to distinguish simulations with entrainment or not')
        log.info('It will be an option developped late sorry')
        log.info('So for the moment, please put simTypeList to null')
        n=input()

    for generalParam in config['GENERAL']:
        if generalParam not in paramNameList:
            print('generalParameter = ', generalParam)
            print(config['GENERAL'][generalParam])
            paramNameList.append(generalParam)
            paramValueList.append(config['GENERAL'][generalParam])

    # Open config file
    outCfg = 'ana1Tests/' + viscName + 'Cfg.ini'
    if os.path.exists(outCfg):
        os.remove(outCfg)
    f = open(outCfg, 'w')
    f.write('[GENERAL]\n')

    index=0
    for paramName in paramNameList:
        f.write(paramName + '=' + paramValueList[index] + '\n')
        index = index + 1

    # Close config file
    f.close()
