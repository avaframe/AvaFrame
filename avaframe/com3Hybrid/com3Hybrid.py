"""
Hybrid model module

This module contains the core function of the hybrid model.
This hybrid model combines the DFA simulation and the alpha beta model
"""
import pathlib
import logging
from configupdater import ConfigUpdater
import numpy as np
import copy
import matplotlib.pyplot as plt

# Local imports
# import config and init tools
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in3Utils import cfgUtils
from avaframe.in1Data import getInput

# import computation modules
from avaframe.com1DFA import com1DFA, particleTools, com1DFATools
from avaframe.com2AB import com2AB

# import analysis tools
from avaframe.out3Plot import outAB
from avaframe.runScripts import runAna3AIMEC

# import plotting tools
from avaframe.out3Plot import outCom3Plots

# create local logger
log = logging.getLogger(__name__)


def maincom3Hybrid(cfgMain, cfgHybrid):
    """This is the core function of the com3Hybrid module

    Here the DFA and AB models are run one after the other to first produce (from com1DFA)
    the avalanche path, then get the friction angle corresponding to the topography (from com2AB)
    and finally run the DFA simultaion one last time with the correct friction angle and get a 3D
    output of the avalanche
    """
    avalancheDir = cfgMain['MAIN']['avalancheDir']
    demOri = getInput.readDEM(avalancheDir)
    # get comDFA configuration path for hybrid model
    hybridModelDFACfg = pathlib.Path('com3Hybrid', 'hybridModel_com1DFACfg.ini')

    # prepare plots
    figDict = outCom3Plots.initializeFigures()

    # get initial mu value
    muArray = np.array([cfgHybrid.getfloat('DFA', 'mu')])
    alphaArray = np.array([np.degrees(np.arctan(muArray[0]))])

    # prepare for iterating
    nIterMax = cfgHybrid.getint('GENERAL', 'nIterMax')
    iteration = 0
    iterate = True
    while iteration < nIterMax and iterate:
        # update the com1DFA mu value
        updater = ConfigUpdater()
        updater.read(hybridModelDFACfg)
        updater['GENERAL']['mu'].value = muArray[-1]
        updater.update_file()

        # Run dense flow with coulomb friction
        dem, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=hybridModelDFACfg)
        particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName='', flagAvaDir=True, comModule='com1DFA')
        # postprocess to extract path and energy line
        avaProfilePart, avaProfileMass, avaProfileKE = particleTools.getCom1DFAPath(particlesList, demOri)
        avaProfileMassExt = com1DFATools.extendCom1DFAPath(cfgHybrid, demOri, particlesList[0], avaProfileMass.copy())
        # save profile as AB profile in Inputs
        pathAB = pathlib.Path(avalancheDir, 'Inputs', 'LINES', 'pathAB_aimec')
        name = 'massAvaPath'
        shpConv.writeLine2SHPfile(avaProfileMassExt, name, pathAB)
        # also save it to work
        pathAB = pathlib.Path(avalancheDir, 'Work', 'com3Hybrid', 'pathAB_aimec' + str(iteration))
        name = 'massAvaPath'
        shpConv.writeLine2SHPfile(avaProfileMassExt, name, pathAB)

        # Run Alpha Beta
        hybridModelABCfg = pathlib.Path('com3Hybrid', 'hybridModel_com2ABCfg.ini')
        cfgAB = cfgUtils.getModuleConfig(com2AB, fileOverride=hybridModelABCfg)
        # take the path extracted from the DFA model as input
        resAB = com2AB.com2ABMain(cfgAB, avalancheDir)
        # make AB plot
        reportDictList = []
        _, plotFile, writeFile = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

        # make custom com3 profile plot
        resAB = outAB.processABresults(resAB, name)
        figDict = outCom3Plots.updateProfilePlot(resAB, name, figDict, iteration)
        figDict = outCom3Plots.updatePathPlot(demOri, avaProfileMassExt, resAB, name, figDict, iteration)

        # update alpha (result from AB model on new profile)
        alpha = resAB[name]['alpha']
        log.info('Alpha is %.2f' % alpha)
        # append new alpha and corresponding mu value to the list
        alphaArray = np.append(muArray, alpha)
        muArray = np.append(muArray, np.tan(np.radians(alpha)))
        # keep iterating if the change in alpha is too big
        iterate = keepIterating(cfgHybrid, alphaArray)
        iteration = iteration + 1

    pathDict, rasterTransfo, resAnalysisDF = runAna3AIMEC.runAna3AIMEC(avalancheDir=avalancheDir)
    simID = simDF.index[0]
    refSimulation = pathDict['refSimulation']

    # fetch fields for desired time step
    timeSteps = [timeStepInfo[-1]]
    fields = com1DFATools.fetchField(timeSteps, 'pta', avalancheDir, inputDir='')
    outCom3Plots.finalizePathPlot(avalancheDir, figDict, resAnalysisDF, refSimulation, dem, demOri, particlesList[-1], fields[0])
    outCom3Plots.finalizeProfilePlot(avalancheDir, figDict, resAnalysisDF, refSimulation)

    hybridModelDFACfg = pathlib.Path('com3Hybrid', 'hybridModel_com1DFACfg.ini')
    cfgDFA = cfgUtils.getModuleConfig(com1DFA, fileOverride=hybridModelDFACfg)
    outCom3Plots.plotEnergyProfile(avalancheDir, cfgDFA, resAB, name, simID, demOri, avaProfileMass)


def keepIterating(cfgHybrid, alphaArray):
    """Check the change in alpha between two iterations

    If this change is bigger than the alphaThreshold, keep iterationg and return True
    """
    alphaThreshold = cfgHybrid.getfloat('GENERAL', 'alphaThreshold')
    iterate = np.abs(alphaArray[-1] - alphaArray[-2]) >= alphaThreshold
    return iterate
