
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
from avaframe.out3Plot import outAna5Plots

# create local logger
log = logging.getLogger(__name__)


def mainAna5Hybrid(cfgMain, cfgHybrid):
    avalancheDir = cfgMain['MAIN']['avalancheDir']
    demOri = getInput.readDEM(avalancheDir)
    # get comDFA configuration path for hybrid model
    hybridModelDFACfg = pathlib.Path('ana5Hybrid', 'hybridModel_com1DFACfg.ini')

    # prepare plots
    figDict = outAna5Plots.initializeFigures()
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
        # ----------------
        # Run dense flow with coulomb friction
        particlesList, fieldsList, _, dem, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=hybridModelDFACfg)
        # postprocess to extract path and energy line
        avaProfilePart, avaProfileMass, avaProfileKE = particleTools.getCom1DFAPath(particlesList, demOri)
        avaProfileMassExt = com1DFATools.extendCom1DFAPath(cfgHybrid, demOri, particlesList[0], avaProfileMass.copy())
        # save profile as AB profile in Inputs
        pathAB = pathlib.Path(avalancheDir, 'Inputs', 'LINES', 'pathAB_aimec')
        name = 'massAvaPath'
        shpConv.writeLine2SHPfile(avaProfileMassExt, name, pathAB)
        # also save it to work
        pathAB = pathlib.Path(avalancheDir, 'Work', 'ana5Hybrid', 'pathAB_aimec' + str(iteration))
        name = 'massAvaPath'
        shpConv.writeLine2SHPfile(avaProfileMassExt, name, pathAB)

        # Run Alpha Beta
        hybridModelABCfg = pathlib.Path('ana5Hybrid', 'hybridModel_com2ABCfg.ini')
        cfgAB = cfgUtils.getModuleConfig(com2AB, fileOverride=hybridModelABCfg)
        # take the path extracted from the DFA model as input
        resAB = com2AB.com2ABMain(cfgAB, avalancheDir)
        # make AB plot
        reportDictList = []
        _, plotFile, writeFile = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

        # make custom ana5 profile plot
        resAB = outAB.processABresults(resAB, name)
        figDict = outAna5Plots.updateProfilePlot(resAB, name, figDict, iteration)
        figDict = outAna5Plots.updatePathPlot(demOri, avaProfileMassExt, resAB, name, figDict, iteration)

        # updat alpha (result from AB model on new profile)
        alpha = resAB[name]['alpha']
        log.info('Alpha is %.2f' % alpha)
        alphaArray = np.append(muArray, alpha)
        muArray = np.append(muArray, np.tan(np.radians(alpha)))
        iterate = keepIterating(cfgHybrid, alphaArray)
        iteration = iteration + 1

    pathDict, rasterTransfo, newRasters, resAnalysis = runAna3AIMEC.runAna3AIMEC(avalancheDir=avalancheDir)
    simID = simDF.index[0]
    indSim = pathDict['simID'].index(simID)
    outAna5Plots.finalizePathPlot(avalancheDir, figDict, resAnalysis, indSim, dem, demOri, particlesList[-1], fieldsList[-1])
    outAna5Plots.finalizeProfilePlot(avalancheDir, figDict, resAnalysis, indSim)

    hybridModelDFACfg = pathlib.Path('ana5Hybrid', 'hybridModel_com1DFACfg.ini')
    cfgDFA = cfgUtils.getModuleConfig(com1DFA, fileOverride=hybridModelDFACfg)
    outAna5Plots.plotHybridRes(avalancheDir, cfgDFA, resAB, name, simID, demOri, avaProfileMass)


def keepIterating(cfgHybrid, alphaArray):
    alphaThreshold = cfgHybrid.getfloat('GENERAL', 'alphaThreshold')
    iterate = np.abs(alphaArray[-1] - alphaArray[-2]) >= alphaThreshold
    return iterate
