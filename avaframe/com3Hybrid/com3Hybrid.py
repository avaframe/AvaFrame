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

# Local imports
# import config and init tools
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in3Utils import cfgUtils
from avaframe.in1Data import getInput

# import computation modules
from avaframe.com1DFA import com1DFA, particleTools
import avaframe.ana5Utils.DFAPathGeneration as DFAPath
from avaframe.com2AB import com2AB

# import analysis tools
from avaframe.out3Plot import outAB

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

    # get initial mu value
    muArray = np.array([cfgHybrid.getfloat('DFA', 'mu')])
    alphaArray = np.array([np.degrees(np.arctan(muArray[0]))])

    # prepare for iterating
    nIterMax = cfgHybrid.getint('GENERAL', 'nIterMax')
    iteration = 0
    iterate = True
    resultsHybrid = {}
    while iteration < nIterMax and iterate:
        # update the com1DFA mu value
        updater = ConfigUpdater()
        updater.read(hybridModelDFACfg)
        updater['GENERAL']['mu'].value = ('%.4f' % muArray[-1])
        updater.update_file()

        # Run dense flow with coulomb friction
        dem, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=hybridModelDFACfg)
        simID = simDF.index[0]
        particlesList, timeStepInfo = particleTools.readPartFromPickle(avalancheDir, simName=simID, flagAvaDir=True,
                                                                       comModule='com1DFA')
        # postprocess to extract path and energy line
        avaProfileMass = DFAPath.getDFAPathFromPart(particlesList, addVelocityInfo=True)
        # make a copy because extendDFAPathKernel might modify avaProfileMass
        avaProfileMassExt = DFAPath.extendDFAPath(cfgHybrid['PATH'], avaProfileMass.copy(), dem, particlesList[0])
        avaProfileMassExtOri = copy.deepcopy(avaProfileMassExt)
        avaProfileMassExtOri['x'] = avaProfileMassExtOri['x'] + demOri['header']['xllcenter']
        avaProfileMassExtOri['y'] = avaProfileMassExtOri['y'] + demOri['header']['yllcenter']
        # also save it to work
        pathAB = pathlib.Path(avalancheDir, 'Outputs', 'com3Hybrid', 'pathAB_aimec' + str(iteration))
        name = 'massAvaPath'
        shpConv.writeLine2SHPfile(avaProfileMassExtOri, name, pathAB)

        # Run Alpha Beta
        hybridModelABCfg = pathlib.Path('com3Hybrid', 'hybridModel_com2ABCfg.ini')
        cfgAB = cfgUtils.getModuleConfig(com2AB, fileOverride=hybridModelABCfg)
        cfgAB['ABSETUP']['path2Line'] = str(pathAB) + '.shp'
        # take the path extracted from the DFA model as input
        pathDict, demAB, splitPoint, eqParams, resAB = com2AB.com2ABMain(cfgAB, avalancheDir)
        # make AB plot
        reportDictList = []
        _, plotFile, writeFile = outAB.writeABpostOut(pathDict, demAB, splitPoint, eqParams, resAB, cfgAB, reportDictList)

        # make custom com3 profile plot
        alpha = resAB[name]['alpha']
        indAlpha = resAB[name]['indAlpha']
        indBetaPoint = resAB[name]['indBetaPoint']
        sBetaPoint = resAB[name]['s'][indBetaPoint]
        xAB = resAB[name]['x'][indAlpha]
        yAB = resAB[name]['y'][indAlpha]
        resultsHybrid['iterration' + str(iteration)] = {'path': avaProfileMassExt, 'alpha': alpha,
                                                        'xAB': xAB, 'yAB': yAB, 'sBetaPoint': sBetaPoint}

        # update alpha (result from AB model on new profile)
        alpha = resAB[name]['alpha']
        log.info('Alpha is %.2f' % alpha)
        # append new alpha and corresponding mu value to the list
        alphaArray = np.append(alphaArray, alpha)
        muArray = np.append(muArray, np.tan(np.radians(alpha)))
        # keep iterating if the change in alpha is too big
        iterate = keepIterating(cfgHybrid, alphaArray)
        iteration = iteration + 1

    # fetch fields for desired time step
    fields, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['pta'], simName=simID,
                                                       flagAvaDir=True, comModule='com1DFA', timeStep=timeStepInfo[-1])
    outCom3Plots.hybridProfilePlot(avalancheDir, resultsHybrid)
    outCom3Plots.hybridPathPlot(avalancheDir, dem, resultsHybrid, fields[0], particlesList[-1], muArray)

    log.debug('Alpha array is %s' % alphaArray)
    log.debug('mu array is %s' % muArray)


def keepIterating(cfgHybrid, alphaArray):
    """Check the change in alpha between two iterations

    If this change is bigger than the alphaThreshold, keep iterationg and return True
    """
    alphaThreshold = cfgHybrid.getfloat('GENERAL', 'alphaThreshold')
    iterate = np.abs(alphaArray[-1] - alphaArray[-2]) >= alphaThreshold
    return iterate
