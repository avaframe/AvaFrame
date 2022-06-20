"""
Rotation test

This module runs a DFA simulation for different grid alignments
and compares the results
"""
import pathlib
import logging

# Local imports
# import config and init tools
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import geoTrans as gT
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import initializeProject as iP
# import computation modules
from avaframe.com1DFA import com1DFA

# import analysis modules
from avaframe.runScripts import runAna3AIMEC
from avaframe.ana3AIMEC import ana3AIMEC
from avaframe.ana1Tests import energyLineTest
# create local logger
log = logging.getLogger(__name__)


def mainRotationTest(cfgMain, cfgAimec, runDFAModule):
    """This is the core function of the Rotation Test module

    This module runs a DFA simulation and analyzes the results
    """
    avalancheDir = cfgMain['MAIN']['avalancheDir']
    avalancheDir = pathlib.Path(avalancheDir)
    # get path to com1DFA configuration file used for the energy line test
    rotationTestCfgFile = pathlib.Path('ana1Tests', 'rotationTest_com1DFACfg.ini')
    # run the com1DFA module or load the results from com1DFA
    dem, simDF, resTypeList = com1DFA.runOrLoadCom1DFA(avalancheDir, cfgMain, runDFAModule=runDFAModule,
                                                       cfgFile=rotationTestCfgFile)

    # generate mass averaged path profile
    energyLineTestCfgFile = pathlib.Path('ana1Tests', 'energyLineTest_com1DFACfg.ini')
    energyLineTestCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=energyLineTestCfgFile, modInfo=True)
    for simHash in simDF.index:
        # rotate results to be able to proceede to the aimec analysis
        simName = simDF.loc[simHash, 'simName']
        prepareInputsRotationTest(avalancheDir, simDF, simHash, resTypeList)

        # make energy line analysis
        _ = energyLineTest.mainEnergyLineTest(avalancheDir, energyLineTestCfg, simName, dem)

    # proceede to aimec analysis
    iP.cleanModuleFiles(avalancheDir, ana3AIMEC)
    # write configuration to file
    cfgUtils.writeCfgFile(avalancheDir, ana3AIMEC, cfgAimec)
    runAna3AIMEC.runAna3AIMEC(avalancheDir, cfgAimec)

    return True


def prepareInputsRotationTest(avalancheDir, simDF, simHash, resTypeList):
    log.info('Rotating simulation: %s' % simHash)
    simName = simDF.loc[simHash, 'simName']
    relName = (simName.split('_'))[0]
    theta = float(relName[3:])
    for resType in resTypeList:
        fName = simDF.loc[simHash, resType]
        rasterDict = IOf.readRaster(fName)
        rotatedRaster = gT.rotateRaster(rasterDict, theta, deg=True)
        outDir = avalancheDir / 'Outputs' / 'com1DFARotated' / 'peakFiles'
        fU.makeADir(outDir)
        outFileName = outDir / fName.name
        IOf.writeResultToAsc(rotatedRaster['header'], rotatedRaster['rasterData'], outFileName, flip=False)
