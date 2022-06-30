"""
Rotation test

This module runs a DFA simulation for different grid alignments
and compares the results
"""
import logging
import shutil
from datetime import datetime
import pathlib

# Local imports
# import config and init tools
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import geoTrans as gT
import avaframe.in2Trans.ascUtils as IOf
from avaframe.version import getVersion
# import computation modules
from avaframe.com1DFA import com1DFA
# import analysis modules
from avaframe.ana1Tests import energyLineTest
# import plotting module
import avaframe.out3Plot.plotUtils as pU
# create local logger
log = logging.getLogger(__name__)


def mainRotationTest(avalancheDir, dem, simDF, resTypeList, flagMass):
    """This is the core function of the Rotation Test module

    This module runs the energy line test and rotates the simulation results for the aimec analysis
    """
    energyLineTestCfgFile = pathlib.Path('ana1Tests', 'energyLineTest_com1DFACfg.ini')
    energyLineTestCfg, modInfo = cfgUtils.getModuleConfig(com1DFA, fileOverride=energyLineTestCfgFile, modInfo=True)
    inDir = avalancheDir / 'Outputs' / 'com1DFA'
    outDir = avalancheDir / 'Outputs' / 'com1DFARotated'
    fU.makeADir(outDir)
    for simHash in simDF.index:
        # rotate results to be able to proceede to the aimec analysis
        simName = simDF.loc[simHash, 'simName']
        if 'ent' in simName:
            flagMass = True
        log.info('Rotating simulation: %s' % simName)
        relName = (simName.split('_'))[0]
        theta = float(relName[3:])
        simDF.loc[simHash, 'relAngle'] = theta
        # copy mass file
        massFile = list(inDir.glob('mass_*' + simName + '*.txt'))
        shutil.copy(massFile[0], outDir)
        for resType in resTypeList:
            fName = simDF.loc[simHash, resType]
            rasterDict = IOf.readRaster(fName)
            rotatedRaster = gT.rotateRaster(rasterDict, theta, deg=True)
            fU.makeADir(outDir / 'peakFiles')
            outFileName = outDir / 'peakFiles' / fName.name
            IOf.writeResultToAsc(rotatedRaster['header'], rotatedRaster['rasterData'], outFileName, flip=False)

        # make energy line analysis
        errorEnergyTest, savePath = energyLineTest.mainEnergyLineTest(avalancheDir, energyLineTestCfg, simName, dem)
        simDF.loc[simHash, 'runOutSError'] = errorEnergyTest['runOutSError']
        simDF.loc[simHash, 'runOutZError'] = errorEnergyTest['runOutZError']
        simDF.loc[simHash, 'runOutAngleError'] = errorEnergyTest['runOutAngleError']
        simDF.loc[simHash, 'rmseVelocityElevation'] = errorEnergyTest['rmseVelocityElevation']
        simDF.loc[simHash, 'energyLinePlotPath'] = savePath

    return simDF, flagMass


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


def initializeRotationTestReport(avalancheDir, resTypeList, comModule):
    """Initialize report for rotation test

    """
    dateTimeInfo = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    reportRotationTest = {'headerLine': {'type': 'title', 'title': 'Rotation Test for DFA Simulation'},
                          'avaName': {'type': 'avaName', 'name': str(avalancheDir)},
                          'time': {'type': 'time', 'time': dateTimeInfo},
                          'Simulation Parameters': {
                          'type': 'list',
                          'Program version': getVersion(),
                          'DFA module': comModule},
                          'Rotation test input simulations': {
                          'type': 'pandasDF', 'column names': {'simName': 'Simulation name',
                                                               'releaseScenario': 'Release name', 'relAngle': 'Angle °'}},
                          'Rotation test Energy line result table': {
                          'type': 'pandasDF', 'column names': {'simName': 'Simulation name', 'relAngle': 'Angle °',
                                                               'runOutSError': 's Error [m]', 'runOutZError': 'z Error [m]',
                                                               'runOutAngleError': 'angle Error [°]',
                                                               'rmseVelocityElevation': 'velocity elevation rsme [m]'}},
                          'Rotation test Energy line figures': {'type': 'image'},
                          'Rotation test AIMEC result table': {
                          'type': 'pandasDF', 'column names': {'simName': 'Simulation name', 'relAngle': 'Angle °',
                                                               'sRunout': 'Runout [m]', 'TP': 'True positive area [-]',
                                                               'FP': 'False positive area [-]',
                                                               'FN': 'False negative area [-]'}},
                          'Rotation test AIMEC figures': {'type': 'image'}}
    resTypeList = list(set(resTypeList).intersection(['ppr', 'pfv', 'pft']))
    for resType in resTypeList:
        colName = 'max%sCrossMax' % resType
        colLabel = 'Maximum %s' % resType
        unit = pU.cfgPlotUtils['unit' + resType]
        colLabel = '%s [%s]' % (colLabel, unit)
        reportRotationTest['Rotation test AIMEC result table']['column names'][colName] = colLabel
    return reportRotationTest


def buildRotationTestReport(avalancheDir, reportRotationTest, simDF, resAnalysisDF, aimecPlotDict, flagMass):
    outDir = avalancheDir / 'Outputs' / 'ana1Tests' / 'rotationTest'
    fU.makeADir(outDir)
    energyLinePlotsDF = simDF.set_index('simName')['energyLinePlotPath']
    energyLinePlotsDict = energyLinePlotsDF.to_dict()
    inPlotsDict = {'energyLinePlots': energyLinePlotsDict}
    inPlotsDict.update(aimecPlotDict['slCompPlot'])
    inPlotsDict.update(aimecPlotDict['areasPlot'])
    if flagMass:
        inPlotsDict.update(aimecPlotDict['massAnalysisPlot'])
    outPlotDict = gR.copyPlots(inPlotsDict, outDir)
    # add plot info to general report Dict
    # reportRotationTest['Simulation Results'] = plotPaths

    simDF = simDF.sort_values(by=['relAngle'], ascending=True)
    resAnalysisDF = resAnalysisDF.sort_values(by=['relAngle'], ascending=True)

    reportRotationTest['Rotation test input simulations']['dataDF'] = simDF
    reportRotationTest['Rotation test Energy line result table']['dataDF'] = simDF
    reportRotationTest['Rotation test Energy line figures'].update(outPlotDict['energyLinePlots'])
    reportRotationTest['Rotation test AIMEC result table']['dataDF'] = resAnalysisDF
    reportRotationTest['Rotation test AIMEC figures'].update(
        inPlotsDict['Aimec comparison of mean and max values along path'])
    reportRotationTest['Rotation test AIMEC figures'].update(outPlotDict['Aimec area analysis'])

    if flagMass:
        reportRotationTest['Rotation test AIMEC figures'].update(outPlotDict['Aimec mass analysis'])
    # add energy line and aimec results table to report
    gR.writeReport(outDir, [reportRotationTest], True, plotDict='', standaloneReport=False,
                   reportName='rotationTestReport.md')
