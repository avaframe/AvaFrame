"""
    Run script for running python DFA kernel
"""

import time
import copy
import os

# Local imports
import avaframe.in3Utils.initialiseDirs as inDirs
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in1Data import getInput as gI
import avaframe.com1DFAPy.com1DFA as com1DFA
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in3Utils.fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def setDEMorigin(demOri):
    """ set origin of DEM to 0,0 """

    dem = copy.deepcopy(demOri)
    dem['header'].xllcenter = 0
    dem['header'].yllcenter = 0
    dem['header'].xllcorner = 0
    dem['header'].yllcorner = 0

    return dem


def runCom1DFAPy(avaDir='', cfgFile='', relTh='', flagAnalysis=True):
    """ run com1DFAPy module """

    # +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
    # log file name; leave empty to use default runLog.log
    logName = 'testKernel'

    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avaDir != '':
        avalancheDir = avaDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']
    # set module name, reqiured as long we are in dev phase
    # - because need to create e.g. Output folder for com1DFAPy to distinguish from
    # current com1DFA
    modName = 'com1DFAPy'

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    # Load configuration
    if cfgFile !='':
        cfg = cfgUtils.getModuleConfig(com1DFA, cfgFile)
    else:
        cfg = cfgUtils.getModuleConfig(com1DFA)
    cfgGen = cfg['GENERAL']
    flagDev = cfg['FLAGS'].getboolean('flagDev')

    # for timing the sims
    startTime = time.time()


    # +++++++++Inputs++++++++++++++++++++++++
    # ------------------------
    # fetch input data - dem, release-, entrainment- and resistance areas
    demFile, relFiles, entFiles, resFile, flagEntRes = gI.getInputData(
        avalancheDir, cfg['FLAGS'], flagDev)
    demOri = IOf.readRaster(demFile)
    # derive line from release area polygon
    releaseLine = shpConv.readLine(relFiles[0], 'release1', demOri)
    dem = setDEMorigin(demOri)

    # -----------------------
    # Initialize mesh
    dem = com1DFA.initializeMesh(dem)
    # ------------------------
    # process release info to get it as a raster
    relRaster = com1DFA.prepareArea(releaseLine, demOri)
    if len(relTh) == 0:
        relTh = cfgGen.getfloat('relTh')

    relRaster = relRaster * relTh

    # ------------------------
    # initialize simulation : create directories
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName)
    # create particles, create resistance and
    # entrainment matrix, initialize fields, get normals and neighbours
    particles, fields, Cres, Ment = com1DFA.initializeSimulation(cfgGen, relRaster, dem, avalancheDir)

    # +++++++++PERFORM SIMULAITON++++++++++++++++++++++
    # ------------------------
    #  Start time step computation
    Tsave, T, U, Z, S, Particles, Fields, Tcpu = com1DFA.DFAIterate(
        cfgGen, particles, fields, dem, Ment, Cres)

    tcpuDFA = time.time() - startTime
    log.info(('cpu time DFA = %s s' % (tcpuDFA)))

    # export particles dictionaries of saving time steps
    outDirData = os.path.join(outDir, 'particles')
    fU.makeADir(outDirData)
    com1DFA.savePartToPickle(Particles, outDirData)

    if flagAnalysis:
        # +++++++++POSTPROCESS++++++++++++++++++++++++
        # -------------------------------
        # Analyse resutls
        outDirPlots = os.path.join(outDir, 'reports')
        fU.makeADir(outDirPlots)
        com1DFA.analysisPlots(Particles, Fields, cfg, demOri, dem, outDirPlots)

        # +++++++++EXPORT RESULTS AND PLOTS++++++++++++++++++++++++
        # Result parameters to be exported
        com1DFA.exportFields(cfgGen, Tsave, Fields, relFiles[0], demOri, outDir)

        # Generate plots for all peakFiles
        plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'], modName)
        # Setup simulation dictionaries for report genereation
        reportDictList = []
        # Create dictionary
        reportNull = {}
        reportNull = {'headerLine': {'type': 'title', 'title': 'com1DFA Simulation'},
        'simName': {'type': 'simName', 'name': 'null'},
            'Simulation Parameters': {
                'type': 'list',
                'Release Area Scenario': os.path.splitext(os.path.basename(relFiles[0]))[0],
                'Release Area': releaseLine['Name'],
                'Entrainment Area': '',
                'Resistance Area': '',
                'Parameter variation on': '',
                'Parameter value': '',
                'Mu': cfgGen['mu'],
                'Release thickness [m]': cfgGen.getfloat('relTh')},
            'Release area': {'type': 'columns', 'Release area scenario': os.path.splitext(os.path.basename(relFiles[0]))[0]}}

        if relTh is not float:
            reportNull['Release thickness [m]'] = 'release function'
        # Add to report dictionary list
        reportDictList.append(reportNull)

        return Particles, Fields, Tsave, dem, plotDict, reportDictList
    else:
        return Particles, Fields, Tsave, dem


if __name__ == "__main__":
    runCom1DFAPy()
