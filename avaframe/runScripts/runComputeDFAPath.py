"""
    Create an avalanche path from DFA simulation results and create a split point
    Configuration should be specified in DFAPathGenerationCfg.ini
    or the local version of it.
    It is possible to generate a path from particles or from fields.
    From particles, you need to save the particles dictionaries at
    multiple time steps (first, some in the middle and last)
    From fields, you need the FT, FM at multiple time steps
"""
import pathlib
import numpy as np
# Local imports
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana5Utils.DFAPathGeneration as DFAPath
from avaframe.in1Data import getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.geoTrans as gT
# import plotting tools
from avaframe.out3Plot import outCom3Plots

# set to true if you want to run com1DFA first to create the DFA input needed for the path generation
runDFAModule = True

# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runComputeDFAPath'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avaList = [cfgMain['MAIN']['avalancheDir']]
# name of avalanche directory as list, multiple possible
# avaList = ['avaAlr', 'avaHit', 'avaKot', 'avaMal', 'avaWog', 'avaGar', 'avaHelix',
#            'avaBowl', 'avaParabola', 'avaHockeySmall', 'avaHockeyChannel']
# avaList = ['data/' + name for name in avaList]

for avaName in avaList:
    # set avaDir
    avalancheDir = pathlib.Path(avaName)
    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)

    if runDFAModule:
        # Clean input directory of old work and output files from module
        initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)
        # create com1DFA configuration
        # read the default one (and only this one, no local is read)
        com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False, onlyDefault=True)
        # and adapt settings. Especially one needs to save multiple intermediate time steps in order to compute the
        # path. One also needs to save the correct particle or fields variables in order to compute the path
        com1DFACfg['GENERAL']['tSteps'] = '0:5'
        com1DFACfg['GENERAL']['resType'] = 'pta|particles'
        com1DFACfg['GENERAL']['simTypeList'] = 'null'
        com1DFACfg['GENERAL']['sphKernelRadiusTimeStepping'] = 'True'
        com1DFACfg['GENERAL']['cMax'] = '0.01'
        com1DFACfg['GENERAL']['sphOption'] = '2'
        com1DFACfg['GENERAL']['massPerParticleDeterminationMethod'] = 'MPPKR'
        com1DFACfg['GENERAL']['explicitFriction'] = '1'
        com1DFACfg['GENERAL']['frictModel'] = 'Coulomb'
        com1DFACfg['GENERAL']['mu'] = '0.42'
        outDir = pathlib.Path(avalancheDir, 'Outputs', 'DFAPath')
        fU.makeADir(outDir)
        # write configuration to file
        com1DFACfgFile = outDir / 'com1DFAPathGenerationCfg.ini'
        with open(com1DFACfgFile, 'w') as configfile:
            com1DFACfg.write(configfile)
        # call com1DFA and perform simulations
        dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=com1DFACfgFile)
    else:
        # read simulation dem
        demOri = gI.readDEM(avalancheDir)
        dem = com1DFA.setDEMoriginToZero(demOri)
        dem['originalHeader'] = demOri['header'].copy()
        # load DFA results (use runCom1DFA to generate these results for example)
        # here is an example with com1DFA but another DFA computational module can be used
        # as long as it produces some particles or FV, FD and FM results
        # create dataFrame of results (read DFA data)
        simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

    # load config for path generation (fron DFAPathGenerationCfg.ini or its local)
    DFAPathCfg = cfgUtils.getModuleConfig(DFAPath)
    # loop over all simulations
    for simName, simDFrow in simDF.iterrows():
        log.info('Computing avalanche path from simulation: %s', simName)
        pathFromPart = DFAPathCfg['PATH'].getboolean('pathFromPart')
        resampleDistance = DFAPathCfg['PATH'].getfloat('nCellsResample') * dem['header']['cellsize']
        # get the mass average path
        avaProfileMass, particlesIni = DFAPath.generateAveragePath(avalancheDir, pathFromPart, simName, dem)
        avaProfileMass, _ = gT.prepareLine(dem, avaProfileMass, distance=resampleDistance, Point=None)
        avaProfileMass['indStartMassAverage'] = 1
        avaProfileMass['indEndMassAverage'] = np.size(avaProfileMass['x'])
        # make the parabolic fit
        parabolicFit = DFAPath.getParabolicFit(DFAPathCfg['PATH'], avaProfileMass, dem)
        # here the avaProfileMass given in input is overwriten and returns only a x, y, z extended profile
        avaProfileMass = DFAPath.extendDFAPath(DFAPathCfg['PATH'], avaProfileMass, dem, particlesIni)
        # resample path and keep track of start and end of mass averaged part
        avaProfileMass = DFAPath.resamplePath(DFAPathCfg['PATH'], dem, avaProfileMass)
        # get split point
        splitPoint = DFAPath.getSplitPoint(DFAPathCfg['PATH'], avaProfileMass, parabolicFit)
        # make analysis and generate plots
        _ = outCom3Plots.generateCom1DFAPathPlot(avalancheDir, DFAPathCfg['PATH'], avaProfileMass, dem,
                                                 parabolicFit, splitPoint, simName)
        # now save the path and split point as shape file
        DFAPath.saveSplitAndPath(avalancheDir, simDFrow, splitPoint, avaProfileMass, dem)
