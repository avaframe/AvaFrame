"""
    Create an avalanch path from DFA simulation results
    Configuration should be specified in DFAPathGenerationCfg.ini
    or the local version of it.
    It is possible to generate a path from particles or from fields.
    From particles, you need to save the particles dictionaries at
    multiple time steps (first, some in the middle and last)
    From fields, you need the FT, FM at multiple time steps
"""
import pathlib
# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana5Utils.DFAPathGeneration as DFAPath
from avaframe.in1Data import getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.com1DFA import com1DFA

# set to true if you want to run com1DFA first to create the DFA input needed for the path generation
runDFAModule = True

# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runComputeDFAPath'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

if runDFAModule:
    # Clean input directory of old work and output files from module
    initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)
    # create com1DFA configuration
    # read the default one (and only this one, no local is read)
    com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', toPrint=False, onlyDefault=False)
    # and adapt setting (especially time steps saved and outpus)
    com1DFACfg['GENERAL']['tSteps'] = '0:5'
    com1DFACfg['GENERAL']['resType'] = 'particles'
    com1DFACfg['GENERAL']['simTypeList'] = 'null'
    outDir = pathlib.Path(avalancheDir, 'Outputs', 'DFAPath')
    fU.makeADir(outDir)
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
    # as long as it produces some particles or FV, FT and FM results
    # create dataFrame of results (read DFA data)
    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)


DFAPathCfgFile = pathlib.Path('ana5Utils', 'DFAPathGenerationCfg.ini')
DFAPathCfg = cfgUtils.getModuleConfig(DFAPath)
for simName, simDFrow in simDF.iterrows():
    log.info('Computing avalanche path from simulation: %s', simName)
    pathFromPart = DFAPathCfg['PATH'].getboolean('pathFromPart')
    avaProfileMass = DFAPath.generateAveragePath(avalancheDir, pathFromPart, simName, dem)
    # here the avaProfileMass given in input is overwriten and returns only a x, y, z extended profile
    avaProfileMass = DFAPath.extendDFAPath(avalancheDir, DFAPathCfg['PATH'], dem, simName, avaProfileMass)
    # put path back in original location
    avaProfileMass['x'] = avaProfileMass['x'] + dem['originalHeader']['xllcenter']
    avaProfileMass['y'] = avaProfileMass['y'] + dem['originalHeader']['yllcenter']
    # save profile in Inputs
    pathAB = pathlib.Path(avalancheDir, 'Outputs', 'DFAPath', 'massAvgPath_%s_AB_aimec' % simName)
    name = 'massAvaPath'
    shpConv.writeLine2SHPfile(avaProfileMass, name, pathAB)
    log.info('Saved path to: %s', pathAB)
