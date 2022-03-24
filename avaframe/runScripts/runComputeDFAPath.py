"""
    Create an avalanch path from DFA simulation results
    Configuration should be specified in DFAPathGenerationCfg.ini
    or the local version of it.
    It is possible to generate a path from particles or from fields.
    From particles, you need to save the particles dictionaries at
    multiple time steps (first, some in the middle and last)
    From fields, you need the FD, FM at multiple time steps
"""
import pathlib
# Local imports
import avaframe.in2Trans.shpConversion as shpConv
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.ana5Utils.DFAPathGeneration as DFAPath
from avaframe.in1Data import getInput as gI
from avaframe.com1DFA import com1DFA


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

# call com1DFA and perform simulations
demOri = gI.readDEM(avalancheDir)
dem = com1DFA.setDEMoriginToZero(demOri)
dem['originalHeader'] = demOri['header'].copy()
# here is an example with com1DFA but another DFA computational module can be used
# as long as it produces some FV, FD and FM results
getDFAPathCfgFile = pathlib.Path('ana5Utils', 'DFAPathGenerationCfg.ini')
getDFAPathCfg = cfgUtils.getModuleConfig(DFAPath)

# create dataFrame of results (read DFA data)
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)
for simName, simDFrow in simDF.iterrows():
    log.info('Computing avalanche path from simulation: %s', simName)
    pathFromPart = getDFAPathCfg['PATH'].getboolean('pathFromPart')
    avaProfileMass = DFAPath.generateAveragePath(avalancheDir, pathFromPart, simName, dem)
    # here the avaProfileMass given in input is overwriten and returns only a x, y, z extended profile
    avaProfileMass = DFAPath.extendDFAPath(avalancheDir, getDFAPathCfg['PATH'], dem, simName, avaProfileMass)
    # put path back in original location
    avaProfileMass['x'] = avaProfileMass['x'] + dem['originalHeader']['xllcenter']
    avaProfileMass['y'] = avaProfileMass['y'] + dem['originalHeader']['yllcenter']
    # save profile in Inputs
    pathAB = pathlib.Path(avalancheDir, 'Outputs', 'DFAPath', 'path_%s_AB_aimec' % simName)
    name = 'massAvaPath'
    shpConv.writeLine2SHPfile(avaProfileMass, name, pathAB)
    log.info('Saved path to: %s', pathAB)
