"""
    Create an avalanche path from DFA simulation results and create a split point
    Configuration should be specified in DFAPathGenerationCfg.ini
    or the local version of it.
    It is possible to generate a path from particles or from fields.
    From particles, you need to save the particles dictionaries at
    multiple time steps (first, some in the middle and last)
    From fields, you need the FT, FM at multiple time steps
"""
# import general python modules
import pathlib
import numpy as np
import time
import argparse

# local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initializeProject as initProj
import avaframe.ana5Utils.DFAPathGeneration as DFAPath
from avaframe.in1Data import getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.com1DFA import com1DFA
import avaframe.in3Utils.geoTrans as gT
from avaframe.in3Utils import cfgHandling
from avaframe.out3Plot import outCom3Plots


def runComputeDFAPath(avalancheDir='', runDFAModule=False):
    """
    Run DFA path generation in the default configuration with only an avalanche directory as input.

    Parameters
    ----------
    avalancheDir: str
        path to avalanche directory (setup e.g. with init scripts)
    runDFAModule: bool
        True if com1DFA should be run to create required DFA results before path generation; overrides ini setting

    Returns
    -------
    massAvgPath: str
        file path to the mass-averaged path result saved as a shapefile.
    splitPoint: str
        file path to the split point result saved as a shapefile.
    """

    # Time the whole routine
    startTime = time.time()

    # log file name; leave empty to use default runLog.log
    logName = 'runComputeDFAPath'

    # Load avalanche directory from general configuration file
    # More information about the configuration can be found here
    # on the Configuration page in the documentation
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info("MAIN SCRIPT")
    log.info('Current avalanche: %s', avalancheDir)

    # load config for path generation (from DFAPathGenerationCfg.ini or its local)
    cfgDFAPath = cfgUtils.getModuleConfig(DFAPath)

    # Clean input directory(ies) of old work and output files
    # If you just created the ``avalancheDir`` this one should be clean but if you
    # already did some calculations you might want to clean it::
    initProj.cleanSingleAvaDir(avalancheDir, deleteOutput=False)

    if runDFAModule:
        # Clean input directory of old work and output files from module
        initProj.cleanModuleFiles(avalancheDir, com1DFA, deleteOutput=True)
        # create and read the default com1DFA config (no local is read)
        com1DFACfg = cfgUtils.getModuleConfig(com1DFA, fileOverride='', modInfo=False, toPrint=False,
                                              onlyDefault=cfgDFAPath['com1DFA_com1DFA_override'].getboolean(
                                                  'defaultConfig'))
        # and override with settings from DFAPath config
        com1DFACfg, cfgDFAPath = cfgHandling.applyCfgOverride(com1DFACfg, cfgDFAPath, com1DFA,
                                                                     addModValues=False)
        outDir = pathlib.Path(avalancheDir, 'Outputs', 'DFAPath')
        fU.makeADir(outDir)
        # write configuration to file
        com1DFACfgFile = outDir / 'com1DFAPathGenerationCfg.ini'
        with open(com1DFACfgFile, 'w') as configfile:
            com1DFACfg.write(configfile)
        # call com1DFA and perform simulations
        dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(cfgMain, cfgInfo=com1DFACfgFile)
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

    #Clean DFAPath output in avalanche directory
    initProj.cleanModuleFiles(avalancheDir, DFAPath, deleteOutput=True)

    for simName, simDFrow in simDF.iterrows():
        log.info('Computing avalanche path from simulation: %s', simName)
        pathFromPart = cfgDFAPath['PATH'].getboolean('pathFromPart')
        resampleDistance = cfgDFAPath['PATH'].getfloat('nCellsResample') * dem['header']['cellsize']
        # get the mass average path
        avaProfileMass, particlesIni = DFAPath.generateAveragePath(avalancheDir, pathFromPart, simName, dem,
                                                                   addVelocityInfo=cfgDFAPath['PATH'].getboolean('addVelocityInfo'))
        avaProfileMass, _ = gT.prepareLine(dem, avaProfileMass, distance=resampleDistance, Point=None)
        avaProfileMass['indStartMassAverage'] = 1
        avaProfileMass['indEndMassAverage'] = np.size(avaProfileMass['x'])
        # make the parabolic fit
        parabolicFit = DFAPath.getParabolicFit(cfgDFAPath['PATH'], avaProfileMass, dem)
        # here the avaProfileMass given in input is overwritten and returns only an x, y, z extended profile
        avaProfileMass = DFAPath.extendDFAPath(cfgDFAPath['PATH'], avaProfileMass, dem, particlesIni)
        # resample path and keep track of start and end of mass averaged part
        avaProfileMass = DFAPath.resamplePath(cfgDFAPath['PATH'], dem, avaProfileMass)
        # get split point
        splitPoint = DFAPath.getSplitPoint(cfgDFAPath['PATH'], avaProfileMass, parabolicFit)
        # make analysis and generate plots
        _ = outCom3Plots.generateCom1DFAPathPlot(avalancheDir, cfgDFAPath['PATH'], avaProfileMass, dem,
                                                 parabolicFit, splitPoint, simName)
        # now save the path and split point as shapefiles
        massAvgPath,splitPoint = DFAPath.saveSplitAndPath(avalancheDir, simDFrow, splitPoint, avaProfileMass, dem)

    log.info("DFA path generation completed")

    # Print time needed
    endTime = time.time()
    log.info('Took %6.1f seconds to calculate.' % (endTime - startTime))

    return str(massAvgPath), str(splitPoint)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run computeDFAPath workflow')
    parser.add_argument('avadir', metavar='avalancheDir', type=str, nargs='?', default='',
                        help='the avalanche directory')
    parser.add_argument('--runDFA', action="store_true",
                        help='set flag to run com1DFA with cfg override to generate usable example results before'
                             ' generating path')

    args = parser.parse_args()
    runComputeDFAPath(str(args.avadir), bool(args.runDFA))
