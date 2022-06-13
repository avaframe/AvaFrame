"""
    Run script for running the standard tests with com1DFA
    in this test all the available tests tagged standardTest are performed
"""

# Load modules
import time
import pathlib
import shutil
from configupdater import ConfigUpdater

# Local imports
from avaframe.version import getVersion
from avaframe.com1DFA import com1DFA
from avaframe.ana1Tests import testUtilities as tU
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# +++++++++REQUIRED+++++++++++++
# Which result types do we want to save in the benchmarks
outputVariable = ['ppr', 'pft', 'pfv']
# ++++++++++++++++++++++++++++++

# log file name; leave empty to use default runLog.log
logName = 'runUpdateBenchmarkTestsCom1DFA'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False)

# filter benchmarks to extract only desired ones
type = 'TAGS'
valuesList = ['standardTest']
testList = tU.filterBenchmarks(testDictList, type, valuesList, condition='and')

# Set directory for full standard test report
outDir = pathlib.Path('..', 'benchmarks')
log = logUtils.initiateLogger(outDir, logName)
fU.makeADir(outDir)

version = getVersion()

log.info('The following benchmark tests will be updated ')
for test in testList:
    log.info('%s' % test['NAME'])

# run Standard Tests sequentially
for test in testList:

    avaDir = test['AVADIR']

    simNameRef = test['simNameRef']
    refDir = pathlib.Path('..', 'benchmarks', test['NAME'])

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir, keep=logName)

    # RunCom1DFA
    # Load input parameters from configuration file for standard tests
    avaName = pathlib.Path(avaDir).name
    standardCfg = refDir / ('%s_com1DFACfg.ini' % test['AVANAME'])
    modName = 'com1DFA'

    # change the ini file to force particle initialization
    updater = ConfigUpdater()
    updater.read(standardCfg)
    updater['GENERAL']['initialiseParticlesFromFile'].value = 'False'
    updater.update_file()
    # Set timing
    startTime = time.time()
    # call com1DFA run
    dem, plotDict, reportDictList, simDF = com1DFA.com1DFAMain(avaDir, cfgMain, cfgFile=standardCfg)
    endTime = time.time()
    timeNeeded = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeeded)))

    # change the ini file to read particles from file
    updater = ConfigUpdater()
    updater.read(standardCfg)
    updater['GENERAL']['initialiseParticlesFromFile'].value = 'False'
    updater['GENERAL']['particleFile'].value = '../benchmarks/' + test['NAME']
    updater.update_file()

    # Update benchmarks
    # copy Simulation Parameters to benchmark dict
    rep = reportDictList[0]
    test['simName'] = rep['simName']
    test['Simulation Parameters'] = rep['Simulation Parameters']
    # here we need to reset the version because after updating the first benchmark,
    # the version will be marked as dirty...
    test['Simulation Parameters']['Program version'] = version
    test['simNameRef'] = rep['simName']['name']
    # get results file names (.asc) add them to the dictionary and copy the files to benchmark
    # first clean the benchmark directory
    ascFiles = list(refDir.glob('*.asc'))
    for file in ascFiles:
        if file.is_file():
            file.unlink()
    partDir = refDir / 'particles'
    if partDir.is_dir():
        shutil.rmtree(partDir)
    # set copy peak results
    resDir = pathlib.Path(avaDir, 'Outputs', modName, 'peakFiles')
    simName = rep['simName']['name']
    files = []
    for suf in outputVariable:
        simFile = resDir / (simName + '_' + suf + '.asc')
        if simFile.is_file():
            # add file name to dict
            files.append(simFile.stem)
            # copy file to benchmark
            destFile = refDir / (simName + '_' + suf + '.asc')
            shutil.copy2(simFile, destFile)
        else:
            log.warning('did not find the file %s' % simFile)
    # set copy mass result
    resDir = pathlib.Path(avaDir, 'Outputs', modName)
    simFile = resDir / ('mass_' + simName + '.txt')
    if simFile.is_file():
        # add file name to dict
        files.append(simFile.stem)
        # copy file to benchmark
        destFile = refDir / ('mass_' + simName + '.txt')
        shutil.copy2(simFile, destFile)
    else:
        log.warning('did not find the file %s' % simFile)
    test['FILES'] = files
    # set copy particles
    resDir = pathlib.Path(avaDir, 'Outputs', modName, 'particles')
    simFile = sorted(list(resDir.glob('*.p')))[0]
    if simFile.is_file():
        # copy file to benchmark
        destFile = refDir / 'particles'
        destFile.mkdir()
        simComponents = rep['simName']['name'].split('_')
        destFile = destFile / (simComponents[0] + '_' + simComponents[2])
        destFile.mkdir()
        destFile = destFile / (simFile.name)
        shutil.copy2(simFile, destFile)
    else:
        log.warning('did not find the file %s' % simFile)
    # write the benchmark dict as JSON
    tU.writeDesDicttoJson(test, test['NAME'], refDir)
