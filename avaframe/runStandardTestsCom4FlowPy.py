"""
    Run script for running the standard tests with com4FlowPy
    in this test all the available tests tagged standardTest are performed
"""

# Load modules
import time
import pathlib
import numpy as np
import rasterio

# Local imports
from avaframe.com4FlowPy import com4FlowPy
from avaframe.runCom4FlowPy import readFlowPyinputs
from avaframe.ana1Tests import testUtilities as tU
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def read_raster(path):
    raster = rasterio.open(path)
    output = raster.read(1)
    return output


def compare(path, pathRef):
    raster = read_raster(path)
    rasterRef = read_raster(pathRef)
    diff = rasterRef - raster
    equal = np.array_equal(rasterRef, raster)
    closeArray = np.isclose(raster, rasterRef, rtol=1e-04, equal_nan=True)
    closePercent = np.count_nonzero(closeArray[np.logical_or(raster > 0, rasterRef > 0)]) / rasterRef[np.logical_or(raster > 0, rasterRef > 0)].size
    return diff, equal, closePercent


# Which result types for comparison plots
outputVariable = ['FP_travel_angle', 'z_delta']

# log file name; leave empty to use default runLog.log
logName = 'runStandardTestsCom4FlowPy'

# Load settings from general configuration file
cfgMain = cfgUtils.getGeneralConfig()

# load all benchmark info as dictionaries from description files
testDictList = tU.readAllBenchmarkDesDicts(info=False, inDir=pathlib.Path('..', 'benchmarksCom4FlowPy'))

# filter benchmarks for tag standardTest
#  filterType = 'TAGS'
#  valuesList = ['resistance']
filterType = 'TAGS'
valuesList = ['standardTest']
testList = tU.filterBenchmarks(testDictList, filterType, valuesList, condition='or')

# Set directory for full standard test report
outDir = pathlib.Path.cwd() / 'tests' / 'reportsCom4FlowPy'
fU.makeADir(outDir)

print(outDir)

# Start writing markdown style report for standard tests
reportFile = outDir / 'standardTestsReportCom4FlowPy.md'
with open(reportFile, 'w') as pfile:

    # Write header
    pfile.write('# Standard Tests Report \n')
    pfile.write('## Compare com4FlowPy simulations to benchmark results \n')

log = logUtils.initiateLogger(outDir, logName)
log.info('The following benchmark tests will be fetched ')
for test in testList:
    log.info('%s' % test['NAME'])

# run Standard Tests sequentially

for test in testList:

    avaDir = test['AVADIR']
    cfgMain['MAIN']['avalancheDir'] = avaDir

    # Fetch benchmark test info
    refDir = pathlib.Path('..', 'benchmarksCom4FlowPy', test['NAME'])

    # Clean input directory(ies) of old work and output files
    initProj.cleanSingleAvaDir(avaDir)

    # Load input parameters from configuration file for standard tests
    standardCfg = refDir / ('%s_com4FlowPyCfg.ini' % test['AVANAME'])
    modName = 'com4FlowPy'
    cfg = cfgUtils.getModuleConfig(com4FlowPy, fileOverride=standardCfg)
    cfgGen = cfg["GENERAL"]
    cfgGen["cpuCount"] = str(cfgUtils.getNumberOfProcesses(cfgMain, 9999))

    compDir = pathlib.Path(avaDir, 'Outputs', modName, 'peakFiles')
    avalancheDir = cfgMain["MAIN"]["avalancheDir"]
    cfgPath = readFlowPyinputs(avalancheDir, cfg, log)
    cfgPath["customDirs"] = False
    cfgPath["resDir"] = compDir
    fU.makeADir(cfgPath["resDir"])
    cfgPath["thalwegDir"] = cfgPath["resDir"] / "thalwegData"
    cfgPath["tempDir"] = cfgPath["workDir"] / "temp"
    fU.makeADir(cfgPath["tempDir"])
    cfgPath["deleteTemp"] = "False"

    # Set timing
    startTime = time.time()
    # call com4FlowPy run
    com4FlowPy.com4FlowPyMain(cfgPath, cfgGen)
    endTime = time.time()
    timeNeeded = endTime - startTime
    log.info(('Took %s seconds to calculate.' % (timeNeeded)))

    for variable in outputVariable:
        pathRasterRef = refDir / ('%s.tif' % variable)
        pathRaster = compDir / ('%s.tif' % variable)
        diff, eq, close = compare(pathRaster, pathRasterRef)

        if eq and np.sum(abs(diff[diff != 0])) == 0:
            message = f'{test['NAME']}: for {variable}: rasters are equal \n'
        else:
            message = f'{test['NAME']}: for {variable}: rasters are *NOT* equal, but {np.round(close, 4) * 100}% \
            of the affected area is close (relative tolerance: 10^-4) \n'
        log.info(message)
        with open(reportFile, 'a') as pfile:
           pfile.write(message)