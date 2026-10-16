#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Run script to execute com4FlowPy benchmarks, utilizing maximum available
    CPU cores to evaluate performance, and reporting realTime and peak RAM usage.
    Produces a Markdown report under tests/reportsCom4FlowPy/
"""

import time
import pathlib
import numpy as np
from datetime import datetime
import tempfile
import os
import threading
import psutil

# Local imports
import avaframe as avaf
from avaframe.com4FlowPy import com4FlowPy
from avaframe.runCom4FlowPy import readFlowPyinputs
from avaframe.ana1Tests import testUtilities as tU
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in2Trans.rasterUtils as rasterUtils


class PeakMemoryTracker:
    """Recursively tracks memory of parent and all child processes to capture peak RSS."""
    def __init__(self, interval=0.01):
        self.interval = interval
        self.peakMem = 0.0
        self.running = False
        self._thread = None

    def _track(self):
        parent = psutil.Process()
        while self.running:
            try:
                # Add memory of parent process
                totalMem = parent.memory_info().rss
                # Add memory of all children recursively
                for child in parent.children(recursive=True):
                    try:
                        totalMem += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                
                totalMemMB = totalMem / (1024 * 1024)
                if totalMemMB > self.peakMem:
                    self.peakMem = totalMemMB
            except Exception:
                pass
            time.sleep(self.interval)

    def start(self):
        self.peakMem = 0.0
        self.running = True
        self._thread = threading.Thread(target=self._track, daemon=True)
        self._thread.start()

    def stop(self):
        self.running = False
        if self._thread:
            self._thread.join(timeout=1.0)
        return self.peakMem


def _checkNumbaInstalled() -> bool:
    try:
        import numba
        return True
    except ImportError:
        return False


def compareRasters(path, pathRef):
    """Compare two rasters and compute the difference between them."""
    rasterDict = rasterUtils.readRaster(path, noDataToNan=False)
    raster = rasterDict["rasterData"]
    rasterRefDict = rasterUtils.readRaster(pathRef, noDataToNan=False)
    rasterRef = rasterRefDict["rasterData"]
    
    diff = rasterRef - raster
    equal = np.array_equal(rasterRef, raster)

    closeArray = np.isclose(raster, rasterRef, rtol=1e-04, equal_nan=True)
    mask = np.logical_or(raster > 0, rasterRef > 0)
    num_close = np.count_nonzero(closeArray[mask])
    total = rasterRef[mask].size
    closePercent = num_close / total if total > 0 else 1.0

    return diff, equal, closePercent


def main():
    # avaframe directory
    _avaframeDir  = pathlib.Path(avaf.__file__).parents[0]
    _benchmarkDir = pathlib.Path(_avaframeDir, '..', 'benchmarks')

    logName = 'runBenchmarksRAMcom4FlowPy'

    # Load settings from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()

    # Load all benchmark info as dictionaries from description files
    testDictList = tU.readAllBenchmarkDesDicts(info=False, inDir=_benchmarkDir)

    # Filter benchmarks for tag standardTest
    filterType = 'TAGS'
    valuesList = ['standardTest', 'com4FlowPy']
    testListAll = tU.filterBenchmarks(testDictList, filterType, valuesList, condition='and')

    valuesList = ['standardTest', 'com4FlowPy', 'numba']
    testListNumba = tU.filterBenchmarks(testDictList, filterType, valuesList, condition='and')

    if _checkNumbaInstalled():
        testList = testListAll
    else:
        testList = [item for item in testListAll if item not in testListNumba]

    # Set directory for benchmark report
    outDir = _avaframeDir / 'tests' / 'reportsCom4FlowPy'
    fU.makeADir(outDir)

    reportFile = outDir / 'benchmarksRAMcom4FlowPyReport.md'

    _startDate = datetime.now()
    maxCores = psutil.cpu_count(logical=True) or 1

    log = logUtils.initiateLogger(outDir, logName)
    log.info('Starting Com4FlowPy RAM Performance Benchmarks')

    benchmark_summary = []

    # Create temporary directory for isolated runs
    with tempfile.TemporaryDirectory(prefix="avaframe_ramBenchmarks_") as tempDir:
        tmpTestsDir = pathlib.Path(tempDir)

        for i, test in enumerate(testList):
            log.info(f"Running benchmark {i+1}/{len(testList)}: {test['NAME']}")

            # Define avaDir relative to _avaframeDir
            avaDir = str(_avaframeDir / pathlib.Path(test['AVADIR']))
            cfgMain['MAIN']['avalancheDir'] = avaDir

            # Fetch benchmark test info
            refDir = pathlib.Path(_avaframeDir, '..', 'benchmarks', test['NAME'])

            # Clean input directory of old work and output files
            initProj.cleanSingleAvaDir(_avaframeDir / avaDir, deleteOutput=False)

            # Load input parameters from configuration file for standard tests
            benchmarkCfg = refDir / ('%s' % test['INI'])
            cfg = cfgUtils.getModuleConfig(com4FlowPy, fileOverride=benchmarkCfg)
            cfgGen = cfg["GENERAL"]
            
            # Explicitly force CPU count to the maximum available logical cores to target multiprocessing
            cfgGen["cpuCount"] = str(maxCores)
            cfgGen["overwriteResults"] = "reRunAndOverwrite"

            avalancheDir = cfgMain["MAIN"]["avalancheDir"]
            cfgPath = readFlowPyinputs(avalancheDir, cfg, log)

            avaDirTempOutPut = test['AVADIR']
            compDir = tmpTestsDir / pathlib.Path(avaDirTempOutPut)
            
            cfgPath["customDirs"] = False
            cfgPath["resDir"] = compDir
            fU.makeADir(cfgPath["resDir"])
            cfgPath["thalwegDir"] = cfgPath["resDir"] / "thalwegData"
            cfgPath["tempDir"] = cfgPath["workDir"] / "temp"
            fU.makeADir(cfgPath["tempDir"])
            cfgPath["deleteTemp"] = "False"
            cfgPath["outputFiles"] = cfg["PATHS"]["outputFiles"]
            cfgPath["outputNoDataValue"] = cfg["PATHS"].getfloat("outputNoDataValue")
            cfgPath["useCompression"] = cfg["PATHS"].getboolean("useCompression")
            cfgPath["uid"] = cfgUtils.cfgHash(cfg)
            cfgPath["timeString"] = datetime.now().strftime("%Y%m%d_%H%M%S")

            # Setup memory tracker and timer
            tracker = PeakMemoryTracker(interval=0.01)
            tracker.start()
            
            startTime = time.perf_counter()
            com4FlowPy.com4FlowPyMain(cfgPath, cfgGen)
            endTime = time.perf_counter()
            
            peakMemUsed = tracker.stop()
            timeNeeded = endTime - startTime

            log.info(f"Test {test['NAME']} finished. Time: {timeNeeded:.2f}s, Peak RAM: {peakMemUsed:.2f} MB")

            # Compare generated files to benchmark references to verify correctness
            outputVariable = test['FILES']
            comparison_results = []

            for variable in outputVariable:
                pathRasterRef = None
                for file in os.listdir(refDir):
                    if file.endswith('%s.tif' % variable):
                        pathRasterRef = refDir / file
                        break
                
                pathRaster = None
                for file in os.listdir(compDir):
                    if file.endswith('%s.tif' % variable) and (str(file).split('_')[1] == cfgPath["uid"]):
                        pathRaster = compDir / file
                        break

                if pathRasterRef and pathRaster and os.path.isfile(pathRasterRef) and os.path.isfile(pathRaster):
                    diff, eq, close = compareRasters(pathRaster, pathRasterRef)
                    if eq and np.sum(abs(diff[diff != 0])) == 0:
                        status = "rasters are equal |&check;"
                    else:
                        status = f"rasters are __NOT(!)__ equal - {np.round(close, 4) * 100}% of the affected area is close |&cross;"
                else:
                    status = "Missing File |&cross;"
                
                comparison_results.append((variable, status))

            benchmark_summary.append({
                "name": test['NAME'],
                "description": test['DESCRIPTION'],
                "reference": test.get('REFERENCE', 'N/A'),
                "benchmarked_version": test.get('BENCHMARKED_AVAFRAME_VERSION', 'N/A'),
                "time": timeNeeded,
                "peak_ram": peakMemUsed,
                "comparisons": comparison_results
            })

    # Write summary table and details to markdown report in the exact style of the standardTests report
    _endDate = datetime.now()
    testAvaFVersion = avaf.version.getVersion()

    with open(reportFile, 'w') as pfile:
        pfile.write('# Com4FlowPy Performance and RAM Benchmark Report \n\n')
        pfile.write('Profiling and comparing __com4FlowPy__ simulations to selected benchmark results with maximum CPU cores.\n\n')
        pfile.write('* * * \n')
        pfile.write(f'**CPU Cores Utilized:** {maxCores} &check;  \n')
        pfile.write(f'**tested AvaFrame Version:** {testAvaFVersion}  \n')
        if _checkNumbaInstalled():
            pfile.write('`numba` __found__: running all benchmarks for python and numba engines &check;\n')
        else:
            pfile.write('`numba` __NOT found__: skipping benchmarks for numba engine &cross;\n')
        pfile.write('* * * \n\n')

        pfile.write('## Performance Summary Table\n\n')
        pfile.write('| Benchmark Name | Elapsed Time (s) | Peak RAM Usage (MB) | Status |\n')
        pfile.write('| :--- | :---: | :---: | :--- |\n')
        for summary in benchmark_summary:
            status_emoji = "&check; Pass"
            for _, stat in summary["comparisons"]:
                if "&cross;" in stat:
                    status_emoji = "&cross; Fail"
                    break
            pfile.write(f'| **{summary["name"]}** | {summary["time"]:.2f}s | {summary["peak_ram"]:.2f} MB | {status_emoji} |\n')
        pfile.write('\n* * *\n\n')

        # Write test details in the exact style of standardTestsReportCom4FlowPy.md
        pfile.write('## Detailed Benchmark Runs\n\n')
        for i, test in enumerate(benchmark_summary):
            pfile.write(f"### Test {i+1}: *{test['name']}*\n")
            for descLine in str(test['description']).split('\n'):
                pfile.write(f"{descLine}{' ' * 2}\n")
            pfile.write(f"__Reference__: {test['reference']}{' ' * 2}\n")
            pfile.write(f"__tested AvaFrame Version__: {testAvaFVersion}{' ' * 2}\n")
            
            benchAvaFVersion = test['benchmarked_version']
            if benchAvaFVersion == testAvaFVersion:
                pfile.write(f"__benchmarked AvaFrame::com4FlowPy version__: {benchAvaFVersion}{' ' * 2}\n")
            else:
                pfile.write(f"__benchmarked AvaFrame::com4FlowPy version__: <span style=\"color:red\"> {benchAvaFVersion} </span>{' ' * 2}\n")

            pfile.write(f"__Elapsed Time__: {test['time']:.2f}s{' ' * 2}\n")
            pfile.write(f"__Peak RAM Usage__: {test['peak_ram']:.2f} MB{' ' * 2}\n\n")

            pfile.write("|Model Output|Result of comparison|status\n")
            pfile.write("|----:|:-----:|:---:|\n")
            for variable, status in test["comparisons"]:
                pfile.write(f"|__{variable}__| {status}\n")
            pfile.write('\n * * * \n\n')

        pfile.write(f"__benchmark(s) finished @__: {_endDate}\n")
        pfile.write(f"__timeDelta__: {_endDate - _startDate}\n")

    log.info(f"Benchmarking complete. Report written to {reportFile}")


if __name__ == "__main__":
    main()
